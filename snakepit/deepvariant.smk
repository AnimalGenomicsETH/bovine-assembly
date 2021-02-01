from pathlib import PurePath

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = 'input'
    elif base == 'output':
        base_dir = 'output_{model}'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{haplotype}_{phase}',**kwargs)
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext)

if Path('config/deepvariant.yaml').exists():
    configfile: 'config/deepvariant.yaml'

wildcard_constraints:
     subset = r'|_child|_parent1|_parent2',
     haplotype = r'asm',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid'

def get_model(wildcards,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    if wildcards['subset'] == '':
        if wildcards['model'] == 'pbmm2':
            return model_location.format('pacbio')
        elif wildcards['model'] == 'hybrid':
            return model_location.format('hybrid_pacbio_illumina')
    else:
        if wildcards['model'] == 'pbmm2':
            if wildcards['subset'] == 'child':
                return model_location.format(f'deeptrio/pacbio/child')
            else:
                return model_location.format(f'deeptrio/pacbio/parent')

include: 'whatshap.smk'

localrules: samtools_faidx

def make_singularity_call(wildcards,extra_args='',tmp_bind='tmp',input_bind=None, output_bind=None, work_bind=None):
    return f'singularity exec {extra_args} -B $TMPDIR:/tmp,{input_bind or get_dir("input",**wildcards)}:/input,{output_bind or get_dir("output",**wildcards)}:/output,{work_bind or get_dir("work",**wildcards)}:/output/intermediate'

for dir in ('input',):
    Path(get_dir(dir)).mkdir(exist_ok=True)

rule all:
    input:
        get_dir('output','asm.trio.merged.unphased.vcf.gz',model='pbmm2'),
        #get_dir('output',f'asm.{config["phased"]}.vcf.gz',model=config['model'])

rule minimap_align:
    input:
        ref = config['reference'],
        reads = config['short_reads']
    output:
        temp(get_dir('input','{haplotype}.unphased.mm2.bam'))
    threads: 24
    resources:
        mem_mb = 5000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax sr -t {threads} {input.ref} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule pbmm2_align:
    input:
        ref = config['reference'],
        reads = config['long_reads']
    output:
        temp(get_dir('input','{haplotype}.unphased.pbmm2.bam'))
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        'pbmm2 align {input.ref} {input.reads} {output} --sort --preset CCS -j {threads}'

rule merge_hybrid:
    input:
        expand(get_dir('input','{{haplotype}}.unphased.{model}.bam'),model=('mm2','pbmm2'))
    output:
        temp(get_dir('input','{haplotype}.unphased.hybrid.bam'))
    threads: 8
    resources:
        mem_mb = 3000
    shell:
        'samtools merge -@ {threads} {output} {input}'

rule samtools_index:
    input:
        '{bam}.bam'
    output:
        '{bam}.bam.bai'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule samtools_faidx:
    output:
        f'{config["reference"]}.fai'
    shell:
        'samtools faidx {input}'

rule deepvariant_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{model}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = f'/output/intermediate/make_examples.tfrecord@{config["shards"]}.gz',
        gvcf = f'/output/intermediate/gvcf.tfrecord@{config["shards"]}.gz',
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--parse_sam_aux_fields={i} --sort_by_haplotypes={i}'.format(i=('true' if wildcards.phase == 'phased' else 'false')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --realign_reads=false --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = lambda wildcards: make_singularity_call(wildcards)
    threads: 1
    resources:
        mem_mb = 4000,
        disk_scratch = 10,
        use_singularity = True
    shell:
        '''
        mkdir -p {params.dir_}
        {params.singularity_call} \
        {config[container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref /{input.ref[0]} \
        --reads /{input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples{{subset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output{subset}.tfrecord.gz'))
    params:
        outfile = lambda wildcards, output: '/output/intermediate/' + PurePath(output[0]).name,
        examples = f'/output/intermediate/make_examples.tfrecord@{config["shards"]}.gz',
        model = lambda wildcards: get_model(wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['container'] if wildcards.subset == '' else config['DT_container']
    threads: 16
    resources:
        mem_mb = 4000,
        disk_scratch = 10,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /bin/bash -c "cd /output; ../opt/deepvariant/bin/call_variants \
        --outfile {params.outfile} \
        --examples {params.examples} \
        --checkpoint {params.model} \
        --use_openvino"
        '''

rule deepvariant_post_postprocess:
    input:
        ref = config['reference'],
        variants = get_dir('work','call_variants_output{subset}.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf{{subset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{haplotype}{subset}.{phase}.vcf.gz'),
        gvcf = get_dir('output','{haplotype}{subset}.{phase}.g.vcf.gz')
    params:
        variants = lambda wildcards, input: '/output/intermediate/' + PurePath(input['variants']).name,
        gvcf = (f'/output/intermediate/gvcf.tfrecord-{N:05}-of-{config["shards"]:05}.gz' for N in range(config['shards'])),
        vcf_out = lambda wildcards, output: '/output/' + PurePath(output['vcf']).name,
        gvcf_out = lambda wildcards, output: '/output/' + PurePath(output['gvcf']).name,
        singularity_call = lambda wildcards: make_singularity_call(wildcards),
        contain = lambda wildcards: config['container'] if wildcards.subset == '' else config['DT_container']
    threads: 1
    resources:
        mem_mb = 30000,
        disk_scratch = 10,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/postprocess_variants \
        --ref /{input.ref} \
        --infile {params.variants} \
        --outfile {params.vcf_out} \
        --gvcf_outfile {params.gvcf_out} \
        --nonvariant_site_tfrecord_path {params.gvcf}
        '''


rule deeptrio_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{model}.bam'),'','.bai'),
        bam_p1 = multiext(get_dir('input','sire.{phase}.{model}.bam'),'','.bai'),
        bam_p2 = multiext(get_dir('input','dam.{phase}.{model}.bam'),'','.bai')
    output:
        example = temp(get_dir('work',f'make_examples{S}.tfrecord-{{N}}-of-{{sharding}}.gz') for S in ('_child','_parent1','_parent2')),
        gvcf = temp(get_dir('work',f'gvcf{S}.tfrecord-{{N}}-of-{{sharding}}.gz') for S in ('_child','_parent1','_parent2'))
    params:
        examples = f'/output/intermediate/make_examples.tfrecord@{config["shards"]}.gz',
        gvcf = f'/output/intermediate/gvcf.tfrecord@{config["shards"]}.gz',
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--{i}parse_sam_aux_fields --{i}sort_by_haplotypes'.format(i=('no' if wildcards.phase == 'unphased' else '')),
        model_args = lambda wildcards: '--alt_aligned_pileup diff_channels --norealign_reads --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = lambda wildcards, output: make_singularity_call(wildcards)
    threads: 1
    resources:
        mem_mb = 4000,
        disk_scratch = 10,
        use_singularity = True
    shell:
        '''
        mkdir -p {params.dir_}
        {params.singularity_call} \
        {config[DT_container]} \
        /opt/deepvariant/bin/deeptrio/make_examples \
        --mode calling \
        --ref /{input.ref[0]} \
        --reads /{input.bam[0]} \
        --reads_parent1 /{input.bam_p1[0]} \
        --reads_parent2 /{input.bam_p2[0]} \
        --sample_name offspring \
        --sample_name_parent1 sire \
        --sample_name_parent2 dam \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --pileup_image_height_child 100 \
        --pileup_image_height_parent 100 \
        --task {wildcards.N}
        '''

rule deeptrio_merge:
    input:
        (get_dir('output',f'{{haplotype}}{S}.{{phase}}.g.vcf.gz') for S in ('_child','_parent1','_parent2'))
    output:
        get_dir('output','{haplotype}.trio.merged.{phase}.vcf.gz')
    params:
        singularity_call = lambda wildcards, output: make_singularity_call(wildcards)
    threads: 12
    resources:
        mem_mb = 30000,
        disk_scratch = 10,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[GL_container]} \
        quay.io/mlin/glnexus:v1.2.7 \
        /usr/local/bin/glnexus_cli \
        --config DeepVariantWGS \
        --threads {threads} \
        {input} \
        | bcftools view - | bgzip -c > {output}
        '''
        #--trim-uncalled-alleles
