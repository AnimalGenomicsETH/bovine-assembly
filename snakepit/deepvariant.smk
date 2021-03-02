from pathlib import PurePath

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = 'input_{animal}'
    elif base == 'output':
        base_dir = 'output_{animal}'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{haplotype}_{phase}_{model}_{ref}',**kwargs)
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

#if Path('config/deepvariant.yaml').exists():
#    configfile: 'config/deepvariant.yaml'

wildcard_constraints:
     subset = r'|_child|_parent1|_parent2',
     haplotype = r'asm|hap1|hap2|parent1|parent2',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid|wgs'

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

localrules: samtools_faidx, rtg_pedigree, rtg_format, rtg_mendelian_concordance

##TODO determine how to pull singularity images
##TODO merge calls with snakemake format

def make_singularity_call(wildcards,extra_args='',tmp_bind='$TMPDIR',input_bind=True, output_bind=True, work_bind=True):
    call = f'singularity exec {extra_args} '
    if input_bind:
        call += f'-B {get_dir("input",**wildcards)}:/input '
    if output_bind:
        call += f'-B {get_dir("output",**wildcards)}:/output '
    if work_bind:
        call += f'-B {get_dir("work",**wildcards)}:/output/intermediate '
    if tmp_bind:
        call += f'-B {tmp_bind}:/tmp '
    return call

for dir in ('input','output','work'):
    Path(get_dir(dir,**config)).mkdir(exist_ok=True)

rule all:
    input:
        get_dir('output','{haplotype}.{phase}.{reference}.{model}.vcf.gz',**config) if not config['trio'] else [],
        multiext(get_dir('output','{haplotype}.trio.merged.{phase}.{reference}.{model}',**config),'.vcf.gz','.mendelian.log') if config['trio'] else []

rule minimap_align:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref],
        reads = lambda wildcards: config['short_reads']['individual' if 'parent' not in wildcards.haplotype else wildcards.haplotype]
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.mm2.bam'))
    threads: 24
    resources:
        mem_mb = 5000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax sr -t {threads} {input.ref} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule pbmm2_align:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref],
        reads = lambda wildcards: config['long_reads']['individual' if 'parent' not in wildcards.haplotype else wildcards.haplotype]
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.pbmm2.bam'))
    threads: 12
    resources:
        mem_mb = 3500,
        walltime = '24:00'
    shell:
        'pbmm2 align {input.ref} {input.reads} {output} --sort --preset CCS -j {threads}'

rule merge_hybrid:
    input:
        (get_dir('input','{haplotype}.unphased.{ref}.{model}.bam',model=X) for X in ('mm2','pbmm2'))
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.hybrid.bam'))
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
        mem_mb = 4000,
        walltime = '30'
    shell:
        'samtools index -@ {threads} {input}'

rule samtools_faidx:
    input:
        '{reference}'
    output:
        '{reference}.fai'
    shell:
        'samtools faidx {input}'

rule deepvariant_make_examples:
    input:
        ref = lambda wildcards: multiext(config['references'][wildcards.ref],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{ref}.{model}.bam'),'','.bai')
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
        mem_mb = 6000,
        disk_scratch = 1,
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
        examples = lambda wildcards: f'/output/intermediate/make_examples{wildcards.subset}.tfrecord@{config["shards"]}.gz',
        model = lambda wildcards: get_model(wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['container'] if wildcards.subset == '' else config['DT_container'],
        vino = lambda wildcards: '--use_openvino' if wildcards.subset == '' else ''
    threads: 24
    resources:
        mem_mb = 1500,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '4:00' if wildcards.subset == '' else '24:00'
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /bin/bash -c "cd /output; ../opt/deepvariant/bin/call_variants \
        --outfile {params.outfile} \
        --examples {params.examples} \
        --checkpoint {params.model} \
        {params.vino}"
        '''

rule deepvariant_postprocess:
    input:
        ref = lambda wildcards: multiext(config['references'][wildcards.ref],'','.fai'),
        variants = get_dir('work','call_variants_output{subset}.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf{{subset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{haplotype}{subset}.{phase}.{ref}.{model}.vcf.gz'),
        gvcf = get_dir('output','{haplotype}{subset}.{phase}.{ref}.{model}.g.vcf.gz')
    params:
        variants = lambda wildcards, input: '/output/intermediate/' + PurePath(input['variants']).name,
        gvcf = lambda wildcards: f'/output/intermediate/gvcf{wildcards.subset}.tfrecord@{config["shards"]}.gz',
        vcf_out = lambda wildcards, output: '/output/' + PurePath(output['vcf']).name,
        gvcf_out = lambda wildcards, output: '/output/' + PurePath(output['gvcf']).name,
        singularity_call = lambda wildcards: make_singularity_call(wildcards),
        contain = lambda wildcards: config['container'] if wildcards.subset == '' else config['DT_container']
    threads: 1
    resources:
        mem_mb = 30000,
        disk_scratch = 1,
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
        ref = lambda wildcards: multiext(config['references'][wildcards.ref],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{ref}.{model}.bam'),'','.bai'),
        bam_p1 = multiext(get_dir('input','parent1.{phase}.{ref}.{model}.bam'),'','.bai'),
        bam_p2 = multiext(get_dir('input','parent2.{phase}.{ref}.{model}.bam'),'','.bai')
    output:
        example = temp(get_dir('work',f'make_examples{S}.tfrecord-{{N}}-of-{{sharding}}.gz') for S in ('_child','_parent1','_parent2')),
        gvcf = temp(get_dir('work',f'gvcf{S}.tfrecord-{{N}}-of-{{sharding}}.gz') for S in ('_child','_parent1','_parent2'))
    params:
        examples = f'/output/intermediate/make_examples.tfrecord@{config["shards"]}.gz',
        gvcf = f'/output/intermediate/gvcf.tfrecord@{config["shards"]}.gz',
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--{i}parse_sam_aux_fields --{i}sort_by_haplotypes'.format(i=('no' if wildcards.phase == 'unphased' else '')),
        model_args = lambda wildcards: '--alt_aligned_pileup diff_channels --norealign_reads --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = lambda wildcards: make_singularity_call(wildcards)
    threads: 1
    resources:
        mem_mb = 4000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = '24:00'
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
        --sample_name_parent1 parent1 \
        --sample_name_parent2 parent2 \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --pileup_image_height_child 100 \
        --pileup_image_height_parent 100 \
        --task {wildcards.N}
        '''

rule deeptrio_GLnexus_merge:
    input:
        (get_dir('output',f'{{haplotype}}{S}.{{phase}}.{{ref}}.{{model}}.g.vcf.gz') for S in ('_child','_parent1','_parent2'))
    output:
        get_dir('output','{haplotype}.trio.merged.{phase}.{ref}.{model}.vcf.gz')
    params:
        gvcfs = lambda wildcards, input: list(f'/output/{PurePath(fpath).name}' for fpath in input),
        out = lambda wildcards, output: f'/output/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        singularity_call = lambda wildcards: make_singularity_call(wildcards)
    threads: 12 #force using 4 threads for bgziping
    resources:
        mem_mb = 5000,
        disk_scratch = 50,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[GL_container]} \
        /usr/local/bin/glnexus_cli \
        --dir {params.DB} \
        --config DeepVariantWGS \
        --threads {threads} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        '''
        #--trim-uncalled-alleles

rule rtg_pedigree:
    output:
        get_dir('input','trio.ped')
    params:
        gender = 1 if config['gender'] == 'male' else 2
    shell:
        '''
        FILE={output}
cat <<EOM >$FILE
#PED format pedigree
#
#fam-id/ind-id/pat-id/mat-id: 0=unknown
#sex: 1=male; 2=female; 0=unknown
#phenotype: -9=missing, 0=missing; 1=unaffected; 2=affected
#
#fam-id ind-id pat-id mat-id sex phen
1 offspring parent1 parent2 {params.gender} 0
1 parent1 0 0 1 0
1 parent2 0 0 2 0
EOM
        '''

rule rtg_format:
    input:
        ref = lambda wildcards: multiext(config['references'][wildcards.ref],'','.fai')
    output:
        sdf = directory(get_dir('input','{ref}.sdf')) #TODO make general to match haplotype?
    params:
        singularity_call = lambda wildcards: make_singularity_call(wildcards,tmp_bind=False,output_bind=False,work_bind=False)
    shell:
        '''
        {params.singularity_call} \
        {config[RTG_container]} \
        rtg format -o {output.sdf} {input.ref[0]}
        '''

rule rtg_mendelian_concordance:
    input:
        sdf = get_dir('input','{ref}.sdf'),
        vcf = get_dir('output','{haplotype}.trio.merged.{phase}.{ref}.{model}.vcf.gz'),
        pedigree = get_dir('input','trio.ped')
    output:
        vcf = get_dir('output','{haplotype}.trio.merged.{phase}.{ref}.{model}.annotated.vcf.gz'),
        log = get_dir('output','{haplotype}.trio.merged.{phase}.{ref}.{model}.mendelian.log')
    params:
        vcf_in = lambda wildcards, input: '/output/' + PurePath(input.vcf).name,
        vcf_annotated = lambda wildcards, output: '/output/' + PurePath(output.vcf).name,
        singularity_call = lambda wildcards: make_singularity_call(wildcards,tmp_bind=False,work_bind=False)
    shell:
        '''
        {params.singularity_call} \
        {config[RTG_container]} \
        rtg mendelian -i {params.vcf_in} -o {params.vcf_annotated} --pedigree={input.pedigree} -t {input.sdf} > {output.log}
        '''
