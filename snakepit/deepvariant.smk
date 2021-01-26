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

include: 'whatshap.smk'

localrules: samtools_faidx

def make_singularity_call(extra_args='',tmp_bind='tmp',input_bind=get_dir('input'),output_bind=get_dir('output')):
    return f'singularity exec --no-home --cleanenv --containall {extra_args} -B {tmp_bind}:/tmp,{input_bind}:/input,{output_bind}:/output'

for dir in ('input',):
    Path(get_dir(dir)).mkdir(exist_ok=True)

rule all:
    input:
        get_dir('output',f'asm.{config["phased"]}.vcf.gz',model=config['model'])

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
        examples = lambda wildcards, output: PurePath(output[0]).with_name('make_examples.tfrecord@{config[shards]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[0]).with_name('gvcf.tfrecord@{config[shards]}.gz'),
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--{i}parse_same_aux_fields --{i}sort_by_haplotypes'.format(i=('no' if wildcards.phase == 'unphased' else '')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --norealign_reads --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = make_singularity_call()
    threads: 1
    resources:
        mem_mb = 4000,
        use_singularity = True
    shell:
        '''
        mkdir -p {params.dir_}
        {params.singularity_call} \
        {config[container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref /{input.ref} \
        --reads /{input.bam[0]} \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output.tfrecord.gz'))
    params:
        examples = get_dir('work','make_examples.tfrecord@{config[shards]}.gz'),
        model = lambda wildcards: 'pacbio' if wildcards.model == 'pbmm2' else 'hybrid_pacbio_illumina',
        singularity_call = lambda wildcards, threads: make_singularity_call('--env OMP_NUM_THREADS={threads}')
    threads: 16
    resources:
        mem_mb = 4000,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[container]} \
        /bin/bash -c "cd /output; ../opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --checkpoint /opt/models/{params.model}/model.ckpt \
        --use_openvino"
        '''

rule deepvaraint_postprocess:
    input:
        ref = config['reference'],
        variants = get_dir('work','call_variants_output.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{haplotype}.{phase}.vcf.gz'),
        gvcf = get_dir('output','{haplotype}.{phase}.g.vcf.gz')
    params:
        lambda wildcards, input: ('/'+record for record in input['gvcf']),
        singularity_call = make_singularity_call()
    threads: 1
    resources:
        mem_mb = 30000,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[container]} \
        /opt/deepvariant/bin/postprocess_variants \
        --ref /input/asm.fasta \
        --infile /{input.variants} \
        --outfile /{output.vcf} \
        --gvcf_outfile /{output.gvcf} \
        --nonvariant_site_tfrecord_path {input.gvcf}
        '''


rule deeptrio_make_examples:
    input:
        ref = multiext(config['reference'],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{model}.bam'),'','.bai'),
        bam_p1 = multiext(get_dir('input','{parent1}.{phase}.{model}.bam'),'','.bai'),
        bam_p2 = multiext(get_dir('input','{parent2}.{phase}.{model}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','DT_make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','DT_gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name('make_examples.tfrecord@{config[shards]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[0]).with_name('gvcf.tfrecord@{config[shards]}.gz'),
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--{i}parse_same_aux_fields --{i}sort_by_haplotypes'.format(i=('no' if wildcards.phase == 'unphased' else '')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --norealign_reads --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = make_singularity_call()
    threads: 1
    resources:
        mem_mb = 4000,
        use_singularity = True
    shell:
        '''
        mkdir -p {params.dir_}
        {params.singularity_call} \
        {config[container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref /{input.ref} \
        --reads /{input.bam} \
        --reads_parent1 /{input.bam_p1} \
        --reads_parent2 /{input.bam_p2} \
        --sample_name offspring \
        --sample_name_parent1 sire \
        --sample_name_parent2 dam \
        --examples {params.examples} \
        --gvcf {params.gvcf} \
        {params.model_args} \
        {params.phase_args} \
        --task {wildcards.N}
        '''
