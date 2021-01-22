from pathlib import PurePath

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = 'input'
    elif base == 'output':
        base_dir = 'output'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results_{haplotype}_{phase}')
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext)

configfile: 'config/deepvariant.yaml'

include: 'whatshap.smk'

def make_singularity_call(extra_args='',tmp_bind='tmp',input_bind=get_dir('input'),output_bind=get_dir('output')):
    return f'singularity exec --no-home --cleanenv --containall {extra_args} -B {tmp_bind}:/tmp,{input_bind}:/input,{output_bind}:/output'

for dir in ('input','output'):
    Path(get_dir(dir)).mkdir(exist_ok=True)

rule all:
    input:
        get_dir('output','asm.phased.vcf.gz')

rule pbmm2_align:
    input:
        ref = config['reference'],
        reads = config['reads']
    output:
        temp(get_dir('input','{haplotype}_hifi_reads.unphased.pbmm2.bam'))
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        'pbmm2 align {input.ref} {input.reads} {output} --sort --preset CCS -j {threads}'

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

rule deepvariant_make_examples:
    input:
        ref = config['reference'],
        bam = get_dir('input','{haplotype}_hifi_reads.{phase}.pbmm2.bam'),
        bai = get_dir('input','{haplotype}_hifi_reads.{phase}.pbmm2.bam.bai')
    output:
        example = temp(get_dir('work','make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name('make_examples.tfrecord@{config[shards]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[0]).with_name('gvcf.tfrecord@{config[shards]}.gz'),
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        phase_args = lambda wildcards: '--{i}parse_same_aux_fields --{i}sort_by_haplotypes'.format(i=('no' if wildcards.phase == 'unphased' else '')),
        singularity_call = make_singularity_call()
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        '''
        mkdir -p {params.dir_}
        {params.singularity_call} \
        {config[container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref /{input.ref} \
        --reads /{input.bam} \
        --examples {params.examples} \
        --add_hp_channel \
        --alt_aligned_pileup diff_channels \
        --gvcf {params.gvcf} \
        --norealign_reads \
        {params.phase_args} \
        --vsc_min_fraction_indels 0.12 \
        --task {wildcards.N}"
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output.tfrecord.gz'))
    params:
        examples = get_dir('work','make_examples.tfrecord@{config[shards]}.gz'),
        singularity_call = lambda wildcards, threads: make_singularity_call('--env OMP_NUM_THREADS={threads}')
    threads: 16
    resources:
        mem_mb = 4000
    shell:
        '''
        {params.singularity_call} \
        {config[container]} \
        /bin/bash -c "cd /output; ../opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --checkpoint /opt/models/pacbio/model.ckpt \
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
