
INPUT='input'
OUTPUT='output'

SHARDING = 16

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work', **kwargs):
    if base == 'input':
        base_dir = 'input'
    elif base == 'output':
        base_dir = 'output'
    elif base == 'work':
        base_dir = 'output/intermediate_results'
    else:
        raise('Base not found')
    return Path(base_dir.format_map(Default(kwargs)))

input_dir = 'input'
output_dir = 'output'
work_dir = 'output/intermediate_results_dir'

configfile: 'config/deepvariant.yaml'

rule all:
    input:
        'asm.output.vcf.gz'

rule pbmm2_align:
    input:
        ref = config['reference'],
        reads = config['reads']
    output:
        temp(input_dir + '{haplotype}_hifi_reads.unphased.pbmm2.bam')
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        'pbmm2 align {input.ref} {input.reads} {output} --sort --preset CCS -j {threads}'

rule deepvariant_make_examples:
    input:
        ref = '',
        bam = '{haplotype}_hifi_reads.{phase}.pbmm2.bam'
    output:
        temp(work_dir / 'make_examples.tfrecord-{N}-of-{sharding}.gz')
    params:
        examples = '/output/intermediate_results_dir/make_examples.tfrecord@{config[shards]}.gz',
        gvcf = '/output/intermediate_results_dir/gvcf.tfrecord@{config[shards]}.gz'
        phase_args = lambda wildcards: '--{i}parse_same_aux_fields --{i}sort_by_haplotypes'.format(i='no' if wildards.phase == 'unphased' else '')
    threads: 1
    resources:
        mem_mb = 4000
    shell:
        '''
        singularity exec --no-home --cleanenv --containall -B tmp:/tmp,input:/input,output:/output \
        deepvariant_1.1.0.sif \
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
        (work_dir / f'make_examples.tfrecord-{N:05}-of-{config["shards"]:05}.gz' for N in range(config['shards']))
    output:
        temp('output/intermediate_results_dir/call_variants_output.tfrecord.gz')
    params:
        examples = '/output/intermediate_results_dir/make_examples.tfrecord@{config[shards]}.gz',
    threads: 16
    resources:
        mem_mb = 4000
    shell:
        '''
        singularity exec --no-home --cleanenv --containall --env OMP_NUM_THREADS={threads} -B tmp:/tmp,input:/input,output:/output \
        deepvariant_1.1.0.sif \
        /bin/bash -c "cd /output; ../opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --checkpoint /opt/models/pacbio/model.ckpt \
        --use_openvino"
        '''

rule deepvaraint_postprocess:
    input:
        ref = config['reference'],
        variants = 'output/intermediate_results_dir/call_variants_output.tfrecord.gz',
        gvcf = (work_dir / f'gvcf.tfrecord-{N:05}-of-{config["shards"]:05}.gz' for N in range(config['shards']))
    output:
        vcf = 'output/asm.output.vcf.gz',
        gvcf = 'test.output.g.vcf.gz'
    params:
        lambda wildcards, input: ('/'+record for record in input['gvcf'])
    shell:
        '''
        singularity exec --no-home --cleanenv --containall -B tmp:/tmp,input:/input,output:/output \
        deepvariant_1.1.0.sif \
        /opt/deepvariant/bin/postprocess_variants \
        --ref /input/asm.fasta \
        --infile /{input.variants} \
        --outfile /{output.vcf} \
        --gvcf_outfile /{output.gvcf} \
        --nonvariant_site_tfrecord_path {input.gvcf}
        '''
