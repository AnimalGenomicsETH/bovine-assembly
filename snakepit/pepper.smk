from pathlib import PurePath
from itertools import product

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

wildcard_constraints:
     subset = r'|_child|_parent1|_parent2',
     haplotype = r'asm|hap1|hap2|parent1|parent2',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid|bwa'

def get_model(wildcards,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    if wildcards['subset'] == '':
        if wildcards['model'] == 'pbmm2':
            return model_location.format('pacbio')
        elif wildcards['model'] == 'hybrid':
            return model_location.format('hybrid_pacbio_illumina')
        elif wildcards['model'] == 'bwa':
            return model_location.format('wgs')
    else:
        if wildcards['model'] == 'pbmm2':
            if wildcards['subset'] == 'child':
                return model_location.format(f'deeptrio/pacbio/child')
            else:
                return model_location.format(f'deeptrio/pacbio/parent')

rule extract_haplotags:
    input:
        ''
    output:
        ''
    shell:
        '''
        zgrep -o ">\S*" {input.hap1} | cut -c 2- | awk '{print $0" H1}' > {output}
        zgrep -o ">\S*" {input.hap2} | cut -c 2- | awk '{print $0" H2}' > {output}
        '''

rule deepvariant_make_examples:
    input:
        ref = lambda wildcards: multiext(config['references'][wildcards.ref],'','.fai'),
        bam = multiext(get_dir('input','{haplotype}.{phase}.{ref}.{model}.bam'),'','.bai'),
        vcf = get_dir('work','PEPPER_HP_OUTPUT_1.vcf.gz')
    output:
        example = temp(get_dir('work','make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name(f'make_examples.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[1]).with_name(f'gvcf.tfrecord@{config["shards"]}.gz'),
        phase_args = lambda wildcards: '--parse_sam_aux_fields={i} --sort_by_haplotypes={i}'.format(i=('true' if wildcards.phase == 'phased' else 'false')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --realign_reads=false --vsc_min_fraction_indels 0.12' if wildcards.model == 'pbmm2' else '',
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}'
    threads: 1
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 1,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        /opt/deepvariant/bin/make_examples \
        --mode calling \
        --ref {params.ref} \
        --reads {input.bam[0]} \
        --examples {params.examples} \
        --proposed_variants {input.vcf} \
        --sample_name {wildcards.haplotype} \
        --noadd_hp_channel \
        --alt_aligned_pileup "rows" \
        --norealign_reads \
        --hp_tag_for_assembly_polishing "1" \
        --parse_sam_aux_fields \
        --sort_by_haplotypes \
        --variant_caller "vcf_candidate_importer" \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('work', f'make_examples{{subset}}.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('work','call_variants_output{subset}.tfrecord.gz'))
    params:
        examples = lambda wildcards: (f'make_examples.tfrecord@{config["shards"]}.gz'),
        model = lambda wildcards: get_model(wildcards),
        dir_ = lambda wildcards: get_dir('work',**wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container']
    threads: 32
    resources:
        mem_mb = 1500,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '4:00' if wildcards.subset == '' else '24:00'
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --model_type=WGS \
        --customized_model=/opt/dv_models/202012_polish_nohp_rows/model.ckpt-28400
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
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container'] if wildcards.subset == '' else config['DT_container']
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
        --ref {params.ref} \
        --infile {input.variants} \
        --outfile {output.vcf} \
        --gvcf_outfile {output.gvcf} \
        --nonvariant_site_tfrecord_path {params.gvcf}
        '''
