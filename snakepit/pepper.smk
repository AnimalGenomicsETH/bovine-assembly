from pathlib import PurePath
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = 'PEPPERv4'
    elif base == 'output':
        base_dir = 'PEPPERv4'
    elif base == 'DV':
        base_dir = get_dir('output','dv_intermediate_outputs_{hap}',**kwargs)
    elif base == 'polished':
        base_dir = 'polished'
    else:
        raise Exception(f'Base: {base} was not found!')
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

def make_singularity_call(wildcards,extra_args='',tmp_bind='$TMPDIR',input_bind=True, output_bind=True, work_bind=True):
    call = f'singularity exec {extra_args} '
    if output_bind:
        call += f'-B {":/".join([get_dir("output",**wildcards)]*2)} '
    if work_bind:
        call += f'-B {":/".join([get_dir("DV",**wildcards)]*2)} '
    if tmp_bind:
        call += f'-B {tmp_bind}:/tmp '
    return call

rule all:
    input:
        expand('output_GxP/intermediate_results_hap1_unphased_pbmm2_ARS/make_examples.tfrecord-{N}-of-{sharding}.gz',N=range(config['shards']),sharding=config['shards'])


#pepper_snp call_variant -b /data/Assembly_ONT_reads.sorted.bam -f /data/Assembly.fasta -t 24 -m /opt/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl -o /data/PEPPERv4/pepper_snp/ -s gaur -w 4 -bs 64 --ont



#pepper_hp call_variant -b /data/asm.tagged.bam -f /data/Assembly.fasta -t 36 -m /opt/pepper_models/PEPPER_HP_R941_ONT_V4.pkl -o /data/PEPPERv5/pepper_hp/ -s gaur -w 4 -bs 64 --ont_asm


rule pepper_make_images:
    input:
        asm = config['assembly'],
        bam = lambda wildcards: multiext(get_dir('input','ONT.sorted.mm2.bam') if wildcard.mode == 'snp' else get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam'),'','.bai')
    output:
        directory(get_dir('output','pepper_{mode}/images'))
    params:
        pepper_mode = lambda wildcards: 'pepper_hp' if wildcard.mode == 'hp' else 'pepper_snp'
    threads: 24
    resources:
        mem_mb = 4000
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        {params.pepper_mode} make_images \
        -b {input.bam} \
        -f {input.asm} \
        -t {threads} \
        -o {output}
        '''

rule pepper_run_inference:
    input:
        get_dir('output','pepper_{mode}/images')
    output:
        directory(get_dir('output','pepper_{mode}/predictions'))
    params:
        pepper_mode = lambda wildcards: 'pepper_hp' if wildcard.mode == 'hp' else 'pepper_snp',
        model = lambda wildcards: f'/opt/pepper_models/PEPPER_{wildcard.mode.upper()}_R941_ONT_V4.pkl'
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        {params.pepper_mode} run_inference \
        -i {input.images} \
        -m {param.model} \
        -o {ouput} \
        -s {wildcards.sample} \
        -w 4 \
        -bs 64 \
        -t {threads}
        '''

rule pepper_find_candidates:
    input:
        predictions = get_dir('output','pepper_{mode}/predictions'),
        bam = get_dir('input','bam'),
        asm = config['assembly']
    output:
        get_dir('output','pepper_{mode}/PEPPER_SNP_OUPUT.vcf.gz')
    params:
        pepper_mode = lambda wildcards: 'pepper_hp' if wildcard.mode == 'hp' else 'pepper_snp',
        mode = lambda wildcards: '--ont_asm' if wildcards.mode == 'hp' else '--ont'
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        {params.pepper_mode} find_candidates \
        -i {input.predictions} \
        -b {input.bam} \
        -f {input.asm} \
        -s {wildcards.sample} \
        -o {output} \
        -t {threads} \
        {param.mode}
        '''

rule pepper_post_vcf:
    input:
        get_dir('output','pepper_snp/PEPPER_SNP_OUPUT.vcf.gz')
    output:
        get_dir('output','PEPPER_SNP_OUPUT.vcf.gz')
    shell:
        '''
        bgzip -@ 4 -c {input} > {output}
        tabix -p vcf {output}
        '''

rule pepper_margin:
    input:
        bam = '',
        asm = config['assembly'],
        vcf = get_dir('output','PEPPER_SNP_OUPUT.vcf.gz')
    output:
        get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam')
    params:
        json = '/opt/margin_dir/params/misc/allParams.ont_haplotag.json'        
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        margin phase {input.bam} {input.asm} {input.vcf} {params.json} -t {threads} -V -o {output}
        '''

rule pepper_margin_post:
    input:
        ''
    output:
        ''
    shell:
        '''
        mv PEPPERv4/*.bam PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam;
        samtools index -@24 PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam
        '''

rule pepper_hp_post:
    input:
        ''
    output:
        ''
    shell:
        '''
        
        bgzip PEPPERv4/PEPPER_HP_OUTPUT_1.vcf;
        bgzip PEPPERv4/PEPPER_HP_OUTPUT_2.vcf;
        tabix -p vcf PEPPERv4/PEPPER_HP_OUTPUT_1.vcf.gz;
        tabix -p vcf PEPPERv4/PEPPER_HP_OUTPUT_2.vcf.gz;
        rm -rf PEPPERv4/pepper_hp/;
        '''

#pepper_hp make_images [-h] -b BAM -f FASTA -t THREADS [-r REGION] [-o OUTPUT_DIR]
#pepper_hp run_inference -i /data/PEPPERv4/pepper_hp/images_04262021_165011 -t 18 -m /opt/pepper_models/PEPPER_HP_R941_ONT_V4.pkl -o /data/PEPPERv4/pepper_hp/predictions_04262021_165011 -s gaur -w 4 -bs 64
#pepper_hp find_candidates -i /data/PEPPERv4/pepper_hp/predictions_04262021_165011 -b /data/Assembly_ONT_reads.sorted.bam -f /data/Assembly.fasta -s gaur -o /data/PEPPERv4/pepper_hp -t 16 --ont_asm



rule pepper_tag_bam:
    input:
        bam = '',
        tags = get_dir('output','reads.haplotags')
    output:
        get_dir('output','haplotagged.bam')
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        tagBam \
        {input.bam} \
        {input.tags} \
        {output}
        '''

#rule extract_haplotags:
#    input:
#        unknown = config['reads']['unknown'],
#        hap1 = config['reads']['hap1'],
#        hap2 = config['reads']['hap2'],
#        reads = config['reads']['all']
#    output:
#        get_dir('output','reads.haplotags')
#    shell:
#        '''
#        zgrep -o ">\S*" {input.unknown} | cut -c 2- | awk '{print $0" H0}' > {output}
#        zgrep -o ">\S*" {input.hap1} | cut -c 2- | awk '{print $0" H1}' >> {output}
#        zgrep -o ">\S*" {input.hap2} | cut -c 2- | awk '{print $0" H2}' >> {output}
#        '''

rule deepvariant_make_examples:
    input:
        ref = lambda wildcards: multiext(config['assembly'],'','.fai'),
        bam = multiext(get_dir('input','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam'),'','.bai'),
        vcf = get_dir('output','PEPPER_HP_OUTPUT_{hap}.vcf.gz')
    output:
        example = temp(get_dir('DV','make_examples.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name(f'make_examples.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B .:/reference'),#{PurePath(input.ref[0]).parent}:/reference/'),
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
        --sample_name {wildcards.hap} \
        --noadd_hp_channel \
        --alt_aligned_pileup "rows" \
        --norealign_reads \
        --hp_tag_for_assembly_polishing "{wildcards.hap}" \
        --parse_sam_aux_fields \
        --sort_by_haplotypes \
        --variant_caller "vcf_candidate_importer" \
        --task {wildcards.N}
        '''

rule deepvariant_call_variants:
    input:
        (get_dir('DV', f'make_examples.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        temp(get_dir('DV','call_variants_output.tfrecord.gz'))
    params:
        examples = f'make_examples.tfrecord@{config["shards"]}.gz',
        model = '/opt/dv_models/202012_polish_nohp_rows/model.ckpt-28400',#lambda wildcards: get_model(wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container'],
        vino = lambda wildcards: '--use_openvino'
    threads: 32
    resources:
        mem_mb = 3000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '24:00'
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /opt/deepvariant/bin/call_variants \
        --outfile {output} \
        --examples {params.examples} \
        --checkpoint {params.model}
        '''

rule deepvariant_postprocess:
    input:
        ref = lambda wildcards: multiext(config['assembly'],'','.fai'),
        variants = get_dir('DV','call_variants_output.tfrecord.gz'),
    output:
        vcf = get_dir('output','PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP{hap}.vcf.gz'),
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
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
        '''

rule bcftools_consensus:
    input:
        asm = config['assembly'],
        vcf = get_dir('output','PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP{hap}.vcf.gz')
    output:
        get_dir('polished','Assembly.pepperv4.hap{hap}.fasta')
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B .:/reference'),
        contain = lambda wildcards: config['DV_container']
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        bcftools consensus \
        -f {input.asm}
        -H 2 \
        -s {wildcards.hap} \
        -o {output} \
        {input.vcf}
        '''
