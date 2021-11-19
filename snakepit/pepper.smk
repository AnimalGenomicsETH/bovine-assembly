from pathlib import PurePath
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = 'PEPPER'
    elif base == 'output':
        base_dir = 'PEPPER'
    elif base == 'DV':
        base_dir = get_dir('output','dv_intermediate_outputs_{hap}',**kwargs)
    else:
        raise Exception(f'Base: {base} was not found!')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

wildcard_constraints:
     subset = r'|_child|_parent1|_parent2',
     haplotype = r'asm|hap1|hap2|parent1|parent2',
     hap = r'1|2|u',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid|bwa',
     ext = r'|_1|_2'

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
        expand(get_dir('output','Assembly.pepperv4.hap{hap}.fasta'),haplotype=config['haplotype'],hap=('1','2')) if not config['trio'] else \
        get_dir('output','Assembly.pepperv4.hap{hap}.fasta',haplotype=config['haplotype'],hap=config['haplotype'])


rule samtools_faidx:
    input:
        '{reference}'
    output:
        '{reference}.fai'
    shell:
        'samtools faidx {input}'

rule map_ONT_reads:
    input:
        reads = lambda wildcards: config['reads']['hap'+wildcards.hap],
        asm = config['assembly']
    output:
        temp(get_dir('input','ONT_reads.hap{hap}.unsorted.mm2.bam'))
    threads: 16
    resources:
        mem_mb = 4500,
        walltime = '24:00'
    shell:
        'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} | samtools view -@ 2 -hb -F 0x904 -o {output}'

rule sort_bam:
    input:
        get_dir('input','ONT_reads.hap{hap}.unsorted.mm2.bam')
    output:
        temp(get_dir('input','ONT_reads.hap{hap}.mm2.bam'))
    threads: 12
    resources:
        mem_mb = 6000,
        disk_scratch = 200
    shell:
        'samtools sort {input} -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule samtools_merge:
    input:
        (get_dir('input','ONT_reads.hap{hap}.mm2.bam',hap=HAP) for HAP in (1,2,'u'))#config['reads'])
    output:
        get_dir('input','ONT_reads.mm2.bam')
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        'samtools merge --threads {threads} {output} {input}'

rule samtools_index_bam:
    input:
        '{bam}.bam'
    output:
        '{bam}.bam.bai'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule pepper_make_images:
    input:
        ref = config['assembly'],
        bam = lambda wildcards: multiext(get_dir('input','ONT_reads.mm2.bam') if wildcards.mode == 'snp' else get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam'),'','.bai')
    output:
        directory(get_dir('output','pepper_{mode}/images'))
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,work_bind=False,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        pepper_mode = lambda wildcards: 'pepper_hp' if wildcards.mode == 'hp' else 'pepper_snp',
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref).name}'
    threads: 24
    resources:
        mem_mb = 12000,
        disk_scratch = 10
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        {params.pepper_mode} make_images \
        -b {input.bam} \
        -f {params.ref} \
        -t {threads} \
        -o {output}
        '''

rule pepper_run_inference:
    input:
        get_dir('output','pepper_{mode}/images')
    output:
        directory(get_dir('output','pepper_{mode}/predictions'))
    params:
        singularity_call = lambda wildcards: make_singularity_call(wildcards,work_bind=False),
        pepper_mode = lambda wildcards: 'pepper_hp' if wildcards.mode == 'hp' else 'pepper_snp',
        model = lambda wildcards: f'/opt/pepper_models/PEPPER_{wildcards.mode.upper()}_R941_ONT_V4.pkl'
    threads: 12
    resources:
        mem_mb = 12000,
        walltime = '24:00',
        disk_scratch = 10
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        {params.pepper_mode} run_inference \
        -i {input} \
        -m {params.model} \
        -o {output} \
        -s {config[sample]} \
        -w 4 \
        -bs 64 \
        -t {threads}
        '''

rule pepper_snp_find_candidates:
    input:
        predictions = get_dir('output','pepper_snp/predictions'),
        bam = get_dir('input','ONT_reads.mm2.bam'),
        ref = config['assembly']
    output:
        get_dir('output','pepper_snp/PEPPER_SNP_OUPUT.vcf.gz')
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,work_bind=False,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref).name}'
    threads: 12
    resources:
        mem_mb = 12000,
        walltime = '24:00',
        disk_scratch = 10
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        pepper_snp find_candidates \
        -i {input.predictions} \
        -b {input.bam} \
        -f {params.ref} \
        -s {config[sample]} \
        -o {output} \
        -t {threads} \
        --ont
        '''

rule pepper_hp_find_candidates:
    input:
        predictions = get_dir('output','pepper_hp/predictions'),
        bam = get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam'),
        ref = config['assembly']
    output:
        (get_dir('output',f'pepper_hp/PEPPER_HP_OUTPUT_{N}.vcf') for N in (1,2))
    params:
        out = lambda wildcards: get_dir('output','pepper_hp',**wildcards),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,work_bind=False,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref).name}'
    threads: 12
    resources:
        mem_mb = 12000,
        walltime = '24:00',
        disk_scratch = 10
    shell:
        '''
        {params.singularity_call} \
        {config[DV_container]} \
        pepper_hp find_candidates \
        -i {input.predictions} \
        -b {input.bam} \
        -f {params.ref} \
        -s {config[sample]} \
        -o {params.out} \
        -t {threads} \
        --ont_asm
        '''

rule pepper_post_vcf:
    input:
        lambda wildcards: get_dir('output',f'pepper_{wildcards.mode.lower()}/PEPPER_{{mode}}_OUTPUT{{ext}}.vcf')
    output:
        get_dir('output','PEPPER_{mode}_OUTPUT{ext}.vcf.gz')
    threads: 4
    resources:
        mem_mb = 4000,
        walltime = '1:00'
    shell:
        '''
        bgzip -@ 4 -c {input} > {output}
        tabix -p vcf {output}
        '''

if not config['trio']:  
    rule pepper_margin:
        input:
            bam = get_dir('input',f'ONT_reads.mm2.bam'),
            asm = config['assembly'],
            vcf = get_dir('output','PEPPER_SNP_OUTPUT.vcf.gz')
        output:
            get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam')
        params:
            singularity_call = lambda wildcards: make_singularity_call(wildcards,work_bind=False,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
            ref = lambda wildcards,input: f'/reference/{PurePath(input.ref).name}',
            json = '/opt/margin_dir/params/misc/allParams.ont_haplotag.json'        
        threads: 24
        resources:
            mem_mb = 5000,
            disk_scratch = 10
        shell:
            '''
            {params.singularity_call} \
            {config[DV_container]} \
            margin phase {input.bam} {input.asm} {input.vcf} {params.json} -t {threads} -V -o {output}
            '''
else:   
    rule pepper_tag_bam:
        input:
            bam = get_dir('input',f'ONT_reads.mm2.bam'),
            tags = get_dir('output','reads.haplotags')
        output:
            get_dir('output','MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam')
        params:
            singularity_call = lambda wildcards: make_singularity_call(wildcards,work_bind=False)
        threads: 1
        resources:
            mem_mb = 10000,
            disk_scratch = 1,
            walltime = '24:00'
        shell:
            '''
            {params.singularity_call} \
            {config[DV_container]} \
            tagBam \
            {input.bam} \
            {input.tags} \
            {output}
            '''

    rule extract_haplotags:
        input:
            hap1 = config['reads']['hap1'],
            hap2 = config['reads']['hap2'],
            hapu = config['reads']['hapu']
        output:
            get_dir('output','reads.haplotags')
        params:
            read_tag = '2' if ('.fa' in PurePath(config['reads']['hap1']).suffixes or '.fasta' in PurePath(config['reads']['hap1']).suffixes) else '4' 
        threads: 1
        resources:
            mem_mb = 5000,
            walltime = '2:00'
        shell:
            '''
            awk 'NR%{params.read_tag}==1 {{print $1"\tH1"}}' <(zcat {input.hap1}) | cut -c 2- >> {output}
            awk 'NR%{params.read_tag}==1 {{print $1"\tH2"}}' <(zcat {input.hap2}) | cut -c 2- >> {output}
            awk 'NR%{params.read_tag}==1 {{print $1"\tnone"}}' <(zcat {input.hapu}) | cut -c 2- >> {output}
            '''

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
        mem_mb = 20000,
        walltime = '24:00',
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
        examples = lambda wildcards, input: PurePath(input[0]).with_name(f'make_examples.tfrecord@{config["shards"]}.gz'),
        model = '/opt/dv_models/202012_polish_nohp_rows/model.ckpt-28400',#lambda wildcards: get_model(wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container'],
        vino = lambda wildcards: '--use_openvino'
    threads: 24
    resources:
        mem_mb = 3000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '120:00'
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
        mem_mb = 90000,
        disk_scratch = 1,
        walltime = '24:00',
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
        ref = lambda wildcards: multiext(config['assembly'],'','.fai'),
        vcf = get_dir('output','PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP{hap}.vcf.gz')
    output:
        get_dir('output','Assembly.pepperv4.hap{hap}.fasta')
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/',tmp_bind=False,work_bind=False),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24:00',
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        bcftools consensus \
        -f {params.ref} \
        -H 2 \
        -s {wildcards.hap} \
        -o {output} \
        {input.vcf}
        '''
