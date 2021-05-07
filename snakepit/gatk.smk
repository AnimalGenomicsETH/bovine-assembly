from pathlib import PurePath
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='work',ext='', **kwargs):
    if base == 'input':
        base_dir = '{animal}'
    elif base == 'output':
        base_dir = 'output_{animal}'
    elif base == 'work':
        base_dir = get_dir('output','intermediate_results',**kwargs)
    elif base == 'main':
        base_dir = ''
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

wildcard_constraints:
     subset = r'|_child|_parent1|_parent2',
     haplotype = r'asm|hap1|hap2|parent1|parent2',
     phase = r'unphased|phased',
     model = r'pbmm2|hybrid|bwa|mm2'

def get_model(wildcards,base='/opt/models',ext='model.ckpt'):
    model_location = f'{base}/{{}}/{ext}'
    if wildcards['subset'] == '':
        if wildcards['model'] == 'pbmm2':
            return model_location.format('pacbio')
        elif wildcards['model'] == 'hybrid':
            return model_location.format('hybrid_pacbio_illumina')
        elif wildcards['model'] == 'bwa':
            return model_location.format('wgs')
        elif wildcards['model'] == 'mm2':
            return model_location.format('wgs')
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
        call += f'-B {":/".join([get_dir("input",**wildcards)]*2)} '
    if output_bind:
        call += f'-B {":/".join([get_dir("output",**wildcards)]*2)} '
    if work_bind:
        call += f'-B {":/".join([get_dir("work",**wildcards)]*2)} '
    if tmp_bind:
        call += f'-B {tmp_bind}:/tmp '
    return call

for dir_ in ('input','output','work'):
    for animal,reference in product(config['animals'],config['references']):
        Path(get_dir(dir_,animal=animal,ref=reference,**config)).mkdir(exist_ok=True)

rule all:
    input:
        get_dir('main','cohort.vcf.gz'),
        expand(get_dir('output','{animal}.bwa.vcf.gz',**config),animal=config['animals'].keys(),ref=config['references'].keys()) if not config['trio'] else [],
        expand(get_dir('output','{haplotype}.trio.merged.{phase}.{ref}.{model}.{ext}',**config),ext=('vcf.gz','mendelian.log'),animal=config['animals'].keys(),ref=config['references'].keys()) if config['trio'] else []

rule minimap_align:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref],
        reads = lambda wildcards: config['animals'][wildcards.animal]['individual' if 'parent' not in wildcards.haplotype else wildcards.haplotype]['short_reads']
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.mm2.bam'))
    threads: 24
    resources:
        mem_mb = 5000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax sr -t {threads} {input.ref} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule bwamem2_index:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref]
    output:
        temp(get_dir('input','{ref}_index.0123'))
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 1
    resources:
        mem_mb = 85000
    shell:
        'bwa-mem2 index -p {params} {input}'

rule bwamem2_align:
    input:
        index = lambda wildcards: str(PurePath(config['references'][wildcards.ref]).parent / f'{wildcards.ref}_index' / f'{wildcards.ref}.0123'),
        reads = lambda wildcards: config['animals'][wildcards.animal]['individual' if 'parent' not in wildcards.haplotype else wildcards.haplotype]['short_reads']
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.bwa.bam'))
    params:
        prefix = lambda wildcards, input: PurePath(input.index).with_suffix('')
    threads: 12
    resources:
        mem_mb = 6000,
        walltime = '24:00',
        disk_scratch = 250
    shell:
        'bwa-mem2 mem -t {threads} {params.prefix} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule pbmm2_align:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref],
        reads = lambda wildcards: config['animals'][wildcards.animal]['individual' if 'parent' not in wildcards.haplotype else wildcards.haplotype]['long_reads']
    output:
        temp(get_dir('input','{haplotype}.unphased.{ref}.pbmm2.bam'))
    threads: 26
    resources:
        mem_mb = 3000,
        walltime = '4:00'
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
        ref = multiext(config['references']['ARS'],'','.fai'),
        bam = multiext(get_dir('input','{animal}.bam'),'','.bai')
    output:
        example = temp(get_dir('work','make_examples.tfrecord-{N}-of-{sharding}.gz')),
        gvcf = temp(get_dir('work','gvcf.tfrecord-{N}-of-{sharding}.gz'))
    params:
        examples = lambda wildcards, output: PurePath(output[0]).with_name(f'make_examples.tfrecord@{config["shards"]}.gz'),
        gvcf = lambda wildcards, output: PurePath(output[1]).with_name(f'gvcf.tfrecord@{config["shards"]}.gz'),
        phase_args = lambda wildcards: '--parse_sam_aux_fields={i} --sort_by_haplotypes={i}'.format(i=('true' if False else 'false')),
        model_args = lambda wildcards: '--add_hp_channel --alt_aligned_pileup diff_channels --realign_reads=false --vsc_min_fraction_indels 0.12' if 'bwa' == 'pbmm2' else '',
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
        examples = lambda wildcards: (f'make_examples{wildcards.subset}.tfrecord@{config["shards"]}.gz'),
        model = lambda wildcards: get_model({'subset':'','model':'bwa'}),
        dir_ = lambda wildcards: get_dir('work',**wildcards),
        singularity_call = lambda wildcards, threads: make_singularity_call(wildcards,f'--env OMP_NUM_THREADS={threads}'),
        contain = lambda wildcards: config['DV_container'] if wildcards.subset == '' else config['DT_container'],
        vino = lambda wildcards: '--use_openvino' if wildcards.subset == '' else ''
    threads: 32
    resources:
        mem_mb = 4000,
        disk_scratch = 1,
        use_singularity = True,
        walltime = lambda wildcards: '4:00' if wildcards.subset == '' else '24:00'
    shell:
        '''
        {params.singularity_call} \
        {params.contain} \
        /bin/bash -c "cd {params.dir_}; /opt/deepvariant/bin/call_variants \
        --outfile /{output} \
        --examples {params.examples} \
        --checkpoint {params.model} \
        {params.vino}"
        '''

rule deepvariant_postprocess:
    input:
        ref = multiext(config['references']['ARS'],'','.fai'),
        variants = get_dir('work','call_variants_output.tfrecord.gz'),
        gvcf = (get_dir('work', f'gvcf.tfrecord-{N:05}-of-{config["shards"]:05}.gz') for N in range(config['shards']))
    output:
        vcf = get_dir('output','{animal}.bwa.vcf.gz'),
        gvcf = get_dir('output','{animal}.bwa.g.vcf.gz')
    params:
        gvcf = lambda wildcards,input: PurePath(input.gvcf[0]).with_suffix('').with_suffix(f'.tfrecord@{config["shards"]}.gz'),
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/'),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
        contain = lambda wildcards: config['DV_container']
    threads: 1
    resources:
        mem_mb = 60000,
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

rule GLnexus_merge:
    input:
        expand(get_dir('output','{animal}.bwa.g.vcf.gz'),animal=config['animals'])
    output:
        get_dir('main','cohort.vcf.gz')
    params:
        gvcfs = lambda wildcards, input: list('/data/' / PurePath(fpath) for fpath in input),
        out = lambda wildcards, output: f'/data/{PurePath(output[0]).name}',
        DB = lambda wildcards, output: f'/tmp/GLnexus.DB',
        singularity_call = lambda wildcards: make_singularity_call(wildcards,'-B .:/data', input_bind=False, output_bind=False, work_bind=False),
        mem = lambda wildcards,threads,resources: threads*resources['mem_mb']/1000
    threads: 24 #force using 4 threads for bgziping
    resources:
        mem_mb = 5000,
        disk_scratch = 100,
        use_singularity = True
    shell:
        '''
        {params.singularity_call} \
        {config[GL_container]} \
        /bin/bash -c " /usr/local/bin/glnexus_cli \
        --dir {params.DB} \
        --config DeepVariantWGS \
        --threads {threads} \
        --mem-gbytes {params.mem} \
        {params.gvcfs} \
        | bcftools view - | bgzip -@ 4 -c > {params.out}"
        '''
        #| bcftools view - | bgzip -@ 4 -c > {params.out}
        #'''
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
        sdf = directory(get_dir('input','{ref}.sdf'))
    params:
        singularity_call = lambda wildcards,input: make_singularity_call(wildcards,extra_args=f'-B {PurePath(input.ref[0]).parent}:/reference/',tmp_bind=False,output_bind=False,work_bind=False),
        ref = lambda wildcards,input: f'/reference/{PurePath(input.ref[0]).name}',
    shell:
        '''
        {params.singularity_call} \
        {config[RTG_container]} \
        rtg format -o {output.sdf} {params.ref}
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
