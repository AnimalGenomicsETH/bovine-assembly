from pathlib import Path,PurePath


localrules: meryl_print_histogram, merfin_lookup

if Path('config/ratatosk.yaml').exists():
    configfile: 'config/ratatosk.yaml'

wildcard_constraints:
    N = r'\d+',
    merfin_op = r'hist|completeness'

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='main',ext='', **kwargs):
    if base == 'main':
        base_dir = 'merfin_{animal}'
    elif base == 'lookup':
        base_dir = get_dir('main','lookup_{haplotype}',**kwargs)
    elif base == 'work':
        base_dir = get_dir('main','polishing_{haplotype}',**kwargs)
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext)

rule all:
    input:
        expand(get_dir('work','{polish}.{op}',haplotype=config['haplotype'],animal=config['animal']),op=('hist','completeness'),polish=('unpolished','polished'))

rule meryl_count_reads:
    input:
        config['reads']
    output:
        directory(get_dir('main','{haplotype}.readmers.meryl'))
    threads: 16
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[kmer]} memory={params.mem} threads={threads} output {output} {input}'

rule meryl_count_asm:
    input:
        lambda wildcards: config['assembly'] if wildcards.polished == 'unpolished' else get_dir('work',config['polished'])
    output:
        directory(get_dir('work','seqmers.{polished}.meryl'))
    threads: 8
    resources:
        mem_mb = 2500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[kmer]} memory={params.mem} threads={threads} output {output} {input}'

rule meryl_print_histogram:
    input:
        get_dir('main','{haplotype}.readmers.meryl')
    output:
        get_dir('lookup','read_hist.tsv')
    shell:
        'meryl histogram {input} > {output}'

rule merfin_lookup:
    input:
        get_dir('lookup','read_hist.tsv')
    output:
        get_dir('lookup','lookup_table.txt'),
        get_dir('lookup','model.txt')
    params:
        out = lambda wildcards, output: PurePath(output[0]).parent,
        ploidy = lambda wildcards: 2 if wildcards.haplotype == 'asm' else 2
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        'Rscript {config[lookup]} -i {input} -k {config[kmer]} -o {params.out} --fitted_hist -p {params.ploidy}'

rule merfin_vmer:
    input:
        fasta = config['assembly'],
        seqmers = get_dir('work','seqmers.unpolished.meryl'),
        readmers = get_dir('main','{haplotype}.readmers.meryl'),
        vcf = config['vcf'],
        lookup = get_dir('lookup','lookup_table.txt'),
        model = get_dir('lookup','model.txt')
    output:
        get_dir('work','merfin.polish.vcf')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_name('merfin'),
        coverage = lambda wildcards, input: float([line for line in open(input.model)][8].split()[1])
    threads: 12
    resources:
        mem_mb = 6000
    shell:
        'merfin -vmer -sequence {input.fasta} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup} -threads {threads} -peak {params.coverage} -vcf {input.vcf} -output {params.out}'

rule merfin_hist:
    input:
        fasta = lambda wildcards: config['assembly'] if wildcards.polished == 'unpolished' else get_dir('work',config['polished']),
        seqmers = get_dir('work','seqmers.{polished}.meryl'),
        readmers = get_dir('main','{haplotype}.readmers.meryl'),
        lookup = get_dir('lookup','lookup_table.txt'),
        model = get_dir('lookup','model.txt')
    output:
        get_dir('work','{polished}.{merfin_op}')
    params:
        coverage = lambda wildcards, input: float([line for line in open(input.model)][8].split()[1])
    threads: 16
    resources:
        mem_mb = 5500
    shell:
        'merfin -{wildcards.merfin_op} -sequence {input.fasta} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup} -threads {threads} -peak {params.coverage} -output {output}'

rule merfin_polish:
    input:
        fasta = config['assembly'],
        vcf = get_dir('work','merfin.polish.vcf')
    output:
        vcf = multiext(get_dir('work','merfin.polish.vcf.gz'),'','.csi'),
        fasta = get_dir('work',config['polished'])
    shell:
        '''
        bcftools view -Oz {input.vcf} > {output.vcf[0]}
        bcftools index {output.vcf[0]}
        bcftools consensus {output.vcf[0]} -f {input.fasta} -H 1 > {output.fasta}
        '''


