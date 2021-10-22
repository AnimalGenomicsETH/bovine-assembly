from pathlib import Path,PurePath
from itertools import product

localrules: meryl_print_histogram, merfin_lookup, samtools_faidx

wildcard_constraints:
    merfin_op = r'hist|completeness|dump',
    read = r'hifi|SR',
    polished_read = r'unpolished|polished_hifi|polished_SR'

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='main',ext='', **kwargs):
    if base == 'main':
        base_dir = 'merfin_{animal}'
    elif base == 'lookup':
        base_dir = get_dir('main','lookup_{haplotype}_{read}',**kwargs)
    elif base == 'work':
        base_dir = get_dir('main','polishing_{haplotype}_{read}',**kwargs)
    else:
        raise Exception(f'Base: "{base}" not found in structure.')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

def capture_logic():
    conf_wildcards = {'animal':config['animal'],'haplotype':config['haplotype']}

    target_files = []
    for read in config['reads']:
        target_files.extend([get_dir('work','{polish}.{op}',read=read,op=OP,polish=POLISH,**conf_wildcards) for OP,POLISH in product(('hist','completeness'),('unpolished','polished'))])
    if 'hifi' in config['reads'] and 'SR' in config['reads']:
        target_files.extend([get_dir('main','{haplotype}.{polish}.correlation.svg',polish=POLISH,**conf_wildcards) for POLISH in ('unpolished','polished')])

    return target_files


rule all:
    input:
        capture_logic()

rule meryl_count_reads:
    input:
        lambda wildcards: config['reads'][wildcards.read]
    output:
        directory(get_dir('main','{haplotype}.readmers.{read}.meryl'))
    threads: 24
    resources:
        mem_mb = 5000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[kmer]} memory={params.mem} threads={threads} output {output} {input}'

rule meryl_count_asm:
    input:
        lambda wildcards: config['assembly'] if wildcards.polished_read == 'unpolished' else get_dir('work','{haplotype}.merfin.fasta',read=wildcards.polished_read.split('_')[1])
    output:
        directory(get_dir('main','{haplotype}.seqmers.{polished_read}.meryl'))
    threads: 8
    resources:
        mem_mb = 2500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[kmer]} memory={params.mem} threads={threads} output {output} {input}'

rule meryl_print_histogram:
    input:
        get_dir('main','{haplotype}.readmers.{read}.meryl')
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
        ploidy = 2#lambda wildcards: 2 if wildcards.haplotype == 'asm' else 2
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        'Rscript {config[lookup]} -i {input} -k {config[kmer]} -o {params.out} --fitted_hist -p {params.ploidy}'

rule merfin_vmer:
    input:
        fasta = config['assembly'],
        seqmers = get_dir('main','{haplotype}.seqmers.unpolished.meryl'),
        readmers = get_dir('main','{haplotype}.readmers.{read}.meryl'),
        vcf = config['vcf'],
        lookup = get_dir('lookup','lookup_table.txt'),
        model = get_dir('lookup','model.txt')
    output:
        temp(get_dir('work','merfin.polish.vcf'))
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_name('merfin'),
        coverage = lambda wildcards, input: float([line for line in open(input.model)][8].split()[1])
    threads: 12
    resources:
        mem_mb = 8000
    shell:
        'merfin -vmer -sequence {input.fasta} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup} -threads {threads} -peak {params.coverage} -vcf {input.vcf} -output {params.out}'

rule merfin_hist:
    input:
        fasta = lambda wildcards: config['assembly'] if wildcards.polished == 'unpolished' else get_dir('work','{haplotype}.merfin.fasta'),
        seqmers = lambda wildcards: get_dir('main','{haplotype}.seqmers.{X}.meryl',X= 'unpolished' if wildcards.polished == 'unpolished' else f'polished_{wildcards.read}'),
        readmers = get_dir('main','{haplotype}.readmers.{read}.meryl'),
        lookup = get_dir('lookup','lookup_table.txt'),
        model = get_dir('lookup','model.txt')
    output:
        get_dir('work','{polished}.{merfin_op}')
    params:
        coverage = lambda wildcards, input: float([line for line in open(input.model)][8].split()[1])
    threads: 16
    resources:
        mem_mb = lambda wildcards: 7500 if wildcards.merfin_op != 'dump' else 15000,
        walltime = '4:00'
    shell:
        'merfin -{wildcards.merfin_op} -sequence {input.fasta} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup} -threads {threads} -peak {params.coverage} -output {output}'

rule bcftools_polish:
    input:
        fasta = config['assembly'],
        vcf = get_dir('work','merfin.polish.vcf')
    output:
        vcf = temp(multiext(get_dir('work','merfin.polish.vcf.gz'),'','.csi')),
        fasta = get_dir('work','{haplotype}.merfin.fasta')
    threads: 1
    resources:
        mem_mb = 2000,
        walltime = '24:00'
    shell:
        '''
        bcftools view -Oz {input.vcf} > {output.vcf[0]}
        bcftools index {output.vcf[0]}
        bcftools consensus {output.vcf[0]} -f {input.fasta} -H 1 > {output.fasta}
        '''

rule merfin_correlation:
    input:
        (get_dir('work','{polished}.dump',read=read) for read in ('hifi','SR'))
    output:
        kstar = temp(get_dir('main','{haplotype}.{polished}.kstar_cor.txt')),
        plot = get_dir('main','{haplotype}.{polished}.correlation.svg')
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    threads: 2
    resources:
        mem_mb = 3000
    shell:
        '''
        cut -f3-5 {input[1]} | paste {input[0]} - | awk '{{if($3==0){{a="NA"}}else{{a=$5}};if($6==0){{b="NA"}}else{{b=$8}}; dups[a"\\t"b]++}} END{{for (num in dups) {{printf dups[num]"\\t"num"\\n"}}}}' | sort -nrk1 -nk2 > {output.kstar}
        Rscript {config[cartesian]} {output.kstar} {output.plot}
        '''

rule samtools_faidx:
    input:
        '{reference}'
    output:
        '{reference}.fai'
    shell:
        'samtools faidx {input}'

rule dump_to_bigWig:
    input:
        dump = get_dir('work','{polished}.dump'),
        fai = lambda wildcards: f"{config['assembly'] if wildcards.polished == 'unpolished' else get_dir('work','{haplotype}.merfin.fasta')}.fai"
    output:
        wig = temp(get_dir('work','{polished}.wig')),
        sizes =  temp(get_dir('work','{polished}.sizes')),
        bw = get_dir('work','{polished}.bw')
    threads: 1
    resources:
        mem_mb = 80000
    shell:
        '''
        awk 'BEGIN{{print "track autoScale=on"}}{{if($1!=chr){{chr=$1; print "variableStep chrom="chr" span=1"}};if($3!=0){{printf $2+1"\t"$5"\n"}}}' {input.dump} > {output.wig}
        cut -f1,2 {input.fai} > {output.sizes}
        wigToBigWig {output.wig} {output.sizes} {output.bw}
        '''
