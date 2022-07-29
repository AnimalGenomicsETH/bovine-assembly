from pathlib import PurePath

rule all:
    input:
        expand('hifiasm/{sample}.dip.{N}.p_ctg.gfa', N = ('hap1','hap2'),sample=config['samples'])


rule hifiasm_primary:
    input:
        lambda wildcards: expand('EXTRACTED_FASTQ/{{sample}}.{cell}.fq.gz',cell=config['samples'][wildcards.sample]['long_reads'])
    output:
        'hifiasm/{sample}.bp.p_ctg.gfa'
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix(''),
        settings = '-r 3 -a 5 -n 5'
    threads: 32
    resources:
        mem_mb = lambda wildcards, input, threads: max(int(input.size_mb*1.75/threads),6000),
        walltime = '24:00'
    shell:
        'hifiasm -o {params.out} -t {threads} {params.settings} {input}'


rule yak_count_SR:
    input:
        lambda wildcards: config['samples'][wildcards.sample]['trio_reads'][wildcards.haplotype]
    output:
        temp('hifiasm/{sample}.{haplotype}.yak')
    threads: 24
    resources:
        mem_mb = 3000
    params:
        kmer = 31
    shell:
        'yak count -k {params.kmer} -b 37 -t {threads} -o {output} <(zcat {input}) <(zcat {input})'

rule hifiasm_trio:
    input:
        parents = expand('hifiasm/{{sample}}.{N}.yak',N=('hap1','hap2')),
        asm = 'hifiasm/{sample}.bp.p_ctg.gfa'
    output:
        expand('hifiasm/{{sample}}.dip.{N}.p_ctg.gfa', N=('hap1','hap2'))
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').with_suffix('').with_suffix(''),
        settings = '-r 3 -a 5 -n 5'
    threads: 32
    resources:
        mem_mb = 4000,
        walltime = '4:00'
    shell:
        'hifiasm -o {params.out} -t {threads} {params.settings} -1 {input.parents[0]} -2 {input.parents[1]} /dev/null'
