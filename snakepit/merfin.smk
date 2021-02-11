from pathlib import Path,PurePath
import shutil

localrules: ratatosk_make_bins, ratatosk_merge_bin1, ratatosk_merge_bin2, ratatosk_finish

if Path('config/ratatosk.yaml').exists():
    configfile: 'config/ratatosk.yaml'

wildcard_constraints:
    N = r'\d+'

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='segments',ext='', **kwargs):
    if base == 'result':
        base_dir = get_dir('main','corrected',**kwargs)
    elif base == 'main':
        base_dir = 'ratatosk_{animal}'
    elif base == 'segments':
        base_dir = get_dir('main','segments',**kwargs)
    else:
        raise Exception('Base not found')


rule meryl_count:
    input:
        get_dir('work','{haplotype}.contigs.fasta')
    output:
        directory(get_dir('work','{haplotype}.contigs.meryl'))
    threads: 8
    resources:
        mem_mb = 2500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
    shell:
        'meryl count k=21 memory={params.mem} threads={threads} output {output} {input}'

rule merfin_lookup:

rule merfin_vmer:
    input:
        fasta = config['fasta']
        seqmers = f'fasta.meryl',
        readmers = 'reads.meryl',
        vcf = ''
        lookup = 'lookup/lookup_table.txt'
    output:
        'merfin.polish.vcf'
    params:
        coverage = 30,
        low_mem = 10,
        high_mem = 60
    threads: 16
    resources:
        mem_mb = 5500
    shell:
        'merfin -vmer -sequence {input.fasta} -memory1 {params.low_mem} -memory2 {params.high_mem} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup}-threads {threads} -peak {params.coverage} -vcf {input.vcf} -output {output}'

rule merfin_hist:
    input:
        fasta = config['fasta']
        seqmers = f'fasta.meryl',
        readmers = 'reads.meryl',
        vcf = ''
        lookup = 'lookup/lookup_table.txt'
    output:
        'merfin.polish.vcf'
    params:
        coverage = 30,
        low_mem = 10,
        high_mem = 60
    threads: 16
    resources:
        mem_mb = 5500
    shell:
        'merfin -hist -sequence {input.fasta} -memory1 {params.low_mem} -memory2 {params.high_mem} -seqmers {input.seqmers} -readmers {input.readmers} -lookup {input.lookup}-threads {threads} -peak {params.coverage} -vcf {input.vcf} -output {output}'


rule merfin_polish:
    input:
        fasta = '{haplotype}.scaffolds.fasta',
        vcf = 'merfin.polish.vcf'
    output:
        vcf = multiext('merfin.polish.vcf','','.csi')
        fasta = 'asm.merfin.fasta'
    shell:
        '''
        bcftools view -Oz {input.vcf} > {output.vcf[0]}
        bcftools index {output.vcf[0]}
        bcftools consensus {output.vcf[0]} -f {input.fasta} -H 1 > {output.fasta}
        '''
