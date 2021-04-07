from pathlib import PurePath

rule all:
    input:
        str(PurePath(config['fasta']).with_suffix('').with_suffix('.fq.gz'))

rule resolve_haplotype_fastq:
    input:
        fastq_data = config['fastq'],
        fasta_reads = config['fasta']
    output:
        str(PurePath(config['fasta']).with_suffix('').with_suffix('.fq.gz'))
    threads: 8
    resources:
        mem_mb = 1500,
        walltime = '12:00'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        zcat {input.fasta_reads} | grep ">" | cut -c 2- | seqtk subseq {input.fastq_data} - | pigz -4 -p {threads} > {output}
        '''

