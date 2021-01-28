from pathlib import PurePath

rule resolve_haplotype_fastq:
    input:
        fastq_data = config['fastq'],
        fasta_reads = config['fasta']
    output:
        PurePath(config['fasta']).with_suffix('').with_suffix('.fq.gz')
    threads: 8
    resources:
        mem_mb = 1500
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        zcat {input.fasta_reads} | grep ">" | cut -c 2- | seqtk subseq {input.fastq_data} - | pigz -4 -p {threads} > {output}
        '''

rule pbsv_align:
    input:
        asm = WORK_PATH + '{haplotype}.scaffolds.fasta',
        reads = 'data/haplotype_N{N}.hifi.fq.gz'
    output:
        WORK_PATH + '{haplotype}_hap{N}.pbmm2.bam'
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        'pbmm2 align {input.asm} {input.reads} {output} --sort --preset CCS -j {threads}'

rule pbsv_discover:
    input:
        WORK_PATH + '{haplotype}_hap{N}.pbmm2.bam'
    output:
        WORK_PATH + '{haplotype}_hap{N}.svsig.gz'
    shell:
        'pbsv discover {input} {output}'

rule pbsv_call:
    input:
        asm = WORK_PATH + '{haplotype}.scaffolds.fasta',
        sig = WORK_PATH + '{haplotype}_hap{N}.svsig.gz'
    output:
        WORK_PATH + '{haplotype}_hap{N}.pbsv.vcf'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'pbsv call --ccs -j {threads} {input.asm} {input.sig} {output}'
