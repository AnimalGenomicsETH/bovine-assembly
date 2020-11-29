localrules: resolve_haplotype_reads

rule resolve_haplotype_reads:
    input:
        data = 'data/offspring.raw.hifi.fq.gz',
        asm = 'canu_100/hap{N}.contigs_all.fa'
    output:
        'data/haplotype_N{N}.hifi.fq.gz'
    params:
        fasta = 'canu_100/trio/haplotype/haplotype-{N}.fasta.gz',
        headers = 'data/haplotype_N{N}.headers'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        grep ">" {params.fasta} | cut -c 2- | seqtk subseq {input.data} - | pigz -8 -p 4 > {output}
        '''

rule pbsv_align:
    input:
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        reads = 'data/haplotype_N{N}.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{haplotype}_hap{N}.pbmm2.bam'
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        'pbmm2 align {input.asm} {input.reads} {output} --sort --preset CCS -j {threads}'

rule pbsv_discover:
    input:
        '{assembler}_{sample}/{haplotype}_hap{N}.pbmm2.bam'
    output:
        '{assembler}_{sample}/{haplotype}_hap{N}.svsig.gz'
    shell:
        'pbsv discover {input} {output}'

rule pbsv_call:
    input:
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        sig = '{assembler}_{sample}/{haplotype}_hap{N}.svsig.gz'
    output:
        '{assembler}_{sample}/{haplotype}_hap{N}.pbsv.vcf'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'pbsv call --ccs -j {threads} {input.asm} {input.sig} {output}'
