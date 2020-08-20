rule count_read_kmers:
    input:
        'data/{animal}.fq.gz'
    output:
        directory('{animal}_hifi.meryl')
    threads: 24
    resources:
        mem_mb = 5000
    params:
        K=21,
        mem_limit = resources.mem_mb * threads * .9
    shell:
        'meryl count k={params.K} memory={params.mem_limit} threads={threads} output {output} {input}'
       
rule count_asm_kmers:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        directory('{assembler}_{animal}_asm.meryl')
    threads: 24
    resources:
        mem_mb = 5000
    params:
        K=21,
        mem_limit = resources.mem_mb * threads * .9
    shell:
        'meryl count k={params.K} memory={params.mem_limit} threads={threads} output {output} {input}'
 
rule spectra_cn:

