rule telomer_content:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{assembler}.telo.txt'
    script:
        '../scripts/count_telomers.py'

rule raw_QC:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        'data/{animal}.QC.txt'
    shell:
        '''
        src/fasterqc {input} {output}
        '''

rule sample_data:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        'data/{animal}SAMP{sample}.hifi.fq.gz'
    shell:
        'seqtk sample 100 {input} {output}'


