rule telomer_content:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{assembler}.telo.txt'
    script:
        '../scripts/count_telomers.py'

rule count_scaffold_gaps:
    input:
        '{assembler}/{animal}.scaffolds.fasta'
    output:
        'results/{animal}_{assembler}.gaps.txt'
    run:
        import re, screed
        with open(output[0],'w') as fout:
            fout.write('scaffold\tgaps\tlengths\n')
            for scaffold in screed.open(input[0]):
                contig_lengths = list(map(len, re.split('[nN]+', str(scaffold.sequence))))
                fout.write('\t'.join((scaffold.name,str(len(contig_lengths)-1),','.join(map(str,contig_lengths))))+'\n')

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


