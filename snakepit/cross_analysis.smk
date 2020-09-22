localrules: telomer_content, count_scaffold_gaps, scaffold_chromosomes

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
        gap_sequence = re.compile(r'[nN]+')
        with open(output[0],'w') as fout:
            fout.write('scaffold\tgaps\tlengths\tglen\n')
            for scaffold in screed.open(input[0]):
                contig_lengths = list(map(len, gap_sequence.split(str(scaffold.sequence))))
                gap_lengths = list(map(len, gap_sequence.findall(str(scaffold.sequence))))
                #for F in (re.split, re.findall):
                #     map(lambda s: str(len(s)), F(scaffold.sequence))
                fout.write('\t'.join((scaffold.name,str(len(contig_lengths)-1),','.join(map(str,contig_lengths)),','.join(map(str,gap_lengths))))+'\n')

#ragtag.py correct, full reads mapping, c_mapping
rule ragtag_correct_alignment:
    input:
        reads = 'data/{animal}.hifi.fq.gz',
        asm = '{assembler}/{animal}.contigs.fasta'
    output:
        '{assembler}/{animal}_reads_aln.sam'
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        'minimap2 -ax asm20 -t {threads} {input.asm} {input.reads} > {output}'
     
rule ragtag_correct:
    input:
        asm = '{assembler}/{animal}.contigs.fasta',
        aln = '{assembler}/{animal}_reads_aln.sam',
        reads = '{data/{animal}.hifi.fq.gz'
    output:
        '{assembler}/{animal}.corrected.fasta'
    shell:
        'ragtag.py correct {config[ref_genome]} {input.asm} -o {wildcards.assembler} -R {input.reads} -T corr'

rule scaffold_alignment:
    input:
        asm = '{assembler}/{animal}.contigs.corrected.fasta'
    output:
        '{assembler}/{animal}_scaffold.paf'
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        'minimap2 -c -x asm5 -t {threads} {config[ref_genomes]} {input.asm} > {output}'

rule scaffold_chromosomes:
    input:
        asm = '{assembler}/{animal}.contigs.fasta',
        aln = '{assembler}/{animal}_scaffold.paf'
    output:
        '{assembler}/{animal}.scaffold.fasta'
    shell:
        'ragtag.py scaffold {contig[ref_genomes]} {input.asm} -o {wildcards.assembler}'

rule remap_reads:
    input:
        reads = 'data/{animal}.hifi.fq.gz',
        asm = '{assembler}/{animal}.scaffolds.fasta'
    output:
        '{assembler}/{animal}_scaffold_read.sam'
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        'minimap2 -ax asm20 -t {threads} {input.asm} {input.reads} > {output}'

rule sam_to_bam:
    input:
        '{assembler}/{file}.sam'
    output:
        '{assembler}/{file}.bam'
    params:
        memory = 4000
    threads: 16
    resources:
        mem_mb = '{params}',
        disk_scratch = 200
    shell:
        'samtools sort {input} -m {params}M -@ {threads} -T \$TMPDIR -o {output}'

rule assembly_coverage:
    input:
        '{assembler}/{file}.bam'
    output:
        'cov'
    shell:
        'samtools coverage {input} > {output}'

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
