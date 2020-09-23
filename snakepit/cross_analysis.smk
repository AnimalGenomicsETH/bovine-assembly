localrules: count_telomers, count_scaffold_gaps, scaffold_chromosomes

rule telomer_asdcontent:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{assembler}.teasdlo.txt'
    script:
        '../scripts/count_telomers.py'

rule count_telomers:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{assembler}.telo.txt'
    run:
        import re, screed
        from scipy.stats import binom
        region = 1000
        telomere = re.compile("TTAGGG", re.IGNORECASE)

        with open(output[0],'w') as fout:
            out.write('name\trepeat_count\tprobability\n')
            for seq in screed.open(snakemake.input[0]):
                c_repeats = len(telomere.findall(seq.sequence[:region]))
                fout.write(f'{seq.name}\t{c_repeats}\t{binom.sf(c_repeats,region,0.25**6):.4f}\n')

rule count_scaffold_gaps:
    input:
        '{assembler}/{animal}.scaffolds.fasta'
    output:
        'results/{animal}_{assembler}.gaps.txt'
    run:
        import re, screed
        gap_sequence = re.compile(r'[nN]+')
        with open(output[0],'w') as fout:
            fout.write('scaffold\tgaps\tlengths\n')
            for scaffold in screed.open(input[0]):
                contig_lengths = list(map(len, gap_sequence.split(str(scaffold.sequence))))
                #gap_lengths = list(map(len, gap_sequence.findall(str(scaffold.sequence))))
                #for F in (re.split, re.findall):
                #     map(lambda s: str(len(s)), F(scaffold.sequence))
                fout.write('\t'.join((scaffold.name,str(len(contig_lengths)-1),','.join(map(str,contig_lengths))))+'\n')

#ragtag.py correct, full reads mapping, c_mapping
rule ragtag_correct:
    input:
        asm = '{assembler}/{animal}.contigs.fasta'
        reads = '{data/{animal}.hifi.fq.gz'
    output:
        '{assembler}/{animal}.contigs.corrected.fasta'
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        'ragtag.py correct {config[ref_genome]} {input.asm} -o {wildcards.assembler} -R {input.reads} -T corr -t {threads}'

rule ragtag_scaffold:
    input:
        '{assembler}/{animal}.contigs.fasta'
        #or corrected
    output:
        '{assembler}/ragtag.scaffolds.fasta'
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        'ragtag.py scaffold {config[ref_genome]} {input} -o {wildcards.assembler} -t {threads}'

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
        '{assembler}/{file}.cov'
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
        'seqtk sample {input} {wildcards.sample} > {output}'
