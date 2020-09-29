localrules: count_telomers, count_scaffold_gaps, prep_window, window_coverage, chromosome_coverage, sample_data, raw_QC

rule count_telomers:
    input:
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{sample}_{assembler}.telo.txt'
    run:
        import re, screed
        from scipy.stats import binom
        region = 1000
        telomere = re.compile("TTAGGG", re.IGNORECASE)

        with open(output[0],'w') as fout:
            fout.write('name\trepeat_count\tprobability\n')
            for seq in screed.open(input[0]):
                c_repeats = len(telomere.findall(seq.sequence[:region]))
                fout.write(f'{seq.name}\t{c_repeats}\t{binom.sf(c_repeats,region,0.25**6):.4f}\n')

rule count_scaffold_gaps:
    input:
        '{assembler}_{sample}/{animal}.scaffolds.fasta'
    output:
        'results/{animal}_{sample}_{assembler}.gaps.txt'
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
        asm = '{assembler}_{sample}/{animal}.contigs.fasta',
        reads = 'data/{animal}.cleaned.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{animal}.contigs.corrected.fasta'
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        '''
        ragtag.py correct {config[ref_genome]} {input.asm} -o {wildcards.assembler}_{wildcards.sample} -R {input.reads} -T corr -t {threads} --mm2-params "-x asm20"
        #mv {wildcards.assembler}_{wildcards.sample}/ragtag.contigs.corrected {output}
        '''

rule ragtag_scaffold:
    input:
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        '{assembler}_{sample}/{animal}.scaffolds.fasta'
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        '''
        ragtag.py scaffold {config[ref_genome]} {input} -o {wildcards.assembler}_{wildcards.sample} -t {threads} --mm2-params "-c -x asm5"
        mv {wildcards.assembler}_{wildcards.sample}/ragtag.scaffolds.fasta {output}
        '''

rule remap_reads:
    input:
        reads = 'data/{animal}.cleaned.hifi.fq.gz',
        asm = '{assembler}_{sample}/{animal}.scaffolds.fasta'
    output:
        '{assembler}_{sample}/{animal}_scaffolds_reads.sam'
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
    threads: 16
    resources:
        mem_mb = 4000,
        disk_scratch = 200
    shell:
        'samtools sort {input} -m 3900M -@ {threads} -T \$TMPDIR -o {output}'

rule prep_window:
    input:
        '{assembler}_{sample}/{animal}.scaffolds.fasta'
    output:
        fai = '{assembler}_{sample}/{animal}.scaffolds.fasta.fai',
        genome = '{assembler}_{sample}/{animal}.genome',
        bed = '{assembler}_{sample}/{animal}.windows.bed'
    params:
        10000
    shell:
        '''
        samtools faidx {input}
        awk -v OFS='\\t' {{'print $1,$2'}} {output.fai} > {output.genome}
        /cluster/work/pausch/alex/software/bedtools2/bin/windowMaker -g {output.genome} -w {params} > {output.bed}
        '''

rule window_coverage:
    input:
        windows = '{assembler}_{sample}/{animal}.windows.bed',
        bam = '{assembler}_{sample}/{animal}_scaffolds_reads.bam'
    output:
        'results/{animal}_{sample}_{assembler}.windows.coverage.txt'
    shell:
        'samtools bedcov {input.windows} {input.bam} > {output}'
        
rule chromosome_coverage:
    input:
        '{assembler}_{sample}/{animal}_scaffolds_reads.bam'
    output:
        'results/{animal}_{sample}_{assembler}.chrm.coverage.txt'
    shell:
        'samtools coverage {input} -o {output}'

rule raw_QC:
    input:
        'data/{animal}.{read_t}.hifi.fq.gz'
    output:
        'data/{animal}.{read_t}.QC.txt'
    shell:
        '''
        {workflow.basedir}/src/fasterqc {input} {output}
        '''

rule sample_data:
    input:
        'data/{animal}.cleaned.hifi.fq.gz'
    output:
        'data/{animal}.{sample}.hifi.fq.gz'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        if [ {wildcards.sample} -eq 100 ]
        then
            ln -s $(pwd)/{input} {output}
        else
            seqtk sample {input} $(bc <<<"scale=2;{wildcards.sample}/100") | pigz -p 4 > {output}
        fi
        '''

rule filter_data:
    input:
        'data/{animal}.raw.hifi.fq.gz'
    output:
        'data/{animal}.cleaned.hifi.fq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input} -o {output} --average_qual {config[filtering][avg_qual]} --length_required {config[filtering][min_length]} --thread {threads} --html data/{wildcards.animal}.html --json /dev/null'

checkpoint split_chromosomes:
    input:
        '{assembler}_{sample}/{animal}.scaffolds.fasta'
    output:
        directory('split_{animal}_{sample}_{assembler}')
    shell:
        '''
        mkdir -p {output}
        awk '$0 ~ "^>" {{ match($1, /^>([^:|\s]+)/, id); filename=id[1]}} {{print >> "{output}/"filename".chrm.fa"}}' {input}
        '''

rule repeat_masker:
    input:
        'split_{animal}_{sample}_{assembler}/{chunk}.chrm.fa'
    output:
        'split_{animal}_{sample}_{assembler}/{chunk}.chrm.fa.masked'
    threads: 8
    resources:
        mem_mb = 400,
        walltime =  '1:00'
    shell:
        'echo $(({threads}/2)); RepeatMasker -qq -xsmall -pa $(({threads}/2)) -species "Bos taurus" {input}'

def aggregate_chrm_input(wildcards):
    checkpoint_output = checkpoints.split_chromosomes.get(**wildcards).output[0]
    return expand(f'split_{wildcards.animal}_{wildcards.sample}_{wildcards.assembler}/{{chunk}}.chrm.fa.masked',chunk=glob_wildcards(os.path.join(checkpoint_output, '{chunk}.chrm.fa')).chunk)

rule merge_repeat_masked:
    input:
        aggregate_chrm_input
    output:
        '{assembler}_{sample}/{animal}.scaffolds.fasta.masked'
    shell:
        'cat {input} > {output}'
