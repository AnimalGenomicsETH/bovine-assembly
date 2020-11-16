localrules: count_telomers, count_scaffold_gaps, prep_window, window_coverage, chromosome_coverage, split_chromosomes, merge_masked_chromosomes

rule count_telomers:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        'results/{haplotype}_{sample}_{assembler}.telo.txt'
    run:
        import re, screed
        from scipy.stats import binom
        region = 1000
        telomere = re.compile("TTAGGG", re.IGNORECASE)

        with open(output[0],'w') as fout:
            fout.write('name\trepeat_count\tprobability\n')
            for seq in screed.open(input[0]):
                c_repeats = len(telomere.findall(seq.sequence[region:]))
                fout.write(f'{seq.name}\t{c_repeats}\t{binom.sf(c_repeats,region,0.25**6):.4f}\n')

rule count_scaffold_gaps:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        'results/{haplotype}_{sample}_{assembler}.gaps.txt'
    run:
        import re, screed
        gap_sequence = re.compile(r'[nN]+')
        with open(output[0],'w') as fout:
            fout.write('scaffold\tgaps\tlengths\twidths\n')
            for scaffold in screed.open(input[0]):
                contig_lengths = list(map(len, gap_sequence.split(str(scaffold.sequence))))
                gap_lengths = list(map(len, gap_sequence.findall(str(scaffold.sequence))))
                #for F in (re.split, re.findall):
                #     map(lambda s: str(len(s)), F(scaffold.sequence))
                fout.write('\t'.join((scaffold.name,str(len(contig_lengths)-1),','.join(map(str,contig_lengths)),','.join(map(str,gap_lengths))))+'\n')

#ragtag.py correct, full reads mapping, c_mapping
rule ragtag_correct:
    input:
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta',
        reads = 'data/offspring.cleaned.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{haplotype}.contigs.corrected.fasta'
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        '''
        ragtag.py correct {config[ref_genome]} {input.asm} -o {wildcards.assembler}_{wildcards.sample} -R {input.reads} -T corr -t {threads} --mm2-params "-c -x asm20"
        #mv {wildcards.assembler}_{wildcards.sample}/ragtag.contigs.corrected {output}
        '''

rule ragtag_scaffold:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    params:
        '{assembler}_{sample}/{haplotype}'
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        '''
        ragtag.py scaffold {config[ref_genome]} {input} -o {params} -t {threads} --mm2-params "-c -x asm5" -r -m 1000000
        mv {params}/ragtag.scaffolds.fasta {output}
        '''

rule prep_window:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        fai = '{assembler}_{sample}/{haplotype}.scaffolds.fasta.fai',
        genome = '{assembler}_{sample}/{haplotype}.genome',
        bed = '{assembler}_{sample}/{haplotype}.windows.bed'
    params:
        10000
    shell:
        '''
        samtools faidx {input}
        awk -v OFS='\\t' {{'print $1,$2'}} {output.fai} > {output.genome}
        bedtools makewindows -g {output.genome} -w {params} > {output.bed}
        '''

rule index_bam:
    input:
        '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam'
    output:
        '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam.bai'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule window_coverage:
    input:
        windows = '{assembler}_{sample}/{haplotype}.windows.bed',
        bam = '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam',
        bai = '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam.bai'
    output:
        'results/{haplotype}_{sample}_{assembler}.windows.{type}.coverage.txt'
    shell:
        'samtools bedcov {input.windows} {input.bam} > {output}'

rule chromosome_coverage:
    input:
        bam = '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam',
        bai = '{assembler}_{sample}/{haplotype}_scaffolds_{type}_reads.bam.bai'
    output:
        'results/{haplotype}_{sample}_{assembler}.chrm.{type}.coverage.txt'
    shell:
        'samtools coverage {input.bam} -o {output}'

checkpoint split_chromosomes:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        directory('{assembler}_{sample}/{haplotype}_split_chrm'),
    params:
        headers = 'headers.temp',
        chrm = 'chrm.temp',
        ur_tigs = 'unplaced_ref_contigs.chrm.fa',
        ua_tigs = 'unplaced_asm_contigs.chrm.fa',
    shell:
        '''
        mkdir -p {output}
        grep ">" {input} | cut -c 2- > {output}/{params.headers}
        grep "{config[ref_chrm]}" {output}/{params.headers} | seqtk subseq {input} - > {output}/{params.chrm}
        grep "{config[ref_tig]}" {output}/{params.headers} | seqtk subseq {input} - > {output}/{params.ur_tigs}
        grep -v -e "{config[ref_chrm]}" -e "{config[ref_tig]}" {output}/{params.headers} | seqtk subseq {input} - > {output}/{params.ua_tigs}
        awk '$0 ~ "^>" {{ match($1, /^>([^:|\s]+)/, id); filename=id[1]}} {{print >> "{output}/"filename".chrm.fa"}}' {output}/{params.chrm}

        for val in ref asm; do
            paste -d " " - - < unplaced_$val_contigs.chrm.fa > {output}/collapsed_$val.txt
            split -a 2 -d -C 50MiB --additional-suffix=.chrm.fa {output}/collapsed_$val.txt {output}/unplaced_$val_contigs_
            find {output}/unplaced_$val_contigs_* -exec bash -c "cat {{}} | tr -s ' ' '\n' > {{}}.temp && mv {{}}.temp {{}}" \;
        done
        #split -a 2 -d -l 50 --additional-suffix=.chrm.fa {output}/{params.ua_tigs} {output}/unplaced_asm_contigs_
        #split -a 2 -d -l 50 --additional-suffix=.chrm.fa {output}/{params.ur_tigs} {output}/unplaced_ref_contigs_
        rm {output}/{params.headers} {output}/{params.chrm} {output}/{params.ua_tigs} {output}/{params.ur_tigs}
        '''

rule repeat_masker:
    input:
        '{assembler}_{sample}/{haplotype}_split_chrm/{chunk}.chrm.fa'
    output:
        #NOTE repeatmasker doesn't output .masked if no masking, so just wrap the plain sequence via seqtk
        '{assembler}_{sample}/{haplotype}_split_chrm/{chunk}.chrm.fa.masked'
    threads: 8
    resources:
        mem_mb = 600,
        walltime =  '2:00'
    shell:
        '''
        RepeatMasker -qq -xsmall -pa $(({threads}/2)) -species "Bos taurus" {input}
        if [ ! -f {output} ]; then
          seqtk seq -l60 {input} > {output}
        fi
        '''

def aggregate_chrm_input(wildcards):
    checkpoint_output = checkpoints.split_chromosomes.get(**wildcards).output[0]
    return expand('{fpath}/{chunk}.chrm.fa.masked',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('{chunk}.chrm.fa')).chunk)

rule merge_masked_chromosomes:
    input:
        aggregate_chrm_input
    output:
        masked = '{assembler}_{sample}/{haplotype}.scaffolds.fasta.masked',
        csv = 'results/{haplotype}_{sample}_{assembler}.repeats.csv'
    shell:
        '''
        cat {input} > {output.masked}
        python {workflow.basedir}/scripts/masker_table.py --haplotype {wildcards.haplotype} --sample {wildcards.sample} --assembler {wildcards.assembler}
        '''

rule TGS_gapcloser:
    input:
        scaffolds = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        reads = lambda wildcards: expand('data/{parent}.hifi.fasta', parent = 'sire' if wildcards.haplotype == 'hap1' else 'dam')
    output:
        '{assembler}_{sample}/{haplotype}.scaff_seq'
    params:
        out = '{assembler}_{sample}/{haplotype}'
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        '''
        {config[tgs_root]}/TGS-GapCloser.sh --scaff {input.scaffolds} --reads {input.reads} --output {params.out} --minmap_arg '-x asm20' --tgstype pb --ne --thread {threads}
        '''

rule polish_scaffolds:
    input:
        scaffolds = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        aln = '{assembler}_{sample}/{haplotype}_scaffolds_reads.sam',
        reads = 'data/offspring.{sample}.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{haplotype}.polished.fasta'
    threads: 16
    resources:
        mem_mb = 28000,
        walltime = '2:30'
    shell:
        'racon -t {threads} {input.reads} {input.aln} {input.scaffolds} > {output}'
