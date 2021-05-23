localrules: count_telomers, count_scaffold_gaps, prep_window, window_coverage, chromosome_coverage, split_chromosomes, merge_masked_chromosomes

rule count_telomers:
    input:
        WORK_PATH + '{haplotype}.scaffolds.fasta'
    output:
        RESULT_PATH + '.telo.txt'
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
        WORK_PATH + '{haplotype}.scaffolds.fasta'
    output:
        RESULT_PATH + '.gaps.txt'
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
        asm = WORK_PATH + '{haplotype}.contigs.fasta',
        reads = 'data/offspring.cleaned.hifi.fq.gz'
    output:
        WORK_PATH + '{haplotype}.contigs.corrected.fasta'
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
        WORK_PATH + '{haplotype}.contigs.fasta'
    output:
        WORK_PATH + '{haplotype}.scaffolds.fasta'
    params:
        WORK_PATH + '{haplotype}_scaf'
    threads: 12
    resources:
        mem_mb = 3000
    shell:
        '''
        ragtag.py scaffold {config[ref_genome]} {input} -o {params} -t {threads} --unimap-params "-c -x asm5 -t {threads}" --aligner unimap -r -m 1000000
        cp {params}/ragtag.scaffold.fasta {output}
        '''

rule prep_window:
    input:
        get_dir('work','{haplotype}.scaffolds.fasta.fai')
    output:
        genome = WORK_PATH + '{haplotype}.genome',
        bed = WORK_PATH + '{haplotype}.windows.bed'
    params:
        10000
    shell:
        '''
        awk -v OFS='\\t' {{'print $1,$2'}} {input} > {output.genome}
        bedtools makewindows -g {output.genome} -w {params} > {output.bed}
        '''

rule window_coverage:
    input:
        windows = WORK_PATH + '{haplotype}.windows.bed',
        bam = WORK_PATH + '{haplotype}_scaffolds_{type}_reads.bam',
        bai = WORK_PATH + '{haplotype}_scaffolds_{type}_reads.bam.bai'
    output:
        RESULT_PATH + '.windows.{type}.coverage.txt'
    shell:
        'samtools bedcov {input.windows} {input.bam} > {output}'

rule chromosome_coverage:
    input:
        bam = WORK_PATH + '{haplotype}_scaffolds_{type}_reads.bam',
        bai = WORK_PATH + '{haplotype}_scaffolds_{type}_reads.bam.bai'
    output:
        RESULT_PATH + '.chrm.{type}.coverage.txt'
    shell:
        'samtools coverage {input.bam} -o {output}'

rule megadepth_coverage:
    input:
        multiext(get_dir('work','{haplotype}.{type}.{mapper}.bam'),'','.bai')
    output:
        get_dir('result','.coverage.{type}.{mapper}.all.bw')
    params:
        out = lambda wilcards, output: PurePath(output[0]).with_suffix('').with_suffix(''),
        opt = lambda wildcards: '--longreads' if wildcards.type == 'hifi' else ''
    threads: 4
    resources:
        mem_mb = 7000
    shell:
        'megadepth {input[0]} --threads {threads} {params.opt} --filter-out 260 --bigwig --prefix {params.out}'
        
checkpoint split_chromosomes:       
    input:
        WORK_PATH + '{haplotype}.scaffolds.fasta'
    output:
        directory(WORK_PATH + '{haplotype}_split_chrm'),
    params:
        asm = lambda wildcards, input: '../' + PurePath(input[0]).name,
        headers = 'headers.temp',
        chrm = 'chrm.temp',
        ur_tigs = 'unplaced_ref_contigs.chrm.fa',
        ua_tigs = 'unplaced_asm_contigs.chrm.fa',
    shell:
        '''
        mkdir -p {output} && cd {output}
        grep ">" {params.asm} | cut -c 2- > {params.headers}
        grep -P "^[YX\d]" {params.headers} | seqtk subseq {params.asm} - > {params.chrm}
        grep "^{config[ref_tig]}" {params.headers} | seqtk subseq {params.asm} - > {params.ur_tigs}
        grep -v -P "^([XY\d]|{config[ref_tig]})" {params.headers} | seqtk subseq {params.asm} - > {params.ua_tigs}
        awk '$0 ~ "^>" {{ match($1, /^>([^:|\s]+)/, id); filename=id[1]}} {{print >> filename".chrm.fa"}}' {params.chrm}

        for val in ref asm; do
            paste -d " " - - < unplaced_${{val}}_contigs.chrm.fa > collapsed_${{val}}.txt
            split -a 2 -d -C 50MiB --additional-suffix=.chrm.fa collapsed_${{val}}.txt unplaced_${{val}}_contigs_
            find unplaced_${{val}}_contigs_* -exec bash -c "cat {{}} | tr -s ' ' '\n' > {{}}.temp && mv {{}}.temp {{}}" \;
        done
        rm {params.headers} {params.chrm} {params.ua_tigs} {params.ur_tigs}
        '''

rule repeat_masker:
    input:
        WORK_PATH + '{haplotype}_split_chrm/{chunk}.chrm.fa'
    output:
        #NOTE repeatmasker doesn't output .masked if no masking, so just wrap the plain sequence via seqtk
        WORK_PATH + '{haplotype}_split_chrm/{chunk}.chrm.fa.masked'
    threads: 6#lambda wildcards, input: 18 if input.size_mb < 100 else 24
    resources:
        mem_mb = 1000,
        walltime = '14:00'
    shell:
        '''
        RepeatMasker -xsmall -pa $(({threads}/2)) -lib {config[repeat_library]} -qq -no_is {input} #-species "Bos taurus"
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
        masked = WORK_PATH + '{haplotype}.scaffolds.fasta.masked',
        csv = RESULT_PATH + '.repeats.csv',
        tbl = RESULT_PATH + '.repeats.tbl'
    params:
        out_list = lambda wildcards, input: [PurePath(fin).with_suffix('.out') for fin in input]
    shell:
        '''
        cat {input} > {output.masked}
        export PERL5LIB={config[perl_lib]}
        cat {params.out_list} | buildSummary.pl - > {output.tbl}
        python {workflow.basedir}/scripts/masker_table.py --haplotype {wildcards.haplotype} --sample {wildcards.sample} --assembler {wildcards.assembler}
        '''

rule TGS_gapcloser:
    input:
        scaffolds = get_dir('work','{haplotype}.scaffolds.fasta'),
        reads = lambda wildcards: f'data/{"sire" if wildcards.haplotype == "hap1" else "dam"}.hifi.fasta'
    output:
        get_dir('work','{haplotype}.filled.fasta')
    params:
        dir_ = get_dir('work','{haplotype}_TGS'),
        scaffolds = lambda wildcards, input: '../' + PurePath(input['scaffolds']).name,
        reads = lambda wildcards, input: '../../' + input['reads']
    threads: 12
    resources:
        mem_mb = 6000,
        walltime = '24:00'
    shell:
        '''
        mkdir -p {params.dir_}
        (cd {params.dir_} && {config[tgs_root]}/TGS-GapCloser.sh --scaff {params.scaffolds} --reads {params.reads} --output {wildcards.haplotype} --minmap_arg '-x asm20' --tgstype pb --ne --thread {threads})
        cp {params.dir_}/{params.out}.scaff_seqs {output}
        '''
