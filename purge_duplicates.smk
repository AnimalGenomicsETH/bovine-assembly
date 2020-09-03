localrules: cut_and_split, purge_dups

rule map_reads:
    input:
        asm = 'canu/{animal}.contigs_raw.fasta',
        reads = 'data/{animal}.hifi.fq.gz'
    output:
        'canu/{animal}_read_aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xmap-pb -t {threads} {input.asm} {input.reads} > {output}'
 
rule cut_and_split:
    input:
        paf = 'canu/{animal}_read_aln.paf',
        contigs = 'canu/{animal}.contigs_raw.fasta'
    output:
        splits = 'canu/{animal}.split.fasta'
    params:
        config['pd_root']
    shell:
        '''
        {params}/bin/pbcstat {input.paf} -O canu
        {params}/bin/calcuts PB.stat > cutoffs -O canu
        {params}/bin/split_fa {input.contigs} > {output.splits}
        '''

rule map_splits:
    input:
        'canu/{animal}.split.fasta'
    output:
        'canu/{animal}_ctg-aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xasm5 -t {threads} -DP {input} {input} > {output}'

rule purge_dups:
    input:
        paf = 'canu/{animal}_ctg-aln.paf',
        contigs = 'canu/{animal}.contigs_raw.fasta'
    output:
        'canu/{animal}.contigs.fasta'
        #rename contigs
        #'canu/{animal}.spectra.png'
    shell:
        '''
        {config["pd_root"]}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > canu/dups.bed
        {config["pd_root"]}/bin/get_seqs canu/dups.bed {input.contigs} -p canu/{wildcards.animal}
        touch {output}
        '''

rule KMC_reads:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        'canu/{animal}_kmer_reads'
    shell:
        '{config["kmc_root"}/bin/kmc -ci0 -cs5000 -t {threads} -m {params.memory} -k {config["kmer"]} {input} {output}'

rule KMC_ref:
    input:
        'canu/{animal}.contigs.fasta'
    output:
        'canu/{animal}_kmer_asm'
    shell:
        '{config["kmc_root"]}/bin/kmc -ci0 -t {threads} -fm -m {params.memory} -k {config["kmer"]} {input} {output}'

rule KMC_analysis:
    input:
        reads = 'canu/{animal}_kmer_reads',
        asm = 'canu/{animal}_kmer_asm'
    output:
        matrix = 'canu/{animal}.matrix',
        plot = 'canu/{animal}.spectra.png'
    shell:
        '''
        {config["kmc_root"]}/bin/kmc_tools analyze {input.reads} {inputs.asm} {output.matrix}
        python3 {config["kmc_root"]}/spectra.py {output.matrix} {output.plot}
        '''
