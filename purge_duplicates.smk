localrules: cut_and_split, purge_dups, assess_purging, KMC_analysis

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
        {params}/bin/calcuts canu/PB.stat > canu/cutoffs
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
        contigs = 'canu/{animal}.purged.fa',
    params:
        config['pd_root']
    shell:
        '''
        {params}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > canu/dups.bed
        {params}/bin/get_seqs canu/dups.bed {input.contigs} -p canu/{wildcards.animal}
        '''

rule assess_purging:
    input:
        'canu/{animal}.matrix'
    output:
        'canu/{animal}.contigs.fasta'
    shell:
        'mv canu/{wildcards.animal}.purged.fa {output}'

rule KMC_reads:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        multiext('canu/{animal}_kmer_reads','','kmc_pre','kmc_suf')
    threads: 12
    resources:
        mem_mb = 3000,
        disk_scratch = 15000
    params:
        out = lambda wildcards: f'canu/{wildcards.animal}_kmer_reads',
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj'],
        root = config['kmc_root'],
        K = config['k-mers']
    shell:
        '''
        {params.root}/bin/kmc -ci0 -cs5000 -t{threads} -fq -m{params.memory} -k{params.K} {input} {params.out} $TMPDIR
        touch canu/{wildcards.animal}_kmer_reads
        '''

rule KMC_ref:
    input:
        'canu/{animal}.purged.fa'
    output:
        multiext('canu/{animal}_kmer_asm','','kmc_pre','kmer_suf')
    threads: 12
    resources:
        mem_mb = 3000,
        disk_scratch = 15000
    params:
        out = lambda wildcards: f'canu/{wildcards.animal}_kmer_asm',
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj'],
        root = config['kmc_root'],
        K = config['k-mers']
    shell:
        '''
        {params.root}/bin/kmc -ci0 -t{threads} -fm -m{params.memory} -k{params.K} {input} {params.out} $TMPDIR
        touch canu/{wildcards.animal}_kmer_asm
        '''

rule KMC_analysis:
    input:
        reads = 'canu/{animal}_kmer_reads',
        asm = 'canu/{animal}_kmer_asm'
    output:
        matrix = 'canu/{animal}.matrix',
        plot = 'canu/{animal}.spectra.png'
    params:
        config['kmc_root']
    shell:
        '''
        {params}/bin/kmc_tools analyze {input.reads} {input.asm} {output.matrix}
        python3 {params}/spectra.py {output.matrix} {output.plot}
        '''
