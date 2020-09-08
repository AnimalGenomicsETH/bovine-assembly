localrules: cut_and_split, purge_dups, assess_purging

rule map_reads:
    input:
        asm = 'canu/{animal}.contigs_raw.fa',
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
        contigs = 'canu/{animal}.contigs_raw.fa'
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
        contigs = 'canu/{animal}.contigs_raw.fa'
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
        expand('canu/{{animal}}.{sample}.matrix',sample=['purged','contigs_raw'])
    output:
        'canu/{animal}.contigs.fasta'
    shell:
        'mv canu/{wildcards.animal}.purged.fa {output}'

rule KMC_reads:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        base = 'canu/{animal}_kmer_reads',
        raw = multiext('canu/{animal}_kmer_reads','','.kmc_pre','.kmc_suf')
    threads: 12
    resources:
        mem_mb = 4000,
        disk_scratch = 150
    params:
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj'],
        root = config['kmc_root'],
        K = config['k-mers']
    shell:
        '''
        {params.root}/bin/kmc -ci0 -cs1023 -t{threads} -fq -m{params.memory} -k{params.K} {input} {output.base} $TMPDIR
        touch {output.base}
        '''

rule KMC_ref:
    input:
        'canu/{animal}.{sample}.fa'
    output:
        base = 'canu/{animal}_kmer_{sample}',
        raw = multiext('canu/{animal}_kmer_{sample}','.kmc_pre','.kmc_suf')
    threads: 12
    resources:
        mem_mb = 3000,
        disk_scratch = 150,
        walltime = '20'
    params:
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj'],
        root = config['kmc_root'],
        K = config['k-mers']
    shell:
        '''
        {params.root}/bin/kmc -ci0 -t{threads} -fm -m{params.memory} -k{params.K} {input} {output.base} $TMPDIR
        touch {output.base}
        '''

rule KMC_analysis:
    input:
        reads = 'canu/{animal}_kmer_reads',
        asm = 'canu/{animal}_kmer_{sample}'
    output:
        matrix = 'canu/{animal}.{sample}.matrix',
        plot = 'canu/{animal}.{sample}.spectra.png'
    threads: 2
    resources:
        mem_mb = 20000
    params:
        config['kmc_root']
    shell:
        '''
        {params}/bin/kmc_tools analyze {input.reads} {input.asm} {output.matrix}
        python3 {params}/spectra.py {output.matrix} {output.plot}
        '''
