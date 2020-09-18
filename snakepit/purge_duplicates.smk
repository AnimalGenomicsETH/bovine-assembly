localrules: cut_and_split, purge_dups, assess_purging

rule map_reads:
    input:
        asm = '{assembler}/{animal}.contigs_raw.fa',
        reads = 'data/{animal}.hifi.fq.gz'
    output:
        '{assembler}/{animal}_read_aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xmap-pb -t {threads} {input.asm} {input.reads} > {output}'

rule cut_and_split:
    input:
        paf = '{assembler}/{animal}_read_aln.paf',
        contigs = '{assembler}/{animal}.contigs_raw.fa'
    output:
        splits = '{assembler}/{animal}.split.fasta'
    params:
        config['pd_root']
    shell:
        '''
        {config[pd_root]}/bin/pbcstat {input.paf} -O {wildcards.assembler}
        {config[pd_root]}/bin/calcuts {wildcards.assembler}/PB.stat > {wildcard.assembler}/cutoffs
        {config[pd_root]}/bin/split_fa {input.contigs} > {output.splits}
        '''

rule map_splits:
    input:
        '{assembler}/{animal}.split.fasta'
    output:
        '{assembler}/{animal}_ctg-aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xasm5 -t {threads} -DP {input} {input} > {output}'

rule purge_dups:
    input:
        paf = '{assembler}/{animal}_ctg-aln.paf',
        contigs = '{assembler}/{animal}.contigs_raw.fa'
    output:
        contigs = '{assembler}/{animal}.purged.fa',
        bed = '{assembler}/{animal}_dups.bed'
    shell:
        '''
        {config[pd_root]}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > {output.bed}
        {config[pd_root]}/bin/get_seqs {output.bed} {input.contigs} -p {wildcards.assembler}/{wildcards.animal}
        '''

rule assess_purging:
    input:
        expand('{{assembler}}/{{animal}}.{sample}.matrix',sample=['purged','contigs_raw'])
    output:
        '{assembler}/{animal}.contigs.fasta'
    shell:
        'mv {wildcards.assembler}/{wildcards.animal}.purged.fa {output}'

rule KMC_reads:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        base = '{assembler}/{animal}_kmer_reads',
        raw = multiext('{assembler}/{animal}_kmer_reads','','.kmc_pre','.kmc_suf')
    threads: 12
    resources:
        mem_mb = 4000,
        disk_scratch = 150
    params:
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']#,
    shell:
        '''
        {config[kmc_root]}/bin/kmc -ci0 -cs1023 -t{threads} -fq -m{params.memory} -k{config[k-mers]} {input} {output.base} $TMPDIR
        touch {output.base}
        '''

rule KMC_ref:
    input:
        '{assembler}/{animal}.{sample}.fa'
    output:
        base = '{assembler}/{animal}_kmer_{sample}',
        raw = multiext('{assembler}/{animal}_kmer_{sample}','.kmc_pre','.kmc_suf')
    threads: 12
    resources:
        mem_mb = 3000,
        disk_scratch = 150,
        walltime = '20'
    params:
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']#,
    shell:
        '''
        {config[kmc_root]}/bin/kmc -ci0 -t{threads} -fm -m{params.memory} -k{config[k-mers]} {input} {output.base} $TMPDIR
        touch {output.base}
        '''

rule KMC_analysis:
    input:
        reads = '{assembler}/{animal}_kmer_reads',
        asm = '{assembler}/{animal}_kmer_{sample}'
    output:
        matrix = '{assembler}/{animal}.{sample}.matrix',
        plot = '{assembler}/{animal}.{sample}.spectra.png'
    threads: 2
    resources:
        mem_mb = 20000
    shell:
        '''
        {config[kmc_root]}/bin/kmc_tools analyze {input.reads} {input.asm} {output.matrix}
        python3 {config[kmc_root]}/spectra.py {output.matrix} {output.plot}
        '''
