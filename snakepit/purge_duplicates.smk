localrules: cut_and_split, purge_dups, assess_purging

rule map_reads:
    input:
        asm = '{assembler}_{sample}/{animal}.{haplotype}.contigs_raw.fa',
        reads = 'data/{animal}.{sample}.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{haplotype}/{animal}_read_aln.paf'
    params:
        '{assembler}_{sample}/{haplotype}'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        '''
        mkdir -p {params}
        minimap2 -I 6G -xmap-pb -t {threads} {input.asm} {input.reads} > {output}
        '''

rule cut_and_split:
    input:
        paf = '{assembler}_{sample}/{haplotype}/{animal}_read_aln.paf',
        contigs = '{assembler}_{sample}/{animal}.{haplotype}.contigs_raw.fa'
    output:
        splits = '{assembler}_{sample}/{haplotype}/{animal}.split.fasta'
    params:
        '{assembler}_{sample}/{haplotype}'
    shell:
        '''
        mkdir -p {params} && cd {params}
        {config[pd_root]}/bin/pbcstat ../../{input.paf}
        {config[pd_root]}/bin/calcuts PB.stat > cutoffs
        {config[pd_root]}/bin/split_fa ../../{input.contigs} > ../../{output.splits}
        '''

rule map_splits:
    input:
        '{assembler}_{sample}/{haplotype}/{animal}.split.fasta'
    output:
        '{assembler}_{sample}/{haplotype}/{animal}_ctg-aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xasm5 -t {threads} -DP {input} {input} > {output}'

rule purge_dups:
    input:
        paf = '{assembler}_{sample}/{haplotype}/{animal}_ctg-aln.paf',
        contigs = '{assembler}_{sample}/{animal}.{haplotype}.contigs_raw.fa'
    output:
        contigs = '{assembler}_{sample}/{haplotype}/{animal}.purged.fa',
        bed = temp('{assembler}_{sample}/{animal}_dups.bed')
    shell:
        '''
        {config[pd_root]}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > {output.bed}
        {config[pd_root]}/bin/get_seqs {output.bed} {input.contigs} -p {wildcards.assembler}_{wildcards.sample}/{wildcards.haplotype}/{wildcards.animal}
        '''

rule assess_purging:
    input:
        expand('{{assembler}}_{{sample}}/{{haplotype}}/{{animal}}.{progress}.matrix',progress=['purged','contigs_raw'])
    output:
        '{assembler}_{sample}/{animal}.{haplotype}.contigs.fasta'
    shell:
        'mv {wildcards.assembler}_{wildcards.sample}/{wildcards.haplotype}/{wildcards.animal}.purged.fa {output}'

rule KMC_reads:
    input:
        'data/{animal}.{sample}.hifi.fq.gz'
    output:
        base = '{assembler}_{sample}/{haplotype}/{animal}_kmer_reads',
        raw = multiext('{assembler}_{sample}/{haplotype}/{animal}_kmer_reads','','.kmc_pre','.kmc_suf')
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
        '{assembler}_{sample}/{haplotype}/{animal}.{progress}.fa'
    output:
        base = '{assembler}_{sample}/{haplotype}/{animal}_kmer_{progress}',
        raw = multiext('{assembler}_{sample}/{haplotype}/{animal}_kmer_{progress}','.kmc_pre','.kmc_suf')
    threads: 12
    resources:
        mem_mb = 3000,
        disk_scratch = 150,
        walltime = '2:00'
    params:
        memory = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']#,
    shell:
        '''
        {config[kmc_root]}/bin/kmc -ci0 -t{threads} -fm -m{params.memory} -k{config[k-mers]} {input} {output.base} $TMPDIR
        touch {output.base}
        '''

rule KMC_analysis:
    input:
        reads = '{assembler}_{sample}/{haplotype}/{animal}_kmer_reads',
        asm = '{assembler}_{sample}/{haplotype}/{animal}_kmer_{progress}'
    output:
        matrix = '{assembler}_{sample}/{haplotype}/{animal}.{progress}.matrix',
        plot = '{assembler}_{sample}/{haplotype}/{animal}.{progress}.spectra.png'
    threads: 2
    resources:
        mem_mb = 20000
    shell:
        '''
        {config[kmc_root]}/bin/kmc_tools analyze {input.reads} {input.asm} {output.matrix}
        python3 {config[kmc_root]}/spectra.py {output.matrix} {output.plot}
        '''
