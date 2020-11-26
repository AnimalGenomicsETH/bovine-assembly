localrules: cut_and_split, purge_dups, assess_purging

rule map_reads:
    input:
        asm = '{assembler}_{sample}/{haplotype}.contigs_raw.fa',
        reads = 'data/offspring.{sample}.hifi.fq.gz'
    output:
        '{assembler}_{sample}/{haplotype}_pd/read_aln.paf'
    params:
        lambda wildcards, output: PurePath(output[0]).parent
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
        paf = '{assembler}_{sample}/{haplotype}_pd/read_aln.paf',
        asm = '{assembler}_{sample}/{haplotype}.contigs_raw.fa'
    output:
        splits = '{assembler}_{sample}/{haplotype}_pd/split.fasta',
        asm = '{assembler}_{sample}/{haplotype}_pd/contigs_raw.fa'
    params:
        lambda wildcards, output: PurePath(output[0]).parent
    shell:
        '''
        cp {input.asm} {output.asm}
        cd {params}
        {config[pd_root]}/bin/pbcstat ../../{input.paf}
        {config[pd_root]}/bin/calcuts PB.stat > cutoffs
        {config[pd_root]}/bin/split_fa ../../{input.asm} > ../../{output.splits}
        '''

rule map_splits:
    input:
        '{assembler}_{sample}/{haplotype}_pd/split.fasta'
    output:
        '{assembler}_{sample}/{haplotype}_pd/ctg-aln.paf'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'minimap2 -I 6G -xasm5 -t {threads} -DP {input} {input} > {output}'

rule purge_dups:
    input:
        paf = '{assembler}_{sample}/{haplotype}_pd/ctg-aln.paf',
        contigs = '{assembler}_{sample}/{haplotype}_pd/contigs_raw.fa'
    output:
        contigs = '{assembler}_{sample}/{haplotype}_pd/purged.fa',
        bed = '{assembler}_{sample}/{haplotype}_pd/dups.bed'
    params:
        dir_ = lambda wildcards, output: PurePath(output['bed']).parent,
        bed = lambda wildcards, output: PurePath(output['bed']).name,
        contigs = lambda wildcards, input: PurePath(input['contigs']).name
    shell:
        '''
        {config[pd_root]}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > {output.bed}
        cd {params.dir_}
        {config[pd_root]}/bin/get_seqs {params.bed} {params.contigs}
        '''

rule assess_purging:
    input:
        matrix = expand('canu_{{sample}}/{{haplotype}}_pd/{progress}.matrix',progress=['purged','contigs_raw']),
        purged = 'canu_{sample}/{haplotype}_pd/purged.fa'
    output:
        'canu_{sample}/{haplotype}.contigs.fasta'
    shell:
        'cp {input.purged} {output}'

rule KMC_reads:
    input:
        'data/offspring.{sample}.hifi.fq.gz'
    output:
        base = 'data/offspring.{sample}.kmer_reads',
        raw = multiext('data/offspring.{sample}.kmer_reads','','.kmc_pre','.kmc_suf')
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
        '{assembler}_{sample}/{haplotype}_pd/{progress}.fa'
    output:
        base = '{assembler}_{sample}/{haplotype}_pd/{haplotype}_kmer_{progress}',
        raw = multiext('{assembler}_{sample}/{haplotype}_pd/{haplotype}_kmer_{progress}','.kmc_pre','.kmc_suf')
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
        reads = 'data/offspring.{sample}.kmer_reads',
        asm = '{assembler}_{sample}/{haplotype}_pd/{haplotype}_kmer_{progress}'
    output:
        matrix = '{assembler}_{sample}/{haplotype}_pd/{progress}.matrix',
        plot = '{assembler}_{sample}/{haplotype}_pd/{progress}.spectra.png'
    params:
        reads = '{assembler}_{sample}/{haplotype}_pd/{haplotype}_kmer_reads'
    threads: 2
    resources:
        mem_mb = 20000
    shell:
        '''
        for ex in pre suf; do
          ln -sfn ../../{input.reads}.kmc_$ex {params.reads}.kmc_$ex
        done
        {config[kmc_root]}/bin/kmc_tools analyze {params.reads} {input.asm} {output.matrix}
        python3 {config[kmc_root]}/spectra.py {output.matrix} {output.plot}
        '''
