localrules: merqury_formatting, merqury_formatting_parental

rule count_asm_kmers:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        directory('{assembler}_{sample}/{haplotype}.contigs.meryl')
    threads: 8
    resources:
        mem_mb = 2500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule count_SR_reads:
    input:
        'data/{individual}.read_R{N}.SR.fq.gz'
    output:
        directory('data/{individual}.read_R{N}.meryl')
    threads: 24
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merge_SR_reads:
    input:
        expand('data/{{individual}}.read_R{N}.meryl', N = (1,2))
    output:
        directory('data/{individual}.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule generate_hapmers:
    input:
        dam = 'data/dam.meryl',
        sire = 'data/sire.meryl',
        child = 'data/offspring.meryl'
    output:
        directory(expand('data/{parent}.hapmer.meryl', parent = ('dam', 'sire')))
    params:
        dam = lambda wildcards,input: PurePath(input['dam']).name,
        sire = lambda wildcards,input: PurePath(input['sire']).name,
        child = lambda wildcards,input: PurePath(input['child']).name
    threads: 12
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd data
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hapmers.sh {params.sire} {params.dam} {params.child}
        '''

rule hap_blob:
    input:
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    output:
        '{assembler}_{sample}/{hap}.hapmers.blob.png'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 8
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hap_blob.sh {params.hapmers} {params.asm} {params.out}
        '''

rule merqury_spectra_cn:
    input:
        read_db = 'data/offspring.meryl',
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        asm_dbs = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.meryl', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        filt = 'data/offspring.filt'
    output:
        multiext('{assembler}_{sample}/{hap}', '.qv', '.completeness.stats')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        base = lambda wildcards, input: PurePath(input['read_db']).stem,
        out = '{hap}',
        read_db = lambda wildcards, input: PurePath(input['read_db']).name,
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        ln -sfn  ../{input.filt}
        ln -sfn ../{input.read_db}
        find ../data/ -name "{params.base}.gt*.meryl" -exec ln -sfn ../data/{{}} . \;
        $MERQURY/eval/spectra-cn.sh {params.read_db} {params.asm} {params.out}
        '''

rule merqury_spectra_hap:
    input:
        reads = 'data/offspring.meryl',
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        stats = '{assembler}_{sample}/{hap}.completeness.stats'
    output:
        '{assembler}_{sample}/{hap}.hap.completeness.stats'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}.hap',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 12
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/spectra-hap.sh ../{input.reads} {params.hapmers} {params.asm} {params.out}
        '''

rule merqury_phase_block:
    input:
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta',
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam'))
    output:
        multiext('{assembler}_{sample}/{haplotype}.100_20000.','phased_block.bed','switch.bed', 'switches.txt')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{haplotype}',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: PurePath(input['asm']).name
    threads: 12
    resources:
        mem_mb = 3500
    envmodules:
       'gcc/8.2.0',
       'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/phase_block.sh {params.asm} {params.hapmers} {params.out}
        '''

rule merqury_block_n_stats:
    input:
        asm_block = multiext('{assembler}_{sample}/{haplotype}','.contigs.fasta','.100_20000.phased_block.bed')
    output:
        multiext('{assembler}_{sample}/{haplotype}.{haplotype}.contigs','.continuity.NG.png','.block.NG.png')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{haplotype}',
        asm_block = lambda wildcards,input: tuple(PurePath(i).name for i in input['asm_block'])
    threads: 6
    resources:
        mem_mb = 2500
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/block_n_stats.sh {params.asm_block} {params.out} $( echo "({config[genome_est]}*1000000000)/1" | bc)
        '''

rule merqury_prep:
    input:
        'data/offspring.meryl'
    output:
        hist = 'data/offspring.hist',
        filt = 'data/offspring.filt',
    threads: 12
    resources:
        mem_mb = 2000
    shell:
        '''
        meryl histogram {input} > {output.hist}
        java -jar -Xmx1g {config[merqury_root]}/eval/kmerHistToPloidyDepth.jar {output.hist} > {output.hist}.ploidy
        filt=`sed -n 2p {output.hist}.ploidy | awk '{{print $NF}}'`
        echo $filt > {output.filt}
        meryl greater-than $filt output data/offspring.gt$filt.meryl {input}
        '''

rule merqury_formatting:
    input:
        '{assembler}_{sample}/{haplotype}.{haplotype}.contigs.continuity.NG.png',
        stats = lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.completeness.stats', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        hap_stats = lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.hap.completeness.stats', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        qv = lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.qv', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        blob = lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.hapmers.blob.png', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        switches = '{assembler}_{sample}/{haplotype}.100_20000.switches.txt'
    output:
        'results/{haplotype}_{sample}_{assembler}.merqury.stats.txt'
    shell:
        '''
        awk '/{wildcards.haplotype}/ {{print "QV "$4}}' {input.qv} > {output}
        awk '/{wildcards.haplotype}/ {{print "completeness "$5}}' {input.stats} >> {output}
        awk '/{wildcards.haplotype}/ && /sire/ {{print "sire "$5}}' {input.hap_stats} >> {output}
        awk '/{wildcards.haplotype}/ && /dam/ {{print "dam "$5}}' {input.hap_stats} >> {output}
        awk '{{print "switches "$14+0}}' {input.switches} >> {output}
        '''

rule merqury_formatting_parental:
    input:
        stats = '{assembler}_{sample}/{parent}.completeness.stats',
        qv = '{assembler}_{sample}/{parent}.qv',
    output:
        'results/{parent}_{sample}_{assembler}.merqury.stats.simple'
    shell:
        '''
        awk '/{wildcards.parent}/ {{print "QV "$4}}' {input.qv} > {output}
        awk '/{wildcards.parent}/ {{print "completeness "$5}}' {input.stats} >> {output}
        '''
