localrules: split_reads, merqury_formatting

rule merge_SR_reads:
    input:
        expand('data/{{parent}}.read_R{N}.SR.meryl', N = (1,2))
    output:
        directory('data/{parent}.SR.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

#NOTE technically should be w.r.t. the read sampling?
rule generate_hapmers:
    input:
        dam = 'data/dam.{read_t}.meryl',
        sire = 'data/sire.{read_t}.meryl',
        child = 'data/reads.100.hifi.meryl'
    output:
        directory(expand('data/{parent}.{{read_t}}.hapmer.meryl', parent = ('dam', 'sire')))
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
        hapmers = expand('data/{parent}.{{read_t}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    output:
        '{assembler}_{sample}/{hap}.{read_t}.hapmers.blob.png'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}.{read_t}',
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
        read_db = 'data/reads.{sample}.hifi.meryl',
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        asm_dbs = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.meryl', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        filt = 'data/reads.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{hap}.{read_t}', '.qv', '.completeness.stats')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}.{read_t}',
        read_db = lambda wildcards, input: '..' / PurePath(input['read_db']),
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2',
        'igv/2.8.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        ln -sfn ../{input.filt}
        find ../data/ -name "*gt*" -exec ln -sfn ../data/{{}} . \;
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} {params.asm} {params.out}
        '''

rule spectra_hap:
    input:
        reads = 'data/reads.{sample}.hifi.meryl',
        hapmers = expand('data/{parent}.{{read_t}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{X}.contigs.fasta', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        stats = multiext('{assembler}_{sample}/{hap}.{read_t}','.qv','.completeness.stats')
    output:
        '{assembler}_{sample}/{hap}.{read_t}.hap.completeness.stats'
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}.{read_t}',
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
        hapmers = expand('data/{parent}.SR.hapmer.meryl',parent=('sire','dam'))
    output:
        multiext('{assembler}_{sample}/{haplotype}.100_20000.','phased_block.bed','switch.bed')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{haplotype}',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: PurePath(input['asm']).name
    threads: 12
    resources:
        mem_mb = 6000
    envmodules:
       'gcc/8.2.0',
       'r/4.0.2',
       'igv/2.8.2'
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

checkpoint split_reads:
    input:
        'data/{data}.{modifier}.{read_t}.fq.gz'
    output:
        directory('split_{data}_{modifier}_{read_t}')
    params:
        lambda wildcards: config["split_size"][wildcards.read_t]
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        mkdir -p {output}
        zcat {input} | split -a 2 -d -l {params} --filter='pigz -p 6 > $FILE.fq.gz' - {output}/chunk_
        '''

rule count_many:
    input:
        'split_{data}_{modifier}_{read_t}/chunk_{chunk}.fq.gz'
    output:
        temp(directory('split_{data}_{modifier}_{read_t}/chunk_{chunk}.meryl'))
    threads: lambda wildcards: 18 if wildcards.read_t == 'hifi' else 4
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

def aggregate_split_input(wildcards):
    checkpoint_output = checkpoints.split_reads.get(**wildcards).output[0]
    return expand('split_{data}_{modifier}_{read_t}/chunk_{chunk}.meryl',data=wildcards.data,modifier=wildcards.modifier,read_t=wildcards.read_t,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('chunk_{chunk}.fq.gz')).chunk)

rule merge_many:
    input:
        aggregate_split_input
    output:
        directory('data/{data}.{modifier}.{read_t}.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule count_asm_kmers:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        directory('{assembler}_{sample}/{haplotype}.contigs.meryl')
    threads: 12
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merqury_prep:
    input:
        'data/reads.{sample}.hifi.meryl'
    output:
        hist = 'data/reads.{sample}.hifi.hist',
        filt = 'data/reads.{sample}.hifi.filt'
    threads: 12
    resources:
        mem_mb = 2000
    shell:
        '''
        meryl histogram {input} > {output.hist}
        java -jar -Xmx1g {config[merqury_root]}/eval/kmerHistToPloidyDepth.jar {output.hist} > {output.hist}.ploidy
        filt=`sed -n 2p {output.hist}.ploidy | awk '{{print $NF}}'`
        echo $filt > {output.filt}
        meryl greater-than $filt output data/{wildcards.sample}.hifi.gt$filt.meryl {input}
        '''

rule merqury_formatting:
    input:
        '{assembler}_{sample}/{haplotype}.{haplotype}.contigs.continuity.NG.png',
        lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.SR.completeness.stats', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.SR.hap.completeness.stats', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
        lambda wildcards: expand('{{assembler}}_{{sample}}/{hap}.SR.hapmers.blob.png', hap = 'asm' if wildcards.haplotype == 'asm' else 'trio'),
    output:
        'results/{haplotype}_{sample}_{assembler}.merqury.stats.txt'
    shell:
        'touch {output}'
