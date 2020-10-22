localrules: split_reads, merqury_formatting

rule count_SR_reads:
    input:
        'data/{parent}_R{N}.fastq.gz'
    output:
        directory('data/{parent}.read_R{N}.meryl')
    threads: 6
    resources:
        mem_mb = 10000,
        walltime = '12:00'
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merge_SR_reads:
    input:
        expand('data/{{parent}}.read_R{N}.meryl',N=(1,2))
    output:
        directory('data/{parent}.meryl')
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
        child = 'data/{animal}.{sample}.hifi.meryl'
    output:
        expand('data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('dam','sire'))
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
        $MERQURY/trio/hapmers.sh ../{input.sire} ../{input.dam} ../{input.child}
        mv dam.hapmer.meryl dam.{wildcards.animal}.{wildcards.sample}.hapmer.meryl
        mv sire.hapmer.meryl sire.{wildcards.animal}.{wildcards.sample}.hapmer.meryl
        '''

rule hap_blob:
    input:
        hapmers = expand('data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    output:
        '{assembler}_{sample}/{animal}.{hap}.hapmers.blob.png'
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{hap}',
        hapmers = expand('../data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{animal}}.{X}.contigs.fasta',X='asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
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
        read_db = 'data/{animal}.{sample}.hifi.meryl',
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        asm_dbs = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.contigs.meryl',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        filt = 'data/{animal}.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{animal}.{hap}','.qv','.completeness.stats')
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{hap}',
        asm = lambda wildcards: expand('{{animal}}.{X}.contigs.fasta',X = ('asm' if wildcards.hap == 'asm' else ('hap1','hap2'))),
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
        reads = 'data/{animal}.{sample}.hifi.meryl',
        hapmers = expand('data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2')),
        stats = multiext('{assembler}_{sample}/{animal}.{hap}','.qv','.completeness.stats')    
    output:
        '{assembler}_{sample}/{animal}.{hap}.completeness.stats'
        #lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{{hap}}.{{animal}}.{X}.contigs.dam.{{animal}}.{{sample}}.hapmer.spectra-cn.ln.png',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{hap}',
        hapmers = expand('../data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards,input: list(PurePath(i).name for i in input['asm'])#('/'expand('{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
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
        asm = '{assembler}_{sample}/{animal}.{haplotype}.contigs.fasta',
        hapmers = expand('data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam'))
    output:
        '{assembler}_{sample}/{animal}.{haplotype}.100_20000.phased_block.bed'
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{haplotype}',
        asm = '{animal}.{haplotype}.contigs.fasta',
        hapmers = expand('../data/{parent}.{{animal}}.{{sample}}.hapmer.meryl',parent=('sire','dam'))
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

rule merqury_block_n_stats_new:
    input:
        asm_block = multiext('{assembler}_{sample}/{animal}.{haplotype}','.contigs.fasta','.100_20000.phased_block.bed')
    output:
        multiext('{assembler}_{sample}/{animal}.{haplotype}.{animal}.{haplotype}.contigs','.continuity.NG.png','.block.NG.png')
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{haplotype}',
        asm_block = multiext('{animal}.{haplotype}','.contigs.fasta','.100_20000.phased_block.bed')
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

rule merqury_block_n_stats:
    input:
        asms_blocks = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.{ext}',ext = ('.contigs.fasta','.100_20000.phased_block.bed'), X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))#,
        #block = lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    output:
        #fake file from merqury
        '{assembler}_{sample}/{animal}.{hap}.scaff.complete'
        #lambda wildcards: expand('{{assembler}}_{{sample}}/{{animal}}.{X}.scaff.sizes', X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    params:
        dir_ = '{assembler}_{sample}',
        out = '{animal}.{hap}',
        asms_blocks = lambda wildcards: expand('{{animal}}.{X}.{ext}',ext = ('.contigs.fasta','.100_20000.phased_block.bed'), X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
        #asm = lambda wildcards: expand('{{animal}}.{X}.contigs.fasta',X = 'asm' if wildcards.hap == 'asm' else ('hap1','hap2'))
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {param.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/block_n_stats.sh {params.asms_blocks} {params.out} {config[genome_est]}
        '''

checkpoint split_reads:
    input:
        'data/{animal}.{sample}.hifi.fq.gz'
    output:
        directory('split_{animal}_{sample}/')
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        #can replace with {output}
        mkdir -p split_{wildcards.animal}_{wildcards.sample}
        zcat {input} | split -a 2 -d -C {config[split_size]}GiB --filter='pigz -p 6 > $FILE.fq.gz' - split_{wildcards.animal}_{wildcards.sample}/chunk_
        '''

rule count_many:
    input:
        'split_{animal}_{sample}/chunk_{chunk}.fq.gz'
    output:
        directory('split_{animal}_{sample}/chunk_{chunk}.meryl')
    threads: 18
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

def aggregate_split_input(wildcards):
    checkpoint_output = checkpoints.split_reads.get(**wildcards).output[0]
    return expand('split_{animal}_{sample}/chunk_{chunk}.meryl',animal=wildcards.animal,sample=wildcards.sample,chunk=glob_wildcards(os.path.join(checkpoint_output, 'chunk_{chunk}.fq.gz')).chunk)

rule merge_many:
    input:
        aggregate_split_input
    output:
        directory('data/{animal}.{sample}.hifi.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule count_asm_kmers:
    input:
        '{assembler}_{sample}/{animal}.{haplotype}.contigs.fasta'
    output:
        directory('{assembler}_{sample}/{animal}.{haplotype}.contigs.meryl')
    threads: 12
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merqury_prep:
    input:
        'data/{animal}.{sample}.hifi.meryl'
    output:
        hist = 'data/{animal}.{sample}.hifi.hist',
        filt = 'data/{animal}.{sample}.hifi.filt'
    threads: 12
    resources:
        mem_mb = 2000
    shell:
        '''
        meryl histogram {input} > {output.hist}
        java -jar -Xmx1g {config[merqury_root]}/eval/kmerHistToPloidyDepth.jar {output.hist} > {output.hist}.ploidy
        filt=`sed -n 2p {output.hist}.ploidy | awk '{{print $NF}}'`
        echo $filt > {output.filt}
        meryl greater-than $filt output data/{wildcards.animal}.{wildcards.sample}.hifi.gt$filt.meryl {input}
        '''

rule merqury_spectra:
    input:
        read_db = 'data/{animal}.{sample}.hifi.meryl',
        asm_db = '{assembler}_{sample}/{animal}.asm.contigs.meryl',
        asm = '{assembler}_{sample}/{animal}.asm.contigs.fasta',
        filt = 'data/{animal}.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{animal}.asm.f','.qv','.completeness.stats')
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2',
        'igv/2.8.2'
    shell:
        '''
        cd {wildcards.assembler}_{wildcards.sample}
        export MERQURY={config[merqury_root]}
        if [ ! -L {wildcards.animal}.{wildcards.sample}.hifi.filt ];
        then
            ln -s ../{input.filt}
            find ../data/ -name "*gt*" -exec ln -s ../data/{{}} . \;
        fi
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} ../{input.asm} {wildcards.animal}.asm
        '''

rule merqury_spectra_trio:
    input:
        read_db = 'data/{animal}.{sample}.hifi.meryl',
        hap_dbs = expand('{{assembler}}_{{sample}}/{{animal}}.hap{N}.contigs.meryl',N=(1,2)),
        hap1 = '{assembler}_{sample}/{animal}.hap1.contigs.fasta',
        hap2 = '{assembler}_{sample}/{animal}.hap2.contigs.fasta',
        filt = 'data/{animal}.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{animal}.trio.f','.qv','.completeness.stats')
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2',
        'igv/2.8.2'
    shell:
        '''
        cd {wildcards.assembler}_{wildcards.sample}
        export MERQURY={config[merqury_root]}
        ln -sfn ../{input.filt}
        find ../data/ -name "*gt*" -exec ln -sfn ../data/{{}} . \;     
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} ../{input.hap1} ../{input.hap2} {wildcards.animal}.trio
        '''

rule merqury_formatting:
    input:
        qv = '{assembler}_{sample}/{animal}.{haplotype}.qv',
        completeness = '{assembler}_{sample}/{animal}.{haplotype}.completeness.stats'
    output:
        'results/{animal}_{haplotype}_{sample}_{assembler}.merqury.stats.txt'
    shell:
        '''
        awk '{{print "QV " $4}}' {input.qv} > {output}
        awk '{{print "CV " $5}}' {input.completeness} >> {output}
        '''
