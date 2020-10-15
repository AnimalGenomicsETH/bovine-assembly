localrules: split_reads, merqury_formatting

rule count_SR_reads:
    input:
        'data/{parent_R{N}.fastq.gz'
    output:
        directory('data/{parent}.read_R{N}.meryl')
    threads: 18
    resources:
        mem_mb = 3000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merge_SR_reads:
    input:
        expand('data/{{parent}}.read_R{N}.meryl',N=(1,2))
    output:
        directory('data/{parent}.meryl')
    threads: 12ß
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
        directory('data/dam.hapmers.meryl'),
        directory('data/sire.hapmers.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hapmers.sh {input.sire} {input.dam} {input.child}
        '''

rule hap_blob:
    input:
        expand('data/{parent}.hapmers.meryl',parent=('sire','dam')),
        expand('{{assembler}}_{{sample}}/{{animal}}.hap{N}.contigs.fasta',N=(1,2))
    output:
        '{assembler}_{sample}/{animal}.blob.hapmers.count'
    params:
        out = '{assembler}_{sample}/{animal}.blob'
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hap_blob.sh {input} {params.out}
        '''

rule spectra_hap:
    input:
        'data/{animal}.{sample}.hifi.meryl',
        expand('data/{parent}.hapmers.meryl',parent=('sire','dam')),
        expand('{{assembler}}_{{sample}}/{{animal}}.hap{N}.contigs.fasta',N=(1,2))
    output:
        '$name.$asm.$read_hap.spectra-cn'
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/spectra-hap.sh {input} {params.out}
        '''

rule merqury_phase_block:
    input:
        asm = '{assembler}_{sample}/{animal}.hap{N}.contigs.fasta',
        hapmers = expand('data/{parent}.hapmers.meryl',parent=('sire','dam'))
    output:
        '{assembler}_{sample}/{animal}.hap{N}.phased_block.bed'
    params:
        '{assembler}_{sample}'
    shell:
        '''
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/phase_block.sh {input.asm} {input.hapmers} {params.out}
        '''

rule merqury_block_n_stats:
    input:
        hap1 = '{assembler}_{sample}/{animal}.hap1.contigs.fasta',
        hap2 = '{assembler}_{sample}/{animal}.hap2.contigs.fasta',
        block1 = '{assembler}_{sample}/{animal}.hap1.phased_block.bed',
        block2 = '{assembler}_{sample}/{animal}.hap2.phased_block.bed'
    output:
        expand('{assembler}_{sample}/{animal}.hap{N}.scaff.sizes',N=(1,2))
    params:
        out = '{assembler}_{sample}/{animal}'
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/block_n_stats.sh {input.hap1} {input.block1} {input.hap2} {input.block2} {params.out} {config[genome_est]}
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
        multiext('{assembler}_{sample}/{animal}.asm','.qv','.completeness.stats')
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {wildcards.assembler}_{wildcards.sample}
        export MERQURY={config[merqury_root]}
        if [ ! -L {wildcards.animal}.{wildcards.sample}.hifi.filt ];
        then
            ln -s ../{input.filt}
            find ../data/ -name "*gt*" -exec ln -s ../data/{{}} . \;
        fi
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} ../{input.asm} {wildcards.animal}_asm
        '''

rule merqury_spectra_trio:
    input:
        read_db = 'data/{animal}.{sample}.hifi.meryl',
        hap_dbs = '{assembler}_{sample}/{animal}.{haplotype}.contigs.meryl',
        hap1 = '{assembler}_{sample}/{animal}.hap1.contigs.fasta',
        hap2 = '{assembler}_{sample}/{animal}.hap2.contigs.fasta',
        filt = 'data/{animal}.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{animal}.{haplotype}','.qv','.completeness.stats')
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
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} ../{input.hap1} ../{input.hap2} {wildcards.animal}_trio
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
