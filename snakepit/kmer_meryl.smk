localrules: split_reads, merqury_submit, merqury_formatting

checkpoint split_reads:
    input:
        'data/{animal}.hifi.fq.gz'
    output:
        directory('split_{animal}/')
    params:
        config['split_size']
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        mkdir -p split_{wildcards.animal}
        zcat {input} | split -a 2 -d -C {params}GiB --filter='pigz -p 6 > $FILE.fq.gz' - split_{wildcards.animal}/chunk_
        '''

rule count_many:
    input:
        'split_{animal}/chunk_{sample}.fq.gz'
    output:
        directory('split_{animal}/chunk_{sample}.meryl')
    threads: 18
    resources:
        mem_mb = 3000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'

def aggregate_split_input(wildcards):
    checkpoint_output = checkpoints.split_reads.get(**wildcards).output[0]
    return expand('split_{animal}/chunk_{sample}.meryl',animal=wildcards.animal,sample=glob_wildcards(os.path.join(checkpoint_output, 'chunk_{sample}.fq.gz')).sample)

rule merge_many:
    input:
        aggregate_split_input
    output:
        directory('data/{animal}.hifi.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'echo {input}; echo "hello"; meryl union-sum k={params.K} memory={params.mem} threads={threads} output {output} {input} '

rule count_asm_kmers:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        directory('{assembler}/{animal}.contigs.meryl')
    threads: 12
    resources:
        mem_mb = 3000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'

rule merqury_prep:
    input:
        'data/{animal}.hifi.meryl'
    output:
        hist = 'data/{animal}.hifi.hist',
        filt = 'data/{animal}.hifi.filt'
    threads: 12
    resources:
        mem_mb = 2000
    shell:
        '''
        meryl histogram {input} > {output.hist}
        java -jar -Xmx1g {config[merqury_root]}/eval/kmerHistToPloidyDepth.jar {output.hist} > {output.hist}.ploidy
        filt=`sed -n 2p {output.hist}.ploidy | awk '{{print $NF}}'`
        echo $filt > {output.filt}
        meryl greater-than $filt output data/{wildcards.animal}.hifi.gt$filt.meryl {input}
        '''

rule merqury_spectra:
    input:
        read_db = 'data/{animal}.hifi.meryl',
        asm_db = '{assembler}/{animal}.contigs.meryl',
        asm = '{assembler}/{animal}.contigs.fasta',
        filt = 'data/{animal}.hifi.filt'
    output:
        multiext('{assembler}/{animal}','.qv','.completeness.stats')
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {wildcards.assembler}
        export MERQURY={config[merqury_root]}
        ln -s ../{input.filt}
        find ../data/ -name "*gt*" -exec ln -s ../data/{{}} . \;
        $MERQURY/eval/spectra-cn.sh ../{input.read_db} {wildcards.animal}.contigs.fasta {wildcards.animal}
        '''

rule merqury_submit:
    input:
        read_db = 'data/{animal}.hifi.meryl',
        asm_db = '{assembler}/{animal}.contigs.meryl',
        asm = '{assembler}/{animal}.contigs.fasta'
    output:
        multiext('{assembler}/{animal}','.dddqv','.completeness.qqqstats')
    params:
        mq_dir = config['merqury_root']
    shell:
        '''
        export MERQURY={params.mq_dir}
        $MERQURY/_submit_merqury.sh {input.read_db} {input.asm} {wildcards.assembler}/{wildcards.animal}
        touch {output}
        '''

rule merqury_formatting:
    input:
        qv = '{assembler}/{animal}.qv',
        completeness = '{assembler}/{animal}.completeness.stats'
    output:
        'results/{animal}_{assembler}.merqury.stats.txt'
    shell:
        '''
        awk '{{print "QV " $4}}' {input.qv} > {output}
        awk '{{print "CV " $5}}' {input.completeness} >> {output}
        '''
