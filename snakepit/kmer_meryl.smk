localrules: split_reads, merqury_formatting

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
        'echo {input}; echo "hello"; meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input} '

rule count_asm_kmers:
    input:
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        directory('{assembler}_{sample}/{animal}.contigs.meryl')
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
        meryl greater-than $filt output data/{wildcards.animal}.hifi.gt$filt.meryl {input}
        '''

rule merqury_spectra:
    input:
        read_db = 'data/{animal}.{sample}.hifi.meryl',
        asm_db = '{assembler}_{sample}/{animal}.contigs.meryl',
        asm = '{assembler}_{sample}/{animal}.contigs.fasta',
        filt = 'data/{animal}.{sample}.hifi.filt'
    output:
        multiext('{assembler}_{sample}/{animal}','.qv','.completeness.stats')
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

rule merqury_formatting:
    input:
        qv = '{assembler}_{sample}/{animal}.qv',
        completeness = '{assembler}_{sample}/{animal}.completeness.stats'
    output:
        'results/{animal}_{sample}_{assembler}.merqury.stats.txt'
    shell:
        '''
        awk '{{print "QV " $4}}' {input.qv} > {output}
        awk '{{print "CV " $5}}' {input.completeness} >> {output}
        '''
