configfile: 'merqury_config.yaml'

localrules: kmer_completeness, QV_stats, populate_hist, asm_CN_only
    
rule all:
    input: 
        'cow_hifiasm.qv.stats',
        'cow_hifiasm.completeness.stats',
        'cow_hifiasm.only.hist',
        'cow_hifiasm.spectra-cn.hist'
        
rule count_read_kmers:
    input:
        'data/{animal}.fq.gz'
    output:
        directory('{animal}.hifi.meryl')
    threads: 24
    resources:
        mem_mb = 5000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'
       
rule count_asm_kmers:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        directory('{animal}_{assembler}.meryl')
    threads: 24
    resources:
        mem_mb = 5000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'


rule difference_reads:
    input:
        asm = '{animal}_{assembler}.meryl',
        reads = '{animal}.hifi.meryl'
    output:
        directory('read.k.{animal}_{assembler}.0.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        '''
        meryl difference threads={threads} memory={params.mem} output {output} {input.reads} {input.asm}
        '''

rule difference_CN:
    input:
        asm = '{animal}_{assembler}.meryl',
        reads = '{animal}.hifi.meryl'
    output:
        directory('{animal}_{assembler}.0.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        '''
        meryl difference threads={threads} memory={params.mem} output {output} {input.asm} {input.reads}
        '''

rule intersect_CN:
    input:
        asm = '{animal}_{assembler}.meryl',
        reads = '{animal}.hifi.meryl'
    output:
        directory('read.k.{animal}_{assembler}.{code}.{i}.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        '''
        meryl intersect threads={threads} memory={params.mem} output {output} {input.reads} {wildcards.code} {wildcards.i} {input.asm}
        '''

rule solid_reads:
    input:
        '{animal}.hifi.meryl'
    output:
        directory('{animal}.solid.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        hist = '{animal}.hist',
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000,
        path = config['merqury_path']
    shell:
        '''
        meryl histogram {input} > {params.hist}
        java -jar -Xmx1g {params.path}/eval/kmerHistToPloidyDepth.jar {params.hist} > {params.hist}.ploidy
        filt=`sed -n 2p {params.hist}.ploidy | awk '{{print $NF}}'`
        meryl greater-than $filt threads={threads} memory={params.mem} output {output} {input}
        '''

rule solid_kmers:
    input:
        asm = '{animal}_{assembler}.meryl',
        solid_read = '{animal}.solid.meryl'
    output:
        directory('{animal}_{assembler}.solid.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        '''
        meryl intersect threads={threads} memory={params.mem} output {output} {input.asm} {input.solid_read}
        '''
        

rule populate_hist:
    input:
        read_only = 'read.k.{animal}_{assembler}.0.meryl',
        cn1 = 'read.k.{animal}_{assembler}.equal-to.1.meryl',
        cn2 = 'read.k.{animal}_{assembler}.equal-to.2.meryl',
        cn3 = 'read.k.{animal}_{assembler}.equal-to.3.meryl',
        cn4 = 'read.k.{animal}_{assembler}.equal-to.4.meryl',
        cngt4 = 'read.k.{animal}_{assembler}.greater-than.4.meryl'
    output:
        '{animal}_{assembler}.spectra-cn.hist'
    shell:
        '''
        meryl histogram {input.read_only} | awk '{{print "read-only\t"$0}}' >> {output}
        meryl histogram {input.cn1} | awk -v cn=1 '{{print cn"\t"$0}}' >> {output}
        meryl histogram {input.cn2} | awk -v cn=2 '{{print cn"\t"$0}}' >> {output}
        meryl histogram {input.cn3} | awk -v cn=3 '{{print cn"\t"$0}}' >> {output}
        meryl histogram {input.cn4} | awk -v cn=4 '{{print cn"\t"$0}}' >> {output}     
        meryl histogram {input.cngt4} | awk -v cn=">4" '{{print cn"\t"$0}}' >> {output}
        '''

rule asm_CN_only:
    input:
        '{animal}_{assembler}.0.meryl'
    output:
        '{animal}_{assembler}.only.hist'
    shell:
        '''
        PRESENT=`meryl statistics {input} | head -n4 | tail -n1 | awk '{{print $2}}'`
	    DISTINCT=`meryl statistics {input} | head -n3 | tail -n1 | awk '{{print $2}}'`
	    MULTI=$(($PRESENT-$DISTINCT))
	    printf "1\t0\t$DISTINCT\n2\t0\t$MULTI" > {output}
        '''

rule QV_stats:
    input:
        asm_0 = '{animal}_{assembler}.0.meryl',
        asm = '{animal}_{assembler}.meryl'
    output:
        '{animal}_{assembler}.qv.stats'
    params:
        K = config['k-mers']
    shell:
        #NOTE set +e required due to "The Infamous SIGPIPE signal"
        '''
        set +e > /dev/null
        ASM_ONLY=$(meryl statistics {input.asm_0}  | head -n 4 | tail -n 1 | awk '{{print $2}}')
        TOTAL=$(meryl statistics {input.asm}  | head -n 4 | tail -n 1 | awk '{{print $2}}')
        ERROR=$(echo "$ASM_ONLY $TOTAL" | awk -v k={params.K} '{{print (1-(1-$1/$2)^(1/k))}}')
        QV=$(echo "$ASM_ONLY $TOTAL" | awk -v k={params.K} '{{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}}')
        printf "{wildcards.animal}_{wildcards.assembler}\\t$ASM_ONLY\\t$TOTAL\\t$QV\\t$ERROR" >> {output}
        '''

rule kmer_completeness:
    input:
        asm_solid = '{animal}_{assembler}.solid.meryl',
        read_solid = '{animal}.solid.meryl'
    output:
        '{animal}_{assembler}.completeness.stats'
    shell:
        #NOTE set +e required due to "The Infamous SIGPIPE signal"
        '''
        set +e > /dev/null
        TOTAL=$(meryl statistics {input.read_solid}  | head -n 3 | tail -n 1 | awk '{{print $2}}')
        ASM=$(meryl statistics {input.asm_solid}  | head -n 3 | tail -n 1 | awk '{{print $2}}')
        printf "{wildcards.animal}_{wildcards.assembler}\\tall\\t$ASM\\t$TOTAL" | awk '{{print $0"\t"((100*$3)/$4)}}' >> {output}
        '''
