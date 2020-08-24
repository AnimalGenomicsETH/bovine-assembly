from glob import glob

configfile: 'merqury_config.yaml'

localrules: split_reads, kmer_completeness, QV_stats, populate_hist, asm_CN_only, make_plots
    
rule all:
    input: 
        'cow_hifiasm.qv.stats',
        'cow_hifiasm.completeness.stats',
        'cow_hifiasm.only.hist',
        'cow_hifiasm.spectra-cn.hist',
        'cow_hifiasm.spectra-asm.ln.png',
        'cow_hifiasm.spectra-cn.ln.png',
        'cow.hifi.meryl'

checkpoint split_reads:
    input:
        'data/{animal}.fq.gz'
    output:
        directory('split_{animal}')
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        mkdir -p split_{wildcards.animal}
        zcat {input} | split -a 2 -d -C 10GiB --filter='pigz -p 6 > $FILE.fq.gz' - split_{wildcards.animal}/chunk_
        '''

rule count_many:
    input:
        'split_{animal}/{sample}.fq.gz'
    output:
        directory('split_{animal}/{sample}.meryl')
    threads: 18
    resources:
        mem_mb = 3000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'

def aggregate_split_input(wildcards):
    checkpoint_output = checkpoints.split_reads.get(**wildcards).output[0]
    print('CPO',checkpoint_output)
    return expand('split_{animal}/chunk_{sample}.meryl',animal=wildcards.animal,sample=glob_wildcards(os.path.join(checkpoint_output, 'chunk_{sample}.fq.gz')).sample)
    
rule merge_many:
    input:
        aggregate_split_input
    output:
        directory('{animal}.hifi.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        K = config['k-mers']
    shell:
        '''
        meryl union-sum k={params.K} threads={threads} output {output} {input}
        '''

#rule count_read_kmers:
#    input:
#        'data/{animal}.fq.gz'
#    output:
#        directory('{animal}.hifi.meryl')
#    threads: 18
#    resources:
#        mem_mb = 3000
#    params:
#        K = config['k-mers'],
#        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
#    run:
#        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'
       
rule count_asm_kmers:
    input:
        '{assembler}/{animal}.contigs.fasta'
    output:
        directory('{animal}_{assembler}.meryl')
    threads: 12
    resources:
        mem_mb = 3000
    params:
        K = config['k-mers'],
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1100
    shell:
        'meryl count k={params.K} memory={params.mem} threads={threads} output {output} {input}'


rule difference_reads:
    input:
        asm = '{animal}_{assembler}.meryl',
        reads = '{animal}.hifi.meryl'
    output:
        reads = directory('read.k.{animal}_{assembler}.0.meryl'),
        asm = directory('{animal}_{assembler}.0.meryl')
    threads: 12
    resources:
        mem_mb = 2000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
    shell:
        '''
        meryl difference threads={threads} memory={params.mem} output {output.reads} {input.reads} {input.asm}
        meryl difference threads={threads} memory={params.mem} output {output.asm} {input.    asm} {input.reads}
        '''

#rule difference_CN:
#    input:
#        asm = '{animal}_{assembler}.meryl',
#        reads = '{animal}.hifi.meryl'
#    output:
#        directory('{animal}_{assembler}.0.meryl')
#    threads: 12
#    resources:
#        mem_mb = 2000
#    params:
#        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/1000
#    shell:
#        '''
#        meryl difference threads={threads} memory={params.mem} output {output} {input.asm} {input.reads}
#        '''

rule intersect_CN:
    input:
        asm = '{animal}_{assembler}.meryl',
        reads = '{animal}.hifi.meryl'
    output:
        directory('read.k.{animal}_{assembler}.{code,.*}.{i,.*}.meryl')
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
        directory('{animal}.hifi.solid.meryl')
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
        solid_read = '{animal}.hifi.solid.meryl'
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
        printf "Copies\tkmer_multiplicity\tCount\n" > {output}
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

rule asm_hist:
    input:
        asm_0 = 'read.k.{animal}_{assembler}.0.meryl',
        asm = 'read.k.{animal}_{assembler}...meryl'
    output:
        '{animal}_{assembler}.spectra-asm.hist'
    shell:
        '''
        printf "Assembly\tkmer_multiplicity\tCount\n" > {output}
        meryl histogram {input.asm_0} | awk '{{print "read-only\t"$0}}' >> {output}
        meryl histogram {input.asm} | awk -v hap="{wildcards.animal}_{wildcards.assembler}" '{{print hap"\t"$0}}' >> {output}
        '''

rule asm_dst:
    input:
        '{animal}_{assembler}.0.meryl'
    output:
        '{animal}_{assembler}.dist_only.hist'
    shell:
        '''
        set +e > /dev/null
        ASM1_ONLY=$(meryl statistics {input} | head -n3 | tail -n1 | awk '{{print $2}}')
		echo -e "{wildcards.animal}_{wildcards.assembler}\t0\t$ASM1_ONLY" > {output}
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
        read_solid = '{animal}.hifi.solid.meryl'
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

rule make_plots:
    input:
        spc_hist = '{animal}_{assembler}.spectra-cn.hist',
        asm_hist = '{animal}_{assembler}.spectra-asm.hist',
        histonly = '{animal}_{assembler}.only.hist',
        distonly = '{animal}_{assembler}.dist_only.hist' 
    output:
        '{animal}_{assembler}.spectra-asm.ln.png',
        '{animal}_{assembler}.spectra-cn.ln.png'
    params:
        path = config['merqury_path'],
        r_lib = config['r_lib']
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        export R_LIBS_USER={params.r_lib}
        Rscript {params.path}/plot/plot_spectra_cn.R -f {input.asm_hist} -o {wildcards.animal}_{wildcards.assembler}.spectra-asm -z {input.distonly}
        Rscript {params.path}/plot/plot_spectra_cn.R -f {input.spc_hist} -o {wildcards.animal}_{wildcards.assembler}.spectra-cn -z {input.histonly}
        '''
