rule trio_yak:
    input:
        expand(f'{config["raw_data"]}/OB_trio_SR/NovaSeq_20200929_NOV475_o7428_DataDelivery/20200929.B-BOV_{{trio}}_R{{read}}.fastq.gz',read=(1,2),trio=lambda wildcards: config['trio'][wildcards.parent])
        #r1 = 'parental.fq.gz',
        #{config["trio"]["sire"]}.
        #r2 = f'{config["raw_data"]}/OB_trio_SR/NovaSeq_20200929_NOV475_o7428_DataDelivery/20200929.B-BOV_'
    output:
        'data/{parent}.yak'
    threads: 16
    resources:
        mem_mb = 2500
    params:
        K = 31 #config['k-mers']
    shell:
        'yak count -k {params.K} -b 37 -t {threads} -o {output} <(cat {input}) <(cat {input})'

rule trio_hifiasm:
    input:
        reads = 'data/{anima}_{sample}.hifi.fq.gz',
        mat = 'data/dam.yak',
        #mat = 'data/dam{config["trio"]["dam"]}.yak',
        pat = 'data/sire.yak',
        force = 'hifiasm_{sample}/{animal}.asm.contigs.fasta'
    output:
        expand('hifiasm_{sample}/{animal}.trio.hap{N}.p_ctg.gfa',N=(1,2))
    threads: 36
    resources:
        mem_mb = 4000,
        walltime = '24:00'
    shell:
        #NOTE /dev/null is used since the overlaps already exist from original hifiasm run
        'hifiasm -o hifiasm_{wildcards.sample}/{wildcards.animal}.trio -t {threads} -1 {input.pat} -2 {input.mat} /dev/null'
#rule trio_canu:
