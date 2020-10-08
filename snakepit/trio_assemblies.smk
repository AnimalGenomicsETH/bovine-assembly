rule trio_yak:
    input:
        #expand(f'{config["raw_data"]}/OB_trio_SR/NovaSeq_20200929_NOV475_o7428_DataDelivery/20200929.B-BOV_{{trio}}_R{{read}}.fastq.gz',read=(1,2),trio=
        lambda wildcards: (f'{config["raw_data"]}/OB_trio_SR/NovaSeq_20200929_NOV475_o7428_DataDelivery/20200929.B-BOV_{config["trio"][wildcards.parent]}_R{R}.fastq.gz' for R in (1,2))
        #r1 = 'parental.fq.gz',
        #{config["trio"]["sire"]}.
        #r2 = f'{config["raw_data"]}/OB_trio_SR/NovaSeq_20200929_NOV475_o7428_DataDelivery/20200929.B-BOV_'
    output:
        'data/{parent}.yak'
    threads: 24
    resources:
        mem_mb = 3000
    params:
        K = 31 #config['k-mers']
    shell:
        'set +o pipefail; yak count -k {params.K} -b 37 -t {threads} -o {output} <(zcat {input}) <(zcat {input})'

rule trio_hifiasm:
    input:
        reads = 'data/{animal}.{sample}.hifi.fq.gz',
        mat = 'data/dam.yak',
        #mat = 'data/dam{config["trio"]["dam"]}.yak',
        pat = 'data/sire.yak',
        force = 'hifiasm_{sample}/{animal}.asm.contigs.fasta'
    output:
        expand('hifiasm_{{sample}}/{{animal}}.hap{N}.p_ctg.gfa',N=(1,2))
    threads: 32
    resources:
        mem_mb = 4000,
        walltime = '2:00'
    shell:
        #NOTE /dev/null is used since the overlaps already exist from original hifiasm run
        '''
        hifiasm -o hifiasm_{wildcards.sample}/{wildcards.animal}.asm -t {threads} -1 {input.pat} -2 {input.mat} /dev/null
        mv hifiasm_{wildcards.sample}/{wildcards.animal}.asm.hap1.p_ctg.gfa hifiasm_{wildcards.sample}/{wildcards.animal}.hap1.p_ctg.gfa
        mv hifiasm_{wildcards.sample}/{wildcards.animal}.asm.hap2.p_ctg.gfa hifiasm_{wildcards.sample}/{wildcards.animal}.hap2.p_ctg.gfa
        '''

#rule trio_canu:
