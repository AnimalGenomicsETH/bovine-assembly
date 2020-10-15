localrules: trio_canu, haplotype_canu

rule trio_yak:
    input:
        expand('data/{{parent}}_R{I}.fastq.gz',I=(1,2))
    output:
        'data/{parent}.yak'
    threads: 24
    resources:
        mem_mb = 3000
    params:
        K = 31
    shell:
        'set +o pipefail; yak count -k {params.K} -b 37 -t {threads} -o {output} <(zcat {input}) <(zcat {input})'

rule trio_hifiasm:
    input:
        reads = 'data/{animal}.{sample}.hifi.fq.gz',
        mat = 'data/dam.yak',
        pat = 'data/sire.yak',
        asm = 'hifiasm_{sample}/{animal}.asm.contigs.fasta'
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

rule trio_canu:
    input:
        reads = 'data/{animal}.{sample}.hifi.fq.gz',
        dam = expand('data/dam_R{N}.fastq.gz',N=(1,2)),
        sire =  expand('data/sire_R{N}.fastq.gz',N=(1,2))
    output:
        expand('canu_{{sample}}/trio/{{animal}}-haplotype{N}.sh',N=(1,2))
    params:

    shell:
        '''
        canu -haplotype -p {wildcards.animal} -d canu_{wildcards.sample}/trio genomesize={config[genome_est]} -haplotype1 {input.sire} -haplotype2 {input.dam} -pacbio-raw {input.reads} -batMemory=60 executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"'
        '''

rule haplotype_canu:
    input:
        'canu_{sample}/trio/{animal}-haplotype{N}.sh'
    output:
        'canu_{sample}/animal.hap{N}.contigs_raw.fa'
    params:
        temp = 'hap{N}.complete',
        dir_ = 'canu_{sample}/trio/{animal}-haplotype{N}'
    shell:
        '''
        sed -i -e {input}
        sed -i '4a\ onSuccess="touch {params.temp}"'
        ./{input}

        while [ ! -e {params.dir_}/{params.temp} ]; do sleep 60; done
        echo "complete file found, ending sleep loop"
        rm {params.dir_}/{params.temp}
        mv {params.dir_}/{wildcards.animal}-haplotype{wildcards.N}.contigs.fasta {output}
        '''
