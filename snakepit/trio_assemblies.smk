localrules: trio_canu, prep_haplotype_canu, haplotype_canu

rule trio_yak:
    input:
        expand('data/{{parent}}.read_R{N}.SR.fq.gz', N = (1,2))
    output:
        'data/{parent}.yak'
    threads: 24
    resources:
        mem_mb = 3000
    params:
        kmer = 31
    shell:
        'set +o pipefail; yak count -k {params.kmer} -b 37 -t {threads} -o {output} <(zcat {input}) <(zcat {input})'

rule trio_hifiasm:
    input:
        reads = 'data/offspring.{sample}.hifi.fq.gz',
        mat = 'data/dam.yak',
        pat = 'data/sire.yak',
        asm = 'hifiasm_{sample}/asm.contigs.fasta'
    output:
        expand('hifiasm_{{sample}}/hap{N}.p_ctg.gfa', N = (1,2))
    params:
        out = 'hifiasm_{sample}/asm',
        old = expand('hifiasm_{{sample}}/asm.hap{N}.p_ctg.gfa', N = (1,2)),
        new = expand('hifiasm_{{sample}}/hap{N}.p_ctg.gfa', N = (1,2))
    threads: 32
    resources:
        mem_mb = 4000,
        walltime = '2:00'
    shell:
        #NOTE /dev/null is used since the overlaps already exist from original hifiasm run
        '''
        hifiasm -o {params.out} -t {threads} -1 {input.pat} -2 {input.mat} /dev/null
        mv {params.old[0]} {params.new[0]} && mv {params.old[1]} {params.new[1]}
        '''

rule trio_canu:
    input:
        reads = 'data/offspring.{sample}.hifi.fq.gz',
        dam = expand('data/dam.read_R{N}.SR.fq.gz', N = (1,2)),
        sire =  expand('data/sire.read_R{N}.SR.fq.gz', N = (1,2))
    output:
        expand('canu_{{sample}}/trio/asm-haplotype{N}.sh', N=(1,2))
    log:
        'logs/assembler_canu/sample-{sample}.partion.out'
    params:
        temp = 'asm.complete',
        full = 'canu_{sample}/trio/asm.complete'
    shell:
        '''
        canu -haplotype -p asm -d canu_{wildcards.sample}/trio genomesize={config[genome_est]}g -haplotype1 {input.sire} -haplotype2 {input.dam} -pacbio-raw {input.reads} -batMemory=60 executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' > {log}
        while [ ! -e {output[1]} ]; do sleep 60; done
        echo "complete file found, ending sleep loop"
        '''

rule prep_haplotype_canu:
    input:
        'canu_{sample}/trio/asm-haplotype{N}.sh'
    output:
        'canu_{sample}/trio/asm-haplotype{N}_edited.sh'
    params:
        temp = 'hap{N}.complete'
    shell:
        '''
        #ONSEC=" 'onSuccess=\\\"touch {params.temp}\\\"'"
        ONSEC='  onSuccess="touch {params.temp}" \\\\'
        ONFAIL='  onFailure="touch {params.temp}" \\\\'
        sed -e '/raw/d' -e 's/-pacbio/-pacbio-hifi/' -e "10 a $ONSEC" -e "11 a $ONFAIL" {input} > {output}
        '''

rule haplotype_canu:
    input:
        'canu_{sample}/trio/asm-haplotype{N}_edited.sh'
    output:
        'canu_{sample}/hap{N}.contigs_all.fa'
    params:
        temp = 'hap{N}.complete',
        dir_ = 'asm-haplotype{N}'
    log:
        'logs/assembler_canu_trio/sample-{sample}.haplotype-{N}.out'
    shell:
        '''
        cd canu_{wildcards.sample}/trio
        bash ../../{input} > ../../{log}

        while [ ! -e {params.dir_}/{params.temp} ]; do sleep 60; done
        echo "complete file found, ending sleep loop"
        rm {params.dir_}/{params.temp}
        mv {params.dir_}/{params.dir_}.contigs.fasta ../../{output}
        '''
