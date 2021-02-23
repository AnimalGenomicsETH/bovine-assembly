localrules: trio_canu, prep_haplotype_canu, haplotype_canu

rule yak_count_SR:
    input:
        expand('data/{{individual}}.read_R{N}.SR.fq.gz', N = (1,2))
    output:
        'data/{individual}.yak'
    threads: 24
    resources:
        mem_mb = 3000
    params:
        kmer = 31
    shell:
        'yak count -k {params.kmer} -b 37 -t {threads} -o {output} <(zcat {input}) <(zcat {input})'

rule trio_hifiasm:
    input:
        reads = 'data/offspring.{sample}.hifi.fq.gz',
        mat = 'data/dam.yak',
        pat = 'data/sire.yak',
        asm = get_dir('work','asm.contigs.fasta',assembler='hifiasm')
    output:
        expand('hifiasm_{{sample}}/hap{N}.p_ctg.gfa', N = (1,2))
    params:
        out = get_dir('work','asm',assembler='hifiasm'),
        old = expand('hifiasm_{{sample}}/asm.hap{N}.p_ctg.gfa', N = (1,2)),
        new = expand('hifiasm_{{sample}}/hap{N}.p_ctg.gfa', N = (1,2)),
        settings = '-r 4 -a 5 -n 5'
    threads: 32
    resources:
        mem_mb = 4000,
        walltime = '2:00'
    shell:
        #NOTE /dev/null is used since the overlaps already exist from original hifiasm run
        '''
        hifiasm -o {params.out} -t {threads} {params.settings} -1 {input.pat} -2 {input.mat} /dev/null
        mv {params.old[0]} {params.new[0]} && mv {params.old[1]} {params.new[1]}
        '''

rule trio_canu:
    input:
        reads = 'data/offspring.{sample}.hifi.fq.gz',
        dam = expand('data/dam.read_R{N}.SR.fq.gz', N = (1,2)),
        sire =  expand('data/sire.read_R{N}.SR.fq.gz', N = (1,2))
    output:
         (get_dir('work',f'trio/asm-haplotype{N}.sh',assembler='canu') for N in (1,2))
    log:
        'logs/assembler_canu/sample-{sample}.partion.out'
    params:
        temp = 'asm.complete',
        full = 'canu_{sample}/trio/asm.complete',
        dir_ = lambda wildcards, output: PurePath(output[0]).parent
    shell:
        '''
        canu -haplotype -p asm -d {params.dir_} genomesize={config[genome_est]}g -haplotype1 {input.sire} -haplotype2 {input.dam} -pacbio-raw {input.reads} executiveThreads=4 executiveMemory=8g -batMemory=60 -minInputCoverage=0 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' > {log}
        while [ ! -e {output[1]} ]; do sleep 60; done
        echo "complete file found, ending sleep loop"
        '''

rule prep_haplotype_canu:
    input:
        get_dir('work','trio/asm-haplotype{N}.sh',assembler='canu')
    output:
        get_dir('work','trio/asm-haplotype{N}_edited.sh',assembler='canu')
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
        get_dir('work','trio/asm-haplotype{N}_edited.sh',assembler='canu')
    output:
        get_dir('work','hap{N}.contigs_all.fa',assembler='canu')
    params:
        temp = 'hap{N}.complete',
        dir_ = 'asm-haplotype{N}',
        dir_in = lambda wildcards, input: PurePath(input[0]).parent,
        #can wait for the haplotype.success file?
    log:
        'logs/assembler_canu_trio/sample-{sample}.haplotype-{N}.out'
    shell:
        '''
        cd {params.dir_in}
        bash ../../{input} > ../../{log}

        while [ ! -e {params.dir_}/{params.temp} ]; do sleep 60; done
        echo "complete file found, ending sleep loop"
        rm {params.dir_}/{params.temp}
        mv {params.dir_}/{params.dir_}.contigs.fasta ../../{output}
        '''
