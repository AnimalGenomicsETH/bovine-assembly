rule assembler_canu_ont:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        'asm'
    shell:
        '''
        canu -p asm -d canu_OBV genomeSize=2.7g -nanopore-raw OBV.pion.fq.gz executiveThreads=4 executiveMemory=8g -batMemory=100 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' -correctedErrorRate=0.144 -corFilter=quick -corPartitions=1500 -minReadLength=10000 -mhapMemory=20 -mhapThreads=12 -MhapBlockSize=1000 -obtOverlapper=mhap -utgOverlapper=mhap
        '''

rule assembler_shasta:
    input:
        'data/offspring.{sample}.pion.fastq'
    output:
        dir_ = directory('shasta_100'),
        asm = 'shasta_100/asm.contigs.fasta'
    threads: 36
    resources:
        mem_mb = 40000,
        walltime = '24:00'
    shell:
        '''
        shasta --input {input} --threads {threads} --assemblyDirectory={output.dir_}
        mv shasta_100/Assembly.fasta {output.asm}
        '''

rule assembler_raven:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        asm = 'raven_100/asm.contigs.fasta',
        gfa = 'raven_100/asm.contigs.gfa'
    threads: 36
    resources:
        mem_mb = lambda wildcards, input, threads: int(input.size_mb*3.5/threads),
        walltime = '24:00'
    shell:
        'raven -t {threads} --graphical-fragment-assembly={output.gfa} {input} > {output.asm}'
