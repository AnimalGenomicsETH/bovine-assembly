rule assembler_canu_ont:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        'asm'
    shell:
        '''
        canu -p asm -d canu_fast genomeSize=2.7g -nanopore raw_data.pion.fq.gz executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' -correctedErrorRate=0.134 -corFilter=quick -corPartitions=1000 -minReadLength=10000 -mhapMemory=30 -mhapThreads=12 -MhapBlockSize=5000 -obtOverlapper=mhap
        '''

rule assembler_shasta:
    input:
        'data/offspring.{sample}.pion.fastq'
    output:
        'asm'
    threads: 36
    resources:
        mem_mb = 40000,
        walltime = '24:00'
    shell:
        'shasta --input {input} --threads {threads} --assemblyDirectory={params.out}'

rule assembler_raven:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        'asm'
    threads: 36
    resources:
        mem_mb = lambda wildcards, input, threads: int(input.size_mb*3.5/threads),
        walltime = '24:00'
    shell:
        'raven --graphical-fragment-assembly -t {threads} {input}'
