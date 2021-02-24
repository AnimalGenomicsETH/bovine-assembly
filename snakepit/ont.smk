rule assembler_canu_ont:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        'asm'
    shell:
        '''
        canu -p asm -d canu_OBV genomeSize=2.7g -nanopore-raw OBV.pion.fq.gz executiveThreads=4 executiveMemory=8g -batMemory=100 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' -correctedErrorRate=0.144 -corFilter=quick -corPartitions=1500 -minReadLength=10000 -mhapMemory=20 -mhapThreads=12 -MhapBlockSize=1000 -obtOverlapper=mhap -utgOverlapper=mhap purgeOverlaps=aggressive
        canu -p asm -d canu_OBV genomeSize=2.7g -fast -trimmed -nanopore canu_OBV/asm.trimmedReads.fasta.gz stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"'  -mhapMemory=16 -mhapThreads=8 -MhapBlockSize=200 purgeOverlaps=aggressive
        '''

rule assembler_shasta:
    input:
        'data/offspring.{sample}.pion.fastq'
    output:
        dir_ = get_dir('work','shasta_{haplotype}',assembler='shasta')
        asm = get_dir('work','shasta_{haplotype}/{haplotype}.contigs.fasta',assembler='shasta')
    params:
        Path('config/shasta_config.conf').resolve() #get full path as required by shasta
    threads: 36
    resources:
        mem_mb = 20000,
        walltime = '4:00'
    shell:
        '''
        shasta --input {input} --threads {threads} --assemblyDirectory={output.dir_} --config {params}
        mv {output.dir_}/Assembly.fasta {output.asm}
        '''

rule assembler_raven:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        asm = get_dir('work','{haplotype}.contigs.fasta',assembler='raven')
        asm = get_dir('work','{haplotype}.gfa',assembler='raven')
        gfa = 'raven_100/asm.contigs.gfa'
    threads: 36
    resources:
        mem_mb = lambda wildcards, input, threads: int(input.size_mb*3.5/threads),
        walltime = '24:00'
    shell:
        'raven -t {threads} -p 0 --graphical-fragment-assembly={output.gfa} {input} > {output.asm}' #no polishing with racon

rule assembler_flye_ont_assemble:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        dir_ = get_dir('work','flye_{haplotype}',assembler='flye')
        asm = get_dir('work','flye_{haplotype}/00-assembly/draft_assembly.fasta',assembler='flye')
    threads: 34
    resources:
        mem_mb = 7500,
        walltime = '24:00'
    shell:
        'flye --nano-raw {input} --threads {threads} --stop-after assembly -o {output.dir_} --genome-size=2.7g'

rule assembler_flye_ont_consensus:
    input:
        get_dir('work','flye_{haplotype}/00-assembly/draft_assembly.fasta',assembler='flye')
    output:
        get_dir('work','flye_{haplotype}/30-contigger/contigs.fasta',assembler='flye')
    params:
        dir_ = lambda wildcards, input: PurePath(input[0]).parents[1]
    threads: 28
    resources:
        mem_mb = 15000,
        walltime = '24:00'
    shell:
        'flye --nano-raw {input} --threads {threads} --resume-from consensus --stop-after contigger -o {params.dir_} --genome-size=2.7g'


rule assembler_flye_ont_scaffold:
    input:
        get_dir('work','flye_{haplotype}/30-contigger/contigs.fasta',assembler='flye')
    output:
        get_dir('work','flye_{haplotype}/asm.flye_scaffold.fasta',assembler='flye')
    params:
        dir_ = lambda wildcards, input: PurePath(input[0]).parents[1]
    threads: 4
    resources:
        mem_mb = 4000,
        walltime = '30'
    shell:
        '''
        flye --nano-raw {input} --threads {threads} --resume-from finalize -o {params.dir_} --genome-size=2.7g
        mv {params.dir}/assembly.fasta {output}
        '''
