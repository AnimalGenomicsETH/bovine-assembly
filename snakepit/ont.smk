localrules: polish_shasta

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

rule sample_hap:
    input:
        config['ont_reads']
    output:
        'data/offspring.{sample}.pion.fasta'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    threads: 4
    resources:
        mem_mb = 3000
    shell:
        '''
        seqtk seq -A -f $(bc <<<"scale=2;{wildcards.sample}/100") {input} | pigz -p {threads} > {output}
        '''

rule assembler_shasta:
    input:
        'data/offspring.{sample}.pion.fasta'
    output:
        dir_ = get_dir('work','shasta_{haplotype}',assembler='shasta')
    params:
        Path('config/shasta_config.conf').resolve() #get full path as required by shasta
    threads: 42
    resources:
        mem_mb = 18000,
        walltime = '4:00'
    shell:
        '''
        shasta --input /cluster/scratch//alleonard/offspring.{wildcards.sample}.pion.fasta --threads {threads} --assemblyDirectory={output.dir_} --config {params}
        '''

rule polish_shasta:
    input:
        dir_ = get_dir('work','shasta_{haplotype}',assembler='shasta')
    output:
        get_dir('work','{haplotype}.contigs.fasta',assembler='shasta')
    params:
        data = lambda wildcards: Path(f'data/offspring.{wildcards.sample}.pion.fasta').resolve(),
        dir_ = lambda wildcards, output: PurePath(output[0]).parent
    shell:
        '''
        cd {input}
        snakemake -s /cluster/work/pausch/alex/assembly/bovine-assembly/snakepit/polishers.smk --config reads={params.data} --profile "lsf_nt" --quiet
        snakemake -s /cluster/work/pausch/alex/assembly/bovine-assembly/snakepit/deepvariant.smk --configfile /cluster/work/pausch/alex/assembly/BSWCHEM120151536851/config/shasta_DV.yaml --profile "lsf_nt" --quiet
        snakemake -s /cluster/work/pausch/alex/assembly/bovine-assembly/snakepit/merfin.smk --configfile /cluster/work/pausch/alex/assembly/BSWCHEM120151536851/config/shasta_merfin.yaml --profile "lsf_nt" --quiet
        cp merfin_NxB/polishing_asm_hifi/asm.merfin.fasta ../../{output}
        '''

rule assembler_raven:
    input:
        'data/offspring.{sample}.pion.fq.gz'
    output:
        asm = get_dir('work','{haplotype}.contigs.fasta',assembler='raven'),
        gfa = get_dir('work','{haplotype}.gfa',assembler='raven'),
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
        dir_ = get_dir('work','flye_{haplotype}',assembler='flye'),
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
