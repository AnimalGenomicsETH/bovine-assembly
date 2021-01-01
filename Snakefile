from pathlib import Path, PurePath
from glob import glob
from itertools import product

configfile: 'snakepit/run_parameters.yaml'

if config['animal'] != 'test':
    workdir: PurePath(config['workdir']).joinpath(config['animal'])

#GLOBAL VAR
raw_long_reads = f'{config["data"][config["animal"]]["long_reads"]["offspring"]}{{read_name}}.ccs.bam'

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: all, analysis_report, plot_dot, generate_reffai, validation_auN, validation_refalign, validation_yak_completeness, validation_asmgene

for _dir in ['data','results']:
    Path(_dir).mkdir(exist_ok=True)

for assembler, sample in product(config['assemblers'],config['sampling']):
    Path(f'{assembler}_{sample}').mkdir(exist_ok=True)

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

WORK_PATH = '{assembler}_{sample}/'
RESULT_PATH = 'results/{haplotype}_{sample}_{assembler}'

include: 'snakepit/kmer_meryl.smk'
include: 'snakepit/cross_analysis.smk'
include: 'snakepit/data_preparation.smk'
include: 'snakepit/trio_assemblies.smk'
include: 'snakepit/purge_duplicates.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/mappers.smk'
include: 'snakepit/capture_logic.smk'

wildcard_constraints:
    assembler = r'[^\W_]+',
    parent = r'dam|sire',
    individual = r'dam|sire|offspring',
    haplotype = r'asm|hap1|hap2|sire|dam|ref',
    hap = r'\w+',
    sample = r'\d+',
    data = r'[^\W_]+',
    modifier = r'\w+',
    read_t  = r'hifi|SR',
    mapper = r'mm2|wm2'

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        f'{config["animal"]}_analysis_report.pdf'

if 'hifiasm' in config['assemblers']:
    rule assembler_hifiasm:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            WORK_PATH.format_map(Default({'assembler':'hifiasm'})) + 'asm.p_ctg.gfa'
        params:
            out = lambda wildcards, output: PurePath(output[0]).with_name('asm'),
            settings = '-r 4 -a 5 -n 5'
        threads: 36
        resources:
            mem_mb = lambda wildcards, input, threads: max(int(input.size_mb*1.75/threads),2000),
            walltime = '24:00'
        shell:
            'hifiasm -o {params.out} -t {threads} {params.settings} {input}'

    ##Requires gfatools installed
    rule assembler_hifiasm_conversion:
        input:
            WORK_PATH.format_map(Default({'assembler':'hifiasm'})) + '{haplotype}.p_ctg.gfa'
        output:
            WORK_PATH.format_map(Default({'assembler':'hifiasm'})) + '{haplotype}.contigs.fasta'
        resources:
            mem_mb = 10000,
            walltime = '60'
        shell: 'gfatools gfa2fa {input} > {output}'

    rule assembler_hifiasm_parental:
        input:
            'data/{parent}.cleaned.hifi.fq.gz'
        output:
            WORK_PATH.format_map(Default({'assembler':'hifiasm','sample':100})) + '{parent}.p_ctg.gfa'
        params:
            out = lambda wildcards, output: PurePath(output[0]).with_name(wildcards.parent),
            settings = '-r 4 -a 5 -n 5'
        threads: 36
        resources:
            mem_mb = lambda wildcards, input, threads: int(input.size_mb*1.75/threads),
            walltime = '24:00'
        shell:
            'hifiasm -o {params.out} -t {threads} {params.settings} {input}'

if 'canu' in config['assemblers']:
    localrules: assembler_canu, strip_canu_bubbles, assembler_canu_parents
    rule assembler_canu:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            WORK_PATH.format_map(Default({'assembler':'canu'})) + 'asm.contigs_all.fa'
        log: 'logs/assembler_canu/sample-{sample}.asm.out'
        params:
            dir_ = lambda wildcards, output: PurePath(output[0]).parent,
            temp = 'asm.complete'
        shell:
            '''
            canu -p asm -d {params.dir_} genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' -minInputCoverage=0 onSuccess="touch {params.temp}" onFailure="touch {params.temp}" > {log}
            while [ ! -e {params.dir_}/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm {params.dir_}/{params.temp}
            mv {params.dir_}/asm.contigs.fasta {output}
            '''

    rule assembler_canu_parents:
        input:
            'data/{parent}.cleaned.hifi.fq.gz'
        output:
            WORK_PATH.format_map(Default({'assembler':'canu','sample':100})) + '{parent}.contigs_all.fa'
        log: 'logs/assembler_canu_parents/parent-{parent}.out'
        params:
            dir_ = WORK_PATH.format_map(Default({'assembler':'canu','sample':100})) + '{parent}_asm',
            temp = '{parent}.complete'
        shell:
            '''
            canu -p {wildcards.parent} -d {params.dir_} genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' onSuccess="touch {params.temp}" onFailure="touch {params.temp}" > {log}
            while [ ! -e {params.dir_}/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm {params.dir_}/{params.temp}
            mv {params.dir_}/{wildcards.parent}.contigs.fasta {output}
            '''

    rule strip_canu_bubbles:
        input:
            WORK_PATH.format_map(Default({'assembler':'canu'})) + '{haplotype}.contigs_all.fa'
        output:
            WORK_PATH.format_map(Default({'assembler':'canu'})) + '{haplotype}.contigs_raw.fa'
        shell:
            'seqtk seq -l0 {input} | grep "suggestBubble=no" -A 1 --no-group-separator > {output}'

if 'flye' in config['assemblers']:
    rule assembler_flye:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            WORK_PATH.format_map(Default({'assembler':'flye'})) + 'asm.contigs.fasta'
        params:
            out = lambda output: PurePath(output[0]).parent
        threads: 36
        resources:
            mem_mb = 6500,
            walltime = '24:00'
        shell:
            '''
            flye --pacbio-hifi {input} -t {threads} --keep-haplotypes -o {params.out}
            mv {params.out}/assembly.fasta {output}
            '''
if 'IPA' in config['assemblers']:
    rule assembler_IPA:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            'IPA_{sample}/{haplotype}.contigs.fasta'
        shell:
            'snakemake --profile "lsf" -n --snakefile /cluster/work/pausch/alex/software/miniconda3/envs/pbipa/etc/ipa.snakefile --configfile config.yaml  -U finish -d {wildcards.haplotype}'

##Requires yak installed
rule count_yak_asm:
    input:
        contigs = WORK_PATH + '{haplotype}.contigs.fasta'
    output:
        yak = temp(WORK_PATH + '{haplotype}.yak')
    params:
        kmer = 31
    threads: 4
    resources:
        mem_mb = 6000,
        walltime = '0:30'
    shell:
        'yak count -k {params.kmer} -b 37 -t {threads} -o {output.yak} {input.contigs}'

rule validation_yak_qv:
    input:
        yak = 'data/offspring.yak',
        contigs = WORK_PATH + '{haplotype}.contigs.fasta'
    output:
        RESULT_PATH + '.yak.qv.txt'
    params:
        kmer = 31
    threads: 12
    resources:
        mem_mb = 5000
    shell:
        'yak qv -t {threads} -l100k -K2g {input.yak} {input.contigs} > {output}'

rule validation_yak_completeness:
    input:
        yak = 'data/offspring.yak',
        asm_yak = WORK_PATH + '/{haplotype}.yak'
    output:
        RESULT_PATH + '.yak.completeness.txt'
    shell:
        'yak inspect {input.yak} {output.asm_yak} > {output}'

rule validation_yak_trio:
    input:
        contigs = WORK_PATH + '{haplotype}.contigs.fasta',
        parents = expand('data/{parent}.yak', parent = ('sire','dam'))
    output:
        RESULT_PATH + '.yak.trio.txt'
    threads: 16
    resources:
        mem_mb = 6000
    shell:
        'yak trioeval -t {threads} {input.parents} {input.contigs} > {output}'

##Requires k8 and calN50.js installed
rule validation_auN:
    input:
        asm = WORK_PATH + '{haplotype}.{sequence}.fasta',
        ref_fai = f'{config["ref_genome"]}.fai'
    output:
        RESULT_PATH + '.{sequence}.auN.txt'
    shell:
        'calN50.js -s 0.01 -f {input.ref_fai} {input.asm} > {output}'

rule validation_refalign:
    input:
        fai = f'{config["ref_genome"]}.fai',
        paf = WORK_PATH + '{haplotype}_contigs_ref.mg.paf'
    output:
        RESULT_PATH + '.NGA50.txt'
    shell:
        'paftools.js asmstat {input.fai} {input.paf} > {output}'

rule map_splice_reads:
    input:
        lambda wildcards: (WORK_PATH + '{haplotype}.contigs.fasta') if wildcards.reference == 'asm' else f'{config["ref_genome"]}'
    output:
        temp(WORK_PATH + '{haplotype}_{reference}_splices.paf')
    threads: 8
    resources:
        mem_mb = 7000
    shell:
        'minimap2 -cxsplice:hq -t {threads} {input} {config[cDNAs]} > {output}'

rule validation_asmgene:
    input:
        asm = WORK_PATH + '{haplotype}_asm_splices.paf',
        ref = WORK_PATH + '{haplotype}_ref_splices.paf'
    output:
        RESULT_PATH + '.asmgene.{opt}.txt'
    params:
        thresh = lambda wildcards: '-i.97' if wildcards.opt == '97' else '-i.99',
        hap_opt = lambda wildcards: '' if wildcards.haplotype == 'asm' else '-a'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        'paftools.js asmgene {params.thresh} {params.hap_opt} {input.ref} {input.asm} > {output}'

##Requires busco (and metaeuk) installed
rule validation_busco:
    input:
        WORK_PATH + '{haplotype}.contigs.fasta'
    output:
        out_dir = directory(WORK_PATH + '{haplotype}_BUSCO'),
        summary = RESULT_PATH + '.BUSCO.txt'
    params:
        tmp_dir = '{haplotype}_{sample}_{assembler}_busco_results'
    threads: 12
    resources:
        mem_mb = 5000,
        walltime = '16:00'
    shell:
        '''
        busco --cpu {threads} -i {input} -o {params.tmp_dir}
        cp {params.tmp_dir}/short_summary*.txt {output.summary}
        mv {params.tmp_dir} {output.out_dir}
        '''

rule plot_dot:
    input:
        WORK_PATH + '{haplotype}_scaffolds_{reference}.{mapper}.paf'
    output:
        RESULT_PATH + '.{reference}.{mapper}.dot.png'
    shell:
        'minidot -L {input} | convert -density 150 - {output}'

rule generate_reffai:
    output:
        f'{config["ref_genome"]}.fai'
    shell:
        'samtools faidx {config[ref_genome]}'

rule analysis_report:
    input:
        capture_logic
    output:
        f'{config["animal"]}_analysis_report.pdf'
    log:
        'logs/analysis_report/unique.out'
    shell:
        'python {workflow.basedir}/scripts/denovo_assembly_statistics.py --animal {config[animal]} --samples {config[sampling]} --input {input} --css {workflow.basedir}/scripts/report.css --outfile {output} > {log}'
