from pathlib import Path, PurePath
from glob import glob
from itertools import product

configfile: 'config/run_parameters.yaml'

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

def get_dir(base,ext='',**kwargs):
    if base == 'work':
        base_dir = '{assembler}_{sample}'
    elif base == 'result':
        base_dir = 'results/{haplotype}_{sample}_{assembler}'
    elif base == 'data':
        base_dir = 'data'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext)

WORK_PATH = '{assembler}_{sample}/'
RESULT_PATH = 'results/{haplotype}_{sample}_{assembler}'

include: 'snakepit/merqury.smk'
include: 'snakepit/cross_analysis.smk'
include: 'snakepit/data_preparation.smk'
include: 'snakepit/trio_assemblies.smk'
include: 'snakepit/purge_duplicates.smk'
include: 'snakepit/variant_calling.smk'
include: 'snakepit/mappers.smk'
include: 'snakepit/capture_logic.smk'
include: 'snakepit/gap_closing.smk'

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
            get_dir('work','asm.p_ctg.gfa',assembler='hifiasm')
        params:
            out = lambda wildcards, output: PurePath(output[0]).with_name('asm'),
            settings = '-r 4 -a 5 -n 5'
        threads: 36
        resources:
            mem_mb = lambda wildcards, input, threads: max(int(input.size_mb*1.75/threads),1500),
            walltime = '24:00'
        shell:
            'hifiasm -o {params.out} -t {threads} {params.settings} {input}'

    ##Requires gfatools installed
    rule assembler_hifiasm_conversion:
        input:
            get_dir('work','{haplotype}.p_ctg.gfa',assembler='hifiasm')
        output:
            get_dir('work','{haplotype}.contigs.fasta',assembler='hifiasm')
        resources:
            mem_mb = 10000,
            walltime = '60'
        shell: 'gfatools gfa2fa {input} > {output}'

    rule assembler_hifiasm_parental:
        input:
            'data/{parent}.cleaned.hifi.fq.gz'
        output:
            get_dir('work','{parent}.p_ctg.gfa',assembler='hifiasm',sample=100)
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
            get_dir('work','asm.contigs_all.fa',assembler='canu')
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
            get_dir('work','{parent}.contigs_all.fa',assembler='canu',sample=100)
        log: 'logs/assembler_canu_parents/parent-{parent}.out'
        params:
            dir_ = lambda wildcards, output: PurePath(output[0]).with_name(f'{wildcards.parent}_asm'),
            temp = lambda wildcards: f'{wildcards.parent}.complete'
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
            get_dir('work','{haplotype}.contigs_all.fa',assembler='canu')
        output:
            get_dir('work','{haplotype}.contigs_raw.fa',assembler='canu')
        shell:
            'seqtk seq -l0 {input} | grep "suggestBubble=no" -A 1 --no-group-separator > {output}'

if 'flye' in config['assemblers']:
    rule assembler_flye:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            get_dir('work','asm.contigs.fasta',assembler='flye')
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
            get_dir('work','{haplotype}.contigs.fasta',assembler='IPA')
        shell:
            'snakemake --profile "lsf" -n --snakefile /cluster/work/pausch/alex/software/miniconda3/envs/pbipa/etc/ipa.snakefile --configfile config.yaml  -U finish -d {wildcards.haplotype}'

if 'peregrine' in config['assemblers']:
    rule assembler_peregrine:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            get_dir('work','{haplotype}.contigs.fasta',assembler='peregrine')
        params:
            dir_ = lambda wildcards, output: PurePath(output[0]).parent,
            procs = lambda wildcards, threads: ' '.join(str(threads)*9)
        shell:
            '''
            rsync -aq {input} $TMPDIR
            cd $TMPDIR
            singularity run -B ${TMPDIR}:/data /cluster/work/pausch/alex/peregrine_latest.sif asm <(echo {input}) {params.procs} --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output data/asm <<< "yes"
            cp -r ${TMPDIR}/asm {params.dir_}
            #mv final name to {output}
            '''

##Requires yak installed
rule count_yak_asm:
    input:
        contigs = get_dir('work','{haplotype}.contigs.fasta')
    output:
        yak = temp(get_dir('work','{haplotype}.yak'))
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
        contigs = get_dir('work','{haplotype}.contigs.fasta')
    output:
        get_dir('result','.yak.qv.txt')
    params:
        kmer = 31
    threads: 12
    resources:
        mem_mb = 5000
    shell:
        'yak qv -t {threads} -l100k -K2g {input.yak} {input.contigs} > {output}'

rule validation_yak_completeness:
    input:
        read_yak = 'data/offspring.yak',
        asm_yak = get_dir('work','{haplotype}.yak')
    output:
        get_dir('result','.yak.completeness.txt')
    shell:
        'yak inspect {input.read_yak} {output.asm_yak} > {output}'

rule validation_yak_trio:
    input:
        contigs = get_dir('work','{haplotype}.contigs.fasta'),
        parents = expand('data/{parent}.yak', parent = ('sire','dam'))
    output:
        get_dir('result','.yak.trio.txt')
    threads: 16
    resources:
        mem_mb = 6000
    shell:
        'yak trioeval -t {threads} {input.parents} {input.contigs} > {output}'

##Requires k8 and calN50.js installed
rule validation_auN:
    input:
        asm = get_dir('work','{haplotype}.{sequence}.fasta'),
        ref_fai = f'{config["ref_genome"]}.fai'
    output:
        get_dir('result','.{sequence}.auN.txt')
    shell:
        'calN50.js -s 0.01 -f {input.ref_fai} {input.asm} > {output}'

rule validation_refalign:
    input:
        fai = f'{config["ref_genome"]}.fai',
        paf = get_dir('work','{haplotype}_contigs_ref.mg.paf')
    output:
        get_dir('result','.NGA50.txt')
    shell:
        'paftools.js asmstat {input.fai} {input.paf} > {output}'

rule validation_asmgene:
    input:
        asm = get_dir('work','{haplotype}_cDNAs_splices.paf'),
        ref = str(PurePath(f'{config["ref_genome"]}').with_name('ref_cDNAs_splices.paf'))
    output:
        get_dir('result','.asmgene.{opt}.txt')
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
        get_dir('work','{haplotype}.contigs.fasta')
    output:
        out_dir = directory(get_dir('work','{haplotype}_BUSCO')),
        summary = get_dir('result','.BUSCO.txt')
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
        get_dir('work','{haplotype}_scaffolds_{reference}.{mapper}.paf')
    output:
        get_dir('result','.{reference}.{mapper}.dot.png')
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
