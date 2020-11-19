from pathlib import Path, PurePath
from glob import glob
from itertools import product

configfile: 'snakepit/run_parameters.yaml'

if config['animal'] != 'test':
    workdir: PurePath(config['workdir']).joinpath(config['animal'])

#GLOBAL VAR
raw_long_reads = f'{config["data"][config["animal"]]["long_reads"]["offspring"]}{{read_name}}.ccs.bam'

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: analysis_report, plot_dot, generate_reffai, validation_auN

for _dir in ['data','results']:
    Path(_dir).mkdir(exist_ok=True)

for assembler, sample in product(config['assemblers'],config['sampling']):
    Path(f'{assembler}_{sample}').mkdir(exist_ok=True)

FOLDER_PATH = '{assembler}_{sample}/'

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
    haplotype = r'\w+',
    hap = r'\w+',
    sample = r'\d+',
    data = r'[^\W_]+',
    modifier = r'\w+',
    read_t  = r'hifi|SR'

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        expand('{animal}_{sample}_analysis_report.pdf',animal=config['animal'],sample=config['sampling']),

if 'hifiasm' in config['assemblers']:
    rule assembler_hifiasm:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            'hifiasm_{sample}/asm.p_ctg.gfa'
        params:
            out = lambda wildcards, output: PurePath(output[0]).parent,
            settings = '-r 4 -a 5 -n 5'
        threads: 36
        resources:
            mem_mb = lambda wildcards, input, threads: int(input.size_mb*1.75/threads),
            walltime = '24:00'
        shell:
            'hifiasm -o {params.out} -t {threads} {params.settings} {input}'

 ##Requires gfatools installed
    rule assembler_hifiasm_conversion:
        input:
            'hifiasm_{sample}/{haplotype}.p_ctg.gfa'
        output:
            'hifiasm_{sample}/{haplotype}.contigs.fasta'
        resources:
            mem_mb = 6000,
            walltime = '0:45'
        shell: 'gfatools gfa2fa {input} > {output}'

    rule assembler_hifiasm_parental:
        input:
            'data/{parent}.cleaned.hifi.fq.gz'
        output:
            'hifiasm_100/{parent}.p_ctg.gfa'
        params:
            out = 'hifiasm_100/{parent}',
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
            'canu_{sample}/asm.contigs_all.fa'
        log: 'logs/assembler_canu/sample-{sample}.asm.out'
        params:
            temp = 'asm.complete'
        shell:
            '''
            canu -p asm -d canu_{wildcards.sample} genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' onSuccess="touch {params.temp}" onFailure="touch {params.temp}" > {log}
            while [ ! -e canu_{wildcards.sample}/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm canu_{wildcards.sample}/{params.temp}
            mv canu_{wildcards.sample}/asm.contigs.fasta {output}
            '''

    rule assembler_canu_parents:
        input:
            'data/{parent}.cleaned.hifi.fq.gz'
        output:
            'canu_100/{parent}.contigs_all.fa'
        log: 'logs/assembler_canu_parents/parent-{parent}.out'
        params:
            temp = '{parent}.complete'
        shell:
            '''
            canu -p {wildcards.parent} -d canu_100 genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' onSuccess="touch {params.temp}" onFailure="touch {params.temp}" > {log}
            while [ ! -e canu_100/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm canu_100/{params.temp}
            mv canu_100/{wildcards.parent}.contigs.fasta {output}
            '''

    rule strip_canu_bubbles:
        input:
            'canu_{sample}/{haplotype}.contigs_all.fa'
        output:
            'canu_{sample}/{haplotype}.contigs_raw.fa'
        shell:
            'seqtk seq -l0 {input} | grep "suggestBubble=no" -A 1 --no-group-separator > {output}'

if 'flye' in config['assemblers']:
    rule assembler_flye:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            'flye_{sample}/asm.contigs.fasta'
        threads: 36
        resources:
            mem_mb = 6500,
            walltime = '24:00'
        shell:
            '''
            flye --pacbio-hifi {input} -t {threads} --keep-haplotypes -o flye_{wildcards.sample}
            mv flye_{wildcards.sample}/assembly.fasta {output}
            '''
if 'IPA' in config['assemblers']:
    rule assembler_IPA:
        input:
            'data/offspring.{sample}.hifi.fq.gz'
        output:
            'IPA_{sample}/{haplotype}.contigs.fasta'
        shell:
            'snakemake --profile "lsf" -n --snakefile /cluster/work/pausch/alex/software/miniconda3/envs/pbipa/etc/ipa.snakefile --configfile config.yaml  -U finish -d RUN'

##Requires yak installed
rule validation_yak:
    input:
        reads = 'data/offspring.{sample}.hifi.fq.gz',
        contigs = '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        yak = temp('{assembler}_{sample}/{haplotype}.yak'),
        qv = 'results/{haplotype}_{sample}_{assembler}.asm-ccs.qv.txt'
    threads: 16
    resources:
        mem_mb = 6000
    shell:
        '''
        yak count -b 37 -t {threads} -o {output.yak} {input.reads}
        yak qv -t {threads} {output.yak} {input.contigs} > {output.qv}
        '''

rule validation_yak_trio:
    input:
        contigs = '{assembler}_{sample}/{haplotype}.contigs.fasta',
        parents = expand('data/{parent}.yak', parent = ('sire','dam'))
    output:
        'results/{haplotype}_{sample}_{assembler}.trioyak.txt'
    threads: 16
    resources:
        mem_mb = 6000
    shell:
        'yak trioeval -t {threads} {input.parents} {input.contigs} > {output}'

##Requires k8 and calN50.js installed
rule validation_auN:
    input:
        '{assembler}_{sample}/{haplotype}.{sequence}.fasta'
    output:
        'results/{haplotype}_{sample}_{assembler}.{sequence}.auN.txt'
    shell:
        'k8 ~/bin/calN50.js -s 0.01 {input} > {output}'

##Requires minigraph and paftools.js installed
rule validation_refalign:
    input:
        ref_fai  = f'{config["ref_genome"]}.fai',#f'{config['ref_genome']}[{{reference}}].fai'
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        paf = temp('{assembler}_{sample}/{haplotype}_asm.paf'),
        NGA = 'results/{haplotype}_{sample}_{assembler}.NGA50.txt'
    threads: 24
    resources:
        mem_mb = 4000,
        walltime = '2:00'
    shell:
        '''
        minigraph -k 21 -xasm  --show-unmap=yes -t {threads} {config[ref_genome]} {input.asm} > {output.paf}
        paftools.js asmstat {input.ref_fai} {output.paf} > {output.NGA}
        '''

rule validation_asmgene:
    input:
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta',
        reads = 'data/offspring.{sample}.hifi.fq.gz'
    output:
        asm_paf = temp('{assembler}_{sample}/{haplotype}_asm_aln.paf'),
        ref_paf = temp('{assembler}_{sample}/{haplotype}_ref_aln.paf'),
        asmgene = 'results/{haplotype}_{sample}_{assembler}.asmgene.txt'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        '''
        minimap2 -cxsplice:hq -t {threads} {input.asm} {config[cDNAs]} > {output.asm_paf}
        minimap2 -cxsplice:hq -t {threads} {input.reads} {config[cDNAs]} > {output.ref_paf}
        paftools.js asmgene -i.97 {output.ref_paf} {output.asm_paf} > {output.asmgene}
        '''

##Requires busco (and metaeuk) installed
rule validation_busco:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        out_dir = directory('{assembler}_{sample}/{haplotype}_BUSCO'),
        summary = 'results/{haplotype}_{sample}_{assembler}.BUSCO.txt'
    params:
        tmp_dir = '{haplotype}_{sample}_{assembler}_busco_results'
    threads: 24
    resources:
        mem_mb = 3500,
        walltime = '6:00'
    shell:
        '''
        busco --cpu {threads} -i {input} -o {params.tmp_dir}
        cp {params.tmp_dir}/short_summary*.txt {output.summary}
        mv {params.tmp_dir} {output.out_dir}
        '''

rule plot_dot:
    input:
        '{assembler}_{sample}/{haplotype}_scaffolds_{reference}.{mapper}.paf'
    output:
        'results/{haplotype}_{sample}_{assembler}.{reference}.{mapper}.dot.png'
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
        f'{config["animal"]}_{{sample}}_analysis_report.pdf'
    log:
        'logs/analysis_report/sample-{sample}.out'
    shell:
        'python {workflow.basedir}/scripts/denovo_assembly_statistics.py --animal {config[animal]} --sample {wildcards.sample} --haplotypes {config[haplotypes]} --assemblers {config[assemblers]} --css {workflow.basedir}/scripts/report.css --outfile {output} > {log}'
