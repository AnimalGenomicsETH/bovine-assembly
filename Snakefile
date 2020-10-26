from pathlib import Path, PurePath
from glob import glob
from itertools import product

configfile: 'snakepit/run_parameters.yaml'
workdir: PurePath(config['workdir']).joinpath(config['animal'])

def glob_purges(wildcards):
    req_files = []
    for _file in glob('**/*.contigs_raw.fa',recursive=True):
        req_files.extend([_file.replace('.fa','.spectra.png'),_file.replace('contigs_raw.fa','purged.spectra.png')])
    return req_files

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: analysis_report, raw_merge_files, plot_dot, paf_variants, generate_reffai, dnadiff

for _dir in ['data','results','intermediates']:
    Path(_dir).mkdir(exist_ok=True)

for assembler, sample in product(config['assemblers'],config['sampling']):
    Path(f'{assembler}_{sample}').mkdir(exist_ok=True)

include: 'snakepit/kmer_meryl.smk'
include: 'snakepit/cross_analysis.smk'
include: 'snakepit/data_preparation.smk'
include: 'snakepit/trio_assemblies.smk'
include: 'snakepit/purge_duplicates.smk'

wildcard_constraints:
    assembler = r'[^\W_]+',
    parent = r'[^\W_]+',
    haplotype = r'\w+',
    hap = r'\w+',
    sample = r'\d+'

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        #f'canu_100/{config["animal"]}.hap1.scaffolds.fasta',
        #f'canu_100/{config["animal"]}.hap2.scaffolds.fasta',
        #f'hifiasm_100/{config["animal"]}.hap1.scaffolds.fasta',
        #f'hifiasm_100/{config["animal"]}.asm.qv',
        #f'hifiasm_100/{config["animal"]}.hap2.contigs.fasta',
        #f'hifiasm_100/{config["animal"]}.trio.completeness.stats',
        #f'hifiasm_100/{config["animal"]}.hap2.{config["animal"]}.hap2.contigs.continuity.NG.png',
        #f'hifiasm_100/{config["animal"]}.asm.{config["animal"]}.asm.contigs.continuity.NG.png',
        #f'hifiasm_100/{config["animal"]}.asm.dnadiff.report',
        #'results/BSWCHEF1201525146361_100_hifiasm.gaps.txt',
        expand('{animal}_{sample}_analysis_report.pdf',animal=config['animal'],sample=config['sampling']),
        #f'hifiasm_100/{config["animal"]}.corrected.scaffolds.fasta'

if 'hifiasm' in config['assemblers']:
    rule assembler_hifiasm:
        input:
            'data/reads.{sample}.hifi.fq.gz'
        output:
            'hifiasm_{sample}/asm.p_ctg.gfa'
        threads: 36
        resources:
            mem_mb = 4000,
            walltime = '24:00'
        params:
            out = 'hifiasm_{sample}/asm',
            settings = '-r 4 -a 5 -n 5'
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

if 'canu' in config['assemblers']:
    localrules: assembler_canu, strip_canu_bubbles
    rule assembler_canu:
        input:
            'data/reads.{sample}.hifi.fq.gz'
        output:
            'canu_{sample}/asm.contigs_all.fa'
        log: 'logs/assembler_canu/sample-{sample}.asm.out'
        params:
            temp = 'asm.complete'
        shell:
            '''
            canu -p asm -d canu_{wildcards.sample} genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' onSuccess="touch {params.temp}" > {log}
            while [ ! -e canu_{wildcards.sample}/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm canu_{wildcards.sample}/{params.temp}
            mv canu_{wildcards.sample}/asm.contigs.fasta {output}
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
            'data/reads.{sample}.hifi.fq.gz'
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

##Requires yak installed
rule validation_yak:
    input:
        reads = 'data/reads.{sample}.hifi.fq.gz',
        contigs = '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        yak = temp('intermediates/{haplotype}_{sample}_{assembler}.yak'),
        qv = 'results/{haplotype}_{sample}_{assembler}.asm-ccs.qv.txt'
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        '''
        yak count -b 37 -t {threads} -o {output.yak} {input.reads}
        yak qv -t {threads} {output.yak} {input.contigs} > {output.qv}
        '''
        #inspect completeness

##Requires k8 and calN50.js installed
rule validation_auN:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        'results/{haplotype}_{sample}_{assembler}.auN.txt'
    shell:
        'k8 ~/bin/calN50.js -s 0.01 {input} > {output}'

##Requires minigraph and paftools.js installed
rule validation_refalign:
    input:
        ref_fai  = f'{config["ref_genome"]}.fai',#f'{config['ref_genome']}[{{reference}}].fai'
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        paf = temp('intermediates/{haplotype}_{sample}_{assembler}_asm.paf'),
        NGA = 'results/{haplotype}_{sample}_{assembler}.NGA50.txt'
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '2:00'
    shell:
        '''
        minigraph -k 21 -xasm  --show-unmap=yes -t {threads} {config[ref_genome]} {input.asm} > {output.paf}
        paftools.js asmstat {input.ref_fai} {output.paf} > {output.NGA}
        '''

rule validation_asmgene:
    input:
        asm = '{assembler}_{sample}/{haplotype}.contigs.fasta',
        reads = 'data/reads.{sample}.hifi.fq.gz'
    output:
        asm_paf = temp('intermediates/{haplotype}_{sample}_{assembler}_asm_aln.paf'),
        ref_paf = temp('intermediates/{haplotype}_{sample}_{assembler}_ref_aln.paf'),
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
        mem_mb = 4000,
        walltime = '16:00'
    shell:
        '''
        busco --cpu {threads} -i {input} -o {params.tmp_dir}
        cp {params.tmp_dir}/short_summary*.txt {output.summary}
        mv {params.tmp_dir} {output.out_dir}
        '''

rule generate_dot_paf:
    input:
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        '{assembler}_{sample}/{haplotype}_ref_scaffolds.paf'
    #params:
    #   ref = lambda wildcards: config['ref_genome'][wildcards.reference]
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        'minimap2 -cx asm5 --cs -t {threads} {config[ref_genome]} {input.asm} > {output}'

rule plot_dot:
    input:
        '{assembler}_{sample}/{haplotype}_ref_scaffolds.paf'
    output:
        'results/{haplotype}_{sample}_{assembler}.dot.png'
    shell:
        'minidot -L {input} | convert -density 150 - {output}'

rule paf_variants:
    input:
        '{assembler}_{sample}/{haplotype}_ref_scaffolds.paf'
    output:
        'results/{haplotype}_{sample}_{assembler}.vcf'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[ref_genome]} - > {output}'

rule nucmer:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        '{assembler}_{sample}/{haplotype}.delta'
    threads: 24
    resources:
        mem_mb = 5500
    shell:
        'nucmer --maxmatch -l 100 -c 500 -t {threads} -p {wildcards.assembler}_{wildcards.sample}/{wildcards.haplotype} {config[ref_genome]} {input}'

rule dnadiff:
    input:
        '{assembler}_{sample}/{haplotype}.delta'
    output:
        '{assembler}_{sample}/{haplotype}.dnadiff.report'
    shell:
        'dnadiff -d {input} -p {wildcards.assembler}_{wildcards.sample}/.{wildcards.haplotype}.dnadiff'

rule generate_reffai:
    output:
        f'{config["ref_genome"]}.fai'
    shell:
        'samtools faidx {config[ref_genome]}'

rule analysis_report:
    input:
        expand('results/{haplotype}_{{sample}}_{assembler}.{ext}',haplotype=config['haplotypes'],assembler=config['assemblers'],ext=config['target_metrics']),
        'data/reads.{sample}.QC.txt'#,
        #glob_purges#,
        #expand('{assembler}_{{sample}}/{{animal}}.scaffolds.fasta.masked',assembler=config['assemblers'])
    output:
        f'{config["animal"]}_{{sample}}_analysis_report.pdf'
    log:
        'logs/analysis_report/sample-{sample}.out'
    shell:
        'python {workflow.basedir}/scripts/denovo_assembly_statistics.py --animal {config[animal]} --sample {wildcards.sample} --haplotypes {config[haplotypes]} --assemblers {config[assemblers]} --css {workflow.basedir}/scripts/report.css --outfile {output} > {log}'

#onsuccess:
#    print('Cleaning up intermediate files')
#    shell('rm -rf intermediates')
