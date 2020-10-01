configfile: 'snakepit/run_parameters.yaml'
workdir: config['workdir']

from pathlib import Path
from glob import glob
from itertools import product

def pair_name_to_infiles():
    # recursively find all *.ccs.bam read files under this animal
    read_path = Path(f'{config["raw_data"]}/{config["animal"]}').glob('PacBioSequel2*/*.ccs.bam')

    # pair read name to infile path using a dictionary
    return {f.name[:-8]:str(f) for f in read_path}

reads_infile_dict = pair_name_to_infiles()

def glob_purges(wildcards):
    req_files = []
    for _file in glob('**/*.contigs_raw.fa',recursive=True):
        req_files.extend([_file.replace('.fa','.spectra.png'),_file.replace('contigs_raw.fa','purged.spectra.png')])
    return req_files

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: analysis_report, raw_merge_files, plot_dot, paf_variants, generate_reffai

for _dir in ['data','results','intermediates']:
    Path(_dir).mkdir(exist_ok=True)

for assembler, sample in product(config['assemblers'],config['sampling']):
    Path(f'{assembler}_{sample}').mkdir(exist_ok=True)

include: 'snakepit/kmer_meryl.smk'
include: 'snakepit/cross_analysis.smk'

wildcard_constraints:
    animal = r'\w+',
    assembler = r'\w+',
    sample = r'\d+'

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        #f'hifiasm_100/{config["animal"]}.scaffolds.fasta.masked'
        expand('{animal}_{sample}_analysis_report.pdf',animal=config['animal'],sample=config['sampling'])
        #multiext(f'results/{config["animal"]}_hifiasm','.telo.txt','.gaps.txt')

rule raw_read_conversion:
    input:
        f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam'
    output:
        temp('data/{read_name}.fastq.gz')
    threads: 8
    resources:
        mem_mb = 3000
    shell: 'samtools fastq -@ {threads} -c 6 -0 {output} {input}'

rule raw_merge_files:
    input:
        expand('data/{read_name}.fastq.gz',read_name=glob_wildcards(f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam').read_name)
    output:
        protected('data/{animal}.raw.hifi.fq.gz')
    shell: 'cat {input} > {output}'

if 'hifiasm' in config['assemblers']:
    rule assembler_hifiasm:
        input:
            'data/{animal}.{sample}.hifi.fq.gz'
        output:
            'hifiasm_{sample}/{animal}.asm.p_ctg.gfa'
        threads: 36
        resources:
            mem_mb = 3100,
            walltime = '16:00'
        shell: 'hifiasm -o hifiasm_{wildcards.sample}/{wildcards.animal}.asm -t {threads} {input}'

 ##Requires gfatools installed
rule assembler_hifiasm_conversion:
    input:
        'hifiasm_{sample}/{animal}.asm.p_ctg.gfa'
    output:
        'hifiasm_{sample}/{animal}.contigs.fasta'
    resources:
        mem_mb = 6000
    shell: 'gfatools gfa2fa {input} > {output}'

if 'canu' in config['assemblers']:
    localrules: assembler_canu
    rule assembler_canu:
        input:
            'data/{animal}.{sample}.hifi.fq.gz'
        output:
            'canu_{sample}/{animal}.contigs_raw.fa'
        log: 'logs/assembler_canu/sample-{sample}.animal-{animal}.out'
        params:
            temp = '{animal}.complete'
        shell:
            '''
            canu -p {wildcards.animal} -d canu_{wildcards.sample} genomeSize={config[genome_est]}g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g -batMemory=50 stageDirectory=\$TMPDIR gridEngineStageOption='-R "rusage[scratch=DISK_SPACE]"' onSuccess="touch {params.temp}" > {log}
            while [ ! -e canu_{wildcards.sample}/{params.temp} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm canu_{wildcards.sample}/{params.temp}
            mv canu_{wildcards.sample}/{wildcards.animal}.contigs.fasta {output}
            '''
    include: 'snakepit/purge_duplicates.smk'

if 'flye' in config['assemblers']:
    rule assembler_flye:
        input:
            'data/{animal}.{sample}.hifi.fq.gz'
        output:
            'flye_{sample}/{animal}.contigs.fasta'
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
        reads = 'data/{animal}.{sample}.hifi.fq.gz',
        contigs = '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        yak = temp('intermediates/{animal}_{sample}_{assembler}.yak'),
        qv = 'results/{animal}_{sample}_{assembler}.asm-ccs.qv.txt'
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
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{sample}_{assembler}.auN.txt'
    shell: 'k8 ~/bin/calN50.js -s 0.01 {input} > {output}'

##Requires minigraph and paftools.js installed
rule validation_refalign:
    input:
        ref_fai  = f'{config["ref_genome"]}.fai',#f'{config['ref_genome']}[{{reference}}].fai'
        asm = '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        paf = temp('intermediates/{animal}_{sample}_{assembler}_asm.paf'),
        NGA = 'results/{animal}_{sample}_{assembler}.NGA50.txt'
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
        asm = '{assembler}_{sample}/{animal}.contigs.fasta',
        reads = 'data/{animal}.{sample}.hifi.fq.gz'
    output:
        asm_paf = temp('intermediates/{animal}_{sample}_{assembler}_asm_aln.paf'),
        ref_paf = temp('intermediates/{animal}_{sample}_{assembler}_ref_aln.paf'),
        asmgene = 'results/{animal}_{sample}_{assembler}.asmgene.txt'
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
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        out_dir = directory('{assembler}_{sample}/{animal}_BUSCO'),
        summary = 'results/{animal}_{sample}_{assembler}.BUSCO.txt'
    params:
        tmp_dir = '{animal}_{sample}_{assembler}_busco_results'
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '16:00'
    shell:
        '''
        busco --cpu {threads} -i {input} -o {params.tmp_dir}
        cp {params.tmp_dir}/short_summary*.txt {output.summary}
        mv {params.tmp_dir} {output.out_dir}
        '''

rule generate_dot_paf:
    input:
        asm = '{assembler}_{sample}/{animal}.scaffolds.fasta'
    output:
        '{assembler}_{sample}/{animal}_ref_scaffolds.paf'
    #params:
    #   ref = lambda wildcards: config['ref_genome'][wildcards.reference]
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        'minimap2 -cx asm5 --cs -t {threads} {config[ref_genome]} {input.asm} > {output}'

rule plot_dot:
    input:
        '{assembler}_{sample}/{animal}_ref_scaffolds.paf'
    output:
        'results/{animal}_{sample}_{assembler}.dot.png'
    shell:
        'minidot -L {input} | convert -density 150 - {output}'

rule paf_variants:
    input:
        '{assembler}_{sample}/{animal}_ref_scaffolds.paf'
    output:
        'results/{animal}_{sample}_{assembler}.vcf'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[ref_genome]} - > {output}'

rule nucmer:
    input:
        '{assembler}_{sample}/{animal}.contigs.fasta'
    output:
        '{assembler}_{sample}/{animal}.delta'
    threads: 4
    resources:
        mem_mb = 15000
    shell:
        'nucmer --maxmatch -l 100 -c 500 -t {threads} -p {wildcards.assembler}_{wildcards.sample}/{wildcards.animal} {config[ref_genome]} {input}'

rule generate_reffai:
    output:
        f'{config["ref_genome"]}.fai'
    shell: 'samtools faidx {config[ref_genome]}'

rule analysis_report:
    input:
        expand('results/{{animal}}_{{sample}}_{assembler}.{ext}',assembler=config['assemblers'],ext=config['target_metrics']),
        'data/{animal}.{sample}.QC.txt',
        glob_purges#,
        #expand('{assembler}_{{sample}}/{{animal}}.scaffolds.fasta.masked',assembler=config['assemblers'])
    output:
        '{animal}_{sample}_analysis_report.pdf'
    log:
        'logs/analysis_report/animal-{animal}_sample-{sample}.out'
    shell: 
        'python {workflow.basedir}/scripts/denovo_assembly_statistics.py --animal {wildcards.animal} --sample {wildcards.sample} --assemblers {config[assemblers]} --css {workflow.basedir}/scripts/report.css --outfile {output} > {log}'

#onsuccess:
#    print('Cleaning up intermediate files')
#    shell('rm -rf intermediates')
