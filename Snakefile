configfile: 'run_parameters.yaml'
workdir: config['workdir']

from pathlib import Path   

def pair_name_to_infiles():
    # recursively find all *.ccs.bam read files under this animal
    read_path = Path(f'{config["raw_data"]}/{config["animal"]}').glob('PacBioSequel2*/*.ccs.bam')

    # pair read name to infile path using a dictionary
    return {f.name[:-8]:str(f) for f in read_path}
 
reads_infile_dict = pair_name_to_infiles()

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: analysis_report, raw_merge_files

for _dir in ['data','results','intermediates']:
    Path(_dir).mkdir(exist_ok=True)

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        expand('{animal}_analysis_report.pdf',animal=config['animal'])

rule raw_read_conversion:
    input:
        f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam'
        #lambda wildcards: reads_infile_dict[wildcards.read_name]
    output:
        'data/{read_name}.fastq.gz'
    threads: 8
    resources:
        mem_mb = 3000
    shell: 'samtools fastq -@ {threads} -c 6 -0 {output} {input}' 

rule raw_merge_files:
    input:
        expand('data/{ID}.fastq.gz',ID=glob_wildcards(f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam').read_name)
        #expand('data/{data_coll_run}.fastq.gz',data_coll_run=reads_infile_dict.keys())
    output:
        'data/{animal}.hifi.fq.gz'
    shell: 'cat {input} > {output}'

if 'hifiasm' in config['assemblers']:
    Path('hifiasm').mkdir(exist_ok=True)
    rule assembler_hifiasm:
        input:
            'data/{animal}.hifi.fq.gz'
        output:
            'hifiasm/{animal}.asm.p_ctg.gfa'
        threads: 36
        resources:
            mem_mb = 3100,
            walltime = '16:00'
        shell: 'hifiasm -o hifiasm/{wildcards.animal}.asm -t {threads} {input}'

if 'canu' in config['assemblers']:
    Path('canu').mkdir(exist_ok=True)
    localrules: assembler_canu
    rule assembler_canu:
        input:
            'data/{animal}.hifi.fq.gz'
        output:
            'canu/{animal}.contigs_raw.fasta'
        log: 'logs/assembler_canu/animal-{animal}.out'
        params:
            '{animal}.complete'
        shell:
            '''
            canu -p {wildcards.animal} -d canu genomeSize=3g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g onSuccess="touch {params}" > {log}
            while [ ! -e canu/{params} ]; do sleep 60; done
            echo "complete file found, ending sleep loop"
            rm canu/{params}
            mv canu/{wildcards.animal}.contigs.fasta {output}
            '''
    include: 'purge_duplicates.smk'

if 'flye' in config['assemblers']:
    Path('flye').mkdir(exist_ok=True)
    rule assembler_flye:
        input:
            'data/{animal}.hifi.fq.gz'
        output:
            'flye/{animal}.contigs.fasta'
        threads: 36
        resources:
            mem_mb = 4500,
            walltime = '20:00'
        shell:
            '''
            flye --pacbio-hifi {input} -t {threads} --keep-haplotypes -o flye
            mv flye assembly.fasta {output}
            '''

##Requires gfatools installed
rule assembler_hifi_conversion:
    input:
        'hifiasm/{animal}.asm.p_ctg.gfa'
    output:
        'hifiasm/{animal}.contigs.fasta'
    resources:
        mem_mb = 6000
    shell: 'gfatools gfa2fa {input} > {output}'

##Requires yak installed
rule validation_yak:
    input: 
        reads = 'data/{animal}.hifi.fq.gz',
        contigs = '{assembler}/{animal}.contigs.fasta'
    output:
        yak = temp('intermediates/{animal}_{assembler}.yak'),
        qv = 'results/{animal}_{assembler}.asm-ccs.qv.txt'
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
        '{assembler}/{animal}.contigs.fasta'
    output:
        'results/{animal}_{assembler}.auN.txt'
    shell: 'k8 ~/bin/calN50.js -s 0.01 {input} > {output}'

##Requires minigraph and paftools.js installed
rule validation_refalign:
    input:
        ref = config['ref_genome'],
        ref_fai  = f'{config["ref_genome"]}.fai',
        asm = '{assembler}/{animal}.contigs.fasta'
    output:
        paf = temp('intermediates/{animal}_{assembler}_asm.paf'),
        NGA = 'results/{animal}_{assembler}.NGA50.txt'
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '2:00'
    shell:
        '''
        minigraph -k 21 -xasm  --show-unmap=yes -t {threads} {input.ref} {input.asm} > {output.paf}    
        paftools.js asmstat {input.ref_fai} {output.paf} > {output.NGA}
        '''

rule validation_asmgene:
    input:
        asm = '{assembler}/{animal}.contigs.fasta',
        reads = 'data/{animal}.hifi.fq.gz',
        cDNAs = config['cDNAs']
    output:
        asm_paf = temp('intermediates/{animal}_{assembler}_asm_aln.paf'),
        ref_paf = temp('intermediates/{animal}_{assembler}_ref_aln.paf'),
        asmgene = 'results/{animal}_{assembler}.asmgene.txt'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        '''
        minimap2 -cxsplice:hq -t {threads} {input.asm} {input.cDNAs} > {output.asm_paf}
        minimap2 -cxsplice:hq -t {threads} {input.reads} {input.cDNAs} > {output.ref_paf}
        paftools.js asmgene -i.97 {output.ref_paf} {output.asm_paf} > {output.asmgene}
        '''

##Requires busco (and metaeuk) installed
rule validation_busco:
    input:
        '{assembler}/{animal}.contigs.fasta'    
    output:
        out_dir = directory('{animal}_{assembler}_busco_results'),
        summary = 'results/{animal}_{assembler}_busco_short_summary.txt'
    threads: 24
    resources:
        mem_mb = 3000,
        walltime = '12:00'
    shell: 
        '''
        busco --cpu {threads} -i {input} -o {output.out_dir}
        cp {output.out_dir}/short_summary*.txt {output.summary}
        '''
    
rule generate_reffai:
    input:
        config['ref_genome']
    output:
        '{input}.fai'   
    shell: 'samtools fqidx {input}'

rule analysis_report:
    input:
        expand('results/{{animal}}_{assembler}.{ext}',assembler=config['assemblers'],ext=config['target_metrics']),
        expand('results/{{animal}}_{assembler}_busco_short_summary.txt',assembler=config['assemblers']),
        assemblers = config['assemblers']
    output:
        '{animal}_analysis_report.pdf'
    params:
        base_dir = workflow.basedir
    log:
        'logs/analysis_report/animal-{animal}.out'
    shell: 'python {params.base_dir}/denovo_assembly_statistics.py --animal {wildcards.animal} --assemblers {input.assemblers} --outfile {output} --css {params.base_dir}/github.css > {log}'

onsuccess:
    print('Cleaning up intermediate files')
    shell('rm -rf intermediates')
