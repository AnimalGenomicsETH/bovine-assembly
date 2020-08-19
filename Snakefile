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
        expand('{animal}_{assembler}_analysis_report.pdf',animal=config['animal'],assembler=config['assemblers'])

rule raw_read_conversion:
    input:
        lambda wildcards: reads_infile_dict[wildcards.read_name]
    output:
        'data/{read_name}.fastq.gz'
    threads: 8
    resources:
        mem_mb = 3000
    shell: 'samtools fastq -@ {threads} -c 6 -0 {output} {input}' 

rule raw_merge_files:
    input:
        expand('data/{data_coll_run}.fastq.gz',data_coll_run=reads_infile_dict.keys())
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
        threads: 24
        resources:
            mem_mb = 6000,
            walltime = "8:00"
        shell: 'hifiasm -o hifiasm/{wildcards.animal}.asm -t {threads} {input}'

if 'canu' in config['assemblers']:
    Path('canu').mkdir(exist_ok=True)
    localrules: assembler_canu
    rule assembler_canu:
        input:
            'data/{animal}.hifi.fq.gz'
        output:
            'canu/{animal}.contigs.fasta'
        shell:
            '''
            canu -p {wildcards.animal} -d canu genomeSize=2.7g -pacbio-hifi {input} executiveThreads=4 executiveMemory=8g onSuccess="touch {wildcards.animal}.complete"
            while  [  ! -f {wildcards.animal}.complete ]
              do
                sleep 60
              done
            echo "complete file found, ending sleep loop"
            rm {wildcards.animal}.complete
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
        'results/{animal}_{assembler}.asm-ccs.qv.txt'
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        '''
        yak count -b 37 -t {threads} -o intermediates/ccs.yak {input.reads}
        yak qv -t {threads} intermediates/ccs.yak {input.contigs} > {output}
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
        'results/{animal}_{assembler}.NGA50.txt'
    threads: 16
    resources:
        mem_mb = 3500,
        walltime = '2:00'
    shell:
        '''
        minigraph -xasm -K1.9g --show-unmap=yes -t {threads} {input.ref} {input.asm} > intermediates/asm.paf    
        paftools.js asmstat {input.ref_fai} intermediates/asm.paf > {output}
        '''

rule validation_asmgene:
    input:
        asm = '{assembler}/{animal}.contigs.fasta',
        reads = 'data/{animal}.hifi.fq.gz',
        cDNAs = config['cDNAs']
    output:
        'results/{animal}_{assembler}.asmgene.txt'
    threads: 24
    resources:
        mem_mb = 2500
    shell:
        '''
        minimap2 -cxsplice:hq -t {threads} {input.asm} {input.cDNAs} > intermediates/{animal}_{assembler}_asm_aln.paf
        minimap2 -cxsplice:hq -t {threads} {input.reads} {input.cDNAs} > intermediates/{animal}_{assembler}_ref_aln.paf
        paftools.js asmgene -i.97 intermediates/{animal}_{assembler}_ref_aln.paf intermediates/{animal}_{assembler}_asm_aln.paf > {output}
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
        mem_mb = 3000
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
        expand('results/{{animal}}_{{assembler}}.{ext}',ext=config['target_metrics']),
        'results/{animal}_{assembler}_busco_short_summary.txt'
    output:
        '{animal}_{assembler}_analysis_report.pdf'
    params:
        base_dir = workflow.basedir
    shell: 'python {params.base_dir}/denovo_assembly_statistics.py --animal {wildcards.animal} --assembler {wildcards.assembler} --outfile {output}'
