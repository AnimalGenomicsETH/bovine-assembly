ANIMAL = "BSWCHEF120152514636"
RAW_DATA_PATH = f"/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_data/PacBio_CCS/{ANIMAL}"
PROCESSED_DATA_PATH = f"/cluster/work/pausch/alex/assembly/{ANIMAL}"
#REF_GENOME = "/cluster/work/pausch/alex/REF_DATA/ARS-UCD1-2.fa.gz"

FINAL_EXTS = ["asm-ccs.qv.txt", "NGA50.txt", "auN.txt"]

from pathlib import Path   

def pair_name_to_infiles():
    # recursively find all *.ccs.bam read files under this animal
    read_path = Path(RAW_DATA_PATH).glob('PacBioSequel2*/*.ccs.bam')

    # pair read name to infile path using a dictionary
    return {f.name[:-8]:str(f) for f in read_path}
 
reads_infile_dict = pair_name_to_infiles()

##DEFINE LOCAL RULES FOR MINIMAL EXECUTION
localrules: analysis_report

##DEFINE CONFIG SETUP
configfile: 'run_parameters.yaml'
workdir: config['workdir']

#------------#
#DEFINE RULES#
#------------#
rule all:
    input:
        "analysis_report.pdf"

rule raw_read_conversion:
    input:
        lambda wildcards: reads_infile_dict[wildcards.read_name]
    output:
        f"{PROCESSED_DATA_PATH}/{{read_name}}.fastq.gz"
    params:
        f"{PROCESSED_DATA_PATH}/{{read_name}}"
    shell:
        """
        pbindex {input}
        bam2fastq -o {params} {input}
        """

rule raw_merge_files:
    input:
        expand(f"{PROCESSED_DATA_PATH}/{{data_coll_run}}.fastq.gz",data_coll_run=reads_infile_dict.keys())
    output:
        f"{PROCESSED_DATA_PATH}/{ANIMAL}.fq.gz"
    shell: "cat {input} > {output}"

rule assembler_hifiasm:
    input:
        f"{PROCESSED_DATA_PATH}/{ANIMAL}.fq.gz"
    output:
        f"{ANIMAL}.asm.p_ctg.gfa"
    resources:
        mem_mb = 6000,
        walltime = 8
    threads: 24
    shell: f"hifiasm -o {ANIMAL}.asm -t {{threads}} {{input}}"

##Requires gfatools installed
rule assembler_conversion:
    input:
        f"{ANIMAL}.asm.p_ctg.gfa"
    output:
        f"{ANIMAL}.contigs.fasta"
    resources:
        mem_mb = 6000
    shell: "gfatools gfa2fa {input} > {output}"

##Requires yak installed
rule validation_yak:
    input: 
        reads = f"{PROCESSED_DATA_PATH}/{ANIMAL}.fq.gz",
        contigs = f"{ANIMAL}.contigs.fasta"
    output:
        f"{ANIMAL}.asm-ccs.qv.txt"
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        """
        yak count -b 37 -t {threads} -o ccs.yak {input.reads}
        yak qv -t {threads} ccs.yak {input.contigs} > {output}
        """
        #inspect completeness

##Requires k8 and calN50.js installed
rule validation_auN:
    input:
        f"{ANIMAL}.contigs.fasta"
    output:
        f"{ANIMAL}.auN.txt"
    shell: "k8 ~/bin/calN50.js -s 0.01 {input} > {output}"

##Requires minigraph and paftools.js installed
rule validation_refalign:
    input:
        ref = config["ref_genome"],
        ref_fai  = "{input.ref}.fai",
        asm = f"{ANIMAL}.contigs.fasta"
    output:
        f"{ANIMAL}.NGA50.txt"
    threads: 16
    resources:
        mem_mb = 3500,
        walltime = "2:00"
    shell:
        """
        minigraph -xasm -K1.9g --show-unmap=yes -t {threads} {input.ref} {input.asm} > asm.paf    
        paftools.js asmstat {input.ref_fai} asm.paf > {output}
        """

##Requires busco (and metaeuk) installed
rule validation_busco:
    input:
        f"{ANIMAL}.contigs.fasta"    
    output:
        out_dir = directory("busco_results"),
        summary = "busco_short_summary.txt"
    threads: 24
    resources:
        mem_mb = 3000
    shell: 
        """
        busco --cpu {threads} -i {input} -o {output.out_dir}
        cp {output.out_dir}/short_summary*.txt {output.summary}
        """
    
rule generate_reffai:
    input:
        config["ref_genome"]
    output:
        "{input}.fai"
    envmodules:
        "gcc/8.2.0",
        "samtools/1.6"   
    shell: "samtools faidx {input}"

rule analysis_report:
    input:
        expand(f"{ANIMAL}.{{ext}}",ext=FINAL_EXTS),
        "busco_short_summary.txt"
    output:
        "analysis_report.pdf"
    shell:
        f"python denovo_assembly_statistics.py --animal {ANIMAL} --assembler hifiasm --jobname ./logs/hifiasm.out --outfile {{output}}"
