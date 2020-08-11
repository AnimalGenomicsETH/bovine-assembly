ANIMAL = "BSWCHEF120152514636"
RAW_DATA_PATH = f"/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_data/PacBio_CCS/{ANIMAL}"
PROCESSED_DATA_PATH = f"/cluster/work/pausch/alex/assembly/{ANIMAL}"
REF_GENOME = "/cluster/work/pausch/alex/REF_DATA/ARS-UCD1-2.fa.gz"

FINAL_EXTS = ["asm-ccs.qv.txt", "NGA50.txt", "asm.auN.txt"]

from pathlib import Path   

def pair_name_to_infiles():
    # recursively find all *.ccs.bam read files under this animal
    read_path = Path(RAW_DATA_PATH).glob('PacBioSequel2*/*.ccs.bam')

    # pair read name to infile path using a dictionary
    return {f.name[:-8]:str(f) for f in read_path}
 
reads_infile_dict = pair_name_to_infiles()

rule all:
    input:
        expand(f"{PROCESSED_DATA_PATH}/{ANIMAL}.{{ext}}",ext=FINAL_EXTS)

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

rule assembler_conversion:
    input:
        f"{ANIMAL}.asm.p_ctg.gfa"
    output:
        f"{ANIMAL}.contigs.fasta"
    shell: "gfatools gfa2fa -s {input} > {output}"

rule validation_yak:
    input: 
        reads = f"{PROCESSED_DATA_PATH}/{ANIMAL}.fq.gz",
        contigs = f"{ANIMAL}.contigs.fasta"
    output:
        f"{PROCESSED_DATA_PATH}/{ANIMAL}.asm-ccs.qv.txt"
    shell:
        """
        yak count -b 37 -t {threads} -o ccs.yak {input.reads}
        yak qv -t {threads} ccs.yak {input.contigs} > {output}
        """
        #inspect completeness

rule validation_auN:
    input:
        f"{ANIMAL}.contigs.fasta"
    output:
        f"{PROCESSED_DATA_PATH}/{ANIMAL}.asm.auN.txt"
    shell: "k8 calN50.js {input} > {output}"

rule validation_refalign:
    input:
        ref = f"{REF_GENOME}",
        ref_fai  = f"{REF_GENOME}.fai",
        asm = f"{ANIMAL}.contigs.fasta"
    output:
        f"{PROCESSED_DATA_PATH}/{ANIMAL}.NGA50.txt"
    shell:
        """
        minigraph -xasm -K1.9g --show-unmap=yes -t {threads} {input.ref} {input.asm} > asm.paf
        paftools.js asmstat {input.ref_fai} asm.paf > {output}
        """

rule generate_reffai:
    input:
        f"{REF_GENOME}"
    output:
        "{input}.fai"
    shell: "samtools faidx {input}"

