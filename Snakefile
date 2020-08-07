ANIMAL = "BSWCHEF120152514636"
RAW_DATA_PATH = f"/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long-read_data/PacBio_CCS/{ANIMAL}"
PROCESSED_DATA_PATH="/cluster/work/pausch/alex/assembly/{ANIMAL}"

FINAL_EXTS = ["asm-ccs.qv.txt", "NGA50.txt", "asm.auN.txt"]



from pathlib import Path   

def pair_name_to_infiles():
    print(f"{{RAW_DATA_PATH}}")
    # recursively find all *.ccs.bam read files under this animal
    read_path = Path(RAW_DATA_PATH).glob('**/*.ccs.bam')

    # pair vcf name to infile path using a dictionary
    return {f.name[:-8]:str(f) for f in read_path}
 
#reads_infile_dict = pair_name_to_infiles()
#print(reads_infile_dict)

rule all:
    input:
        expand(f"{{PROCESSED_DATA_PATH}}/{{ANIMAL}}.{ext}",ext=FINAL_EXTS)

rule raw_read_conversion:
    input:
        lambda wildcards: reads_infile_dict[wildcards.read_name]
    output:
        f"{{PROCESSED_DATA_PATH}}.fq.gz"
    shell: "bam2fastx -o {output} {input}"

rule raw_merge_files:
    input:
        expand("{PROCESSED_DATA_PATH}/{data_coll_run}.fq.gz")
    output:
        "{PROCESSED_DATA_PATH}/{animal}.fq.gz"
    shell: "cat *.fq.gz > {animal}.fq.gz"

rule assembler_hifiasm:
    input:
        "{dataset}/{animal}.fq.gz"
    output:
        "{animal}.asm.p_ctg.gfa"
    shell:
        """
        hifiasm -o {animal}.asm -t {threads} {input}
        """

rule assembler_conversion:
    input:
        "{animal}.p_ctg.gfa"
    output:
        "{animal}.contigs.fasta"
    shell: "gfatools gfa2fa -s {input} > {output}"

rule validation_yak:
    input: 
        reads = "{animal}.long-reads.fq.gz",
        contigs = "{animal}.contigs.fasta"
    output:
        "{animal}_asm-ccs.qv.txt"
    shell:
        """
        yak count -b 37 -t {threads} -o ccs.yak {input.reads}
        yak qv -t {threads} ccs.yak {input.contigs} > {output}
        """
        #inspect completeness

rule validation_auN:
    input:
        "{animal}.contigs.fasta"
    output:
        "{animal}_asm.auN.txt"
    shell: "k8 calN50.js {input} > {output}"

rule validation_refalign:
    input:
        ref = "/cluster/work/pausch/alex/REF_DATA/ARC-USD1-2.fa",
        ref_fai  = "{input.ref}.fai",
        asm = "{animal}.contigs.fasta"
    output:
        "{animal}_NGA50.txt"
    shell:
        """
        minigraph -xasm -K1.9g --show-unmap=yes -t {threads} {input.ref} {input.asm} > {asm.paf}
        paftools.js asmstat {input.ref_fair} {asm.paf} > {output}
        """

rule generate_reffai:
    input:
        "/cluster/work/pausch/alex/REF_DATA/ARC-USD1-2.fa"
    output:
        "{input}.fai"
    shell: "samtools faidx {input}"

