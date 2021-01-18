#minimap2 -ax map-ont -t 12 Assembly.fasta ../OBV.pion.fq.gz | samtools view -hb -F 0x904 > long_usrt.bam
#samtools sort long_usrt.bam -@ 8 -m 3000M -T $TMPDIR -o sorted.bam; samtools index -@ 8 sorted.bam

MODEL = '/cluster/work/pausch/alex/software/pepper/models/PEPPER_polish_haploid_guppy360.pkl'

rule pepper_make_images:
    input:
        fasta = '{haplotype}.unpolished.fasta',
        bam = '{haplotype}_ONT_reads.sorted.bam'
    output:
        temp(directory('pepper_images_{haplotype}'))
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        'pepper make_images -b {input.bam} -f {input.asm} -o {output} -t {threads}'

rule pepper_call_consensus:
    input:
        images = 'pepper_images_{haplotype}'
    output:
        temp(directory('pepper_predictions_{haplotype}'))
    params:
        batch_size = 256,
        workers = 4,
        model = MODEL
    threads: 24
    resources:
        mem_mb = 3500
    shell:
        'pepper call_consensus -i {input.images} -bs {params.batch_size} -w {params.workers} -m {params.model} -t {threads} -o {output}'

rule pepper_stitch:
    input:
        'pepper_predictions_{haplotype}'
    output:
        '{haplotype}.contigs.fasta'
    threads: 6
    resources:
        mem_mb = 12000,
        walltime = '12:00'
    shell:
        'pepper stitch -i {input} -o {output} -t {threads}'

rule racon_polish:
    input:
        scaffolds = WORK_PATH + '{haplotype}.scaffolds.fasta',
        aln = WORK_PATH + '{haplotype}_scaffolds_reads.sam',
        reads = 'data/offspring.{sample}.hifi.fq.gz'
    output:
        WORK_PATH + '{haplotype}.polished.fasta'
    threads: 16
    resources:
        mem_mb = 28000,
        walltime = '2:30'
    shell:
        'racon -t {threads} {input.reads} {input.aln} {input.scaffolds} > {output}'
