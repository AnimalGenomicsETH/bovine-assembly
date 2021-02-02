
MODEL = '/cluster/work/pausch/alex/software/pepper/models/PEPPER_polish_haploid_guppy360.pkl'

rule all:
    input:
        'Assembly.contigs.fasta'

rule map_ONT_reads:
    input:
        reads = config['reads'],
        asm = '{haplotype}.fasta'
    output:
        temp('{haplotype}_ONT_reads.unsorted.bam')
    threads: 16
    resources:
        mem_mb = 4500,
        walltime = '24:00'
    shell:
        'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} | samtools view -hb -F 0x904 -o {output}'

rule sort_bam:
    input:
        '{haplotype}_ONT_reads.unsorted.bam'
    output:
        temp('{haplotype}_ONT_reads.sorted.bam')
    threads: 8
    resources:
        mem_mb = 6000,
        disk_scratch = 200
    shell:
        'samtools sort {input} -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule index_bam:
    input:
        '{haplotype}_ONT_reads.sorted.bam'
    output:
        temp('{haplotype}_ONT_reads.sorted.bam.bai')
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule pepper_make_images:
    input:
        asm = '{haplotype}.fasta',
        bam = '{haplotype}_ONT_reads.sorted.bam',
        bai = '{haplotype}_ONT_reads.sorted.bam.bai'
    output:
        temp(directory('pepper_images_{haplotype}'))
    threads: 12
    resources:
        mem_mb = 6000,
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
    threads: 12
    resources:
        mem_mb = 7500,
        walltime = '24:00'
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

rule merfin:
    input:
        fasta = '',
        seqmers = '',
        readmers = '',
        vcf = ''#should have refcall removed
    output:
        'merfin.vcf'
    params:
        coverage = 30, #double for haploid?
        low_mem = 10,
        high_mem = 60
    threads: 16
    resources:
        mem_mb = 5500
    shell:
        'merfin -vmer -sequence {input.fasta} -memory1 {params.low_mem} -memory2 {params.high_mem} -seqmers {input.seqmers} -readmers {input.readmers} -threads {threads} -peak {params.coverage} -vcf {input.vcf} -output {output}'

rule merfin_polish:
    input:
        fasta = '',
        vcf = 'merfin.vcf'
    output:
        vcf = multiext('','','.csi')
        fasta = 'polished.fa'
    shell:
        '''
        bcftools view -Oz {input.vcf} > {output.vcf[0]}
        bcftools index {output.vcf[0]}
        bcftools consensus {output.vcf[0]} -f {input.fasta} -H 2 > {output.fasta}
        '''

# rule racon_polish:
#     input:
#         scaffolds = WORK_PATH + '{haplotype}.scaffolds.fasta',
#         aln = WORK_PATH + '{haplotype}_scaffolds_reads.sam',
#         reads = 'data/offspring.{sample}.hifi.fq.gz'
#     output:
#         WORK_PATH + '{haplotype}.polished.fasta'
#     threads: 16
#     resources:
#         mem_mb = 28000,
#         walltime = '2:30'
#     shell:
#         'racon -t {threads} {input.reads} {input.aln} {input.scaffolds} > {output}'
