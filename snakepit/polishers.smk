rule polish_pepper:
    input:
        fasta = '.fasta',
        bam = '.bam'
    output:
        directory('out')
    params:
     'm'
    threads: 24
    resources:
        mem_mb = 4000
    shell:
        '''
        pepper polish --bam {input.bam} --fasta {input.fasta} --model_path --o {output} --threads {threads}
        '''

rule pepper_make_images:
    input:
        bam = '.bam',
        asm = '.fasta'
    output:
        directory('pepper_images')
    threads: 12
    resources:
        mem_mb = 3000,
        walltime = '24:00'
    shell:
        'pepper make_images -b {input.bam} -f {input.asm} -o {output} -t {threads}'

rule pepper_call_consensus:
    input:
        images = 'pepper_images'
    output:
        directory('predictions')
    params:
        batch_size = 256,
        workers = 4,
        model = '/cluster/work/pausch/alex/software/pepper/models/PEPPER_polish_haploid_guppy360.pkl'
    shell:
        'pepper call_consensus -i {input.images} -bs {params.batch_size} -w {params.workers} -m {params.model}-t {threads} -o {output}'

rule pepper_stitch:
    input:
        'predictions'
    output:
        'polished_genome.fasta'
    threads: 12
    resources:
        mem_mb = 3000
    shell:
        'pepper stitch -i {input} -o {output} -t {threads}'

rule polish_scaffolds:
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
