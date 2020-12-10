ruleorder: convert_sam_to_paf > map_asm_ref

rule generate_winnow_meryl:
    input:
        asm = lambda wildcards: WORK_PATH + '{haplotype}.scaffolds.fasta' if wildcards.haplotype !='ref' else f'{config["ref_genome"]}'
    output:
        db = temp(directory(WORK_PATH + '{haplotype}_winnow_k{K,\d+}.meryl')),
        rep = temp(WORK_PATH + '{haplotype}_repetitive_k{K,\d+}.txt')
    threads: 6
    resources:
        mem_mb = 3000,
        walltime = '1:00'
    shell:
        '''
        meryl count k={wildcards.K} output {output.db} {input.asm}
        meryl print greater-than distinct={config[winnow_threshold]} {output.db} > {output.rep}
        '''

rule map_asm_ref:
    input:
        asm = lambda wildcards: WORK_PATH + '{haplotype}.{seq_type}.fasta' if wildcards.haplotype != 'ref' else f'{config["ref_genome"]}',
        rep = lambda wildcards: WORK_PATH + '{haplotype}_repetitive_k19.txt' if wildcards.mapper == 'wm2' else [],
        ref = lambda wildcards: WORK_PATH + '{reference}.{seq_type}.fasta' if wildcards.reference != 'ref' else f'{config["ref_genome"]}'
    output:
        WORK_PATH + '{haplotype}_{seq_type}_{reference}.{mapper}.{ext,sam|paf}'
    params:
       mapper = lambda wildcards, input: f'winnowmap -W {input.rep}' if wildcards.mapper == 'wm2' else 'minimap2',
       out_t = lambda wildcards: '-c' if wildcards.ext == 'paf' else '-a'
    threads: 24
    resources:
        mem_mb = lambda wildcards: 2000 if wildcards.mapper == 'mm2' else 3500
    shell:
        '{params.mapper} {params.out_t} -x asm5 --cs -t {threads} {input.ref} {input.asm} -o {output}'

rule convert_sam_to_paf:
    input:
        WORK_PATH + '{haplotype}_{seq_type}_{reference}.{mapper}.sam'
    output:
        WORK_PATH + '{haplotype}_{seq_type}_{reference}.{mapper}.paf'
    shell:
        'paftools.js sam2paf {input} > {output}'

rule map_hifi_cell:
    input:
        reads = lambda wildcards: f'data/{wildcards.haplotype if wildcards.haplotype in ("dam","sire") else "offspring"}_{{read_name}}.temp.fastq.gz',
        asm = WORK_PATH + '{haplotype}.scaffolds.fasta',
        rep = WORK_PATH + '{haplotype}_repetitive_k15.txt'
    output:
        WORK_PATH + '{haplotype}_scaffolds_hifi_reads_{read_name}.bam'
    threads: 24
    resources:
        mem_mb = 4000,
        walltime = '20:00'
    shell:
        'winnowmap -t {threads} -W {input.rep} -ax map-pb {input.asm} {input.reads} | samtools view -b -o {output} -'

rule merge_sort_map_cells:
    input:
        lambda wildcards: expand('{{assembler}}_{{sample}}/{{haplotype}}_scaffolds_hifi_reads_{read_name}.bam',read_name=glob_wildcards(f'data/{wildcards.haplotype if wildcards.haplotype in ("dam","sire") else "offspring"}_{{read_name}}.temp.fastq.gz').read_name)
    output:
        sam = WORK_PATH + '{haplotype}_scaffolds_hifi_reads.bam'
    threads: 16
    resources:
        mem_mb = 6000,
        disk_scratch = lambda wildcards, input: int(input.size_mb/750)
    shell:
        'samtools cat --threads {threads} {input} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule map_SR_reads:
    input:
        reads = expand('data/offspring.read_R{N}.SR.fq.gz', N = (1,2)),
        asm = WORK_PATH + '{haplotype}.scaffolds.fasta'
    output:
        WORK_PATH + '{haplotype}_scaffolds_SR_reads.bam'
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax sr -t {threads} {input.asm} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'
