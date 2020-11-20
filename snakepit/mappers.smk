rule generate_winnow_meryl:
    input:
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        db = temp(directory('{assembler}_{sample}/{haplotype}_winnow_k{K,\d+}.meryl')),
        rep = temp('{assembler}_{sample}/{haplotype}_repetitive_k{K,\d+}.txt')
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
        asm = '{assembler}_{sample}/{haplotype}.{seq_type}.fasta',
        rep = lambda wildcards: '{assembler}_{sample}/{haplotype}_repetitive_k19.txt' if wildcards.mapper == 'wm2' else []
    output:
        '{assembler}_{sample}/{haplotype}_{seq_type}_{reference}.{mapper}.paf'
    params:
       ref = lambda wildcards: config['ref_genome'] if wildcards.reference == 'ref' else '{assembler}_{sample}/{reference}.scaffolds.fasta',
       mapper = lambda wildcards, input: f'winnowmap -W {input.rep}' if wildcards.mapper == 'wm2' else 'minimap2'
    threads: 24
    resources:
        mem_mb = 2000
    shell:
        '{params.mapper} -cx asm5 --cs -t {threads} {params.ref} {input.asm} -o {output}'

rule map_hifi_cell:
    input:
        reads = lambda wildcards: f'data/{wildcards.haplotype if wildcards.haplotype in ("dam","sire") else "offspring"}_{{read_name}}.temp.fastq.gz',
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        rep = '{assembler}_{sample}/{haplotype}_repetitive_k15.txt'
    output:
        '{assembler}_{sample}/{haplotype}_scaffolds_hifi_reads_{read_name}.bam'
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
        sam = '{assembler}_{sample}/{haplotype}_scaffolds_hifi_reads.bam'
    threads: 16
    resources:
        mem_mb = 6000,
        disk_scratch = lambda wildcards, input: int(input.size_mb/750)
    shell:
        'samtools cat --threads {threads} {input} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule map_SR_reads:
    input:
        reads = expand('data/offspring.read_R{N}.SR.fq.gz', N = (1,2)),
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        '{assembler}_{sample}/{haplotype}_scaffolds_SR_reads.bam'
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 200
    shell:
        'minimap2 -ax sr -t {threads} {input.asm} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'
