
rule unimap_align:
    input:
        get_dir('work','{haplotype}.scaffolds.fasta')
    output:
        get_dir('work','{haplotype}.ARS.um.bam')
    threads: 16
    resources:
        mem_mb = 4000
    shell:
        'unimap --paf-no-hit -ax asm5 -r2k --cs -t {threads} {config[ref_genome]} {input} | samtools sort - -@ 4 -T $TMPDIR -o {output}'


rule svim_asm_diploid:
    input:
        (multiext(get_dir('work','{hap}.ARS.um.bam',hap=hap),'','.bai') for hap in ('hap1','hap2'))
    output:
        directory(get_dir('work','svim-asm'))
    threads: 1
    resources:
        mem_mb = 20000
    shell:
        '''
        mkdir -p {output}
        svim-asm diploid {output} {input[0]} {input[2]} {config[ref_genome]}
        '''
