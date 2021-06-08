#ragtag.py correct, full reads mapping, c_mapping
rule ragtag_correct:
    input:
        asm = get_dir('work','{haplotype}.contigs.fasta'),
        reads = 'data/offspring.cleaned.hifi.fq.gz'
    output:
        get_dir('work','{haplotype}.contigs.corrected.fasta')
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        '''
        ragtag.py correct {config[ref_genome]} {input.asm} -o {wildcards.assembler}_{wildcards.sample} -R {input.reads} -T corr -t {threads} --mm2-params "-c -x asm20"
        #mv {wildcards.assembler}_{wildcards.sample}/ragtag.contigs.corrected {output}
        '''
#ragtag.py patch hifiasmv15_100/asm.scaffolds.fasta shastamerfin_100/asm.scaffolds.fasta --aligner minimap2 --mm2-params "-cx asm5 -t 8" --fill-only -o asm_patch3

rule ragtag_scaffold:
    input:
        get_dir('work','{haplotype}.contigs.fasta')
    output:
        get_dir('work','{haplotype}.scaffolds.fasta')
    params:
        lambda wildcards,output: PurePath(output[0]).parent / f'{wildcards.haplotype}_scaffolding'
    threads: 12
    resources:
        mem_mb = 3000
    shell:
        '''
        ragtag.py scaffold {config[ref_genome]} {input} -o {params} -u --mm2-params "-cx asm5 -t {threads}" -r -m 1000000
        sed 's/_RagTag//g' {params}/ragtag.scaffold.fasta > {output}
        '''

rule ragtag_patch:
    input:
        asm1 = lambda wildcards: get_dir('work','{haplotype}.{tigs}.fasta',assembler=wildcards.asm1),
        asm2 = lambda wildcards: get_dir('work','{haplotype}.{tigs}.fasta',assembler=wildcards.asm2)
    output:
        '{asm1}_{asm2}_{sample}_patched_{tigs}_{haplotype}/{haplotype}.scaffold.fasta'
    params:
        lambda wildcards,output: PurePath(output[0]).parent
    threads: 12
    resources:
        mem_mb = 2500,
        walltime = '2:00'
    shell:
        '''
        ragtag.py patch {input.asm1} {input.asm2} --aligner minimap2 --mm2-params "-cx asm5 -t {threads}" -o {params}
        ragtag.py scaffold {config[ref_genome]} {params}/ragtag.patch.fasta -o {params} -u --mm2-params "-cx asm5 -t {threads}" -r -m 1000000
        sed 's/_RagTag//g' {params}/ragtag.scaffold.fasta > {output}
        '''



