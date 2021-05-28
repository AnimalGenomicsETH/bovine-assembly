DIR='synteny/'

rule jupiter_prepare:
    input:
        ref = config['ref'],
        fasta = config['fasta']
    output:
        ref = DIR + '{hap}_reference.fa',
        fasta = DIR + '{hap}-agp.fa'
    shell:
        '''
        jupiter name={wildcards.hap} ref={input.ref} fa={input.fasta} ng=95 m=1000000 g=100 gScaff=100 {wildcards.hap}.agp {wildcards.hap}_reference.karyotype 2> /dev/null
        '''

rule minimap_align:
    input:
        ref = DIR + '{hap}_reference.fa',
        fasta = DIR + '{hap}-agp.fa'
    output:
        DIR + '{hap}-agp.sam'
    threads: 8
    resources:
        mem_mb = 4000,
        walltime = '1:00'
    shell:
        'minimap2 -ax asm5 -t {threads} {input.ref} {input.fasta} > {output}'

rule jupiter_plot:
    input:
        DIR + '{hap}-agp.sam'
    output:
        '{hap}.png'
    shell:
        '''
        jupiter name={wildcards.hap} ref=config[ref] fa=config[fasta] ng=95 m=1000000 g=100 gScaff=100 2> /dev/null
        '''

rule minimap_align:
    input:
        ref = config['ref'],
        fasta = config['fasta']
    output:
        DIR + '{hap}_ref.paf'
    threads: 8
    resources:
        mem_mb = 4000,
        walltime = '1:00'
    shell:
        'minimap2 -cx asm5 -t {threads} {input.ref} {input.fasta} > {output}'


rule plot_dot:
    input:
        DIR + '{hap}_ref.paf'
    output:
        DIR + '{hap}_ref.png'
    shell:
        'minidot -L {input} | convert -density 150 - {output}'

