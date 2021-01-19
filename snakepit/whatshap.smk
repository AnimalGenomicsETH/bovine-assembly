rule whatshap_phase:
    input:
        vcf = '',
        bam = '',
        ref = ''
    output:
        'phased.vcf'
    shell:
        '''
        whatshap phase \
        --output {output} \
        --reference {input.ref} \
        {inpu.vcf} \
        {input.bam}
        '''
rule whatshap_tabix:
    input:
        'whatshap/deepvariant1.phased.vcf.gz'
    output:
        ''
    shell:
        'tabix -p vcf {input}'

rule whatshap_haplotag:
    input:
        vcf = 'phased.vcf'
    output:
        'haplotagged'
    shell:
        '''
        whatshap haplotag \
        --output {output} \
        --reference {input.ref} \
        {input.vcf} \
        {input.bam}
        '''

rule whatshap_index:
    input:
        ''
    output:
        ''
    shell:
        'samtools index whatshap/HG003.GRCh38.chr20.haplotagged.bam'
