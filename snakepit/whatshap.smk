rule whatshap_phase:
    input:
        vcf = '',
        bam = '',
        ref = ''
    output:
        'phased.vcf.gz'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24:00'
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
        'whatshap/deepvariant1.phased.vcf.gz.tbi'
    shell:
        'tabix -p vcf {input}'

rule whatshap_haplotag:
    input:
        vcf = 'phased.vcf.gz',
        bam = '',
        ref = ''
    output:
        'haplotagged.bam'
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24:00'
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
