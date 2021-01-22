rule whatshap_phase:
    input:
        vcf = get_dir('output','{haplotype}.unphased.vcf.gz'),
        bam = get_dir('input','{haplotype}_hifi_reads.unphased.pbmm2.bam'),
        ref = config['reference']
    output:
        temp(get_dir('input','{haplotype}.phasing.vcf.gz'))
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''
        whatshap phase \
        --output {output} \
        --reference {input.ref} \
        {input.vcf} \
        {input.bam}
        '''

rule whatshap_tabix:
    input:
        get_dir('input','{haplotype}.phasing.vcf.gz')
    output:
        temp(get_dir('input','{haplotype}.phasing.vcf.gz.tbi'))
    shell:
        'tabix -p vcf {input}'

rule whatshap_haplotag:
    input:
        vcf = get_dir('input','{haplotype}.phasing.vcf.gz'),
        tbi = get_dir('input','{haplotype}.phasing.vcf.gz.tbi'),
        bam = get_dir('input','{haplotype}_hifi_reads.unphased.pbmm2.bam'),
        ref = config['reference']
    output:
        temp(get_dir('input','{haplotype}_hifi_reads.phased.pbmm2.bam')),
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
