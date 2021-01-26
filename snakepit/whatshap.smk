localrules: whatshap_tabix

rule whatshap_phase:
    input:
        vcf = get_dir('output','{haplotype}.unphased.vcf.gz'),
        bam = get_dir('input','{haplotype}.unphased.{model}.bam'),
        ref = config['reference']
    output:
        temp(get_dir('input','{haplotype}.phasing.{model}.vcf.gz'))
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
        get_dir('input','{haplotype}.phasing.{model}.vcf.gz')
    output:
        temp(get_dir('input','{haplotype}.phasing.{model}.vcf.gz.tbi'))
    shell:
        'tabix -p vcf {input}'

rule whatshap_haplotag:
    input:
        vcf = multiext(get_dir('input','{haplotype}.phasing.{model}.vcf.gz'),'','.tbi'),
        bam = get_dir('input','{haplotype}.unphased.{model}.bam'),
        ref = config['reference']
    output:
        temp(get_dir('input','{haplotype}.phased.{model}.bam')),
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''
        whatshap haplotag \
        --output {output} \
        --reference {input.ref} \
        {input.vcf[0]} \
        {input.bam}
        '''
