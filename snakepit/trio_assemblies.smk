rule trio_yak:
    input:
        'parental.fq.gz'
    output:
        'parental.yak'
    threads: 16
    resources:
        mem_mb = 2500
    params:
        K = config['k-mers']
    shell:
        'yak count -k {params.K} -b 37 -t {threads} -o {output} {input}'

rule trio_hifiasm:
    input:
        reads = '',
        mat = '',
        pat = '',
        force = '.asm'
    output:
        ''
    threads: 36
    resources:
        mem_mb = 2000,
        walltime = '16:00'
    shell:
        #NOTE /dev/null is used since the overlaps already exist from original hifiasm run
        'hifiasm -o {params.out} -t {threads} -1 {input.pat} -2 {input.mat} /dev/null'


rule trio_canu:

