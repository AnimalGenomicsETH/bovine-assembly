rule all:
    input:
        expand('{asm}.path.{L}.bed',asm=config['assemblies'],L=config['L'])

rule make_unique_names:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        temp(get_dir('SV','{asm}.fasta'))
    shell:
        'sed "/^>/ s/$/_{wildcards.asm}}/" {input} > {output}'

rule minigraph_ggs:
    input:
        (get_dir('SV',f'{asm}.fasta') for asm in config['assemblies'])
    output:
        get_dir('SV','graph.L{L}.gfa')
    threads: 16
    resources:
        mem_mb = 4500,
        walltime = '4:00'
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {input} > {output}'

rule minigraph_call:
    input:
        gfa = get_dir('SV','graph.L{L}.gfa'),
        fasta = get_dir('SV','{asm}.fasta')
    output:
        get_dir('SV','{asm}.path.{L}.bed')
    threads: 8
    resources:
        mem_mb = 7500,
        walltime = '2:00'
    shell:
        'minigraph --call -xasm -t {threads} {input.gfa} {input.fa} > {output}'
