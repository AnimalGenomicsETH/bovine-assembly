rule all:
    input:
        expand('{asm}.path.{L}.bed',asm=config['asm'],L=config['L'])

print(config['asm'])
rule minigraph_ggs:
    input:
        config['asm'].values()
    output:
        'graph.L{L}.gfa'
    threads: 16
    resources:
        mem_mb = 4500,
        walltime = '2:00'
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {input} > {output}'

rule minigraph_call:
    input:
        gfa = 'graph.L{L}.gfa',
        fa = lambda wildcards: config['asm'][wildcards.asm]
    output:
        '{asm}.path.{L}.bed'
    threads: 8
    resources:
        mem_mb = 7500,
        walltime = '1:00'
    shell:
        'minigraph --call -xasm -t {threads} {input.gfa} {input.fa} > {output}'
