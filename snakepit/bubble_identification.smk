rule minigraph_ggs:
    input:
    output:
        'graph.L{L}.gfa'
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {input} > {output}'

rule minigraph_call:
    input:
        gfa = 'graph.L{L}.fa',
        fa = '{asm}.fa'
    output:
        '{asm}.path.{L}.bed'
    threads: 8
    resources:
        mem_mb = 7500
        walltime = '1:00'
    shell:
        'minigraph --call -xasm -t {threads} {input.gfa} {input.fa} > {output}
