from pathlib import PurePath

localrules: generate_genomes

rule all:
    input:
        'output_plot.png'

rule samtools_faidx:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        'syri_test/{asm}.fa'
    params:
        regions = ' '.join(list(map(str,range(1,30)))+['X1','X2'])
    shell:
        '''
        samtools faidx {input} {params.regions} | seqtk rename - chr > {output}
        '''

rule minimap2_align:
    input:
        A = 'syri_test/{A}.fa',
        B = 'syri_test/{B}.fa'
    output:
        'syri_test/{A}_{B}.paf'
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        '''
        minimap2 -t {threads} -cxasm10 --cs --eqx {input.A} {input.B} > {output}
        '''

rule syri_annotate:
    input:
        A = 'syri_test/{A}.fa',
        B = 'syri_test/{B}.fa',
        paf = 'syri_test/{A}_{B}.paf'
    output:
        'syri_test/{A}_{B}.syri.out'
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').name + '.',
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 1
    resources:
        mem_mb = 25000
    shell:
        '''
        syri -c {input.paf} -F P -r {input.A} -q {input.B} --dir {params._dir} --prefix {params.prefix}
        '''

rule generate_genomes:
    output:
        'syri_test/genomes.txt'
    run:
        with open(output[0],'w') as fout:
            fout.write('\t'.join(('#file','name','tags'))+'\n')
            for name,asm in config['assemblies'].items():
                fout.write('\t'.join((f'syri_test/{name}.fa',name,'lw:1.5'))+'\n')

def chain_input():
    assemblies = list(config['assemblies'].keys())
    for A,B in zip(assemblies,assemblies[1:]):
        yield f'syri_test/{A}_{B}.syri.out'

rule plotsr:
    input:
        outs = chain_input(),
        genomes = 'syri_test/genomes.txt'
    output:
        'output_plot.png'
    params:
        SRs = lambda wildcards, input: ' '.join((f'--sr {out}' for out in input.outs)) 
    shell:
        '''
        plotsr \
        {params.SRs} \
        --genomes {input.genomes} \
        -o {output}
        '''
