localrules: plot_dot

from pathlib import PurePath, Path

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base,ext='',**kwargs):
    if base == 'synteny':
        base_dir = 'synteny_{ref}'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))


rule all:
    input:
        (get_dir('synteny','{hap}.jupiter.png',ref=REF,hap=HAP) for REF in config['pairs'] for HAP in config['pairs'][REF])

rule jupiter_plot:
    input:
        ref = lambda wildcards: config['references'][wildcards.ref],
        fasta = lambda wildcards: config['pairs'][wildcards.ref][wildcards.hap]
    output:
        paf = get_dir('synteny','{hap}-agp.paf'),
        png = get_dir('synteny','{hap}.jupiter.png')
    params:
        lambda wildcards,output: PurePath(output.png).parent
    threads: 8
    resources:
        mem_mb = 4000,
        walltime = '2:00'
    shell:
        '''
        mkdir -p {params} && cd {params}
        jupiter name={wildcards.hap} ref=../{input.ref} fa=../{input.fasta} ng=100 maxScaff=30 m=20000000 g=99 gScaff=99 t=8
        mv {wildcards.hap}.png ../{output.png}
        '''

rule plot_dot:
    input:
        get_dir('synteny','{hap}-agp.paf')
    output:
        get_dir('synteny','{hap}.dot.png')
    shell:
        'minidot -L {input} | convert -density 150 - {output}'

