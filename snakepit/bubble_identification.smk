localrules: make_unique_names, extract_chromosome, extract_unmapped_regions, samtools_faidx

from pathlib import PurePath, Path

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base,ext='',**kwargs):
    if base == 'SV':
        base_dir = 'SVs'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

def valid_chromosomes(chr_target):
    all_chrs = [*map(str,range(1,30)),'X']
    if chr_target == 'all':
        return all_chrs
    else:
        return chr_target

rule all:
    input:
        expand(get_dir('SV','{asm}.path.{chr}.{L}.unmapped.bed'),L=config['L'],asm=list(config['assemblies'].keys())[1:],chr=valid_chromosomes(config['chromosomes']))
        #(get_dir('SV','{asm}.path.7.{L}.unmapped.bed',L=config['L'],asm=ASM) for ASM in list(config['assemblies'].keys())[1:])

rule make_unique_names:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        temp(get_dir('SV','{asm}.fasta'))
    shell:
        'sed "/^>/ s/$/_{wildcards.asm}/" {input} > {output}'

rule extract_chromosome:
    input:
        multiext(get_dir('SV','{asm}.fasta'),'','.fai')
    output:
        temp(get_dir('SV','{asm}.{chr}.fasta'))
    shell:
        'samtools faidx {input[0]} {wildcards.chr}_{wildcards.asm} > {output}'

rule minigraph_ggs:
    input:
        (get_dir('SV',f'{asm}.{{chr}}.fasta') for asm in config['assemblies'])
    output:
        get_dir('SV','{chr}.L{L}.gfa')
    threads: lambda wildcards: 16 if wildcards.chr == 'all' else 2
    resources:
        mem_mb = 5000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '30'
    shell:
        'minigraph -xggs -t {threads} -l 300k -L {wildcards.L} {input} > {output}'

rule minigraph_align:
    input:
        gfa = get_dir('SV','{chr}.L{L}.gfa'),
        fasta = get_dir('SV','{asm}.{chr}.fasta')
    output:
        get_dir('SV','{asm}.path.{chr}.{L,\d+}.{ext,gaf|bed}')
    params:
        lambda wildcards: '--call' if wildcards.ext == 'bed' else ''
    threads: lambda wildcards: 12 if wildcards.chr == 'all' else 2
    resources:
        mem_mb = 6000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '30'
    shell:
        'minigraph {params} -xasm -t {threads} -l 300k --show-unmap=yes {input.gfa} {input.fasta} > {output}'

rule extract_unmapped_regions:
    input:
        fai = get_dir('SV','{asm}.{chr}.fasta.fai'),
        gaf = get_dir('SV','{asm}.path.{chr}.{L}.gaf'),
        bed = get_dir('SV','{asm}.path.{chr}.{L}.bed')
    output:
        genome = temp(get_dir('SV','{asm}.path.{chr}.{L}.genome')),
        gaf_bed = temp(get_dir('SV','{asm}.path.{chr}.{L}.aligned.bed')),
        unmapped = temp(get_dir('SV','{asm}.path.{chr}.{L}.unmapped.bed'))
    shell:
        '''
        awk '{{print $1"\t"$2}}' {input.fai} > {output.genome}
        awk '/^7_hifiasm/ {{print $1"\t"$3"\t"$4}}' {input.gaf} | sort -n -k 2,2 > {output.gaf_bed}
        #also bed file regions that aren't "."
        #merge?
        bedtools complement -i {output.gaf_bed} -g {output.genome} > unmapped.bed
        '''


rule samtools_faidx:
    input:
        '{fasta}.{fa_ext}'
    output:
        '{fasta}.{fa_ext,fa|fasta}.fai'
    shell:
        'samtools faidx {input}'
