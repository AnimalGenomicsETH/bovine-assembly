localrules: make_unique_names, extract_chromosome, extract_unmapped_regions, samtools_faidx

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

rule all:
    input:
        (get_dir('SV','{asm}.path.7.{L}.unmapped.bed',L=config['L'],asm=ASM) for ASM in config['assemblies'][1:])

rule make_unique_names:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        temp(get_dir('SV','{asm}.fasta'))
    shell:
        'sed "/^>/ s/$/_{wildcards.asm}}/" {input} > {output}'

rule extract_chromosome:
    input:
        multiext(get_dir('SV','{asm}.fasta'),'','.fai')
    output:
        temp(get_dir('SV','{asm}.{chr}.fasta'))
    shell:
        'samtools faidx {input} {wildcards.chr}_{wildcards.asm} > {output}'

rule minigraph_ggs:
    input:
        (get_dir('SV',f'{asm}.{{chr}}.fasta') for asm in config['assemblies'])
    output:
        get_dir('SV','{chr}.L{L}.gfa')
    threads: 16
    resources:
        mem_mb = 6000,
        walltime = '4:00'
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {input} > {output}'

rule minigraph_align:
    input:
        gfa = get_dir('SV','{chr}.L{L}.gfa'),
        fasta = get_dir('SV','{asm}.{chr}.fasta')
    output:
        get_dir('SV','{asm}.path.{chr}.{L}.{ext,gaf|bed}')
    params:
        lambda wildcards: '--call' if wildcards.ext == 'bed' else ''
    threads: 12
    resources:
        mem_mb = 7000,
        walltime = '4:00'
    shell:
        'minigraph {params} -xasm -t {threads} {input.gfa} {input.fa} > {output}'

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
