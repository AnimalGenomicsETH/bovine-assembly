localrules: determine_mash_ordering, mash_dist, make_unique_names, extract_chromosome, extract_unmapped_regions, samtools_faidx

from pathlib import PurePath, Path

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base,ext='',**kwargs):
    if base == 'SV':
        base_dir = 'SVs'
    elif base == 'mash':
        base_dir = 'mash'
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
        expand(get_dir('SV','{asm}.path.{chr}.L{L}.unmapped.bed'),L=config['L'],asm=list(config['assemblies'].keys())[1:],chr=valid_chromosomes(config['chromosomes']))
        #(get_dir('SV','{asm}.path.7.{L}.unmapped.bed',L=config['L'],asm=ASM) for ASM in list(config['assemblies'].keys())[1:])

rule mash_sketch:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        get_dir('mash','{asm}.msh')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 4
    resources:
        mem_mb = 3000
    shell:
        'mash sketch -p {threads} {input} -o {params.out}'

rule mash_dist:
    input:
        ref = get_dir('mash','{ref}.msh'),
        asm = get_dir('mash','{asm}.msh')
    output:
        get_dir('mash','{asm}.{ref}.dist')
    threads: 1
    resources:
        mem_mb = 3000
    shell:
        '''
        mash dist -p {threads} {input} | awk '{{print "{wildcards.asm} "$3" "$5}}' > {output}
        '''

rule determine_mash_ordering:
    input:
        expand(get_dir('mash','{asm}.{ref}.dist'),ref=list(config['assemblies'].keys())[0],asm=list(config['assemblies'].keys())[1:])
    output:
        get_dir('mash','order.txt')
    shell:
        'cat {input} | sort -n -k 2,2 > {output}'


rule make_unique_names:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        temp(get_dir('SV','{asm}.fasta'))
    shell:
        '''
        sed 's/_RagTag//g' {input} > {output}
        sed -i '/^>/ s/$/_{wildcards.asm}/' {output}
        '''

rule extract_chromosome:
    input:
        multiext(get_dir('SV','{asm}.fasta'),'','.fai')
    output:
        temp(get_dir('SV','{asm}.{chr}.fasta'))
    threads: 1
    resources:
        mem_mb = 1500,
        walltime = '5'
    shell:
        'samtools faidx {input[0]} {wildcards.chr}_{wildcards.asm} > {output}'

rule minigraph_ggs:
    input:
        fastas = (get_dir('SV',f'{asm}.{{chr}}.fasta') for asm in config['assemblies']),
        distances = get_dir('mash','order.txt')
    output:
        get_dir('SV','{chr}.L{L}.gfa')
    params:
        backbone = lambda wildcards: get_dir('SV',f'{list(config["assemblies"].keys())[0]}.{wildcards.chr}.fasta'),
        ordered_input = lambda wildcards, input: [get_dir('SV',f'{line.split()[0]}.{wildcards.chr}.fasta') for line in open(input.distances)],
        asm = '-g250k -r250k -j0.05 -l250k'
    threads: lambda wildcards: 16 if wildcards.chr == 'all' else 1
    resources:
        mem_mb = 15000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '4:00'
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {params.backbone} {params.ordered_input} > {output}'

rule gfatools_bubble:
    input:
        get_dir('SV','{chr}.L{L}.gfa')
    output:
        get_dir('SV','bubble.L{L}.bed')
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '10'
    shell:
        'gfatools bubble {input} > {output}'

rule minigraph_align:
    input:
        gfa = get_dir('SV','{chr}.L{L}.gfa'),
        fasta = get_dir('SV','{asm}.{chr}.fasta')
    output:
        get_dir('SV','{asm}.path.{chr}.L{L,\d+}.{ext,gaf|bed}')
    params:
        lambda wildcards: '--call' if wildcards.ext == 'bed' else '--show-unmap=yes'
    threads: lambda wildcards: 12 if wildcards.chr == 'all' else 1
    resources:
        mem_mb = 20000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '30'
    shell:
        'minigraph {params} -xasm -t {threads} -l 300k  {input.gfa} {input.fasta} > {output}'

rule extract_unmapped_regions:
    input:
        fai = get_dir('SV','{asm}.{chr}.fasta.fai'),
        gaf = get_dir('SV','{asm}.path.{chr}.L{L}.gaf'),
        bed = get_dir('SV','{asm}.path.{chr}.L{L}.bed')
    output:
        genome = temp(get_dir('SV','{asm}.path.{chr}.L{L}.genome')),
        gaf_bed = temp(get_dir('SV','{asm}.path.{chr}.L{L}.aligned.bed')),
        unmapped = temp(get_dir('SV','{asm}.path.{chr}.L{L}.unmapped.bed'))
    shell:
        '''
        awk '{{print $1"\t"$2}}' {input.fai} > {output.genome}
        awk '/^{wildcards.chr}_{wildcards.asm}/ {{print $1"\t"$3"\t"$4}}' {input.gaf} | sort -n -k 2,2 > {output.gaf_bed}
        awk '$6=="." {{print "{wildcards.chr}_{wildcards.asm}\t"$2"\t"$3}}' {input.bed} >> {output.gaf_bed}

        #also bed file regions that aren't "."
        #merge?
        bedtools complement -i {output.gaf_bed} -g {output.genome} > {output.unmapped}
        '''


def get_bubble_links(bed,bubbles,mode):

    name,ext = get_name(bed,mode)

    with open(bed,'r') as fin:
        for line in fin:
            source,sink,links = line.rstrip().split()[3:6]

            if links[0] == '.':
                continue
            elif links[0] == '*':
                bubbles[name].append(f'{source[1:]}_DEL{ext}')
            else:
                for edge in links.split(':')[0].split('>')[1:]:
                    bubbles[name].append(f'{edge}{ext}')

def get_name(full_name,mode='both'):
    name = PurePath(full_name).name.split('.')[0]
    if mode == 'both':
        return name,''
    else:
        breed,read = name.split('_')
        if mode == 'breed':
            return breed,read
        else:
            return read,breed

rule colour_shared_bubbles:
    input:
        (get_dir('SV','{asm}.path.{chr}.L{L}.bed',asm=ASM) for ASM in list(config['assemblies'].keys())[1:])
    output:
        get_dir('SV','{chr}.L{L}.bubbles.{mode}.png')
    run:
        from upsetplot import from_contents, plot
        from matplotlib import pyplot as plt
        from collections import defaultdict
        bubbles = defaultdict(list)

        for infile in input:
            get_bubble_links(infile,bubbles,wildcards.mode)
        df = from_contents(bubbles)
        plot(df,sort_by='cardinality')
        plt.savefig(output[0])
        df.to_csv(f'{wildcards.chr}.df')

#df = pd.read_csv('/Users/alexleonard/Documents/Tiergenomik/DATA/test.df',index_col=[0,1,2,3,4])
#us = upsetplot.UpSet(df)
#dist['N_O']=us.intersections[:,True,False,:,:].sum() + us.intersections[:,False,True,:,:].sum()
# z=hierarchy.linkage(dist.values(), 'single')
#dn1 = hierarchy.dendrogram(z, above_threshold_color='y',orientation='top',labels=['PM','N','O','BSW','G'])
#
rule samtools_faidx:
    input:
        '{fasta}.{fa_ext}'
    output:
        '{fasta}.{fa_ext,fa|fasta}.fai'
    threads: 1
    resources:
        mem_mb = 1500,
        walltime = '5'
    shell:
        'samtools faidx {input}'
