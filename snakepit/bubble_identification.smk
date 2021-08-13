localrules: determine_ordering, make_unique_names, colour_shared_bubbles, plot_dendrogram, merge_unmapped

from pathlib import PurePath, Path
from collections import defaultdict
from itertools import product

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base,ext='',**kwargs):
    if base == 'SV':
        base_dir = 'SV_run_{run}_' + list(config['assemblies'].keys())[0]
    elif base == 'repeat':
        base_dir = get_base('SV','bubble_repeats',**kwargs)
    elif base == 'mash':
        base_dir = 'mash'
    elif base == 'fasta':
        base_dir = 'fasta'
    else:
        raise Exception('Base not found')
    if ext and ext[0] == '.':
        return f'{base_dir}{ext}'.format_map(Default(kwargs))
    return str(PurePath(base_dir.format_map(Default(kwargs))) / ext.format_map(Default(kwargs)))

def valid_chromosomes(chr_target):
    all_chrs = [*map(str,range(1,30))]
    if chr_target == 'all':
        return all_chrs
    else:
        return chr_target

order_functions = dict()
for r_range, func in config['runs'].items():
    #unique function
    if '-' not in r_range:
        order_functions[r_range] = func
    #repeat function
    else:
        l,h = (int(i) for i in r_range.split('-'))
        for n in range(l,h):
            order_functions[str(n)] = func

rule all:
    input:
        #expand(get_dir('SV','{asm}.path.{chr}.L{L}.unmapped.bed'),L=config['L'],asm=list(config['assemblies'].keys())[1:],chr=valid_chromosomes(config['chromosomes'])),
        #get_dir('SV','dendrogram.L{L}.{mode}.png',run='static',L=config['L'],mode='breed'),
        #(get_dir('SV','{chr}.L{L}.{mode}.nw',chr=c,run='static',L=config['L'],mode='breed') for c in range(1,30)),
        #(get_dir('SV','{chr}.L{L}.{mode}.nw',chr=c,run=n,L=config['L'],mode='breed') for (c,n) in product(range(1,30),target_names)),
        #(get_dir('SV','{chr}.L{L}.{mode}.nw',chr=c,run=n,L=config['L'],mode='both') for (c,n) in product(range(1,30),order_functions.keys())),
        (get_dir('SV','all.L{L}.{mode}.nw',run=n,L=config['L'],mode='both') for n in order_functions.keys())

rule mash_sketch:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        get_dir('mash','{asm}.msh')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 2
    resources:
        mem_mb = 1000,
        walltime = '30'
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
        mem_mb = 100,
        walltime = '5'
    shell:
        '''
        mash dist -p {threads} {input} | awk '{{print "{wildcards.asm} "$3" "$5}}' > {output}
        '''

import random
checkpoint determine_ordering:
    input:
        expand(get_dir('mash','{asm}.{ref}.dist'),ref=list(config['assemblies'].keys())[0],asm=list(config['assemblies'].keys())[1:])
    output:
        get_dir('SV','order.txt')
    run:
        print("creating path", Path(output[0]).parent)
        Path(output[0]).parent.mkdir(exist_ok=True,parents=True)
        if wildcards.run == 'static':
            shell('sort -k2,2n {input} > {output}')
        else:
            idx = list(range(len(input)))
            #idx = rng.integers(0,2,5)*5+[1,2,3,4,5]#np.random.randint(0,2,5)+range(1,10,2) #offset random samples to pick either HiFi or ONT
            choices = eval(order_functions[wildcards.run])

            assemblies = list(config['assemblies'].keys())[1:]
            with open(output[0],'w') as fout:
                fout.write('\n'.join([assemblies[i] for i in choices]))

rule make_unique_names:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        temp(get_dir('fasta','{asm}.fasta'))
    shell:
        '''
        sed 's/_RagTag//g' {input} > {output}
        sed -i '/^>/ s/$/_{wildcards.asm}/' {output}
        '''

rule extract_chromosome:
    input:
        multiext(get_dir('fasta','{asm}.fasta'),'','.fai')
    output:
        temp(get_dir('fasta','{asm}.{chr}.fasta',chr=CHR) for CHR in valid_chromosomes(config['chromosomes']))
    threads: 1
    resources:
        mem_mb = 5000,
        walltime = '5'
    run:
        for chromosome in valid_chromosomes(config['chromosomes']):
            out_file = output[int(chromosome)-1]
            shell(f'samtools faidx {{input[0]}} {chromosome}_{wildcards.asm} > {out_file}')

rule minigraph_ggs:
    input:
        fastas = (get_dir('fasta','{asm}.{chr}.fasta',asm=asm) for asm in config['assemblies']),
        order = get_dir('SV','order.txt')
    output:
        gfa = get_dir('SV','{chr}.L{L}.gfa')
    params:
        backbone = lambda wildcards: get_dir('fasta',f'{list(config["assemblies"].keys())[0]}.{wildcards.chr}.fasta',**wildcards),
        input_order = lambda wildcards, input: [get_dir('fasta',f'{line.split()[0]}.{{chr}}.fasta',**wildcards) for line in open(input.order)],
        asm = ' --gg-match-pen 5' #NOTE not currently set in the shell #'-g250k -r250k -j0.05 -l250k'
    threads: lambda wildcards: 16 if wildcards.chr == 'all' else 1
    resources:
        mem_mb = 10000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '4:00'
    shell:
        'minigraph -xggs -t {threads} -L {wildcards.L} {params.backbone} {params.input_order} > {output}'

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
        fasta = get_dir('fasta','{asm}.{chr}.fasta')
    output:
        get_dir('SV','{asm}.path.{chr}.L{L,\d+}.{ext,gaf|bed}')
    params:
        lambda wildcards: '--call' if wildcards.ext == 'bed' else '--show-unmap=yes'
    threads: lambda wildcards: 12 if wildcards.chr == 'all' else 1
    resources:
        mem_mb = 8000,
        walltime = lambda wildcards: '4:00' if wildcards.chr == 'all' else '4:00'
    shell:
        'minigraph {params} -xasm -t {threads} -l 300k  {input.gfa} {input.fasta} > {output}'

rule extract_unmapped_regions:
    input:
        fai = get_dir('fasta','{asm}.fasta.fai'),
        gaf = get_dir('SV','{asm}.path.{chr}.L{L}.gaf'),
        bed = get_dir('SV','{asm}.path.{chr}.L{L}.bed')
    output:
        genome = temp(get_dir('SV','{asm}.path.{chr}.L{L}.genome')),
        gaf_bed = temp(get_dir('SV','{asm}.path.{chr}.L{L}.aligned.bed')),
        unmapped = temp(get_dir('SV','{asm}.path.{chr}.L{L}.unmapped.bed'))
    threads: 1
    resources:
        walltime = '5'
    shell:
        '''
        awk '/^{wildcards.chr}_{wildcards.asm}/ {{print "{wildcards.chr}\\t"$2}}' {input.fai} > {output.genome}
        awk '/^{wildcards.chr}_{wildcards.asm}/ {{print "{wildcards.chr}\\t"$3"\\t"$4}}' {input.gaf} | sort -k2,2n | bedtools complement -i - -g {output.genome} | awk '{{print $0"\\tunaligned"}}' > {output.gaf_bed}
        awk '$6=="." {{print "{wildcards.chr}\\t"$2"\\t"$3}}' {input.bed} | awk '{{print $0"\\tuncalled"}}' >> {output.gaf_bed}
        sort -k2,2n {output.gaf_bed} | bedtools merge -d 10000 -c 4 -o distinct > {output.unmapped}
        '''

rule merge_unmapped:
    input:
        (get_dir('SV','{asm}.path.{chr}.L{L}.unmapped.bed',chr=CHR) for CHR in range(1,30))
    output:
        get_dir('SV','{asm}.path.L{L}.unmapped.bed')
    shell:
        'cat {input} > {output}'

import re
def get_bubble_links(bed,bubbles,mode):
    name,ext = get_name(bed,mode)
    with open(bed,'r') as fin:
        for line in fin:
            source,sink,links = line.rstrip().split()[3:6]

            if links[0] == '.':
                bubbles[name].append(f'{source[1:]}_MISDEL{ext}')
            elif links[0] == '*':
                if (int(sink[2:])-int(source[2:])) == 1:
                    bubbles[name].append(f'{source[1:]}_REFDEL{ext}')
                else:
                    bubbles[name].append(f's{int(source[2:])+1}_NOVDEL{ext}')
                #note, if source and sink are continuous numbers, then it is a "non-ref" deletion
                #if they are +2, then it is a deletion unique to the additional assemblies
                #bubbles[name].append(f'{source[1:]}_DEL{ext}')
            else:
                linkage = links.split(':')[0]
                nodes = re.split('>|<',linkage)[1:]
                inversions = re.findall('>|<',linkage)
                for edge, inv in zip(nodes,inversions):
                    inv_ext = '_INV' if inv == '<' else ''
                    bubbles[name].append(f'{edge}{ext}{inv_ext}')

def get_name(full_name,mode='both'):
    (name,_,chromosome,*_) = PurePath(full_name).name.split('.')
    if mode == 'both':
        return name,f'_{chromosome}'
    else:
        breed,read = name.split('_')
        if mode == 'breed':
            return breed,f'{read}_{chromosome}'
        else:
            return read,breed

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return f'{leaf_names[node.id]}:{parentdist-node.dist:.2f}{newick}'
    if len(newick) > 0:
        newick = f'):{parentdist-node.dist:.2f}{newick}'
    else:
        newick = ');'
    newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
    newick = getNewick(node.get_right(), f',{newick}', node.dist, leaf_names)
    return '(' + newick

def get_order(wildcards):
    checkpoint_output = PurePath(checkpoints.determine_ordering.get(**wildcards).output[0])
    orders = []
    for ASM in open(get_dir('SV','order.txt',**wildcards),'r'):
        asm = ASM.split()[0]
        for chrom in (range(1,30) if wildcards.chr == 'all' else [wildcards.chr,]):
            orders.append(get_dir('SV','{asm}.path.{chr}.L{L}.bed',asm=asm,chr=chrom))
    return orders

def add_bubble_lengths(df,wildcards):
    size_map = dict()

    for chrom in (range(1,30) if wildcards.chr == 'all' else [wildcards.chr,]):
        for line in open(get_dir('SV','{chrom}.L{L}.gfa',**wildcards,chrom=chrom),'r'):
            if line[0] != 'S':
                continue
            (_,node,_,length_code,*_) = line.split()
            identifier = f'{node}_{chrom}' #if wildcards.chr == 'all' else node
            size_map[identifier] = int(length_code.split(':')[-1])
    
    df['length'] = [size_map[n.replace('_INV','')] if 'DEL' not in n else 0 for n in df['id']]
    return df

rule generate_newick_tree:
    input:
        order = get_dir('SV','order.txt'),
        beds = lambda wildcards: get_order(wildcards)
    output:
        newick = get_dir('SV','{chr}.L{L}.{mode}.nw'),
        df = get_dir('SV','{chr}.L{L}.{mode}.df')
    resources:
        walltime = '5'
    run:
        from upsetplot import from_contents

        bubbles = defaultdict(list)

        for infile in input.beds:
            get_bubble_links(infile,bubbles,wildcards.mode)
        df = from_contents(bubbles)
        df = add_bubble_lengths(df,wildcards)
        df.to_csv(output.df)
        
        df_pruned = df[~df.id.str.contains('REFDEL')]
        df_pruned = df_pruned[~df_pruned.id.str.contains('MISDEL')]
        newick = form_tree(df_pruned)
        with open(output['newick'],'w') as fout:
            fout.write(newick)

def form_tree(data,no_plot=False):
    from scipy.cluster import hierarchy
    names, dist = combine(data)
    z = hierarchy.linkage([float(i) for i in dist.values()], 'average')
    dn1 = hierarchy.dendrogram(z, above_threshold_color='y',orientation='top',labels=names,no_plot=no_plot)
    tree = hierarchy.to_tree(z,False)
    newick = getNewick(tree, "", tree.dist, names)
    return newick

def get_pairs(tree):
    inner_pair = tree.split('(')[-1].split(')')[0]
    return (inner_pair.split(':')[0],inner_pair.split(',')[-1].split(':')[0])

def combine(data):
    import upsetplot
    names = data.index.names
    intersections = upsetplot.UpSet(data).intersections
    condensed_dist = dict()
    for i in range(len(names)-1):
        for j in range(i+1,len(names)):
            try:
                condensed_dist[f'{names[i]}_{names[j]}'] = intersections.xs(True,level=names[i]).xs(False,level=names[j]).sum()+intersections.xs(False,level=names[i]).xs(True,level=names[j]).sum()
            except: #combination doesn't exist, so trivially 0
                condensed_dist[f'{names[i]}_{names[j]}'] = 0
    return names, condensed_dist

from collections import defaultdict
def join_breed_haplotypes(dist):
    joined_breeds = defaultdict(int)
    for k,v in dist.items():
        (p1,_,p2,_) = k.split('_')
        if p1 == p2:
            continue
    joined_breeds[tuple(sorted((p1,p2)))] += v

#R plot
#library(phylogram)
#library(dendextend)
# x <- read.dendrogram(text = "(gaur:0.389885,(nellore:0.118569,(bsw:0.0403885,(obv:0.0394899,pied:0.0394899):0.000898574):0.0781807):0.271316);")
#y
#dy <- ladder(y)
#> dndlist <- dendextend::dendlist(dx,dy) %>% set("highlight_branches_col";)
#tanglegram(dndlist, sort = TRUE, common_subtrees_color_lines = TRUE, highlight_distinct_edges  = TRUE,)
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
