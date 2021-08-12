
import subprocess
from itertools import product
import pandas as pd

taur_h = list('_'.join(p) for p in product(('O','B','P'),('peregrine','hicanu','hifiasm')))
taur_s = list('_'.join(p) for p in product(('O','B','P'),('shasta','flye','raven')))

run_n = 11
root = f'SV_run_{run_n}_ARS/'

g=pd.read_csv(root+'all.L50.both.df',index_col=list(range(10)))
gc = g[~g.id.str.contains('REFDEL')]
gc = gc[~gc.id.str.contains('MISDEL')]
gc = gc[~gc.id.str.contains('INV')]

slicer_h = 'gc'
slicer_s = 'gc'
for i in taur_h[2::3]:
    slicer_h += f'.xs(True,level="{i}",drop_level=False)'
    slicer_s += f'.xs(False,level="{i}",drop_level=False)'
for i in taur_s[::3]:
    slicer_h += f'.xs(False,level="{i}",drop_level=False)'
    slicer_s += f'.xs(True,level="{i}",drop_level=False)'


def shellx(cmd):
    #print('STDERR: running ',cmd)
    subprocess.run(cmd, shell=True,stderr=subprocess.DEVNULL,stdout=subprocess.DEVNULL)

def find_repeats(bubbles,dtype):
    repeat_rows = []
    for i in bubbles:
        split = i.split('_')
        chrom, node = split[-1], split[0]
        #fname = get_dir('SV','{chr}.L{L}.gfa',chr=chrom,**wildcards)
        fname = f'{root}{chrom}.L50.gfa'
        search_string = f'\'$1=="S" && $2=="{node}" {{print ">{node}\\n"$3}}\''
        #fout = get_dir('repeat',f'{chrom}.{node}.fa',**wildcards)
        fout = f'{root}repeat_nodes/{chrom}.{node}.fa'

        shellx(f'awk {search_string} {fname} > {fout}')
        shellx(f'RepeatMasker -xsmall -lib /cluster/work/pausch/alex/Libraries/BosTau9_repeat_library.fasta -qq -no_is {fout}')
        
        try:
            row = extract_repeats(fout+'.tbl')
        except FileNotFoundError:
            row = dict()

        row.update({'data':dtype,'chr':chrom,'node':node})
        repeat_rows.append(row)
    return repeat_rows

def extract_repeats(fname):
    elements = ('SINEs','LINEs','LTR elements', 'Satellites', 'Unclassified')
    repeats = dict()

    with open(fname,'r') as file_:
        for line in file_:
            if 'total length' in line:
                repeats['length'] = int(line.split('(')[1].split()[0])
            elif 'GC level' in line:
                repeats['GC'] = float(line.split()[2])
            elif 'bases masked' in line:
                repeats['masked'] = float(line.split()[5])
            for element in elements:
                if element in line:
                    repeats[element] = float(line.split('bp')[1].split()[0])
    return repeats

h_only = eval(slicer_h)
s_only = eval(slicer_s)
rs = []
rs.extend(find_repeats(h_only[h_only['length']>1000]['id'],'HiFi'))
rs.extend(find_repeats(s_only[s_only['length']>1000]['id'],'ONT'))

df = pd.DataFrame(rs)
df.to_csv(f'bubbles_repeat.{run_n}.df',index=False)
