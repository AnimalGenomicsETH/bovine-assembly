import pandas
import numpy as np
from collections import defaultdict, Counter
import matplotlib.pyplot as plt

def load_bed(fname):
    data = []
    with open(fname,'r') as fin:
        for line in fin:
            row = dict()
            parts = line.rstrip().split()
            row['ref_contig'] = parts[0].split('.')[0]
            row['chromosome'] = row['ref_contig'] in (*map(str,range(1,31)),'X')
            row['ref_start']= int(parts[1])
            row['ref_end'] = int(parts[2])
            row['source'] = int(parts[3][2:])
            row['sink'] = int(parts[4][2:])
            if parts[5][0] != '.':

                if parts[5][0] == '*':
                    if row['sink'] == (row['source']+1):
                        row['link'] = 'match'
                    else:
                        row['link'] = 'alt'
                else:
                    row['edge'] = parts[5][2:].split(':')[0]
                    try:
                        link = int(parts[5][2:].split(':')[0])
                        if link > row['source'] and link < row['sink']:
                            row['link'] = 'match'
                        else:
                            row['link'] = 'alt'
                    except ValueError:
                        link = parts[5][2:].split(':')[0]
                        row['link'] = 'multi'


                row['size'] = int(parts[5].split(':')[1])
                row['map_contig'] = parts[5].split(':')[3].split('.')[0]
                if '_' in row['map_contig']:
                    row['map_contig'] = row['map_contig'].split('_')[0]
                row['map_start']= int(parts[5].split(':')[4])
                row['map_end'] = int(parts[5].split(':')[5])
            data.append(row)
    return pandas.DataFrame(data)


hifiasm = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/hifiasm_graph.bed')
canu = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/canu_graph.bed')
ARS = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/ARS_graph.bed')

print(len(hifiasm),len(canu),len(ARS))

print(sum(hifiasm['edge']==canu['edge']))
print(sum(hifiasm['edge']==ARS['edge']))
print(sum(ARS['edge']==canu['edge']))

x = (hifiasm['edge']==canu['edge']) & (hifiasm['edge']!=ARS['edge'])
print(sum(x))

print(hifiasm[hifiasm['edge']==ARS['edge']].head())
print(ARS[ARS['link']=='alt'])

def extract_overlap_regions(reference,hifiasm,canu):
    BSW_divergences = (hifiasm['edge']==canu['edge']) & (hifiasm['edge']!=ARS['edge'])
    BSW_hifiasm = hifiasm[BSW_divergences]
    BSW_canu = canu[BSW_divergences]
    hifiasm_exclusive = hifiasm[(ARS['link']=='match') & (hifiasm['edge']!=ARS['edge']) & (hifiasm['edge']!=canu['edge'])]
    canu_exclusive = canu[(ARS['link']=='match') & (canu['edge']!=ARS['edge']) & (hifiasm['edge']!=canu['edge'])]
    return hifiasm_exclusive, BSW_hifiasm

def filter_end_regions(df,max_distance=1e6,chr_only=True):
    df_filt = df[(df['ref_contig']==df['map_contig']) & (abs(df['map_start']-df['ref_start']) < max_distance)]
    if chr_only:
        return df_filt[df_filt['chromosome']==True]
    else:
        return df_filt

bc = extract_overlap_regions(ARS,hifiasm,canu)[1]

print(filter_end_regions(bc))

print(extract_overlap_regions(ARS,hifiasm,canu)[1].dropna()[100:110])


variations = np.zeros((len(truth),2))

for f in ('canu','hifiasm','ARS'):
    df_ont = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/{f}_graph.bed')
    print(df_ont.head())
    matches = df_ont[df_ont['ref_contig']==df_ont['map_contig']]
    print(len(matches))
    print(matches.head(n=1))
    variations[1000][1]+=1
    for row in matches.itertuples():
        if row.link == 'match':
            variations[row.Index][0]+=1
        else:
            variations[row.Index][1]+=1




    print(Counter(matches['link']))
    #variations[df['source']].append()
    #print(len(matches))
    #print(np.mean(matches['ref_start']-matches['map_start']))
    #print(np.std(matches['ref_end']-matches['ref_start']-matches['size']))

print(np.sum(variations,axis=0))
plt.figure()
plt.plot(np.diff(variations,axis=1),lw=0,marker='o')
plt.show()
print(Counter([i[0] for i in np.diff(variations,axis=1).tolist()]))
print(np.where(np.diff(variations,axis=1)>3.5))
