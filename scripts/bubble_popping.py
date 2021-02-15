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
                    link = int(parts[5][2:].split(':')[0])
                    if link > row['source'] and link < row['sink']:
                        row['link'] = 'match'
                    else:
                        row['link'] = 'alt'
                row['size'] = int(parts[5].split(':')[1])
                row['map_contig'] = parts[5].split(':')[3].split('.')[0]
                row['map_start']= int(parts[5].split(':')[4])
                row['map_end'] = int(parts[5].split(':')[5])
            data.append(row)
    return pandas.DataFrame(data)


truth = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/sample_hifi.bed')
variations = np.zeros((len(truth),2))

for f in ('3','4','sire','dam'):
    df_ont = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/sample_{f}.bed')
    print(len(df_ont))
    matches = df_ont[df_ont['ref_contig']==df_ont['map_contig']]
    print(matches.head())
    print(len(matches))
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
