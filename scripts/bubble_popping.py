import pandas
import numpy as np

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
                    row['link'] = 'deletion'
                else:
                    row['link'] = parts[5][2:].split(':')[0]
                row['size'] = int(parts[5].split(':')[1])
                row['map_contig'] = parts[5].split(':')[3].split('.')[0]
                row['map_start']= int(parts[5].split(':')[4])
                row['map_end'] = int(parts[5].split(':')[5])
            data.append(row)
    return pandas.DataFrame(data)

for f in ('ont','3','4'):
    df_ont = load_bed(f'/Users/alexleonard/Documents/Tiergenomik/DATA/sample_{f}.bed')
    matches = df_ont[df_ont['ref_contig']==df_ont['map_contig']]
    print(len(matches))
    print(np.mean(matches['ref_start']-matches['map_start']))
    print(np.std(matches['ref_end']-matches['ref_start']-matches['size']))
