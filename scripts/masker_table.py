import pandas
import glob
import argparse

def build_repeat_dataframe(dir_path):
    rows = []
    for f in glob.glob(f'{dir_path}/*.tbl'):
        rows.append(extract_repeats(f))
    return pandas.DataFrame(rows)

def extract_repeats(fname):
    elements = ('Total interspersed repeats', 'Small RNA', 'Satellites', 'Simple repeats', 'Low complexity')
    repeats = dict()

    with open(fname,'r') as file_:
        for line in file_:
            if 'file name' in line:
                repeats['chrm'] = line.split(':')[1].rstrip()[:-8]
            elif 'total length' in line:
                repeats['length'] = int(line.split('(')[1].split()[0])
            for element in elements:
                if element in line:
                    repeats[element] = float(line.split('bp')[1].split()[0])
    return repeats

def main(direct_input=None):
    parser = argparse.ArgumentParser(description='masker.')
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument('--assembler', type=str, required=True)
    parser.add_argument('--haplotype', type=str, required=True)

    args = parser.parse_args(direct_input)

    df = build_repeat_dataframe(f'{args.assembler}_{args.sample}/{args.haplotype}_split_chrm')
    df.to_csv(f'results/{args.haplotype}_{args.sample}_{args.assembler}.repeats.csv',index=False)

def centromere_locations():
    rows = [line.rstrip().split() for line in open('centr_locs.txt')]
    df = pd.DataFrame(rows,columns=['SW','Dv','db','q','name','start','stop','rem','com','repeat','rep_class','rom','rim','nom','ids'])
    g = df.groupby('name')
    print(g.mean()[g.count()['SW']>10][:30])

if __name__ == '__main__':
    main()
