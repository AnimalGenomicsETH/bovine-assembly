import pandas
import glob

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

if __name__ == '__main__':
    df = build_repeat_dataframe('../third_test/split_BSWCHEF120152514636_100_hifiasm')
    df.to_csv('table.csv',index=False)
