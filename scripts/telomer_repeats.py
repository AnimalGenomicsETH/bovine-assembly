import re
import screed

from scipy.stats import binom

def main():
    region = 1000
    telomeric_repeat_for = re.compile("TTAGGG", re.IGNORECASE)
    telomeric_repeat_back = re.compile("GGGATT", re.IGNORECASE)

    with open(snakemake.output[0],'w') as fout:
        fout.write('name\tcount_for\tp_for\tcount_back\tp_back')

        for seq in screed.open(snakemake.input[0]):
            repeats_for_s = len(telomeric_repeat_for.findall(seq.sequence[:region]))
            repeats_for_e = len(telomeric_repeat_for.findall(seq.sequence[-1 * region:]))
            repeats_back_s = len(telomeric_repeat_back.findall(seq.sequence[:region]))
            repeats_back_e len(telomeric_repeat_back.findall(seq.sequence[-1 * region:])

            fout.write(f'{seq.name}' + \
                       '\t'.join(f'{repeats}\t{binom.sf(repeats,region,0.25**6):.4f}' for repeat in (repeats_for_s,repeats_for_e,repeats_back_s,repeats_back_s))
                       )
