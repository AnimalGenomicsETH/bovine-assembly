

def is_chromosomal(name):
    valid_chromosomes = [*range(1,30),'X','Y','MT']
    return name in valid_chromosomes

rule extract_chromosome_lengths:
    input:
        (get_dir('work','{haplotype}.scaffolds.fasta.fai',sample=100,assembler=ASM) for ASM in config['assemblers'])
    output:
        get_dir('summary','{haplotype}.chromosome_lengths.csv')
    run:
        with open(output[0],'w') as fout:
            for fai, assembler in zip(input,config['assemblers']):
                with open(fai,'r') as fin:
                    for line in fin:
                        parts = line.rstrip().split()
                        if is_chromosomal(parts[0]):
                            fout.write(f'{asm},{parts[1]},{parts[2]}\n')

rule plot_chromosome_lengths:
    input:
        get_dir('summary','{haplotype}.chromosome_lengths.csv')
    output:
        get_dir('summary','{haplotype}.chromosome_lengths.png')
    run:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd

        df =pd.read_csv(input[0])
        #df2= df[~df['chromosome'].str.contains("NKL")]
        #df2['chromosome'] = [i.strip('_RagTag') for i in df2['chromosome']]
        sns.catplot(data=df,kind='bar',x='chromosome',y='length',hue='assembly')
        plt.savefig(output[0])
