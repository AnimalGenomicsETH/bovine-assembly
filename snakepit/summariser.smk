

def is_chromosomal(name):
    valid_chromosomes = [*map(str,range(1,30)),'X','Y','MT']
    return name.replace('_RagTag','') in valid_chromosomes

def length_against_reference(df):
    lens = []
    for _,row in df.iterrows():
        lens.append(int(row['length']-df[(df['assembly']=='Ref')&(df['chromosome']==row['chromosome'])]['length']))
    return lens

rule extract_chromosome_lengths:
    input:
        (get_dir('work','{haplotype}.scaffolds.fasta.fai',sample=100,assembler=ASM) for ASM in config['assemblers']),
        config['ref_genome'] + '.fai'
    output:
        get_dir('summary','{haplotype}.chromosome_lengths.csv')
    run:
        with open(output[0],'w') as fout:
            fout.write('assembly,chromosome,length\n')
            for fai, assembler in zip(input,config['assemblers']+['Ref',]):
                with open(fai,'r') as fin:
                    for line in fin:
                        parts = line.rstrip().split()
                        if is_chromosomal(parts[0]):
                            fout.write(f'{assembler},{parts[0].rstrip("_RagTag")},{parts[1]}\n')

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
        f, axes = plt.subplots(2,1)
        sns.barplot(data=df,x='chromosome',y='length',hue='assembly',ax=axes[0])
        df['residues'] = length_against_reference(df)
        df_no_ref = df[df['assembly']!='Ref']
        sns.barplot(data=df_no_ref,x='chromosome',y='residues',hue='assembly',ax=axes[1])
        plt.savefig(output[0])
