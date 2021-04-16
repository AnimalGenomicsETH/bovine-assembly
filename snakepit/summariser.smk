localrules: extract_chromosome_lengths, plot_chromosome_lengths, extract_centr_telo_features, bedtools_merge_centr_telo

VALID_CHR = [*map(str,range(1,30)),'X','Y','MT']
PRIMARY_CHR = VALID_CHR[:-2]

def is_chromosomal(name):
    #TODO remove ragtag
    return name.replace('_RagTag','') in VALID_CHR

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


rule extract_centr_telo_features:
    input:
        (get_dir('work','{haplotype}_split_chrm/{chr}.chrm.fa.out',chr=CHR) for CHR in PRIMARY_CHR)
    output:
        temp(get_dir('work','{haplotype}.centr_telo.unmerged.bed'))
    params:
        regex = r'BTSAT|OSSAT|\(TTAGGG\)n'
    run:
        import re
        def get_colour(ele):
            if ele == '(TTAGGG)n':
                return '44,174,186'
            elif 'BTSAT' in ele:
                return '186,56,44'
            else: #nominally OSSAT
                return '140,50,200'
            
        pattern = re.compile(params.regex)
        with open(output[0],'w') as fout:
            fout.write('track name="centro telo" visibility=2 itemRgb="On"\n')
            for chr in input:
                with open(chr,'r') as fin:
                    for i, line in enumerate(fin):
                        if i < 3:
                            continue
                        parts = line.rstrip().split()
                        if pattern.match(parts[9]):
                            fout.write('\t'.join(map(str,(parts[4],int(parts[5])-1,parts[6],parts[9],'.','-' if parts[8] == 'C' else '+',int(parts[5])-1,parts[6],get_colour(parts[9]))))+'\n')

rule bedtools_merge_centr_telo:
    input:
        get_dir('work','{haplotype}.centr_telo.unmerged.bed')
    output:
        get_dir('result','.centr_telo.bed')
    shell:
        'bedtools merge -i {input} -header -c 4,5,6,7,8,9 -o distinct,first,first,min,max,first > {output}'

