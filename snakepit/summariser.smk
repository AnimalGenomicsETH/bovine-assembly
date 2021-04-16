
VALID_CHR = [*range(1,30),'X','Y','MT']
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


rule extract_centr_telo_features:
    input:
        get_dir('work','{haplotype}_split_chrm/{chr}',chr=CHR) for CHR in VALID_CHR)
    output:
        temp(get_dir('work','{haplotype}.centr_telo.unmerged.bed'))
    params:
        regex = r'BTSAT|\(TTAGGG\)n'
    run:
        import re
        colours = {'BTSAT':'186,56,44','(TTAGGG)n':'44,174,186'}
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
                            fout.write(f'{parts[4]},{int(parts[5])-1},{parts[6]},{parts[9]},.,{"-" if parts[8] =="C" else "+"},{int(parts[5])-1},{parts[6]},{colours[parts[9]]}\n')

rule bedtools_merge_centr_telo:
    input:
        get_dir('work','{haplotype}.centr_telo.unmerged.bed')
    output:
        get_dir('result','.centr_telo.bed')
    shell:
        'bedtools merge -i {input} -c 4,5,6,7,8,9 -o distinct,first,first,min,max,first > {output}'



#(TTAGGG)n
#track name="ItemRGBDemo" description="Item RGB demonstration" visibility=2 itemRgb="On"
#chr7    127471196  127472363  Pos1  0  +  127471196  127472363  255,0,0
#awk 'BEGIN{OFS="\t"}{if(NR>3 && $10 ~ /BTSAT/) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' hap1_split_chrm/2_RagTag.chrm.fa.out > chr2.bed
#bedtools merge -i chr2.bed -c 4 -o distinct > chr2.merge.bed
