import pandas
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def plot_chromosomes(primary_chrm='CM',unplaced_chrm='VLPJ'):
    df = pandas.read_csv('BSWCHEF120152514636_100_hifiasm.chrm.coverage.txt',delimiter='\t')
    df['chromosome'] = ['C' if primary_chrm in ch else ('U' if unplaced_chrm in ch else 'P') for ch in df['#rname']]

    sns.scatterplot(data=df,x='meandepth',y='meanmapq',hue='chromosome')
    plt.show(block=False)

    df2 = pandas.read_csv('BSWCHEF120152514636_100_hifiasm.windows.coverage.txt',delimiter='\t',header=0,names=['chrm','start','stop','depth'])
    df2['avg_depth'] = df2.apply(lambda x: x['depth']/(x['stop']-x['start']), axis=1)
    dx = df2[[primary_chrm in f for f in df2['chrm']]]
    g = sns.relplot(data=dx,x='start',y='avg_depth',col='chrm',col_wrap=6,kind='line',facet_kws={'sharex':False})
    plt.show(block=False)
    mean_stds = dx.groupby('chrm').agg({'avg_depth':['mean','std']})
    print(mean_stds)

    #sns.lmplot(data=dx,x='start',y='avg_depth',col='chrm',col_wrap=6,scatter_kws={'alpha':.001},line_kws={'color':''})
    plt.show(block=False)

def load_chrom_coverage():
    pass

def load_window_coverage(run,read,primary_chrm=None):
    df = pandas.read_csv(f'{run}.windows.{read}.coverage.txt',delimiter='\t',header=0,names=['chrm','start','stop','depth'])
    df['avg_depth'] = df.apply(lambda x: x['depth']/(x['stop']-x['start']), axis=1)

    if primary_chrm:
        return df[[primary_chrm in f for f in df['chrm']]]
    return df

plot_chromosomes()
