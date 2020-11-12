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

def filter_chromosomes(df,primary_chrm=None):
    return df[[primary_chrm in f for f in df['chrm']]] if primary_chrm else df

def load_chrm_coverage(run,read,primary_chrm=None):
    df = pandas.read_csv(f'{run}.chrm.{read}.coverage.txt',delimiter='\t',header=0,names=['chrm','start','stop','reads','bases','coverage','depth','baseQ','mapQ'])
    return filter_chromosomes(df,primary_chrm)

def load_window_coverage(run,read,primary_chrm=None):
    df = pandas.read_csv(f'{run}.windows.{read}.coverage.txt',delimiter='\t',header=0,names=['chrm','start','stop','depth'])
    df['avg_depth'] = df.apply(lambda x: x['depth']/(x['stop']-x['start']), axis=1)
    return filter_chromosomes(df,primary_chrm)

def plot_chromosomes_new(fname,read_t,primary_chrm):
    df = load_window_coverage(fname,read_t,primary_chrm)
    g = sns.relplot(data=df,x='start',y='avg_depth',col='chrm',col_wrap=6,kind='line',facet_kws={'sharex':False})
    g.map_dataframe(add_std,'start','avg_depth')
    plt.show(block=False)

def add_std(x,y,data,**kwargs):
    mean = np.mean(data['y'])
    std = np.std(data['y'])
    ax = plt.gca()
    ax.axhline(mean,c='r')
    print(f'mean {mean:.2f}, std {std:.2f}')
    ax.set_ylim([0,100])
    for sign in (1,-1):
        ax.axhline(mean+sign*std,c='k',ls='--')

plot_chromosomes_new('../data/asm_100_hifiasm','hifi','NC')


plot_chromosomes_new('../data/asm_100_hifiasm','minimap2','NC')
