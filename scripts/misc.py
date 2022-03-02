df= pd.read_csv('hap2.tbl')
df.melt(id_vars=["chrm", "length
 mllt['chr'] = ['RagTag' in c for c in mllt['chrm']]
 srt = mllt.sort_values(['chrm']).reset_index(drop=True)
sns.scatterplot(data=srt,x='chrm',y='value',hue='variable',style='chr')


f = pd.read_csv('../DATA/summary.df',delimiter='\t')
df_no_missing = df[(df['truth']!='MISSING')&(df['call']!='MISSING')]
g= sns.catplot(data=df_no_missing,x='call',y='count',col='truth',hue='caller',sharey=False)
for ax in g.axes.flatten():
     ...:     ax.set_yscale('log')
df2['conc']=df2['truth']==df2['call']
df2.groupby(['caller','conc']).sum()


counts = df_no_missing.groupby(['chip','caller','conc']).sum()
dy = counts.unstack(level=2).reset_index()
dy.columns=['chip','caller','conc','disconc']
dy['gc']=dy['disconc']/dy['conc']




import numpy as np
pLIs = np.array([float(i.split()[19]) for i in open('Documents/Tiergenomik/bovine-assembly/random_interval_shuffle/pLI.txt')])
num_samples = 100000
idx = np.random.randint(0,len(b),size=(num_samples,295))
samples = pLIs[idx]
x=np.sum(samples>0.9,axis=1)
plt.hist(x,bins=np.linspace(0,100,101))
