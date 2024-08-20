import glob
import Colors
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr
from scipy.interpolate import interp1d
import itertools as it
import statsmodels.stats.multitest as multi
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
import sys
import numpy as np
import matplotlib

# Default sans-serif font is HELVETICA
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
#Always use sans-serif fonts
matplotlib.rcParams['font.family'] = "sans-serif"


path= '/Users/kaikretzschmar/Dropbox/_TomoSeq/'

np.bool = np.bool_
np.random.seed(174612)


#Read data
mdfDT = pd.read_csv(path + '/DT_zscore.tsv', sep = '\t', index_col = 0)
mdfHP = pd.read_csv(path + '/HP_zscore.tsv', sep = '\t', index_col = 0)

mndfDT = pd.read_csv(path + '/DT_norm.tsv', sep = '\t', index_col = 0)
mndfHP = pd.read_csv(path + '/HP_norm.tsv', sep = '\t', index_col = 0)

mdfsDT = [mdfDT[[idx for idx in mdfDT.columns if 'DT1' in idx]],
           mdfDT[[idx for idx in mdfDT.columns if 'DT2' in idx]]]

mdfsHP = [mdfHP[[idx for idx in mdfHP.columns if 'HP1' in idx]],
          mdfHP[[idx for idx in mdfHP.columns if 'HP2' in idx]]]

mndfsDT = [mndfDT[[idx for idx in mndfDT.columns if 'DT1' in idx]],
           mndfDT[[idx for idx in mndfDT.columns if 'DT2' in idx]]]

mndfsHP = [mndfHP[[idx for idx in mndfHP.columns if 'HP1' in idx]],
          mndfHP[[idx for idx in mndfHP.columns if 'HP2' in idx]]]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf

mdf_all = mergeDfs([mdfDT,mdfHP])

mndfDT = mergeDfs(mndfsDT)
mndfHP = mergeDfs(mndfsHP)

mndf_all = mergeDfs([mndfDT,mndfHP])

mdfHP = [mdf_all[[idx for idx in mdf_all.columns if 'HP1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'HP2' in idx]]]

mdfDT = [mdf_all[[idx for idx in mdf_all.columns if 'DT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'DT2' in idx]]]

mdfsHP = [mdf_all[[idx for idx in mdf_all.columns if 'HP1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'HP2' in idx]]]

mdfsDT = [mdf_all[[idx for idx in mdf_all.columns if 'DT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'DT2' in idx]]]

mndfsHP = [mndf_all[[idx for idx in mndf_all.columns if 'HP1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'HP2' in idx]]]

mndfsDT = [mndf_all[[idx for idx in mndf_all.columns if 'DT1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'DT2' in idx]]]

0.5*sum([df.sum().sum() for df in mndfsDT]), 0.5*sum([df.sum().sum() for df in mndfsHP])

mndfsDT = [26399688.0*df/df.sum().sum() for df in mndfsDT]
mndfsHP = [29401432.5*df/df.sum().sum() for df in mndfsHP]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfDT = mergeDfs(mdfDT)
mdfHP = mergeDfs(mdfHP)

mndfDT = mergeDfs(mndfsDT)
mndfHP = mergeDfs(mndfsHP)

mndf_all = mergeDfs([mndfDT,mndfHP])


def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

corrGenes_DT = pd.read_table(path + '/corr_genes_DT.tsv', sep = '\t', index_col=0)

pvth = 0.05

gl_DT_dict = {k: corrGenes_DT[k][corrGenes_DT['adj-pv'] < pvth] for k in corrGenes_DT}
genelist_DT = []
for k in gl_DT_dict:
    genelist_DT  += list(gl_DT_dict[k].index)

from collections import Counter

gl_DT_cnt = Counter(genelist_DT)
gl_DT = [g for g in gl_DT_cnt]

corrGenes_HP = pd.read_table(path + '/corr_genes_HP.tsv', sep = '\t', index_col=0)

pvth = 0.05

gl_HP_dict = {k: corrGenes_HP[k][corrGenes_HP['adj-pv'] < pvth] for k in corrGenes_HP}
genelist_HP = []
for k in gl_HP_dict:
    genelist_HP  += list(gl_HP_dict[k].index)


gl_HP_cnt = Counter(genelist_HP)
gl_HP = [g for g in gl_HP_cnt]


import venn

labels = venn.get_labels([set(gl_DT),
                          set(gl_HP)],
                         fill=['number', 'percent'])

for k in labels:
    if labels[k][0:1] == '0':
        labels[k] = ''

oralcolors = [
            # r, g, b, a
            [213, 62, 135, 0.4],
            [75, 155, 121, 0.4],
        ]

oralcolors = [
            [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
            for i in oralcolors
        ]

fig, ax = venn.venn2(labels, names=['DT', 'HP'], colors = oralcolors)
plt.savefig(path + '/venn_genes_p_val_0.05_DT_vs_HP.pdf', bbox_inches='tight')
plt.show()


def averageGastruloid(mzdf,n):
    mdf = mzdf.T
    glds = set([c.rsplit('-')[-1] for c in mdf.index])
    print(glds)
    mdf['x'] = 0
    for gld in glds:
        idxs = [c for c in mdf.index if gld == c.rsplit('-')[-1]]
        xra = np.linspace(0,1,len(idxs))
        for idx, x in zip(idxs, xra):
            mdf.loc[idx,'x'] = x
    mdf = mdf.sort_values(by = 'x', ascending = True)
    xdf = pd.DataFrame(columns = mdf.columns, index = range(n))
    xra = np.linspace(0,1,n+1)
    for i in range(len(xra)-1):
        idxs = [c for c in mdf.index if xra[i]<=mdf.loc[c,'x']<=xra[i+1]]
        xdf.loc[i] = mdf.loc[idxs].mean()
    del xdf['x']
    return xdf.T.fillna(0)


mdfDT_av = averageGastruloid(mdfDT, 19).astype(float)
mnfDT_av = averageGastruloid(mndfDT, 19).astype(float)

mdfHP_av = averageGastruloid(mdfHP, 19).astype(float)
mnfHP_av = averageGastruloid(mndfHP, 19).astype(float)


#Genes to be annotated
g_annotate = ['Sprr2b','Myog','Myh4','Krt81','Krt36','Lft','Itga8','Cxcl2','Krt10',
              'Cnn1','Vegfb','Runx1','Fabp4','Myh11','Krt76']

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

from neupy import algorithms

def SOMclustering(df, gridheight = 10, gridwidth = 10, seed = 54321):
    np.random.seed = seed
    GRID_HEIGHT = gridheight
    GRID_WIDTH = gridwidth
    somfun = algorithms.SOFM(
            n_inputs=len(df.columns),
            features_grid=(GRID_HEIGHT, GRID_WIDTH),
            learning_radius=4,
            reduce_radius_after=50,
            step=0.5,
            std=1,
            shuffle_data=True,
            verbose=False)
    somfun.train(np.array(df), epochs = 300)
    clustersr = somfun.predict(np.array(df)).argmax(axis=1)
    cl_genelist_df = pd.Series(clustersr, index = df.index)
    return cl_genelist_df

def hierarchicalClustering(df, cth = 500, plot = False, method = 'ward', metric = 'euclidean', nolabels = 'True'):
    Z = linkage(df, method=method, metric = metric)
    dg = dendrogram(Z, no_labels=nolabels, color_threshold=cth, no_plot = np.invert(plot))
    plt.show()
    return Z, dg

def getClusterByColor(dg, labels):
    kk = []
    ii = 0
    cluster = 0
    color = dg['color_list'][0]
    clusters = {cluster: []}
    for i in range(len(dg['icoord'])):
        v = dg['icoord'][i]
        for j in [0,2]:
            vj = int(round((v[j]-5.)/10))
            if (v[j]-5.)/10 == vj and vj not in kk:
                kk.append(vj)
                if dg['color_list'][i] == color:
                    clusters[cluster].append(labels[dg['leaves'][vj]])
                else:
                    color = dg['color_list'][i]
                    cluster += 1
                    clusters[cluster] = [labels[dg['leaves'][vj]]]
    return clusters

def comboCluster(df, genes, gridheight = 10, gridwidth = 10, cth = 500, plot = True, method = 'ward', metric = 'euclidean', seed = 12345, nolabels = 'True'):
    somcl = SOMclustering(df.loc[genes], gridheight = gridheight, gridwidth = gridwidth, seed = seed)
    mdf_som = pd.DataFrame({i: df.loc[somcl[somcl==i].index].mean() for i in set(somcl) })
    Z, dg = hierarchicalClustering(mdf_som.T, cth = cth, plot = plot, method = method, metric = metric, nolabels = nolabels)
    if plot:
        plt.show()
    ghcl = getClusterByColor(dg, list(mdf_som.columns))
    cl_df = pd.DataFrame(somcl, columns = ['som'])
    cl_df['hcl'] = [k for s in cl_df['som'] for k in ghcl.keys() if s in ghcl[k]]
    return cl_df


#Extract shared reproducible genes
cnt = Counter(list(gl_DT) + list(gl_HP))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)
#len(gl_both), 5*np.sqrt(len(gl_both))

gl_combined = list(set(list(gl_DT) + list(gl_HP)))
5*np.sqrt(len(gl_combined))

#Shared clustering
somhcl_all_df = comboCluster(mdf_all, gl_combined, cth=10)


plt.figure(figsize=(10,10))

hcldf = somhcl_all_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

Npl = len(mndfsDT + mndfsHP)+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(hcldf['hcl']):
    y = 1.0*len(hcldf[hcldf['hcl'] == i])/len(hcldf)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(mdfsDT+mdfsHP):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(', '.join(df.columns[0].rsplit('-')[1:]))
jidx = 4
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)
plt.yticks([i for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/reprGenes_DT_HP_replicates.pdf', bbox_inches = 'tight')
hcldf.sort_values(by=['hcl','som']).to_csv(path + '/reprGenes_DT_HP.tsv', sep = '\t')
plt.show()

plt.figure(figsize=(8,10))

gene_order = somhcl_all_df.sort_values(by=['hcl','som']).index
Npl = 2+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 3.5
for i in set(somhcl_all_df['hcl']):
    y = 1.0*len(somhcl_all_df[somhcl_all_df['hcl'] == i])/len(somhcl_all_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
xb = 4; y0 = 0; y = 1.0/len(somhcl_all_df)
for g in gene_order:
    if g in gl_DT and g in gl_HP:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_DT and g not in gl_HP:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#E7298A')
    if g not in gl_DT and g in gl_HP:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#1B9E77')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['DT average','HP average']
for i, df in enumerate([mdfDT_av] + [mdfHP_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_DT_HP_averages.pdf', bbox_inches = 'tight')
plt.show()


n10nfdfsDT = 1e5*mnfDT_av/mnfDT_av.sum().sum()
n10nfdfsHP = 1e5*mnfHP_av/mnfHP_av.sum().sum()


dex_DTHP = pd.DataFrame({
    'total-DT': n10nfdfsDT.sum(axis=1),
    'total-HP': n10nfdfsHP.sum(axis=1),
    'total-mean': {idx: np.mean(list(n10nfdfsDT.loc[idx].values) + list(n10nfdfsHP.loc[idx].values)) for idx in set(list(n10nfdfsDT.index)+list(n10nfdfsHP.index)) if idx in n10nfdfsDT.index and idx in n10nfdfsHP.index}
})

dex_DTHP = dex_DTHP.loc[[idx for idx in dex_DTHP.index if idx in mndfDT.index and idx in mndfHP.index]]


dex_DTHP['log2FC'] = np.log2((dex_DTHP['total-DT']+1e-10)/(dex_DTHP['total-HP']+1e-10))
dex_DTHP['log2mu'] = np.log2(0.5*(dex_DTHP['total-DT']+dex_DTHP['total-HP'])) # np.log2(dex_Chi['total-mean'])

dex_DTHP['prob-DT'] = dex_DTHP['total-DT']/dex_DTHP['total-DT'].sum()
dex_DTHP['prob-HP'] = dex_DTHP['total-HP']/dex_DTHP['total-HP'].sum()

dex_DTHP['N-DT'] = dex_DTHP['total-DT'].sum()
dex_DTHP['N-HP'] = dex_DTHP['total-HP'].sum()

fdex_DTHP = dex_DTHP.loc[[g for g in dex_DTHP.index if g in gl_combined]]

from scipy.stats import binom

def dex_Pval(expression, total, prob):
    expression, total = int(expression), int(total)
    pv = binom.cdf(expression, total, prob)
    if pv > 0.5:
        pv = 1-pv
    return pv

fdex_DTHP['PvalDT'] = fdex_DTHP.apply(lambda x: dex_Pval(x['total-DT'],x['N-DT'],x['prob-DT']), axis = 1)
fdex_DTHP['PvalHP'] = fdex_DTHP.apply(lambda x: dex_Pval(x['total-HP'],x['N-HP'],x['prob-HP']), axis = 1)

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_DTHP['log2mu'], fdex_DTHP['log2FC'], s = 5, c = np.log10(fdex_DTHP['PvalDT']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_DTHP[np.log10(fdex_DTHP['PvalDT'])<-2.25].index:
    x = fdex_DTHP.loc[g,'log2mu']
    y = fdex_DTHP.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('log2(mean)')
plt.ylabel('log2(FC)')

fdex_DTHP.to_csv(path + '/dex_DTHP_sameNormalization.tsv', sep = '\t')
plt.savefig(path + '/dex_DTHP_sameNormalization.pdf', bbox_inches = 'tight')

plt.show()

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_DTHP['total-mean'], fdex_DTHP['log2FC'], s = 5, c = np.log10(fdex_DTHP['PvalDT']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_DTHP[np.log10(fdex_DTHP['PvalDT'])<-3].index:
    x = fdex_DTHP.loc[g,'total-mean']
    y = fdex_DTHP.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('mean')
plt.ylabel('log2(FC)')

plt.savefig(path + '/dex_DTHP_mean.pdf', bbox_inches = 'tight')

plt.show()


fdex_DTHP.loc[[g for g in fdex_DTHP[np.log10(fdex_DTHP['PvalDT'])<-2].index]]

[idx for idx in mnfDT_av.index if 'Mog' in idx]

g = 'ENSMUSG00000023084_Lrrc71_ProteinCoding__3'

plt.figure(figsize=(6*1.6,3))
plt.subplot(121)
plt.plot(np.linspace(0,1,len(mndfsDT[0].columns)), mndfsDT[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsDT[1].columns)), mndfsDT[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[g])
plt.legend(['DT1','DT2','DT_av'])

plt.subplot(122)
plt.plot(np.linspace(0,1,len(mndfsHP[0].columns)), mndfsHP[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsHP[1].columns)), mndfsHP[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[g])
plt.legend(['HP1','HP2','HP_av'])

plt.show()


#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'
#oral tissue 50% brighter colour BM, HP, VT, DT: '#fdac6e','#74e7c5','#bab7d9','#f294c4'


xgs = ['ENSMUSG00000026208_Des_ProteinCoding__1','ENSMUSG00000050092_Sprr2b_ProteinCoding__3',
       'ENSMUSG00000024962_Vegfb_ProteinCoding__19','ENSMUSG00000062515_Fabp4_ProteinCoding__3',
       'ENSMUSG00000020913_Krt24_ProteinCoding__11','ENSMUSG00000019761_Krt10_ProteinCoding__11',
       'ENSMUSG00000032068_Plet1_ProteinCoding__9','ENSMUSG00000075402_Krt76_ProteinCoding__15',
       'ENSMUSG00000022952_Runx1_ProteinCoding__16','ENSMUSG00000005836_Gata6_ProteinCoding__18',
       'ENSMUSG00000018830_Myh11_ProteinCoding__16']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_HP_shared_genes_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000057723_Krt33b_ProteinCoding__11','ENSMUSG00000067615_Krt81_ProteinCoding__15',
       'ENSMUSG00000026459_Myog_ProteinCoding__1','ENSMUSG00000020916_Krt36_ProteinCoding__11',
       'ENSMUSG00000048013_Krt35_ProteinCoding__11','ENSMUSG00000057003_Myh4_ProteinCoding__11',
       'ENSMUSG00000050092_Sprr2b_ProteinCoding__3']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_HP_unique_genes_DT_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

xgs = ['ENSMUSG00000030017_Reg3g_ProteinCoding__6','ENSMUSG00000058427_Cxcl2_ProteinCoding__5',
       'ENSMUSG00000001349_Cnn1_ProteinCoding__9','ENSMUSG00000026768_Itga8_ProteinCoding__2',
       'ENSMUSG00000035638_Muc20_ProteinCoding__16','ENSMUSG00000032496_Ltf_ProteinCoding__9']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_HP_unique_genes_HP_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()