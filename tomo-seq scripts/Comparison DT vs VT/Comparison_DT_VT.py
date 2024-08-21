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


path= '/Documents/TomoSeq_project/Tomo/'

np.bool = np.bool_
np.random.seed(174612)


#Read data
mdfDT = pd.read_csv(path + '/DT_zscore.tsv', sep = '\t', index_col = 0)
mdfVT = pd.read_csv(path + '/VT_zscore.tsv', sep = '\t', index_col = 0)

mndfDT = pd.read_csv(path + '/DT_norm.tsv', sep = '\t', index_col = 0)
mndfVT = pd.read_csv(path + '/VT_norm.tsv', sep = '\t', index_col = 0)

mdfsDT = [mdfDT[[idx for idx in mdfDT.columns if 'DT1' in idx]],
           mdfDT[[idx for idx in mdfDT.columns if 'DT2' in idx]]]

mdfsVT = [mdfVT[[idx for idx in mdfVT.columns if 'VT1' in idx]],
          mdfVT[[idx for idx in mdfVT.columns if 'VT2' in idx]]]

mndfsDT = [mndfDT[[idx for idx in mndfDT.columns if 'DT1' in idx]],
           mndfDT[[idx for idx in mndfDT.columns if 'DT2' in idx]]]

mndfsVT = [mndfVT[[idx for idx in mndfVT.columns if 'VT1' in idx]],
          mndfVT[[idx for idx in mndfVT.columns if 'VT2' in idx]]]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf

mdf_all = mergeDfs([mdfDT,mdfVT])

mndfDT = mergeDfs(mndfsDT)
mndfVT = mergeDfs(mndfsVT)

mndf_all = mergeDfs([mndfDT,mndfVT])

mdfVT = [mdf_all[[idx for idx in mdf_all.columns if 'VT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'VT2' in idx]]]

mdfDT = [mdf_all[[idx for idx in mdf_all.columns if 'DT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'DT2' in idx]]]

mdfsVT = [mdf_all[[idx for idx in mdf_all.columns if 'VT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'VT2' in idx]]]

mdfsDT = [mdf_all[[idx for idx in mdf_all.columns if 'DT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'DT2' in idx]]]

mndfsVT = [mndf_all[[idx for idx in mndf_all.columns if 'VT1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'VT2' in idx]]]

mndfsDT = [mndf_all[[idx for idx in mndf_all.columns if 'DT1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'DT2' in idx]]]

0.5*sum([df.sum().sum() for df in mndfsDT]), 0.5*sum([df.sum().sum() for df in mndfsVT])

mndfsDT = [26399688.0*df/df.sum().sum() for df in mndfsDT]
mndfsVT = [22424833.499999993*df/df.sum().sum() for df in mndfsVT]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfDT = mergeDfs(mdfDT)
mdfVT = mergeDfs(mdfVT)

mndfDT = mergeDfs(mndfsDT)
mndfVT = mergeDfs(mndfsVT)

mndf_all = mergeDfs([mndfDT,mndfVT])


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

corrGenes_VT = pd.read_table(path + '/corr_genes_VT.tsv', sep = '\t', index_col=0)

pvth = 0.05

gl_VT_dict = {k: corrGenes_VT[k][corrGenes_VT['adj-pv'] < pvth] for k in corrGenes_VT}
genelist_VT = []
for k in gl_VT_dict:
    genelist_VT  += list(gl_VT_dict[k].index)


gl_VT_cnt = Counter(genelist_VT)
gl_VT = [g for g in gl_VT_cnt]


import venn

labels = venn.get_labels([set(gl_DT),
                          set(gl_VT)],
                         fill=['number', 'percent'])

for k in labels:
    if labels[k][0:1] == '0':
        labels[k] = ''

oralcolors = [
            # r, g, b, a
            [213, 62, 135, 0.4],
            [116, 112, 174, 0.4],
        ]

oralcolors = [
            [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
            for i in oralcolors
        ]

fig, ax = venn.venn2(labels, names=['DT', 'VT'], colors = oralcolors)
plt.savefig(path + '/venn_genes_p_val_0.05_DT_vs_VT.pdf', bbox_inches='tight')
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

mdfVT_av = averageGastruloid(mdfVT, 20).astype(float)
mnfVT_av = averageGastruloid(mndfVT, 20).astype(float)


#Genes to be annotated
#g_annotate = ['Sprr2b','Spink6','Krt33a','Bambi','Dkkl1','Sox21','Grem1','Fgf10','Fgf7','Fgf1']
g_annotate = ['Sprr2b','Krt34','Krt86','Sox21','Grem1','Cdh4','Lrig1','Notch4','Krt4','Nog']

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
cnt = Counter(list(gl_DT) + list(gl_VT))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)
#len(gl_both), 5*np.sqrt(len(gl_both))

gl_combined = list(set(list(gl_DT) + list(gl_VT)))

5*np.sqrt(len(gl_combined))

gl_VT_only = [g for g in gl_VT if g not in gl_DT]

5*np.sqrt(len(gl_VT_only))

#Shared clustering
somhcl_all_df = comboCluster(mdf_all, gl_combined, cth=10)


plt.figure(figsize=(10,10))

hcldf = somhcl_all_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

Npl = len(mndfsDT + mndfsVT)+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(hcldf['hcl']):
    y = 1.0*len(hcldf[hcldf['hcl'] == i])/len(hcldf)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(mdfsDT+mdfsVT):
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

plt.savefig(path + '/reprGenes_DT_VT_replicates.pdf', bbox_inches = 'tight')
hcldf.sort_values(by=['hcl','som']).to_csv(path + '/reprGenes_DT_VT.tsv', sep = '\t')
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
    if g in gl_DT and g in gl_VT:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_DT and g not in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#E7298A')
    if g not in gl_DT and g in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#7570B3')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['DT average','VT average']
for i, df in enumerate([mdfDT_av] + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_DT_VT_averages.pdf', bbox_inches = 'tight')
plt.show()


gl_VT_only = [g for g in gl_VT if g not in gl_DT]

5*np.sqrt(len(gl_VT_only))

#Shared clustering
somhcl_VT_df = comboCluster(mdf_all, gl_VT_only, cth=11)


plt.figure(figsize=(10,10))

hcldf = somhcl_VT_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

plt.figure(figsize=(8,10))

gene_order = somhcl_VT_df.sort_values(by=['hcl','som']).index
Npl = 2+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 3.5
for i in set(somhcl_VT_df['hcl']):
    y = 1.0*len(somhcl_VT_df[somhcl_VT_df['hcl'] == i])/len(somhcl_VT_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
xb = 4; y0 = 0; y = 1.0/len(somhcl_VT_df)
for g in gene_order:
    if g in gl_DT and g in gl_VT:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_DT and g not in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#E7298A')
    if g not in gl_DT and g in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#7570B3')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['DT average','VT average']
for i, df in enumerate([mdfDT_av] + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_DT_VT_VT_only_averages.pdf', bbox_inches = 'tight')
plt.show()



gl_DT_only = [g for g in gl_DT if g not in gl_VT]

5*np.sqrt(len(gl_DT_only))

#Shared clustering
somhcl_DT_df = comboCluster(mdf_all, gl_DT_only, cth=11)

plt.figure(figsize=(10,10))

hcldf = somhcl_DT_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

plt.figure(figsize=(8,10))

gene_order = somhcl_DT_df.sort_values(by=['hcl','som']).index
Npl = 2+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 3.5
for i in set(somhcl_DT_df['hcl']):
    y = 1.0*len(somhcl_DT_df[somhcl_DT_df['hcl'] == i])/len(somhcl_DT_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
xb = 4; y0 = 0; y = 1.0/len(somhcl_DT_df)
for g in gene_order:
    if g in gl_DT and g in gl_VT:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_DT and g not in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#E7298A')
    if g not in gl_DT and g in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#7570B3')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['DT average','VT average']
for i, df in enumerate([mdfDT_av] + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_DT_VT_DT_only_averages.pdf', bbox_inches = 'tight')
plt.show()




5*np.sqrt(len(gl_both))

#Shared clustering
somhcl_both_df = comboCluster(mdf_all, gl_both, cth=11)

plt.figure(figsize=(10,10))

hcldf = somhcl_both_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

plt.figure(figsize=(8,10))

gene_order = somhcl_both_df.sort_values(by=['hcl','som']).index
Npl = 2+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 3.5
for i in set(somhcl_both_df['hcl']):
    y = 1.0*len(somhcl_both_df[somhcl_both_df['hcl'] == i])/len(somhcl_both_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
xb = 4; y0 = 0; y = 1.0/len(somhcl_both_df)
for g in gene_order:
    if g in gl_DT and g in gl_VT:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_DT and g not in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#E7298A')
    if g not in gl_DT and g in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#7570B3')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['DT average','VT average']
for i, df in enumerate([mdfDT_av] + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_DT_VT_both_averages.pdf', bbox_inches = 'tight')
plt.show()



n10nfdfsDT = 1e5*mnfDT_av/mnfDT_av.sum().sum()
n10nfdfsVT = 1e5*mnfVT_av/mnfVT_av.sum().sum()


dex_DTVT = pd.DataFrame({
    'total-DT': n10nfdfsDT.sum(axis=1),
    'total-VT': n10nfdfsVT.sum(axis=1),
    'total-mean': {idx: np.mean(list(n10nfdfsDT.loc[idx].values) + list(n10nfdfsVT.loc[idx].values)) for idx in set(list(n10nfdfsDT.index)+list(n10nfdfsVT.index)) if idx in n10nfdfsDT.index and idx in n10nfdfsVT.index}
})

dex_DTVT = dex_DTVT.loc[[idx for idx in dex_DTVT.index if idx in mndfDT.index and idx in mndfVT.index]]


dex_DTVT['log2FC'] = np.log2((dex_DTVT['total-DT']+1e-10)/(dex_DTVT['total-VT']+1e-10))
dex_DTVT['log2mu'] = np.log2(0.5*(dex_DTVT['total-DT']+dex_DTVT['total-VT'])) # np.log2(dex_Chi['total-mean'])

dex_DTVT['prob-DT'] = dex_DTVT['total-DT']/dex_DTVT['total-DT'].sum()
dex_DTVT['prob-VT'] = dex_DTVT['total-VT']/dex_DTVT['total-VT'].sum()

dex_DTVT['N-DT'] = dex_DTVT['total-DT'].sum()
dex_DTVT['N-VT'] = dex_DTVT['total-VT'].sum()

fdex_DTVT = dex_DTVT.loc[[g for g in dex_DTVT.index if g in gl_combined]]

from scipy.stats import binom

def dex_Pval(expression, total, prob):
    expression, total = int(expression), int(total)
    pv = binom.cdf(expression, total, prob)
    if pv > 0.5:
        pv = 1-pv
    return pv

fdex_DTVT['PvalDT'] = fdex_DTVT.apply(lambda x: dex_Pval(x['total-DT'],x['N-DT'],x['prob-DT']), axis = 1)
fdex_DTVT['PvalVT'] = fdex_DTVT.apply(lambda x: dex_Pval(x['total-VT'],x['N-VT'],x['prob-VT']), axis = 1)

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_DTVT['log2mu'], fdex_DTVT['log2FC'], s = 5, c = np.log10(fdex_DTVT['PvalDT']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_DTVT[np.log10(fdex_DTVT['PvalDT'])<-2.25].index:
    x = fdex_DTVT.loc[g,'log2mu']
    y = fdex_DTVT.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('log2(mean)')
plt.ylabel('log2(FC)')

fdex_DTVT.to_csv(path + '/dex_DTVT_sameNormalization.tsv', sep = '\t')
plt.savefig(path + '/dex_DTVT_sameNormalization.pdf', bbox_inches = 'tight')

plt.show()

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_DTVT['total-mean'], fdex_DTVT['log2FC'], s = 5, c = np.log10(fdex_DTVT['PvalDT']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_DTVT[np.log10(fdex_DTVT['PvalDT'])<-3].index:
    x = fdex_DTVT.loc[g,'total-mean']
    y = fdex_DTVT.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('mean')
plt.ylabel('log2(FC)')

plt.savefig(path + '/dex_DTVT_mean.pdf', bbox_inches = 'tight')

plt.show()


fdex_DTVT.loc[[g for g in fdex_DTVT[np.log10(fdex_DTVT['PvalDT'])<-2].index]]

[idx for idx in mnfDT_av.index if 'Mog' in idx]

g = 'ENSMUSG00000023084_Lrrc71_ProteinCoding__3'

plt.figure(figsize=(6*1.6,3))
plt.subplot(121)
plt.plot(np.linspace(0,1,len(mndfsDT[0].columns)), mndfsDT[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsDT[1].columns)), mndfsDT[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[g])
plt.legend(['DT1','DT2','DT_av'])

plt.subplot(122)
plt.plot(np.linspace(0,1,len(mndfsVT[0].columns)), mndfsVT[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsVT[1].columns)), mndfsVT[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[g])
plt.legend(['VT1','VT2','VT_av'])

plt.show()


#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'
#oral tissue 50% brighter colour BM, HP, VT, DT: '#fdac6e','#74e7c5','#bab7d9','#f294c4'

xgs = ['ENSMUSG00000021194_Chga_ProteinCoding__12','ENSMUSG00000015468_Notch4_ProteinCoding__17',
       'ENSMUSG00000030029_Lrig1_ProteinCoding__6','ENSMUSG00000048616_Nog_ProteinCoding__11',
       'ENSMUSG00000059668_Krt4_ProteinCoding__15','ENSMUSG00000054252_Fgfr3_ProteinCoding__5',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000031565_Fgfr1_ProteinCoding__8']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_VT_shared_genes_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000057723_Krt33b_ProteinCoding__11','ENSMUSG00000067614_Krt86_ProteinCoding__15',
       'ENSMUSG00000078131_Krtap1.3_ProteinCoding__11','ENSMUSG00000044649_Krtap4.2_ProteinCoding__11',
       'ENSMUSG00000063251_Krtap4.1_ProteinCoding__11','ENSMUSG00000043485_Krt34_ProteinCoding__11',
       'ENSMUSG00000046259_Sprr2h_ProteinCoding__3','ENSMUSG00000064165_Krt39_ProteinCoding__11',
       'ENSMUSG00000035805_Mlc1_ProteinCoding__15']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_VT_unique_genes_DT_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

xgs = ['ENSMUSG00000023078_Cxcl13_ProteinCoding__5','ENSMUSG00000021381_Barx1_ProteinCoding__13',
       'ENSMUSG00000000305_Cdh4_ProteinCoding__2']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_DT_VT_unique_genes_VT_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()
