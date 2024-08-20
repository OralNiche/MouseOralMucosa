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
mdfBM = pd.read_csv(path + '/BM_zscore.tsv', sep = '\t', index_col = 0)
mdfVT = pd.read_csv(path + '/VT_zscore.tsv', sep = '\t', index_col = 0)

mndfBM = pd.read_csv(path + '/BM_norm.tsv', sep = '\t', index_col = 0)
mndfVT = pd.read_csv(path + '/VT_norm.tsv', sep = '\t', index_col = 0)

mdfsBM = [mdfBM[[idx for idx in mdfBM.columns if 'BM1' in idx]],
           mdfBM[[idx for idx in mdfBM.columns if 'BM2' in idx]]]

mdfsVT = [mdfVT[[idx for idx in mdfVT.columns if 'VT1' in idx]],
          mdfVT[[idx for idx in mdfVT.columns if 'VT2' in idx]]]

mndfsBM = [mndfBM[[idx for idx in mndfBM.columns if 'BM1' in idx]],
           mndfBM[[idx for idx in mndfBM.columns if 'BM2' in idx]]]

mndfsVT = [mndfVT[[idx for idx in mndfVT.columns if 'VT1' in idx]],
          mndfVT[[idx for idx in mndfVT.columns if 'VT2' in idx]]]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf

mdf_all = mergeDfs([mdfBM,mdfVT])

mndfBM = mergeDfs(mndfsBM)
mndfVT = mergeDfs(mndfsVT)

mndf_all = mergeDfs([mndfBM,mndfVT])

mdfVT = [mdf_all[[idx for idx in mdf_all.columns if 'VT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'VT2' in idx]]]

mdfBM = [mdf_all[[idx for idx in mdf_all.columns if 'BM1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'BM2' in idx]]]

mdfsVT = [mdf_all[[idx for idx in mdf_all.columns if 'VT1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'VT2' in idx]]]

mdfsBM = [mdf_all[[idx for idx in mdf_all.columns if 'BM1' in idx]],
          mdf_all[[idx for idx in mdf_all.columns if 'BM2' in idx]]]

mndfsVT = [mndf_all[[idx for idx in mndf_all.columns if 'VT1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'VT2' in idx]]]

mndfsBM = [mndf_all[[idx for idx in mndf_all.columns if 'BM1' in idx]],
          mndf_all[[idx for idx in mndf_all.columns if 'BM2' in idx]]]

0.5*sum([df.sum().sum() for df in mndfsBM]), 0.5*sum([df.sum().sum() for df in mndfsVT])

mndfsBM = [14306891.999999998*df/df.sum().sum() for df in mndfsBM]
mndfsVT = [22424833.499999993*df/df.sum().sum() for df in mndfsVT]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
        mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfBM = mergeDfs(mdfBM)
mdfVT = mergeDfs(mdfVT)

mndfBM = mergeDfs(mndfsBM)
mndfVT = mergeDfs(mndfsVT)

mndf_all = mergeDfs([mndfBM,mndfVT])


def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

corrGenes_BM = pd.read_table(path + '/corr_genes_BM.tsv', sep = '\t', index_col=0)

pvth = 0.05

gl_BM_dict = {k: corrGenes_BM[k][corrGenes_BM['adj-pv'] < pvth] for k in corrGenes_BM}
genelist_BM = []
for k in gl_BM_dict:
    genelist_BM  += list(gl_BM_dict[k].index)

from collections import Counter

gl_BM_cnt = Counter(genelist_BM)
gl_BM = [g for g in gl_BM_cnt]

corrGenes_VT = pd.read_table(path + '/corr_genes_VT.tsv', sep = '\t', index_col=0)

pvth = 0.05

gl_VT_dict = {k: corrGenes_VT[k][corrGenes_VT['adj-pv'] < pvth] for k in corrGenes_VT}
genelist_VT = []
for k in gl_VT_dict:
    genelist_VT  += list(gl_VT_dict[k].index)


gl_VT_cnt = Counter(genelist_VT)
gl_VT = [g for g in gl_VT_cnt]


import venn

labels = venn.get_labels([set(gl_BM),
                          set(gl_VT)],
                         fill=['number', 'percent'])

for k in labels:
    if labels[k][0:1] == '0':
        labels[k] = ''

oralcolors = [
            # r, g, b, a
            [202, 101, 40, 0.4],
            [116, 112, 174, 0.4],
        ]

oralcolors = [
            [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
            for i in oralcolors
        ]

fig, ax = venn.venn2(labels, names=['BM', 'VT'], colors = oralcolors)
plt.savefig(path + '/venn_genes_p_val_0.05_BM_vs_VT.pdf', bbox_inches='tight')
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


mdfBM_av = averageGastruloid(mdfBM, 17).astype(float)
mnfBM_av = averageGastruloid(mndfBM, 17).astype(float)

mdfVT_av = averageGastruloid(mdfVT, 20).astype(float)
mnfVT_av = averageGastruloid(mndfVT, 20).astype(float)


#Genes to be annotated
g_annotate = ['Zic1','Mhy11','Idi2','Acta2','Tagln','Sox15','Sprr2k','Krtap13',
              'Krt81','Krt34','Fgf1','Fgfr2','Grem1','Vegfb','Des']

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
cnt = Counter(list(gl_BM) + list(gl_VT))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)
#len(gl_both), 5*np.sqrt(len(gl_both))

gl_combined = list(set(list(gl_BM) + list(gl_VT)))
5*np.sqrt(len(gl_combined))

#Shared clustering
somhcl_all_df = comboCluster(mdf_all, gl_combined, cth=9)


plt.figure(figsize=(10,10))

hcldf = somhcl_all_df
gene_order = hcldf.sort_values(by=['hcl','som']).index

Npl = len(mndfsBM + mndfsVT)+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(hcldf['hcl']):
    y = 1.0*len(hcldf[hcldf['hcl'] == i])/len(hcldf)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(mdfsBM+mdfsVT):
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

plt.savefig(path + '/reprGenes_BM_VT_replicates.pdf', bbox_inches = 'tight')
hcldf.sort_values(by=['hcl','som']).to_csv(path + '/reprGenes_BM_VT.tsv', sep = '\t')
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
    if g in gl_BM and g in gl_VT:
       plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_BM and g not in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#D95F02')
    if g not in gl_BM and g in gl_VT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = '#7570B3')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['BM average','VT average']
for i, df in enumerate([mdfBM_av] + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])


plt.savefig(path + '/reprGenes_BM_VT_averages.pdf', bbox_inches = 'tight')
plt.show()


n10nfdfsBM = 1e5*mnfBM_av/mnfBM_av.sum().sum()
n10nfdfsVT = 1e5*mnfVT_av/mnfVT_av.sum().sum()


dex_BMVT = pd.DataFrame({
    'total-BM': n10nfdfsBM.sum(axis=1),
    'total-VT': n10nfdfsVT.sum(axis=1),
    'total-mean': {idx: np.mean(list(n10nfdfsBM.loc[idx].values) + list(n10nfdfsVT.loc[idx].values)) for idx in set(list(n10nfdfsBM.index)+list(n10nfdfsVT.index)) if idx in n10nfdfsBM.index and idx in n10nfdfsVT.index}
})

dex_BMVT = dex_BMVT.loc[[idx for idx in dex_BMVT.index if idx in mndfBM.index and idx in mndfVT.index]]


dex_BMVT['log2FC'] = np.log2((dex_BMVT['total-BM']+1e-10)/(dex_BMVT['total-VT']+1e-10))
dex_BMVT['log2mu'] = np.log2(0.5*(dex_BMVT['total-BM']+dex_BMVT['total-VT'])) # np.log2(dex_Chi['total-mean'])

dex_BMVT['prob-BM'] = dex_BMVT['total-BM']/dex_BMVT['total-BM'].sum()
dex_BMVT['prob-VT'] = dex_BMVT['total-VT']/dex_BMVT['total-VT'].sum()

dex_BMVT['N-BM'] = dex_BMVT['total-BM'].sum()
dex_BMVT['N-VT'] = dex_BMVT['total-VT'].sum()

fdex_BMVT = dex_BMVT.loc[[g for g in dex_BMVT.index if g in gl_combined]]

from scipy.stats import binom

def dex_Pval(expression, total, prob):
    expression, total = int(expression), int(total)
    pv = binom.cdf(expression, total, prob)
    if pv > 0.5:
        pv = 1-pv
    return pv

fdex_BMVT['PvalBM'] = fdex_BMVT.apply(lambda x: dex_Pval(x['total-BM'],x['N-BM'],x['prob-BM']), axis = 1)
fdex_BMVT['PvalVT'] = fdex_BMVT.apply(lambda x: dex_Pval(x['total-VT'],x['N-VT'],x['prob-VT']), axis = 1)

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_BMVT['log2mu'], fdex_BMVT['log2FC'], s = 5, c = np.log10(fdex_BMVT['PvalBM']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_BMVT[np.log10(fdex_BMVT['PvalBM'])<-2.25].index:
    x = fdex_BMVT.loc[g,'log2mu']
    y = fdex_BMVT.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('log2(mean)')
plt.ylabel('log2(FC)')

fdex_BMVT.to_csv(path + '/dex_BMVT_sameNormalization.tsv', sep = '\t')
plt.savefig(path + '/dex_BMVT_sameNormalization.pdf', bbox_inches = 'tight')

plt.show()

# Differential gene expression analysis
plt.figure(figsize=(5*1.6,5))

plt.subplot(111)
plt.scatter(fdex_BMVT['total-mean'], fdex_BMVT['log2FC'], s = 5, c = np.log10(fdex_BMVT['PvalBM']+1e-10))
plt.colorbar(label = 'log10(Pval)')
for g in fdex_BMVT[np.log10(fdex_BMVT['PvalBM'])<-3].index:
    x = fdex_BMVT.loc[g,'total-mean']
    y = fdex_BMVT.loc[g,'log2FC']
    plt.text(x, y, g.rsplit('_')[1], fontsize = 8)
plt.xlabel('mean')
plt.ylabel('log2(FC)')

plt.savefig(path + '/dex_BMVT_mean.pdf', bbox_inches = 'tight')

plt.show()


fdex_BMVT.loc[[g for g in fdex_BMVT[np.log10(fdex_BMVT['PvalBM'])<-2].index]]

[idx for idx in mnfBM_av.index if 'Mog' in idx]

g = 'ENSMUSG00000023084_Lrrc71_ProteinCoding__3'

plt.figure(figsize=(6*1.6,3))
plt.subplot(121)
plt.plot(np.linspace(0,1,len(mndfsBM[0].columns)), mndfsBM[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsBM[1].columns)), mndfsBM[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[g])
plt.legend(['BM1','BM2','BM_av'])

plt.subplot(122)
plt.plot(np.linspace(0,1,len(mndfsVT[0].columns)), mndfsVT[0].loc[g])
plt.plot(np.linspace(0,1,len(mndfsVT[1].columns)), mndfsVT[1].loc[g])
plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[g])
plt.legend(['VT1','VT2','VT_av'])

plt.show()

#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'
#oral tissue 50% brighter colour BM, HP, VT, DT: '#fdac6e','#74e7c5','#bab7d9','#f294c4'


xgs = ['ENSMUSG00000026208_Des_ProteinCoding__1','ENSMUSG00000024962_Vegfb_ProteinCoding__19',
       'ENSMUSG00000031137_Fgf13_ProteinCoding__X','ENSMUSG00000008090_Fgfrl1_ProteinCoding__5',
       'ENSMUSG00000036585_Fgf1_ProteinCoding__18','ENSMUSG00000021732_Fgf10_ProteinCoding__13',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000020327_Fgf22_ProteinCoding__10',
       'ENSMUSG00000039004_Bmp6_ProteinCoding__13','ENSMUSG00000008999_Bmp7_ProteinCoding__2',
       'ENSMUSG00000067336_Bmpr2_ProteinCoding__1','ENSMUSG00000057967_Fgf18_ProteinCoding__11',
       'ENSMUSG00000074934_Grem1_ProteinCoding__2','ENSMUSG00000041415_Dicer1_ProteinCoding__12',
       'ENSMUSG00000026739_Bmi1_ProteinCoding__2','ENSMUSG00000026167_Wnt10a_ProteinCoding__1',
       'ENSMUSG00000033227_Wnt6_ProteinCoding__1']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_BM_VT_shared_genes_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000043485_Krt34_ProteinCoding__11','ENSMUSG00000067615_Krt81_ProteinCoding__15',
       'ENSMUSG00000041287_Sox15_ProteinCoding__11','ENSMUSG00000054215_Sprr2k_ProteinCoding__3',
       'ENSMUSG00000057967_Fgf18_ProteinCoding__11','ENSMUSG00000074934_Grem1_ProteinCoding__2',
       'ENSMUSG00000050224_Krtap13_ProteinCoding__16','ENSMUSG00000051879_Krt71_ProteinCoding__15',
       'ENSMUSG00000022887_Masp1_ProteinCoding__16','ENSMUSG00000053522_Lgals7_ProteinCoding__7',
       'ENSMUSG00000042045_Sln_ProteinCoding__9','ENSMUSG00000038402_Foxf2_ProteinCoding__13',
       'ENSMUSG00000044430_Klk12_ProteinCoding__7','ENSMUSG00000059230_Defb4_ProteinCoding__8',
       'ENSMUSG00000023041_Krt6b_ProteinCoding__15']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_BM_VT_unique_genes_VT_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

xgs = ['ENSMUSG00000032368_Zic1_ProteinCoding__9','ENSMUSG00000033520_Idi2_ProteinCoding__13',
       'ENSMUSG00000018830_Myh11_ProteinCoding__16','ENSMUSG00000035783_Acta2_ProteinCoding__19',
       'ENSMUSG00000038583_Pln_ProteinCoding__10','ENSMUSG00000032085_Tagln_ProteinCoding__9']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0, 1, len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Comparison_BM_VT_unique_genes_BM_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()