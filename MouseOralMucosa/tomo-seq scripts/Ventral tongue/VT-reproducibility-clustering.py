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
filesVT = sorted(glob.glob(path + '/VT*.coutt.tsv'))
print(filesVT)

dfsVT = [pd.read_csv(f, sep = '\t', index_col=0).astype(int) for f in filesVT]

dfsVT[0].columns = [str(i).zfill(2) + '-VT1' for i in range(len(dfsVT[0].columns))]
dfsVT[1].columns = [str(i).zfill(2) + '-VT2' for i in range(len(dfsVT[1].columns))]


#Genes to be annotated
g_annotate = ['Cdh1','Dicer1','Evpl','Fabp5','Cd44','Grem1','Jag1','Klf4','Notch1','Ovol2','Bmi1','Nog']



#Select genes present in all samples
def findUbiGenes(dfs):
    genes = []
    for df in dfs:
        genes += list(df.loc[df.sum(axis=1)>0].index)
    cnt = Counter(genes)
    gs = [g for g in cnt if cnt[g]==len(dfs)]
    print(len(cnt), len(gs))
    return gs

gVT = findUbiGenes(dfsVT)

dfsVT = [df.loc[gVT] for df in dfsVT]


#Number of reads and genes per section
plt.figure(figsize=(3*len(dfsVT),2.4))
for i, df in enumerate(dfsVT):
    plt.subplot(1,len(dfsVT),i+1)
    cnt = df.sum()
    plt.bar(cnt.keys(), cnt)
    plt.yscale('log')
plt.show()


#Normalise and prepare data
def removeGenes(df):
    genes = ['ERCC','__MT']
    exclude = [idx for idx in df.index for g in genes if g in idx]
    return df.loc[[idx for idx in df.index if idx not in exclude]]

def removeGmGenes(df):
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    numbers = [str(i) for i in range(10)]
    gmhumangenes = []
    for idx in df.index:
        g = idx.rsplit('_')[1]
        if len(g)>1 and all([x in alphabet for x in g[:2]]) and all([x in numbers for x in g[2:-2]]) and g[-2]=='.' and g[-1] in numbers:
            gmhumangenes.append(idx)
    return df.loc[[idx for idx in df.index if idx not in gmhumangenes]]

def removePeakyGenes(df):
    genes = []
    for idx in df.index:
        if sum(df.loc[idx]/df.loc[idx].sum()>=0.9)>=1:
            genes.append(idx)
    return df.loc[[idx for idx in df.index if idx not in genes]]

def filter_norm_dfs(dfs):
    fdfs = [removePeakyGenes(removeGmGenes(removeGenes(df))) for df in dfs]
    fdfs = [df.loc[df.index[df.sum(axis=1)>0]] for df in fdfs]
    nfdfs = [df.sum().mean()*df /df.sum() for df in fdfs]
    return nfdfs

nfdfsVT = filter_norm_dfs(dfsVT)

for i, df in enumerate(nfdfsVT):
    df.to_csv(path + 'VT' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

znfdfsVT = [zscore(df) for df in nfdfsVT]
sznfdfsVT = [scalezscore(df) for df in znfdfsVT]

[df.shape for df in znfdfsVT]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
       mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfVT = mergeDfs(znfdfsVT)
mdfVT.to_csv(path + '/VT_zscore.tsv', sep = '\t')


mdfsVT = [mdfVT[[idx for idx in mdfVT.columns if '-VT1' in idx]],
         mdfVT[[idx for idx in mdfVT.columns if '-VT2' in idx]]]


mnfVT = mergeDfs(nfdfsVT)
mnfVT.to_csv(path + '/VT_norm.tsv', sep = '\t')



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


mdfVT_av = averageGastruloid(mdfVT, min([len(df.columns) for df in znfdfsVT]))
mdfVT_av = mdfVT_av.astype(float)
mdfVT_av.to_csv(path + '/VT_average_zscore.tsv', sep = '\t')
mdfVT.to_csv(path + '/VT_zscore_both_libraries.tsv', sep = '\t')

mnfVT_av = averageGastruloid(mnfVT, min([len(df.columns) for df in nfdfsVT]))
mnfVT_av = mnfVT_av.astype(float)
mnfVT_av.to_csv(path + '/VT_average_norm.tsv', sep = '\t')


#Identify reproducible genes

#Background model
#def simulationCorrelations(n1, n2, N=5000, th=0.01):
def simulationCorrelations(n1, n2, N=5000):
    x1 = np.linspace(0,1,n1)
    x2 = np.linspace(0,1,n2)
    ps = [0 for i in range(N)]
    for i in range(N):
        d1 = np.random.rand(n1)
        d2 = np.random.rand(n2)

        id2 = interp1d(x2, d2, kind='quadratic')
        ps[i] = pearsonr(d1, id2(x1))[0]
    psdf = pd.DataFrame({'corr': sorted(ps), 'cum': np.linspace(0,1, N)})
    psdf['pval'] = 1.-psdf['cum']
    return psdf

from scipy.interpolate import interp1d

from collections import Counter

#Correlation between genes
def geneCorrelation(df1, df2):
    x1 = 1.0 * np.arange(len(df1.columns)) / (len(df1.columns) - 1)
    x2 = 1.0 * np.arange(len(df2.columns)) / (len(df2.columns) - 1)
    cnt = Counter(list(df1.index) + list(df2.index))
    genes = [g for g in cnt if cnt[g] == 2]

    cdf = pd.Series(0, index=genes)
    for g in genes:
        idf2 = interp1d(x2, df2.loc[g].values, kind='quadratic')
        cdf.loc[g] = pearsonr(df1.loc[g], idf2(x1))[0]
    return cdf

def assignPval2corr(cdf, pdf):
    cdf2 = pd.DataFrame(cdf, columns = ['corr'])
    cdf2['pv'] = 0
    for g in cdf.index:
        pv = cdf.loc[g]
        if (pdf['corr']>=pv).sum() > 0:
            p0 = pdf[pdf['corr']>=pv].iloc[0]['pval']
        else:
            p0 = 0.
        if (pdf['corr'] <pv).sum() > 0:
            p1 = pdf[pdf['corr'] <pv].iloc[-1]['pval']
        else:
            p1 = 1.
        cdf2.loc[g,'pv'] = 0.5*(p0+p1)
    cdf2['adj-pv'] = multi.fdrcorrection(cdf2['pv'],  method='indep')[1]
    return cdf2

def findCorrGenes(znfdfs):
    pv_sim = {}
    corrGenes = {}
    for i1,i2 in it.combinations(range(len(znfdfs)), 2):
        pv_sim[(i1,i2)] = simulationCorrelations(n1 = len(znfdfs[i1].columns), n2 = len(znfdfs[i2].columns))
        corrGenes[(i1,i2)] = geneCorrelation(znfdfs[i1], znfdfs[i2])
        corrGenes[(i1,i2)] = assignPval2corr(corrGenes[(i1,i2)], pv_sim[(i1,i2)])
    return corrGenes

corrGenes_VT = findCorrGenes(znfdfsVT)

pd.DataFrame(corrGenes_VT[(0, 1)]).to_csv(path + '/corr_genes_VT.tsv', sep = '\t')

corrGenes_VT[(0,1)].head()

#Define p-value threshold
pvth = 0.05

gl_VT_dict = {k: corrGenes_VT[k][corrGenes_VT[k]['adj-pv'] < pvth] for k in corrGenes_VT}
gl_VT = []
for k in gl_VT_dict:
    gl_VT += list(gl_VT_dict[k].index)

from collections import Counter

#Check reproducibility
gl_VT_cnt = Counter(gl_VT)
gl_VT = [g for g in gl_VT_cnt if gl_VT_cnt[g]==len(corrGenes_VT)]

reprGenes = set(list(gl_VT))
len(reprGenes)

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


5*np.sqrt(len(gl_VT))
#5*np.sqrt(1000)

#Clustering
somhcl_VT_df = comboCluster(mdfVT, gl_VT, cth=8)
somhcl_VT_df.sort_values(by=['hcl','som']).to_csv(path + '/clusters_VT_reprGenes.tsv', sep = '\t')

plt.figure(figsize=(10,10))
gene_order = somhcl_VT_df.sort_values(by=['hcl','som']).index

Npl = len(sznfdfsVT)+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4

titles = ['VT replicate 1','VT replicate 2']
for i in set(somhcl_VT_df['hcl']):
    y = 1.0*len(somhcl_VT_df[somhcl_VT_df['hcl'] == i])/len(somhcl_VT_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(sznfdfsVT):
    plt.subplot(1,Npl,i+2)
    plt.imshow(df.loc[gene_order][::-1], aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright='on', length = 2)
plt.yticks(range(len(gene_order)), [g.rsplit('_')[1] if g.rsplit('_')[1] in g_annotate else '' for g in gene_order[::-1]])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
#plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/heatmap_VT_reprGenes_all_replicates.pdf', bbox_inches = 'tight')
plt.show()


plt.figure(figsize=(10,7.5))
gene_order = somhcl_VT_df.sort_values(by=['hcl','som']).index

Npl = len(mdfsVT) + 1 + 2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(somhcl_VT_df['hcl']):
    y = 1.0*len(somhcl_VT_df[somhcl_VT_df['hcl'] == i])/len(somhcl_VT_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

print(mdfVT_av.head())

titles = ['VT replicate 1','VT replicate 2','VT average']
for i, df in enumerate(mdfsVT + [mdfVT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright='on', length = 2)
plt.yticks(range(len(gene_order)), [g.rsplit('_')[1] if g.rsplit('_')[1] in g_annotate else '' for g in gene_order[::-1]])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
#plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/heatmap_VT_reprGenes_all_replicates_and_average.pdf')
plt.show()

from matplotlib.patches import Rectangle, Circle

gene_order = somhcl_VT_df.sort_values(by=['hcl','som']).index

xzdfVT = mdfVT_av.astype(float)

fig, ax = plt.subplots(ncols=1, nrows=1, figsize = (1.5,10))
df = scalezscore(xzdfVT)
im = ax.imshow(df.loc[gene_order][::-1], aspect = 'auto')
ax.set_yticks([])
ax.tick_params('y', length = 2, labelleft=False, labelright=True, left = False, right=True)
ax.set_yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])
ax.set_yticklabels([g.rsplit('_')[1] for g in gene_order[::-1] if g.rsplit('_')[1] in g_annotate])

cbaxes = fig.add_axes([0.25, 0.05,0.5, 0.015])
fig.colorbar(im, label = 'Scaled z-score', orientation = 'horizontal', cax = cbaxes)

dx = 0.1*len(df.columns)
x0 = -2.5*dx
y0 = len(somhcl_VT_df)-0.5
for cl in set(somhcl_VT_df['hcl']):
    c = 'silver' if np.mod(cl,2)==0 else 'k'
    n = (somhcl_VT_df['hcl']==cl).sum()
    rect = Rectangle((x0,y0-n), dx, n, lw = 2, edgecolor = None, facecolor = c, clip_on=False)
    ax.add_patch(rect)
    y0 -= n

plt.savefig(path + '/heatmap_VT_reprGenes_average.pdf', bbox_inches = 'tight')
plt.show()


set(somhcl_VT_df['hcl'])

fig = plt.figure(figsize=(9, 9))
for cl in set(somhcl_VT_df['hcl']):
    ax = fig.add_subplot(3, 3, cl+1)
    genes = somhcl_VT_df[somhcl_VT_df['hcl'] == cl].index
    for i, df in enumerate(sznfdfsVT):
        xra = np.linspace(0, 1, len(df.columns))
        mu = df.loc[genes].mean()
        mup = mu + df.loc[genes].std()
        mup = np.array(mup.tolist(), dtype=float)
        mum = mu - df.loc[genes].std()
        mum = np.array(mum.tolist(), dtype=float)
        ax.plot(xra, mu, c=["#504b8e", "#9793c6"][i])
        ax.fill_between(xra, mup, mum, facecolor=["#504b8e", "#9793c6"][i], interpolate=True, alpha=0.2)

    ax.set_title('cl ' + str(cl) + ', ' + str(len(genes)) + ' genes')

    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1]);
    ax.set_xticklabels([])
    if cl+1 in [1, 4, 7]:
        ax.set_ylabel('Scaled expression (z-score)')
    if cl+1 >= 7:
        plt.xlabel('Relative sections (Epi->LP)')
        ax.set_xticklabels([0, 0.25, 0.50, 0.75, 1])

plt.legend(['VT1-mean','VT1-stdev','VT2-mean','VT2-stdev'], loc=2, bbox_to_anchor=(1, 0.6))
plt.savefig(path + '/reprGeneClusterVT.pdf', bbox_inches = 'tight')
plt.show()



xgs = ['ENSMUSG00000026728_Vim_ProteinCoding__2','ENSMUSG00000045545_Krt14_ProteinCoding__11',
       'ENSMUSG00000025809_Itgb1_ProteinCoding__8','ENSMUSG00000027111_Itga6_ProteinCoding__2']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000033227_Wnt6_ProteinCoding__1','ENSMUSG00000032204_Aqp9_ProteinCoding__9']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_0_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000023041_Krt6b_ProteinCoding__15','ENSMUSG00000030792_Dkkl1_ProteinCoding__7',
       'ENSMUSG00000017485_Top2b_ProteinCoding__14','ENSMUSG00000030265_Kras_ProteinCoding__6',
       'ENSMUSG00000024922_Ovol1_ProteinCoding__19','ENSMUSG00000008999_Bmp7_ProteinCoding__2',
       'ENSMUSG00000041415_Dicer1_ProteinCoding__12','ENSMUSG00000058354_Krt6a_ProteinCoding__15',
       'ENSMUSG00000027533_Fabp5_ProteinCoding__3','ENSMUSG00000036856_Wnt4_ProteinCoding__4',
       'ENSMUSG00000000303_Cdh1_ProteinCoding__8','ENSMUSG00000046095_Krt32_ProteinCoding__11',
       'ENSMUSG00000022587_Ly6e_ProteinCoding__15','ENSMUSG00000006777_Krt23_ProteinCoding__11',
       'ENSMUSG00000054146_Krt15_ProteinCoding__11','ENSMUSG00000027111_Itga6_ProteinCoding__2',
       'ENSMUSG00000074199_Krtdap_ProteinCoding__7','ENSMUSG00000030772_Dkk3_ProteinCoding__7',
       'ENSMUSG00000031757_Mt4_ProteinCoding__8','ENSMUSG00000031762_Mt2_ProteinCoding__8',
       'ENSMUSG00000020758_Itgb4_ProteinCoding__11','ENSMUSG00000001497_Pax9_ProteinCoding__12',
       'ENSMUSG00000034282_Evpl_ProteinCoding__11','ENSMUSG00000059325_Hopx_ProteinCoding__5',
       'ENSMUSG00000061527_Krt5_ProteinCoding__15','ENSMUSG00000020913_Krt24_ProteinCoding__11',
       'ENSMUSG00000063904_Dpp3_ProteinCoding__19','ENSMUSG00000038146_Notch3_ProteinCoding__17',
       'ENSMUSG00000026167_Wnt10a_ProteinCoding__1','ENSMUSG00000070047_Fat1_ProteinCoding__8',
       'ENSMUSG00000020105_Lrig3_ProteinCoding__10','ENSMUSG00000050359_Sprr1a_ProteinCoding__3',
       'ENSMUSG00000048455_Sprr1b_ProteinCoding__3','ENSMUSG00000054252_Fgfr3_ProteinCoding__5',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_1_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000067336_Bmpr2_ProteinCoding__1','ENSMUSG00000020644_Id2_ProteinCoding__12']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_2_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000026923_Notch1_ProteinCoding__2','ENSMUSG00000000078_Klf6_ProteinCoding__13',
       'ENSMUSG00000076431_Sox4_ProteinCoding__13','ENSMUSG00000074934_Grem1_ProteinCoding__2',
       'ENSMUSG00000005148_Klf5_ProteinCoding__14','ENSMUSG00000037279_Ovol2_ProteinCoding__2',
       'ENSMUSG00000005087_Cd44_ProteinCoding__2','ENSMUSG00000003032_Klf4_ProteinCoding__4',
       'ENSMUSG00000025880_Smad7_ProteinCoding__18','ENSMUSG00000024098_Twsg1_ProteinCoding__17',
       'ENSMUSG00000027276_Jag1_ProteinCoding__2','ENSMUSG00000020327_Fgf22_ProteinCoding__10']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_3_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000043485_Krt34_ProteinCoding__11','ENSMUSG00000042589_Cux2_ProteinCoding__5',
       'ENSMUSG00000048616_Nog_ProteinCoding__11','ENSMUSG00000032402_Smad3_ProteinCoding__9',
       'ENSMUSG00000036585_Fgf1_ProteinCoding__18','ENSMUSG00000026739_Bmi1_ProteinCoding__2',
       'ENSMUSG00000031565_Fgfr1_ProteinCoding__8','ENSMUSG00000031137_Fgf13_ProteinCoding__X']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_6_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000021732_Fgf10_ProteinCoding__13']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_7_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000035783_Acta2_ProteinCoding__19','ENSMUSG00000020717_Pecam1_ProteinCoding__11','ENSMUSG00000031871_Cdh5_ProteinCoding__8']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsVT[0].columns)), mdfsVT[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsVT[1].columns)), mdfsVT[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_2_marker_genes_line_graph_all_av_no_legend_VT_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()