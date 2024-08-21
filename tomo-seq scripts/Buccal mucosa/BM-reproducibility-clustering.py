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
filesBM = sorted(glob.glob(path + '/BM*.coutt.tsv'))
print(filesBM)

dfsBM = [pd.read_csv(f, sep = '\t', index_col=0).astype(int) for f in filesBM]

dfsBM[0].columns = [str(i).zfill(2) + '-BM1' for i in range(len(dfsBM[0].columns))]
dfsBM[1].columns = [str(i).zfill(2) + '-BM2' for i in range(len(dfsBM[1].columns))]


#Genes to be annotated
g_annotate = ['Des','Cdh13','Alox8','Hopx','Ly6e']


#Select genes present in all samples
def findUbiGenes(dfs):
    genes = []
    for df in dfs:
        genes += list(df.loc[df.sum(axis=1)>0].index)
    cnt = Counter(genes)
    gs = [g for g in cnt if cnt[g]==len(dfs)]
    print(len(cnt), len(gs))
    return gs

gBM = findUbiGenes(dfsBM)

dfsBM = [df.loc[gBM] for df in dfsBM]


#Number of reads and genes per section
plt.figure(figsize=(3*len(dfsBM),2.4))
for i, df in enumerate(dfsBM):
    plt.subplot(1,len(dfsBM),i+1)
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

nfdfsBM = filter_norm_dfs(dfsBM)

for i, df in enumerate(nfdfsBM):
    df.to_csv(path + 'BM' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

znfdfsBM = [zscore(df) for df in nfdfsBM]
sznfdfsBM = [scalezscore(df) for df in znfdfsBM]

[df.shape for df in znfdfsBM]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
       mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfBM = mergeDfs(znfdfsBM)
mdfBM.to_csv(path + '/BM_zscore.tsv', sep = '\t')


mdfsBM = [mdfBM[[idx for idx in mdfBM.columns if '-BM1' in idx]],
         mdfBM[[idx for idx in mdfBM.columns if '-BM2' in idx]]]


mnfBM = mergeDfs(nfdfsBM)
mnfBM.to_csv(path + '/BM_norm.tsv', sep = '\t')



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


mdfBM_av = averageGastruloid(mdfBM, min([len(df.columns) for df in znfdfsBM]))
mdfBM_av = mdfBM_av.astype(float)
mdfBM_av.to_csv(path + '/BM_average_zscore.tsv', sep = '\t')
mdfBM.to_csv(path + '/BM_zscore_both_libraries.tsv', sep = '\t')

mnfBM_av = averageGastruloid(mnfBM, min([len(df.columns) for df in nfdfsBM]))
mnfBM_av = mnfBM_av.astype(float)
mnfBM_av.to_csv(path + '/BM_average_norm.tsv', sep = '\t')


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

corrGenes_BM = findCorrGenes(znfdfsBM)

pd.DataFrame(corrGenes_BM[(0, 1)]).to_csv(path + '/corr_genes_BM.tsv', sep = '\t')

corrGenes_BM[(0,1)].head()

#Define p-value threshold
pvth = 0.05

gl_BM_dict = {k: corrGenes_BM[k][corrGenes_BM[k]['adj-pv'] < pvth] for k in corrGenes_BM}
gl_BM = []
for k in gl_BM_dict:
    gl_BM += list(gl_BM_dict[k].index)


from collections import Counter

#Check reproducibility
gl_BM_cnt = Counter(gl_BM)
gl_BM = [g for g in gl_BM_cnt if gl_BM_cnt[g]==len(corrGenes_BM)]

reprGenes = set(list(gl_BM))
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


5*np.sqrt(len(gl_BM))
#5*np.sqrt(1000)


#pd.DataFrame(gl_BM.to_csv(path + '/reprGenesBM.tsv', sep = '\t'))

#Clustering
somhcl_BM_df = comboCluster(mdfBM, gl_BM, cth=10)
somhcl_BM_df.sort_values(by=['hcl','som']).to_csv(path + '/clusters_BM_reprGenes.tsv', sep = '\t')

plt.figure(figsize=(10,10))
gene_order = somhcl_BM_df.sort_values(by=['hcl','som']).index

Npl = len(sznfdfsBM)+2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4

titles = ['BM replicate 1','BM replicate 2']
for i in set(somhcl_BM_df['hcl']):
    y = 1.0*len(somhcl_BM_df[somhcl_BM_df['hcl'] == i])/len(somhcl_BM_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(sznfdfsBM):
    plt.subplot(1,Npl,i+2)
    plt.imshow(df.loc[gene_order][::-1], aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright='on', length = 0)
plt.yticks(range(len(gene_order)), [g.rsplit('_')[1] if g.rsplit('_')[1] in g_annotate else '' for g in gene_order[::-1]])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
#plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/heatmap_BM_reprGenes_all_replicates.pdf', bbox_inches = 'tight')
plt.show()


plt.figure(figsize=(10,7.5))
gene_order = somhcl_BM_df.sort_values(by=['hcl','som']).index

Npl = len(mdfsBM) + 1 + 2
plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(somhcl_BM_df['hcl']):
    y = 1.0*len(somhcl_BM_df[somhcl_BM_df['hcl'] == i])/len(somhcl_BM_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

print(mdfBM_av.head())

titles = ['BM replicate 1','BM replicate 2','BM average']
for i, df in enumerate(mdfsBM + [mdfBM_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright='on', length = 0)
plt.yticks(range(len(gene_order)), [g.rsplit('_')[1] if g.rsplit('_')[1] in g_annotate else '' for g in gene_order[::-1]])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
#plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/heatmap_BM_reprGenes_all_replicates_and_average.pdf')
plt.show()

from matplotlib.patches import Rectangle, Circle

gene_order = somhcl_BM_df.sort_values(by=['hcl','som']).index

xzdfBM = mdfBM_av.astype(float)

fig, ax = plt.subplots(ncols=1, nrows=1, figsize = (1.5,10))
df = scalezscore(xzdfBM)
im = ax.imshow(df.loc[gene_order][::-1], aspect = 'auto')
ax.set_yticks([])
ax.tick_params('y', length = 2, labelleft=False, labelright=True, left = False, right=True)
ax.set_yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])
ax.set_yticklabels([g.rsplit('_')[1] for g in gene_order[::-1] if g.rsplit('_')[1] in g_annotate])

cbaxes = fig.add_axes([0.25, 0.05,0.5, 0.015])
fig.colorbar(im, label = 'Scaled z-score', orientation = 'horizontal', cax = cbaxes)

dx = 0.1*len(df.columns)
x0 = -2.5*dx
y0 = len(somhcl_BM_df)-0.5
for cl in set(somhcl_BM_df['hcl']):
    c = 'silver' if np.mod(cl,2)==0 else 'k'
    n = (somhcl_BM_df['hcl']==cl).sum()
    rect = Rectangle((x0,y0-n), dx, n, lw = 2, edgecolor = None, facecolor = c, clip_on=False)
    ax.add_patch(rect)
    y0 -= n

plt.savefig(path + '/heatmap_BM_reprGenes_average.pdf', bbox_inches = 'tight')
plt.show()


set(somhcl_BM_df['hcl'])

fig = plt.figure(figsize=(9, 9))
for cl in set(somhcl_BM_df['hcl']):
    ax = fig.add_subplot(3, 3, cl+1)
    genes = somhcl_BM_df[somhcl_BM_df['hcl'] == cl].index
    for i, df in enumerate(sznfdfsBM):
        xra = np.linspace(0, 1, len(df.columns))
        mu = df.loc[genes].mean()
        mup = mu + df.loc[genes].std()
        mup = np.array(mup.tolist(), dtype=float)
        mum = mu - df.loc[genes].std()
        mum = np.array(mum.tolist(), dtype=float)
        ax.plot(xra, mu, c=["#a24701", "#fd8326"][i])
        ax.fill_between(xra, mup, mum, facecolor=["#a24701", "#fd8326"][i], interpolate=True, alpha=0.2)

    ax.set_title('cl ' + str(cl) + ', ' + str(len(genes)) + ' genes')

    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1]);
    ax.set_xticklabels([])
    if cl+1 in [1, 4, 7]:
        ax.set_ylabel('Scaled expression (z-score)')
    if cl+1 >= 3:
        plt.xlabel('Relative sections (Epi->LP)')
        ax.set_xticklabels([0, 0.25, 0.50, 0.75, 1])

plt.legend(['BM1-mean','BM1-stdev','BM2-mean','BM2-stdev'], loc=2, bbox_to_anchor=(1, 0.6))
plt.savefig(path + '/reprGeneClusterBM.pdf', bbox_inches = 'tight')
plt.show()



#Plot interesting marker genes and reproducible genes found in clusters

xgs = ['ENSMUSG00000026728_Vim_ProteinCoding__2','ENSMUSG00000045545_Krt14_ProteinCoding__11',
       'ENSMUSG00000035493_Tgfbi_ProteinCoding__13','ENSMUSG00000028950_Tas1r1_ProteinCoding__4',
       'ENSMUSG00000025746_Il6_ProteinCoding__5','ENSMUSG00000031230_Fgf16_ProteinCoding__X',
       'ENSMUSG00000025809_Itgb1_ProteinCoding__8','ENSMUSG00000027111_Itga6_ProteinCoding__2']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000028773_Fabp3_ProteinCoding__4','ENSMUSG00000026208_Des_ProteinCoding__1',
       'ENSMUSG00000031137_Fgf13_ProteinCoding__X', 'ENSMUSG00000031841_Cdh13_ProteinCoding__8']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_2_marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000057967_Fgf18_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_4_marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000035557_Krt17_ProteinCoding__11','ENSMUSG00000059325_Hopx_ProteinCoding__5',
       'ENSMUSG00000020891_Alox8_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_5_marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000022587_Ly6e_ProteinCoding__15','ENSMUSG00000038156_Spon1_ProteinCoding__7',
       'ENSMUSG00000039323_Igfbp2_ProteinCoding__1','ENSMUSG00000048373_Fgfbp1_ProteinCoding__5']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_6_marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000020122_Egfr_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2.5, 2.5))
    plt.plot(np.linspace(0,1,len(mdfsBM[0].columns)), mdfsBM[0].loc[x], color='#8e8e8e')
    plt.plot(np.linspace(0,1,len(mdfsBM[1].columns)), mdfsBM[1].loc[x], color='#bebebe')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Cluster_7_marker_genes_line_graph_all_av_no_legend_BM_only_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()
