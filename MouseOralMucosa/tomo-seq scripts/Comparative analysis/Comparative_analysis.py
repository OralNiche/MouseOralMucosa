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

#for i, df in enumerate(nfdfsBM):
#    df.to_csv(path + 'BM' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

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
#mdfBM.to_csv(path + '/BM_zscore.tsv', sep = '\t')


mdfsBM = [mdfBM[[idx for idx in mdfBM.columns if '-BM1' in idx]],
         mdfBM[[idx for idx in mdfBM.columns if '-BM2' in idx]]]


mnfBM = mergeDfs(nfdfsBM)
#mnfBM.to_csv(path + '/BM_norm.tsv', sep = '\t')



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

mnfBM_av = averageGastruloid(mnfBM, min([len(df.columns) for df in nfdfsBM]))
mnfBM_av = mnfBM_av.astype(float)


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

def hierarchicalClustering(df, cth = 100, plot = False, method = 'ward', metric = 'euclidean', nolabels = 'True'):
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

def comboCluster(df, genes, gridheight = 10, gridwidth = 10, cth = 100, plot = True, method = 'ward', metric = 'euclidean', seed = 12345, nolabels = 'True'):
    somcl = SOMclustering(df.loc[genes], gridheight = gridheight, gridwidth = gridwidth, seed = seed)
    mdf_som = pd.DataFrame({i: df.loc[somcl[somcl==i].index].mean() for i in set(somcl) })
    Z, dg = hierarchicalClustering(mdf_som.T, cth = cth, plot = plot, method = method, metric = metric, nolabels = nolabels)
    if plot:
        plt.show()
    ghcl = getClusterByColor(dg, list(mdf_som.columns))
    cl_df = pd.DataFrame(somcl, columns = ['som'])
    cl_df['hcl'] = [k for s in cl_df['som'] for k in ghcl.keys() if s in ghcl[k]]
    return cl_df



np.random.seed(174612)

#Read data
filesHP = sorted(glob.glob(path + '/HP*.coutt.tsv'))
print(filesHP)

dfsHP = [pd.read_csv(f, sep = '\t', index_col=0).astype(int) for f in filesHP]

dfsHP[0].columns = [str(i).zfill(2) + '-HP1' for i in range(len(dfsHP[0].columns))]
dfsHP[1].columns = [str(i).zfill(2) + '-HP2' for i in range(len(dfsHP[1].columns))]


#Genes to be annotated
g_annotate = ['Myh11','Pdgfa','Pdgfrb','Prrx1','Runx1','Krt10','Ovol1','Ovol2','Slurp1','Slurp2']


#Select genes present in all samples
def findUbiGenes(dfs):
    genes = []
    for df in dfs:
        genes += list(df.loc[df.sum(axis=1)>0].index)
    cnt = Counter(genes)
    gs = [g for g in cnt if cnt[g]==len(dfs)]
    print(len(cnt), len(gs))
    return gs

gHP = findUbiGenes(dfsHP)

dfsHP = [df.loc[gHP] for df in dfsHP]


#Number of reads and genes per section
plt.figure(figsize=(3*len(dfsHP),2.4))
for i, df in enumerate(dfsHP):
    plt.subplot(1,len(dfsHP),i+1)
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

nfdfsHP = filter_norm_dfs(dfsHP)

#for i, df in enumerate(nfdfsHP):
#    df.to_csv(path + 'HP' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

znfdfsHP = [zscore(df) for df in nfdfsHP]
sznfdfsHP = [scalezscore(df) for df in znfdfsHP]

[df.shape for df in znfdfsHP]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
       mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfHP = mergeDfs(znfdfsHP)


mdfsHP = [mdfHP[[idx for idx in mdfHP.columns if '-HP1' in idx]],
         mdfHP[[idx for idx in mdfHP.columns if '-HP2' in idx]]]


mnfHP = mergeDfs(nfdfsHP)


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


mdfHP_av = averageGastruloid(mdfHP, min([len(df.columns) for df in znfdfsHP]))
mdfHP_av = mdfHP_av.astype(float)

mnfHP_av = averageGastruloid(mnfHP, min([len(df.columns) for df in nfdfsHP]))
mnfHP_av = mnfHP_av.astype(float)

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

corrGenes_HP = findCorrGenes(znfdfsHP)
corrGenes_HP[(0,1)].head()

#Define p-value threshold
pvth = 0.05

gl_HP_dict = {k: corrGenes_HP[k][corrGenes_HP[k]['adj-pv'] < pvth] for k in corrGenes_HP}
gl_HP = []
for k in gl_HP_dict:
    gl_HP += list(gl_HP_dict[k].index)

from collections import Counter

#Check reproducibility
gl_HP_cnt = Counter(gl_HP)
gl_HP = [g for g in gl_HP_cnt if gl_HP_cnt[g]==len(corrGenes_HP)]

reprGenes = set(list(gl_HP))
len(reprGenes)

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

#for i, df in enumerate(nfdfsVT):
#    df.to_csv(path + 'VT' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

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



mdfsVT = [mdfVT[[idx for idx in mdfVT.columns if '-VT1' in idx]],
         mdfVT[[idx for idx in mdfVT.columns if '-VT2' in idx]]]


mnfVT = mergeDfs(nfdfsVT)




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


mnfVT_av = averageGastruloid(mnfVT, min([len(df.columns) for df in nfdfsVT]))
mnfVT_av = mnfVT_av.astype(float)



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

np.random.seed(174612)

#Read data
filesDT = sorted(glob.glob(path + '/DT*.coutt.tsv'))
print(filesDT)

dfsDT = [pd.read_csv(f, sep = '\t', index_col=0).astype(int) for f in filesDT]

dfsDT[0].columns = [str(i).zfill(2) + '-DT1' for i in range(len(dfsDT[0].columns))]
dfsDT[1].columns = [str(i).zfill(2) + '-DT2' for i in range(len(dfsDT[1].columns))]


#Genes to be annotated

g_annotate = ['Bmi1','Cdh5','Nog','Pdgfa','Pecam1',
              'Krt10','Krt14','Sprr2b','Sprr1b','Krt16',
              'Lor','Slurp1','Slurp2','Itga3','Ppara','Emb']


#Select genes present in all samples
def findUbiGenes(dfs):
    genes = []
    for df in dfs:
        genes += list(df.loc[df.sum(axis=1)>0].index)
    cnt = Counter(genes)
    gs = [g for g in cnt if cnt[g]==len(dfs)]
    print(len(cnt), len(gs))
    return gs

gDT = findUbiGenes(dfsDT)

dfsDT = [df.loc[gDT] for df in dfsDT]


#Number of reads and genes per section
plt.figure(figsize=(3*len(dfsDT),2.4))
for i, df in enumerate(dfsDT):
    plt.subplot(1,len(dfsDT),i+1)
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

nfdfsDT = filter_norm_dfs(dfsDT)

#for i, df in enumerate(nfdfsDT):
#    df.to_csv(path + 'DT' + str(i+1) +'.filtnorm.coutt.tsv', sep = '\t')

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df

znfdfsDT = [zscore(df) for df in nfdfsDT]
sznfdfsDT = [scalezscore(df) for df in znfdfsDT]

[df.shape for df in znfdfsDT]

def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
       mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf


mdfDT = mergeDfs(znfdfsDT)



mdfsDT = [mdfDT[[idx for idx in mdfDT.columns if '-DT1' in idx]],
         mdfDT[[idx for idx in mdfDT.columns if '-DT2' in idx]]]


mnfDT = mergeDfs(nfdfsDT)




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


mdfDT_av = averageGastruloid(mdfDT, min([len(df.columns) for df in znfdfsDT]))
mdfDT_av = mdfDT_av.astype(float)


mnfDT_av = averageGastruloid(mnfDT, min([len(df.columns) for df in nfdfsDT]))
mnfDT_av = mnfDT_av.astype(float)



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

corrGenes_DT = findCorrGenes(znfdfsDT)
corrGenes_DT[(0,1)].head()

#Define p-value threshold
pvth = 0.05

gl_DT_dict = {k: corrGenes_DT[k][corrGenes_DT[k]['adj-pv'] < pvth] for k in corrGenes_DT}
gl_DT = []
for k in gl_DT_dict:
    gl_DT += list(gl_DT_dict[k].index)

from collections import Counter

#Check reproducibility
gl_DT_cnt = Counter(gl_DT)
gl_DT = [g for g in gl_DT_cnt if gl_DT_cnt[g]==len(corrGenes_DT)]

reprGenes = set(list(gl_DT))
len(reprGenes)


mdf = mergeDfs([mdfBM, mdfHP, mdfDT, mdfVT])
gl_all = list(set(list(gl_BM) + list(gl_HP) + list(gl_DT) + list(gl_VT)))

#gl_list = pd.DataFrame(gl_all)

#gl_list.to_csv(path + '/sectioned_genes.tsv', sep = '\t')

import venn

labels = venn.get_labels([set(gl_BM),
                       set(gl_HP),
                     set(gl_DT),
                      set(gl_VT)],
                        fill=['number', 'percent'])

#labels = venn.generate_petal_labels([set(gl_BM),
 #                         set(gl_HP),
  #                      set(gl_DT),
    #                    set(gl_VT)],
     #                    fill=['number', 'percent'])


for k in labels:
    if labels[k][0:1] == '0':
        labels[k] = ''

oralcolors = [
            # r, g, b, a
            [202, 101, 40, 0.4],
            [75, 155, 121, 0.4],
            [213, 62, 135, 0.4],
            [116, 112, 174, 0.4],
        ]


oralcolors = [
    [i[0] / 255.0, i[1] / 255.0, i[2] / 255.0, i[3]]
    for i in oralcolors
]


#fig, ax = venn.venn4(labels, names=['BM', 'HP','DT','VT'])
#fig, ax = venn.venn4(labels, names=['BM', 'HP','DT','VT'])
#fig, ax = venn.venn4(labels, names=['BM', 'HP','DT','VT'], colors= ['#D95F02', '#1B9E77','#E7298A','#7570B3')
fig, ax = venn.venn4(labels, names = ['BM','HP','DT','VT'], colors = oralcolors)

#ax.get_patch_by_id('BM').set_color('#D95F02')
#ax.get_patch_by_id('HP').set_color('#1B9E77')
#ax.get_patch_by_id('VT').set_color('#7570B3')
#ax.get_patch_by_id('DT').set_color('#E7298A')

plt.savefig(path + '/venn_genes_p_val_0.05_all.pdf', bbox_inches='tight')
plt.show()



cnt = Counter(list(gl_DT) + list(gl_VT))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)

len(gl_both), 5*np.sqrt(len(gl_both))

somhcl_all_df = comboCluster(mdf, gl_both, cth=7, gridheight=10, gridwidth=10)

plt.figure(figsize=(10,13))

gene_order = somhcl_all_df.sort_values(by=['hcl','som']).index
Npl = len(sznfdfsDT)+len(sznfdfsVT)+2
titles = ['DT-1', 'DT-2', 'VT-1', 'VT-2']

plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(somhcl_all_df['hcl']):
    y = 1.0*len(somhcl_all_df[somhcl_all_df['hcl'] == i])/len(somhcl_all_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(sznfdfsDT+sznfdfsVT):
    plt.subplot(1,Npl,i+2)
    plt.imshow(df.loc[gene_order][::-1], aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
#plt.tick_params(axis='y', labelleft=False, labelright=True, length = 4, left = False, right = True)
#plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])
jidx = 4
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)
plt.yticks([i for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/reprGenes_DT_VT.pdf', bbox_inches = 'tight')
somhcl_all_df.sort_values(by=['hcl','som']).to_csv(path + '/reprGenes_DT_VT.tsv', sep = '\t')

plt.show()

nndfs_DT = [1e6*df/df.sum() for df in nfdfsDT]
nndfs_VT = [1e6*df/df.sum() for df in nfdfsVT]
nndfs_all = nndfs_DT + nndfs_VT

nndfs_DT = [df.loc[gl_both] for df in nndfs_DT]
nndfs_VT = [df.loc[gl_both] for df in nndfs_VT]
nndfs_all = [df.loc[gl_both] for df in nndfs_all]

dge_df = pd.DataFrame({g: {
    'meanDT': np.mean([df.loc[g].sum() for df in nndfs_DT]),
    'meanVT': np.mean([df.loc[g].sum() for df in nndfs_VT])
} for g in gl_both}).T

dge_df['mean'] = dge_df[['meanDT','meanVT']].mean(axis=1)
dge_df['FC'] = dge_df['meanDT']/dge_df['meanVT']

dge_df.head()

dge_df.to_csv(path + '/reprGenes_DT_VT_expression_levels.tsv', sep = '\t')


fig, ax = plt.subplots()
ax.scatter(dge_df['mean'], np.log2(dge_df['FC']), s = 5, color = 'grey')
ax.axhline(0, c = 'k')

for g in [dge_df[np.log2(dge_df['FC'])>2].index[0]]+list(dge_df[(np.log2(dge_df['FC'])>2)&(dge_df['mean']>10000)].index)+list(dge_df[(np.log2(dge_df['FC'])>1.1)&(dge_df['mean']>40000)].index)+list(dge_df[(np.log2(dge_df['FC'])<-2)].index):
    x = dge_df.loc[g,'mean']
    y = np.log2(dge_df.loc[g,'FC'])
    ax.annotate(g.rsplit('_')[1], (x,  y))

ax.set_ylabel('log2(FC_DT/VT)')
ax.set_xlabel('mean')

plt.savefig(path + '/reprGenes_DT_VT_expression_levels.pdf', bbox_inches = 'tight')
plt.show()


cnt = Counter(list(gl_VT) + list(gl_DT))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)

len(gl_both), 5*np.sqrt(len(gl_both))

somhcl_all_df = comboCluster(mdf, gl_both, cth=11, gridheight=10, gridwidth=10)

plt.figure(figsize=(10,10))

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
    if g in gl_VT and g in gl_DT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'silver')
    if g in gl_VT and g not in gl_DT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'r')
    if g not in gl_VT and g in gl_DT:
        plt.bar([xb], [y], bottom = y0, width = 0.3, color = 'b')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(3,5)

titles = ['VT','DT']
for i, df in enumerate([mdfVT_av] + [mdfDT_av]):
    plt.subplot(1,Npl,i+2)
    plt.imshow(scalezscore(df.loc[gene_order][::-1]), aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
#plt.tick_params(axis='y', labelleft=False, labelright='on', length = 0)
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)

plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
#plt.colorbar(label = 'Scaled z-score')

plt.savefig(path + '/reprGenes_VT_DT_localisation.pdf', bbox_inches = 'tight')
plt.show()




# n10nfdfsChi = mndfChi_av
# n10nfdfsSB = mndfSB_av





cnt = Counter(list(gl_BM) + list(gl_VT))
gl_both = []
for g in cnt:
    if cnt[g]==2:
        gl_both.append(g)

gl_both = pd.Series(gl_both)

len(gl_both), 5*np.sqrt(len(gl_both))

somhcl_all_df = comboCluster(mdf, gl_both, cth=9, gridheight=10, gridwidth=10)

plt.figure(figsize=(10,13))

gene_order = somhcl_all_df.sort_values(by=['hcl','som']).index
Npl = len(sznfdfsBM)+len(sznfdfsVT)+2
titles = ['BM-1', 'BM-2', 'VT-1', 'VT-2']

plt.subplot(1,Npl,1, frameon = False)

y0 = 0; z0 = 0; xb = 4 # 3
for i in set(somhcl_all_df['hcl']):
    y = 1.0*len(somhcl_all_df[somhcl_all_df['hcl'] == i])/len(somhcl_all_df)
    plt.bar([xb], [y], bottom = y0, width = 0.3, color = Colors.bg[i])
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2.5,4.5)

for i, df in enumerate(sznfdfsBM+sznfdfsVT):
    plt.subplot(1,Npl,i+2)
    plt.imshow(df.loc[gene_order][::-1], aspect = 'auto')
    plt.yticks([])
    plt.title(titles[i])
#plt.tick_params(axis='y', labelleft=False, labelright=True, length = 4, left = False, right = True)
#plt.yticks([i for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if g.rsplit('_')[1] in g_annotate])
jidx = 4
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)
plt.yticks([i for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx], [g.rsplit('_')[1] for i, g in enumerate(gene_order[::-1]) if np.mod(i,8)==jidx])

plt.subplot(1,Npl, Npl, frameon = False)
plt.xticks([]); plt.yticks([])
plt.colorbar(label = 'scaled z-score')

plt.savefig(path + '/reprGenes_BM_VT.pdf', bbox_inches = 'tight')
somhcl_all_df.sort_values(by=['hcl','som']).to_csv(path + '/reprGenes_BM_VT.tsv', sep = '\t')

plt.show()

nndfs_BM = [1e6*df/df.sum() for df in nfdfsBM]
nndfs_VT = [1e6*df/df.sum() for df in nfdfsVT]
nndfs_all = nndfs_BM + nndfs_VT

nndfs_BM = [df.loc[gl_both] for df in nndfs_BM]
nndfs_VT = [df.loc[gl_both] for df in nndfs_VT]
nndfs_all = [df.loc[gl_both] for df in nndfs_all]

dge_df = pd.DataFrame({g: {
    'meanBM': np.mean([df.loc[g].sum() for df in nndfs_BM]),
    'meanVT': np.mean([df.loc[g].sum() for df in nndfs_VT])
} for g in gl_both}).T

dge_df['mean'] = dge_df[['meanBM','meanVT']].mean(axis=1)
dge_df['FC'] = dge_df['meanBM']/dge_df['meanVT']
dge_df.head()

dge_df.to_csv(path + '/reprGenes_BM_VT_expression_levels.tsv', sep = '\t')

fig, ax = plt.subplots()
#ax.scatter(dge_df['mean'], np.log2(dge_df['FC'])>0, s = 5, color = '#D95F02')
#ax.scatter(dge_df['mean'], np.log2(dge_df['FC'])==0, s = 5, color = 'black')
ax.scatter(dge_df['mean'], np.log2(dge_df['FC']), s = 5, color = 'grey')
#ax.scatter(dge_df['mean'], np.log2(dge_df['FC'])<0, s = 5, color = '#7570B3')
ax.axhline(0, c = 'k')

#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'

for g in [dge_df[np.log2(dge_df['FC'])>2].index[0]]+list(dge_df[(np.log2(dge_df['FC'])>2)&(dge_df['mean']>10000)].index)+list(dge_df[(np.log2(dge_df['FC'])>1.1)&(dge_df['mean']>40000)].index)+list(dge_df[(np.log2(dge_df['FC'])<-2)].index):
    x = dge_df.loc[g,'mean']
    y = np.log2(dge_df.loc[g,'FC'])
    ax.annotate(g.rsplit('_')[1], (x,  y))

ax.set_ylabel('log2(FC_BM/VT)')
ax.set_xlabel('mean')

plt.savefig(path + '/reprGenes_BM_VT_expression_levels.pdf', bbox_inches = 'tight')
plt.show()

for cl in set(somhcl_all_df['hcl']):
    print(cl, '&', ', '.join(sorted([g.rsplit('_')[1] for g in somhcl_all_df[somhcl_all_df['hcl']==cl].index])), '\\\\')


germlayer_genes = pd.read_excel(path + 'MarkerGenes.xlsx', sheet_name='ClassMarkerGenes')
f_germlayer_genes = germlayer_genes.loc[[idx for idx in germlayer_genes.index if germlayer_genes.loc[idx, 'Cell class'] in ['Epithelial','Immune','Fibroblast','Endothelial','Neural','Muscle']]]
f_germlayer_genes.head()

#f_germlayer_genes['gene name'] = [f_germlayer_genes.loc[idx,'gene name'].rsplit('__')[0] for idx in f_germlayer_genes.index]
f_germlayer_genes.head()

mdfAll = mergeDfs([mdfHP]+[mdfBM]+[mdfVT]+[mdfDT])

hcl_germlayer = linkage(mdfAll.loc[f_germlayer_genes['gene name']], method='ward', metric='euclidean')
dg_germlayer = dendrogram(hcl_germlayer, no_labels=True, color_threshold = 25) # 5
plt.show()

hcl_germlayer = getClusterByColor(dg_germlayer, list(f_germlayer_genes['gene name']))

hcl_germlayer = pd.Series({g: hcl  for hcl in hcl_germlayer for g in hcl_germlayer[hcl]})

go = hcl_germlayer.sort_values().index
N = len(znfdfsHP) + len(znfdfsBM) + len(znfdfsVT) + len(znfdfsDT) + 2

plt.figure(figsize=(13,13))
plt.subplot(1,N,1, frameon = False)
d = {'Epithelial': 0, 'Immune': 1, 'Fibroblast': 2, 'Endothelial': 3, 'Neural': 4, 'Muscle': 5}
cnt = Counter()
y0 = 0; z0 = 0
for g in go[::-1]:
    y = 1.0/len(go)
    c = f_germlayer_genes.set_index('gene name').loc[g, 'Cell class']
    cnt.update([c])
    if cnt[c] == 1:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]])
    else:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]], label='_nolegend_')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2,4)
plt.legend(['Muscle','Neural','Endothelial','Fibroblast','Immune','Epithelial'], loc = 1, bbox_to_anchor = (0,1))

titles = ['HP1','HP2','BM1','BM2','VT1','VT2','DT1','DT2',]
for i, df in enumerate(znfdfsHP+znfdfsBM+znfdfsVT+znfdfsDT):
    plt.subplot(1,N,i+2)
    plt.imshow(scalezscore(df.loc[go]), aspect = 'auto')
    if i == 0:
        plt.ylabel('genes');
    plt.xlabel('Sections (Epi -> LP)')
    plt.yticks([])
    plt.title(titles[i])
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 4, left = False, right = True)
plt.yticks(range(len(go)), [x.rsplit('_')[1]  for i, x in enumerate(go)])
plt.xlabel('Sections (Epi -> LP)')

plt.subplot(1,N,N, frameon= False)
plt.colorbar(label = 'z-score (scaled)', ticks = [-1,0,1], shrink = 0.5)
plt.xticks([]); plt.yticks([])

plt.savefig(path + 'MarkerGeneHeatmap_samples.pdf', bbox_inches ='tight')

plt.show()



col = 'Cell class'

mznfdfchis = mdfAll

hcl_epithelial = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Epithelial']['gene name']], method='ward', metric='euclidean')
hcl_immune = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Immune']['gene name']], method='ward', metric='euclidean')
hcl_fibroblast = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Fibroblast']['gene name']], method='ward', metric='euclidean')
hcl_endo = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Endothelial']['gene name']], method='ward', metric='euclidean')
hcl_neural = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Neural']['gene name']], method='ward', metric='euclidean')
hcl_muscle = linkage(mznfdfchis.loc[f_germlayer_genes[f_germlayer_genes[col]=='Muscle']['gene name']], method='ward', metric='euclidean')

dg_epithelial = dendrogram(hcl_epithelial, no_labels=True, color_threshold = 7, no_plot=True)
dg_immune = dendrogram(hcl_immune, no_labels=True, color_threshold = 7, no_plot=True)
dg_fibroblast = dendrogram(hcl_fibroblast, no_labels=True, color_threshold = 7, no_plot=True)
dg_endo = dendrogram(hcl_endo, no_labels=True, color_threshold = 7, no_plot=True)
dg_neural = dendrogram(hcl_neural, no_labels=True, color_threshold = 7, no_plot=True)
dg_muscle = dendrogram(hcl_muscle, no_labels=True, color_threshold = 7, no_plot=True)

hcl_epithelial = getClusterByColor(dg_epithelial, list(f_germlayer_genes[f_germlayer_genes[col]=='Epithelial']['gene name']))
hcl_immune = getClusterByColor(dg_immune, list(f_germlayer_genes[f_germlayer_genes[col]=='Immune']['gene name']))
hcl_fibroblast = getClusterByColor(dg_fibroblast, list(f_germlayer_genes[f_germlayer_genes[col]=='Fibroblast']['gene name']))
hcl_endo = getClusterByColor(dg_endo, list(f_germlayer_genes[f_germlayer_genes[col]=='Endothelial']['gene name']))
hcl_neural = getClusterByColor(dg_neural, list(f_germlayer_genes[f_germlayer_genes[col]=='Neural']['gene name']))
hcl_muscle = getClusterByColor(dg_muscle, list(f_germlayer_genes[f_germlayer_genes[col]=='Muscle']['gene name']))

hcl_epithelial = pd.Series({g: hcl  for hcl in hcl_epithelial for g in hcl_epithelial[hcl]})
hcl_immune = pd.Series({g: hcl  for hcl in hcl_immune for g in hcl_immune[hcl]})
hcl_fibroblast = pd.Series({g: hcl  for hcl in hcl_fibroblast for g in hcl_fibroblast[hcl]})
hcl_endo = pd.Series({g: hcl  for hcl in hcl_endo for g in hcl_endo[hcl]})
hcl_neural = pd.Series({g: hcl  for hcl in hcl_neural for g in hcl_neural[hcl]})
hcl_muscle = pd.Series({g: hcl  for hcl in hcl_muscle for g in hcl_muscle[hcl]})

go = list(hcl_epithelial.sort_values().index)[::-1] + list(hcl_immune.sort_values().index) + list(hcl_fibroblast.sort_values().index) + list(hcl_endo.sort_values().index) + list(hcl_neural.sort_values().index) + list(hcl_muscle.sort_values().index)

plt.figure(figsize=(4.5,4))
plt.subplot(151, frameon = False)
d = {'Epithelial': 0, 'Immune': 1, 'Fibroblast': 2, 'Endothelial': 3, 'Neural': 4, 'Muscle': 5}
cnt = Counter()
y0 = 0; z0 = 0
for g in go[::-1]:
    y = 1.0/len(go)
    c = f_germlayer_genes.set_index('gene name').loc[g, col]
    cnt.update([c])
    if cnt[c] == 1:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]])
    else:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]], label='_nolegend_')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2,4)
plt.legend(['Muscle','Neural','Endothelial','Fibroblast','Immune','Epithelial'], loc = 1, bbox_to_anchor = (0,1))

plt.subplot(152)
plt.imshow(scalezscore(mdfHP_av.loc[go]), aspect = 'auto')
plt.ylabel('genes'); plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('HP')

plt.subplot(153)
plt.imshow(scalezscore(mdfBM_av.loc[go]), aspect = 'auto')
plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('BM')

plt.subplot(154)
plt.imshow(scalezscore(mdfVT_av.loc[go]), aspect = 'auto')
plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('VT')

plt.subplot(155)
plt.imshow(scalezscore(mdfDT_av.loc[go]), aspect = 'auto')
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)
plt.yticks([i for i, x in enumerate(go) ], [x.rsplit('_')[1] for i, x in enumerate(go)], fontsize = 7)
plt.xlabel('Sections (Epi -> LP)')
#plt.colorbar(label = 'z-score (scaled)', ticks = [-1,0,1], shrink = 0.5)
plt.title('DT')

plt.savefig(path + 'MarkerGeneHeatmap_average_no_legend_with_genes.pdf', bbox_inches ='tight')
plt.show()


plt.figure(figsize=(6,4))
plt.subplot(161, frameon = False)
d = {'Epithelial': 0, 'Immune': 1, 'Fibroblast': 2, 'Endothelial': 3, 'Neural': 4, 'Muscle': 5}
cnt = Counter()
y0 = 0; z0 = 0
for g in go[::-1]:
    y = 1.0/len(go)
    c = f_germlayer_genes.set_index('gene name').loc[g, col]
    cnt.update([c])
    if cnt[c] == 1:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]])
    else:
        plt.bar([3], [y], bottom = y0, width = 0.3, color = Colors.cc[d[c]], label='_nolegend_')
    y0 += y
plt.yticks([]); plt.xticks([])
plt.ylim(0,1)
plt.xlim(2,4)
plt.legend(['Muscle','Neural','Endothelial','Fibroblast','Immune','Epithelial'], loc = 1, bbox_to_anchor = (0,1))

plt.subplot(162)
plt.imshow(scalezscore(mdfHP_av.loc[go]), aspect = 'auto')
plt.ylabel('genes'); plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('HP')

plt.subplot(163)
plt.imshow(scalezscore(mdfBM_av.loc[go]), aspect = 'auto')
plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('BM')

plt.subplot(164)
plt.imshow(scalezscore(mdfVT_av.loc[go]), aspect = 'auto')
plt.xlabel('Sections (Epi -> LP)')
plt.yticks([])
plt.title('VT')

plt.subplot(165)
plt.imshow(scalezscore(mdfDT_av.loc[go]), aspect = 'auto')
plt.tick_params(axis='y', labelleft=False, labelright=True, length = 5, left = False, right = True)
plt.yticks([])
plt.xlabel('Sections (Epi -> LP)')
#plt.colorbar(label = 'z-score (scaled)', ticks = [-1,0,1], shrink = 0.5)
plt.title('DT')

plt.subplot(166, frameon= False)
plt.colorbar(label = 'z-score (scaled)', ticks = [-1,0,1], shrink = 0.5, location = 'right')
plt.xticks([]); plt.yticks([])

plt.savefig(path + 'MarkerGeneHeatmap_average_with_legend_no_genes.pdf', bbox_inches ='tight')
plt.show()



#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'
#oral tissue 50% brighter colour BM, HP, VT, DT: '#fdac6e','#74e7c5','#bab7d9','#f294c4'

g = 'ENSMUSG00000026728_Vim_ProteinCoding__2'

plt.figure(figsize=(3,3))
plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[g], color = '#1B9E77')
plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[g], color = '#7570B3')
plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[g], color = '#E7298A')
plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[g], color = '#D95F02',linewidth=3)
plt.legend(['BM-av','HP-av','VT-av','DT-av'])
plt.xlabel('Relative sections (Epi -> LP)')
plt.ylabel('Normalised expression (z-score)')
plt.title('Vim')
plt.savefig(path + '/Vim_line_graph_all_av_legend.pdf', bbox_inches = 'tight')
plt.show()

plt.figure(figsize=(2,2))
plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[g], color = '#74e7c5')
plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[g], color = '#bab7d9')
plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[g], color = '#f294c4')
plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[g], color = '#D95F02',linewidth=3)
plt.xlabel('Relative sections (Epi -> LP)')
plt.ylabel('z-score (scaled)')
plt.title('Vim')
plt.savefig(path + '/Vim_line_graph_all_av_no_legend.pdf', bbox_inches = 'tight')
plt.show()

xgs = ['ENSMUSG00000022587_Ly6e_ProteinCoding__15','ENSMUSG00000044678_Ly6k_ProteinCoding__15',
       'ENSMUSG00000063522_Ly6m_ProteinCoding__15','ENSMUSG00000075602_Ly6a_ProteinCoding__15']

for x in xgs:
    plt.figure(figsize=(3, 3))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#E7298A')
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Gene expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Marker_genes_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

xgs = ['ENSMUSG00000000303_Cdh1_ProteinCoding__8','ENSMUSG00000045394_Epcam_ProteinCoding__17',
        'ENSMUSG00000026728_Vim_ProteinCoding__2','ENSMUSG00000045545_Krt14_ProteinCoding__11',
        'ENSMUSG00000038151_Prdm1_ProteinCoding__10','ENSMUSG00000035000_Dpp4_ProteinCoding__2',
       'ENSMUSG00000031871_Cdh5_ProteinCoding__8','ENSMUSG00000020717_Pecam1_ProteinCoding__11',
       'ENSMUSG00000028583_Pdpn_ProteinCoding__4','ENSMUSG00000029231_Pdgfra_ProteinCoding__5',
        'ENSMUSG00000026395_Ptprc_ProteinCoding__1','ENSMUSG00000018830_Myh11_ProteinCoding__16',
       'ENSMUSG00000039542_Ncam1_ProteinCoding__9','ENSMUSG00000056569_Mpz_ProteinCoding__1',
       'ENSMUSG00000017300_Tnnc2_ProteinCoding__2','ENSMUSG00000030787_Lyve1_ProteinCoding__7']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#E7298A')
    #plt.xlabel('Sections (Epi -> LP)')
    #plt.ylabel('Gene expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Marker_genes_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


xgs = ['ENSMUSG00000026728_Vim_ProteinCoding__2','ENSMUSG00000045545_Krt14_ProteinCoding__11',
       'ENSMUSG00000035493_Tgfbi_ProteinCoding__13','ENSMUSG00000028950_Tas1r1_ProteinCoding__4',
       'ENSMUSG00000036585_Fgf1_ProteinCoding__18', 'ENSMUSG00000037225_Fgf2_ProteinCoding__3',
       'ENSMUSG00000027208_Fgf7_ProteinCoding__2', 'ENSMUSG00000021732_Fgf10_ProteinCoding__13',
       'ENSMUSG00000020327_Fgf22_ProteinCoding__10',
       'ENSMUSG00000031565_Fgfr1_ProteinCoding__8', 'ENSMUSG00000030849_Fgfr2_ProteinCoding__7',
       'ENSMUSG00000054252_Fgfr3_ProteinCoding__5', 'ENSMUSG00000005320_Fgfr4_ProteinCoding__13']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_all_av_no_legend_BM_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


