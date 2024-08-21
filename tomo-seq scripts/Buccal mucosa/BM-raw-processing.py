#Basic settings

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# Default sans-serif font is HELVETICA
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
#Always use sans-serif fonts
matplotlib.rcParams['font.family'] = "sans-serif"


#Important functions
def readData(filename):
    dfc = pd.read_csv(filename + '.coutc.tsv', sep = ',', index_col=0)
    dfb = pd.read_csv(filename + '.coutb.tsv', sep = ',', index_col=0)
    dft = pd.read_csv(filename + '.coutt.tsv', sep = ',', index_col=0)
    dfc.columns = [c.zfill(3) for c in dfc.columns]
    dfb.columns = [c.zfill(3) for c in dfb.columns]
    dft.columns = [c.zfill(3) for c in dft.columns]
    dfc = dfc[sorted(dfc.columns)]
    dfb = dfb[sorted(dfb.columns)]
    dft = dft[sorted(dft.columns)]
    dfc = dfc[dfc.columns[:20]]
    dfb = dfb[dfb.columns[:20]]
    dft = dft[dft.columns[:20]]
    return dfc, dfb, dft

def getPlotTotalReads(df):
    mt = [idx for idx in df.index if '_chrM' in idx]
    if len(mt) == 0:
        mt = [idx for idx in df.index if '__MT' in idx]
    genes = [idx for idx in df.index if idx and idx not in mt]
    plt.bar(range(len(df.columns)), df.loc[genes].sum()+1, log = 'y')
    plt.bar(range(len(df.columns)), df.loc[mt].sum()+1, log = 'y')
    plt.xticks(range(len(df.columns)), df.columns)
    plt.ylabel('total')
    plt.legend(['genes', 'mt'])
    return plt

def getCBT(dfs):
    dfc, dfb, dft = dfs
    fig = plt.figure(figsize=(15,4))
    plt.subplot(131)
    plt1 = getPlotTotalReads(dfc)
    plt1.title('counts')
    plt.subplot(132)
    plt2 = getPlotTotalReads(dfb)
    plt2.title('barcodes')
    plt.subplot(133)
    plt3 = getPlotTotalReads(dft)
    plt3.title('transcripts')
    return plt.show()

def splitGst(dfs, col):
    return [df[col] for df in dfs]

#Read data
f = 'TomoSeqBM1_total'
df1s = readData(f)
f = 'TomoSeqBM2_total'
df2s = readData(f)

getCBT(df1s)
getCBT(df2s)

print([df.shape for dfs in [df1s, df2s] for df in dfs])

def removeGenes(df):
    genes = ['ERCC','__MT']
    exclude = [idx for idx in df.index for g in genes if g in idx]
    df = df.loc[df.index[df.sum(axis=1)>0]]
    return df.loc[[idx for idx in df.index if idx not in exclude]]

def removeGenesDFS(dfs):
    return [removeGenes(df) for df in dfs]

rdfs1 = removeGenesDFS(df1s)
rdfs2 = removeGenesDFS(df2s)

print(rdfs1[-1].columns[rdfs1[-1].sum()>1e4])
print(rdfs2[-1].columns[rdfs2[-1].sum()>1e4])

dfs_g1 = rdfs1
dfs_g2 = rdfs2

#Pathway for saving data
path= '/Documents/TomoSeq_project/Tomo/'


#QC fpr uniques transcripts and genes per replicate
plt.figure(figsize=(15,5))
titles = ['BM replicate 1','BM replicate 2']
thresholdstranscripts = [1e4,1e4]
thresholdsgenes = [5e3,5e3]

for i, df in enumerate([dfs_g1,dfs_g2]):
    plt.subplot(2,4,i+1)
    plt.bar(range(len(df[-1].columns)), df[-1].sum()/1e6, color = '#D95F02')
    plt.bar(range(len(df[-1].columns)), ((df[-1]>0).sum()< thresholdstranscripts[i])*(df[-1].sum()/1e6), color = 'gray')

    if i == 0:
        plt.ylabel('Unique transcripts (x10^6)')
    plt.title(titles[i])

for i, df in enumerate([dfs_g1,dfs_g2]):
    plt.subplot(2, 4, i + 5)
    plt.bar(range(len(df[-1].columns)), (df[-1] > 0).sum()/1e3, color = '#D95F02')
    plt.bar(range(len(df[-1].columns)), ((df[-1] > 0).sum() < thresholdsgenes[i]) * (df[-1] > 0).sum() / 1e3,
                color='gray')

    if i == 0:
        plt.ylabel('Genes (x10^3)')
    plt.xlabel('Sections (Epi->LP)')

plt.savefig(path + 'QC_BM.pdf', bbox_inches = 'tight')
plt.show()

writer = pd.ExcelWriter(path + '/QC_BM.xlsx')
for i, df in enumerate([dfs_g1,dfs_g2]):
    (df[-1]>0).sum().to_excel(writer, sheet_name = 'num_genes, BM '+str(i+1))
    df[-1].sum().to_excel(writer, sheet_name = 'num_trans, BM ' + str(i+1))
writer.close()

#Select sections based on thresholds for unique transcripts
fdfs_g1 = splitGst(dfs_g1, dfs_g1[-1].columns[(dfs_g1[-1]>0).sum()>thresholdstranscripts[0]])
fdfs_g2 = splitGst(dfs_g2, dfs_g2[-1].columns[(dfs_g2[-1]>0).sum()>thresholdstranscripts[1]])

#Normalisation of data
nfdfs = [df.sum().mean()*df/df.sum() for df in [fdfs_g1[-1], fdfs_g2[-1]]]

#Plot maker genes to check orientation
g = 'ENSMUSG00000000303_Cdh1_ProteinCoding__8' # [idx for idx in dfu_g1[0].index if 'BRCA1' in idx]
plt.figure(figsize=(14,3))
plt.subplot(141)
plt.bar(range(len(nfdfs[0].columns)), nfdfs[0].loc[g], color = '#D95F02'); plt.title('BM1, Cdh1')
plt.subplot(142)
plt.bar(range(len(nfdfs[1].columns)), nfdfs[1].loc[g], color = '#D95F02'); plt.title('BM2, Cdh1')
plt.savefig(path + 'Cdh1_BM.pdf', bbox_inches = 'tight')
plt.show()

g = 'ENSMUSG00000045545_Krt14_ProteinCoding__11' # [idx for idx in dfu_g1[0].index if 'BRCA1' in idx]
plt.figure(figsize=(14,3))
plt.subplot(141)
plt.bar(range(len(nfdfs[0].columns)), nfdfs[0].loc[g], color = '#D95F02'); plt.title('BM1, Krt14')
plt.subplot(142)
plt.bar(range(len(nfdfs[1].columns)), nfdfs[1].loc[g], color = '#D95F02'); plt.title('BM2, Krt14')
plt.savefig(path + 'Krt14_BM.pdf', bbox_inches = 'tight')
plt.show()

g = 'ENSMUSG00000026728_Vim_ProteinCoding__2' # [idx for idx in dfu_g1[0].index if 'BRCA1' in idx]
plt.figure(figsize=(14,3))
plt.subplot(141)
plt.bar(range(len(nfdfs[0].columns)), nfdfs[0].loc[g], color = '#D95F02'); plt.title('BM1, Vim')
plt.subplot(142)
plt.bar(range(len(nfdfs[1].columns)), nfdfs[1].loc[g], color = '#D95F02'); plt.title('BM2, Vim')
plt.savefig(path + 'Vim_BM.pdf', bbox_inches = 'tight')
plt.show()

g = 'ENSMUSG00000018830_Myh11_ProteinCoding__16' # [idx for idx in dfu_g1[0].index if 'BRCA1' in idx]
plt.figure(figsize=(14,3))
plt.subplot(141)
plt.bar(range(len(nfdfs[0].columns)), nfdfs[0].loc[g], color = '#D95F02'); plt.title('BM1, Myh11')
plt.subplot(142)
plt.bar(range(len(nfdfs[1].columns)), nfdfs[1].loc[g], color = '#D95F02'); plt.title('BM2, Myh11')
plt.savefig(path + 'Myh11_BM.pdf', bbox_inches = 'tight')
plt.show()

#Store data
def savedfs(dfs, name):
    dfc, dfb, dft = dfs
    dfc.to_csv(name + '.coutc.tsv', sep = '\t')
    dfb.to_csv(name + '.coutb.tsv', sep = '\t')
    dft.to_csv(name + '.coutt.tsv', sep = '\t')
    return

savedfs(fdfs_g1, path + '/BM1')
savedfs(fdfs_g2, path + '/BM2')

