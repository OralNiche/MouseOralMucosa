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


#Read in all data
mdfHP_av = pd.read_table(path + '/HP_average_zscore.tsv', sep = '\t', index_col=0)
mdfBM_av = pd.read_table(path + '/BM_average_zscore.tsv', sep = '\t', index_col=0)
mdfDT_av = pd.read_table(path + '/DT_average_zscore.tsv', sep = '\t', index_col=0)
mdfVT_av = pd.read_table(path + '/VT_average_zscore.tsv', sep = '\t', index_col=0)

mdfHP = pd.read_table(path + '/HP_zscore_both_libraries.tsv', sep = '\t', index_col=0)
mdfBM = pd.read_table(path + '/BM_zscore_both_libraries.tsv', sep = '\t', index_col=0)
mdfDT = pd.read_table(path + '/DT_zscore_both_libraries.tsv', sep = '\t', index_col=0)
mdfVT = pd.read_table(path + '/VT_zscore_both_libraries.tsv', sep = '\t', index_col=0)

mnfHP_av = pd.read_table(path + '/HP_average_norm.tsv', sep = '\t', index_col=0)
mnfBM_av = pd.read_table(path + '/BM_average_norm.tsv', sep = '\t', index_col=0)
mnfDT_av = pd.read_table(path + '/DT_average_norm.tsv', sep = '\t', index_col=0)
mnfVT_av = pd.read_table(path + '/VT_average_norm.tsv', sep = '\t', index_col=0)


def mergeDfs(dfs):
    mdf = dfs[0]
    for i in range(1, len(dfs)):
       mdf = mdf.merge(dfs[i], how = 'outer', left_index = True, right_index = True)
    mdf = mdf.fillna(0)
    mdf = mdf.loc[mdf.sum(axis=1).sort_values(ascending=False).index]
    return mdf

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

def zscore(df):
    df = df.loc[df.index[df.sum(axis=1)>0]]
    zdf = df.T
    zdf = (zdf-zdf.mean())/zdf.std()
    zdf = zdf.T
    return zdf

def scalezscore(zdf):
    df  = ((zdf.T>=0)*zdf.T/zdf.max(axis=1) + (zdf.T<0)*zdf.T/abs(zdf.min(axis=1))).T
    return df


mdfAll = mergeDfs([mdfHP]+[mdfBM]+[mdfVT]+[mdfDT])

germlayer_genes = pd.read_excel(path + 'MarkerGenes.xlsx', sheet_name='ClassMarkerGenes')
f_germlayer_genes = germlayer_genes.loc[[idx for idx in germlayer_genes.index if germlayer_genes.loc[idx, 'Cell class'] in ['Epithelial','Immune','Fibroblast','Endothelial','Neural','Muscle']]]
f_germlayer_genes.head()
f_germlayer_genes.head()


hcl_germlayer = linkage(mdfAll.loc[f_germlayer_genes['gene name']], method='ward', metric='euclidean')
dg_germlayer = dendrogram(hcl_germlayer, no_labels=True, color_threshold = 25) # 5
plt.show()

hcl_germlayer = getClusterByColor(dg_germlayer, list(f_germlayer_genes['gene name']))

hcl_germlayer = pd.Series({g: hcl  for hcl in hcl_germlayer for g in hcl_germlayer[hcl]})

go = hcl_germlayer.sort_values().index


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
plt.title('DT')

plt.subplot(166, frameon= False)
plt.colorbar(label = 'z-score (scaled)', ticks = [-1,0,1], shrink = 0.5, location = 'right')
plt.xticks([]); plt.yticks([])

plt.savefig(path + 'MarkerGeneHeatmap_average_with_legend_no_genes.pdf', bbox_inches ='tight')
plt.show()

#Normalised z-score expression
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
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Marker_genes_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#Normalised expression
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
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Marker_genes_line_graph_all_norm_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()



#Normalised z-score expression
xgs = ['ENSMUSG00000061527_Krt5_ProteinCoding__15','ENSMUSG00000061527_Krt5_ProteinCoding__15',
       'ENSMUSG00000054146_Krt15_ProteinCoding__11','ENSMUSG00000022510_Trp63_ProteinCoding__16',
       'ENSMUSG00000025064_Col17a1_ProteinCoding__19','ENSMUSG00000020913_Krt24_ProteinCoding__11',
       'ENSMUSG00000038415_Foxq1_ProteinCoding__13','ENSMUSG00000030029_Lrig1_ProteinCoding__6',
       'ENSMUSG00000026739_Bmi1_ProteinCoding__2','ENSMUSG00000025407_Gli1_ProteinCoding__10',
       'ENSMUSG00000027111_Itga6_ProteinCoding__2', 'ENSMUSG00000025809_Itgb1_ProteinCoding__8',
       'ENSMUSG00000020758_Itgb4_ProteinCoding__11','ENSMUSG00000000120_Ngfr_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Basal_epithelial_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#Normalised expression
xgs = ['ENSMUSG00000061527_Krt5_ProteinCoding__15','ENSMUSG00000061527_Krt5_ProteinCoding__15',
       'ENSMUSG00000054146_Krt15_ProteinCoding__11','ENSMUSG00000022510_Trp63_ProteinCoding__16',
       'ENSMUSG00000025064_Col17a1_ProteinCoding__19','ENSMUSG00000020913_Krt24_ProteinCoding__11',
       'ENSMUSG00000038415_Foxq1_ProteinCoding__13','ENSMUSG00000030029_Lrig1_ProteinCoding__6',
       'ENSMUSG00000026739_Bmi1_ProteinCoding__2','ENSMUSG00000025407_Gli1_ProteinCoding__10',
       'ENSMUSG00000027111_Itga6_ProteinCoding__2', 'ENSMUSG00000025809_Itgb1_ProteinCoding__8',
       'ENSMUSG00000020758_Itgb4_ProteinCoding__11','ENSMUSG00000000120_Ngfr_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    #plt.xlabel('Relative sections (Epi -> LP)')
    #plt.ylabel('Normalised expression')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Basal_epithelial_genes_line_graph_all_norm_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


#'ENSMUSG00000063661_Krt73_ProteinCoding__15','ENSMUSG00000035831_Krt25_ProteinCoding__11',
#'ENSMUSG00000055937_Krt28_ProteinCoding__11',

#Normalised z-score expression
xgs = ['ENSMUSG00000006777_Krt23_ProteinCoding__11','ENSMUSG00000023043_Krt18_ProteinCoding__15',
       'ENSMUSG00000044041_Krt13_ProteinCoding__11','ENSMUSG00000020911_Krt19_ProteinCoding__11',
       'ENSMUSG00000061527_Krt5_ProteinCoding__15','ENSMUSG00000049382_Krt8_ProteinCoding__15',
       'ENSMUSG00000053797_Krt16_ProteinCoding__11','ENSMUSG00000020912_Krt12_ProteinCoding__11',
       'ENSMUSG00000059668_Krt4_ProteinCoding__15','ENSMUSG00000037185_Krt80_ProteinCoding__15',
       'ENSMUSG00000061397_Krt79_ProteinCoding__15',
       'ENSMUSG00000048013_Krt35_ProteinCoding__11','ENSMUSG00000048699_Krt90_ProteinCoding__15',
       'ENSMUSG00000067613_Krt83_ProteinCoding__15','ENSMUSG00000035849_Krt222_ProteinCoding__11',
       'ENSMUSG00000044294_Krt84_ProteinCoding__15','ENSMUSG00000035592_Krt33a_ProteinCoding__11',
       'ENSMUSG00000035775_Krt20_ProteinCoding__11','ENSMUSG00000048981_Krt31_ProteinCoding__11',
       'ENSMUSG00000067614_Krt86_ProteinCoding__15',
       'ENSMUSG00000022986_Krt75_ProteinCoding__15','ENSMUSG00000019761_Krt10_ProteinCoding__11',
       'ENSMUSG00000017588_Krt27_ProteinCoding__11',
       'ENSMUSG00000058354_Krt6a_ProteinCoding__15','ENSMUSG00000047641_Krt87_ProteinCoding__15',
       'ENSMUSG00000035557_Krt17_ProteinCoding__11','ENSMUSG00000075402_Krt76_ProteinCoding__15',
       'ENSMUSG00000020916_Krt36_ProteinCoding__11','ENSMUSG00000020913_Krt24_ProteinCoding__11',
       'ENSMUSG00000023039_Krt7_ProteinCoding__15','ENSMUSG00000050463_Krt78_ProteinCoding__15',
       'ENSMUSG00000023041_Krt6b_ProteinCoding__15','ENSMUSG00000054146_Krt15_ProteinCoding__11',
       'ENSMUSG00000046095_Krt32_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Keratin_gene_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#Normalised expression
xgs = ['ENSMUSG00000006777_Krt23_ProteinCoding__11','ENSMUSG00000023043_Krt18_ProteinCoding__15',
       'ENSMUSG00000044041_Krt13_ProteinCoding__11','ENSMUSG00000020911_Krt19_ProteinCoding__11',
       'ENSMUSG00000061527_Krt5_ProteinCoding__15','ENSMUSG00000049382_Krt8_ProteinCoding__15',
       'ENSMUSG00000053797_Krt16_ProteinCoding__11','ENSMUSG00000020912_Krt12_ProteinCoding__11',
       'ENSMUSG00000059668_Krt4_ProteinCoding__15','ENSMUSG00000037185_Krt80_ProteinCoding__15',
       'ENSMUSG00000061397_Krt79_ProteinCoding__15',
       'ENSMUSG00000048013_Krt35_ProteinCoding__11','ENSMUSG00000048699_Krt90_ProteinCoding__15',
       'ENSMUSG00000067613_Krt83_ProteinCoding__15','ENSMUSG00000035849_Krt222_ProteinCoding__11',
       'ENSMUSG00000044294_Krt84_ProteinCoding__15','ENSMUSG00000035592_Krt33a_ProteinCoding__11',
       'ENSMUSG00000035775_Krt20_ProteinCoding__11','ENSMUSG00000048981_Krt31_ProteinCoding__11',
       'ENSMUSG00000067614_Krt86_ProteinCoding__15',
       'ENSMUSG00000022986_Krt75_ProteinCoding__15','ENSMUSG00000019761_Krt10_ProteinCoding__11',
       'ENSMUSG00000017588_Krt27_ProteinCoding__11',
       'ENSMUSG00000058354_Krt6a_ProteinCoding__15','ENSMUSG00000047641_Krt87_ProteinCoding__15',
       'ENSMUSG00000035557_Krt17_ProteinCoding__11','ENSMUSG00000075402_Krt76_ProteinCoding__15',
       'ENSMUSG00000020916_Krt36_ProteinCoding__11','ENSMUSG00000020913_Krt24_ProteinCoding__11',
       'ENSMUSG00000023039_Krt7_ProteinCoding__15','ENSMUSG00000050463_Krt78_ProteinCoding__15',
       'ENSMUSG00000023041_Krt6b_ProteinCoding__15','ENSMUSG00000054146_Krt15_ProteinCoding__11',
       'ENSMUSG00000046095_Krt32_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    #plt.xlabel('Relative sections (Epi -> LP)')
    #plt.ylabel('Normalised expression')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Keratin_genes_line_graph_all_norm_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#BM: 'ENSMUSG00000032093_Cd3e_ProteinCoding__9','ENSMUSG00000002033_Cd3g_ProteinCoding__9','ENSMUSG00000034783_Cd207_ProteinCoding__6',
#VT: 'ENSMUSG00000003379_Cd79a_ProteinCoding__7',
#DT: 'ENSMUSG00000028859_Csf3r_ProteinCoding__4',
#HP: 'ENSMUSG00000001865_Cpa3_ProteinCoding__31',

#Normalised z-score expression
xgs = ['ENSMUSG00000032094_Cd3d_ProteinCoding__9','ENSMUSG00000000409_Lck_ProteinCoding__4',
       'ENSMUSG00000024610_Cd74_ProteinCoding__18','ENSMUSG00000024677_Ms4a6b_ProteinCoding__19',
       'ENSMUSG00000018774_Cd68_ProteinCoding__11','ENSMUSG00000024397_Aif1_ProteinCoding__17',
       'ENSMUSG00000030787_Lyve1_ProteinCoding__7',
       'ENSMUSG00000005947_Itgae_ProteinCoding__11', 'ENSMUSG00000030789_Itgax_ProteinCoding__7',
    'ENSMUSG00000024673_Ms4a1_ProteinCoding__19',
       'ENSMUSG00000037922_Bank1_ProteinCoding__3', 'ENSMUSG00000014030_Pax5_ProteinCoding__4',
       'ENSMUSG00000056054_S100a8_ProteinCoding__3', 'ENSMUSG00000056071_S100a9_ProteinCoding__3',
        'ENSMUSG00000030403_Vasp_ProteinCoding__7',
        'ENSMUSG00000015053_Gata2_ProteinCoding__6',
       'ENSMUSG00000024680_Ms4a2_ProteinCoding__19', 'ENSMUSG00000005672_Kit_ProteinCoding__5']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Immune_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#Normalised expression
xgs = ['ENSMUSG00000032094_Cd3d_ProteinCoding__9','ENSMUSG00000000409_Lck_ProteinCoding__4',
       'ENSMUSG00000024610_Cd74_ProteinCoding__18','ENSMUSG00000024677_Ms4a6b_ProteinCoding__19',
       'ENSMUSG00000018774_Cd68_ProteinCoding__11','ENSMUSG00000024397_Aif1_ProteinCoding__17',
       'ENSMUSG00000030787_Lyve1_ProteinCoding__7',
       'ENSMUSG00000005947_Itgae_ProteinCoding__11', 'ENSMUSG00000030789_Itgax_ProteinCoding__7',
       'ENSMUSG00000024673_Ms4a1_ProteinCoding__19',
       'ENSMUSG00000037922_Bank1_ProteinCoding__3', 'ENSMUSG00000014030_Pax5_ProteinCoding__4',
       'ENSMUSG00000056054_S100a8_ProteinCoding__3', 'ENSMUSG00000056071_S100a9_ProteinCoding__3',
        'ENSMUSG00000030403_Vasp_ProteinCoding__7',
       'ENSMUSG00000015053_Gata2_ProteinCoding__6',
       'ENSMUSG00000024680_Ms4a2_ProteinCoding__19', 'ENSMUSG00000005672_Kit_ProteinCoding__5']


for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    #plt.xlabel('Relative sections (Epi -> LP)')
    #plt.ylabel('Normalised expression')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Immune_genes_line_graph_all_norm_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()




#Normalised z-score expression
xgs = ['ENSMUSG00000031004_Mki67_ProteinCoding__7','ENSMUSG00000027342_Pcna_ProteinCoding__2',
       'ENSMUSG00000019942_Cdk1_ProteinCoding__10','ENSMUSG00000020914_Top2a_ProteinCoding__11',
       'ENSMUSG00000027496_Aurka_ProteinCoding__2','ENSMUSG00000070348_Ccnd1_ProteinCoding__7']

for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#E7298A')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Proliferation_line_graph_all_av_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#Normalised expression
xgs = ['ENSMUSG00000031004_Mki67_ProteinCoding__7','ENSMUSG00000027342_Pcna_ProteinCoding__2',
       'ENSMUSG00000019942_Cdk1_ProteinCoding__10','ENSMUSG00000020914_Top2a_ProteinCoding__11',
       'ENSMUSG00000027496_Aurka_ProteinCoding__2','ENSMUSG00000070348_Ccnd1_ProteinCoding__7']


for x in xgs:
    plt.figure(figsize=(1, 1))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A')
    #plt.xlabel('Relative sections (Epi -> LP)')
    #plt.ylabel('Normalised expression')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Proliferation_genes_line_graph_all_norm_no_legend_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()