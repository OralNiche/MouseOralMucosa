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


#Read in all data
mdfHP_av = pd.read_table(path + '/HP_average_zscore.tsv', sep = '\t', index_col=0)
mdfBM_av = pd.read_table(path + '/BM_average_zscore.tsv', sep = '\t', index_col=0)
mdfDT_av = pd.read_table(path + '/DT_average_zscore.tsv', sep = '\t', index_col=0)
mdfVT_av = pd.read_table(path + '/VT_average_zscore.tsv', sep = '\t', index_col=0)

mnfHP_av = pd.read_table(path + '/HP_average_norm.tsv', sep = '\t', index_col=0)
mnfBM_av = pd.read_table(path + '/BM_average_norm.tsv', sep = '\t', index_col=0)
mnfDT_av = pd.read_table(path + '/DT_average_norm.tsv', sep = '\t', index_col=0)
mnfVT_av = pd.read_table(path + '/VT_average_norm.tsv', sep = '\t', index_col=0)

#oral tissue normal colour BM, HP, VT, DT: '#D95F02','#1B9E77','#7570B3','#E7298A'
#oral tissue 50% brighter colour BM, HP, VT, DT: '#fdac6e','#74e7c5','#bab7d9','#f294c4'


#BM-specific spatially reproducible genes
xgs = ['ENSMUSG00000035557_Krt17_ProteinCoding__11','ENSMUSG00000057967_Fgf18_ProteinCoding__11',
       'ENSMUSG00000021242_Npc2_ProteinCoding__12','ENSMUSG00000074743_Thbd_ProteinCoding__2',
       'ENSMUSG00000026141_Col19a1_ProteinCoding__1','ENSMUSG00000040452_Cdh12_ProteinCoding__15',
       'ENSMUSG00000034177_Rnf43_ProteinCoding__11','ENSMUSG00000020891_Alox8_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_BM_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


#HP-specific spatially reproducible genes
xgs = ['ENSMUSG00000026586_Prrx1_ProteinCoding__1','ENSMUSG00000075602_Ly6a_ProteinCoding__15',
       'ENSMUSG00000002602_Axl_ProteinCoding__7','ENSMUSG00000025809_Itgb1_ProteinCoding__8',
       'ENSMUSG00000074344_Tmigd3_ProteinCoding__3','ENSMUSG00000000562_Adora3_ProteinCoding__3']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#1B9E77', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_HP_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#VT-specific spatially reproducible genes
xgs = ['ENSMUSG00000005087_Cd44_ProteinCoding__2','ENSMUSG00000037279_Ovol2_ProteinCoding__2',
       'ENSMUSG00000059668_Krt4_ProteinCoding__15','ENSMUSG00000070047_Fat1_ProteinCoding__8',
       'ENSMUSG00000003032_Klf4_ProteinCoding__4','ENSMUSG00000008999_Bmp7_ProteinCoding__2',
       'ENSMUSG00000041415_Dicer1_ProteinCoding__12','ENSMUSG00000023041_Krt6b_ProteinCoding__15',
       'ENSMUSG00000031757_Mt4_ProteinCoding__8','ENSMUSG00000054252_Fgfr3_ProteinCoding__5',
       'ENSMUSG00000001497_Pax9_ProteinCoding__12','ENSMUSG00000058354_Krt6a_ProteinCoding__15',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000021732_Fgf10_ProteinCoding__13',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000032125_Robo4_ProteinCoding__9',
       'ENSMUSG00000015468_Notch4_ProteinCoding__17','ENSMUSG00000027496_Aurka_ProteinCoding__2',
       'ENSMUSG00000020105_Lrig3_ProteinCoding__10','ENSMUSG00000000078_Klf6_ProteinCoding__13',
       'ENSMUSG00000070047_Fat1_ProteinCoding__8','ENSMUSG00000068196_Col8a1_ProteinCoding__16',
       'ENSMUSG00000013663_Pten_ProteinCoding__19','ENSMUSG00000029999_Tgfa_ProteinCoding__6',
       'ENSMUSG00000055333_Fat2_ProteinCoding__11','ENSMUSG00000027878_Notch2_ProteinCoding__3']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_VT_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


#DT-specific spatially reproducible genes
xgs = ['ENSMUSG00000030029_Lrig1_ProteinCoding__6','ENSMUSG00000032179_Bmp5_ProteinCoding__9',
       'ENSMUSG00000020916_Krt36_ProteinCoding__11','ENSMUSG00000053797_Krt16_ProteinCoding__11',
       'ENSMUSG00000049382_Krt8_ProteinCoding__15','ENSMUSG00000042284_Itga1_ProteinCoding__13',
       'ENSMUSG00000062352_Itgb1bp1_ProteinCoding__12','ENSMUSG00000032179_Bmp5_ProteinCoding__9',
       'ENSMUSG00000000126_Wnt9a_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mdfHP_av.columns)), mdfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mdfBM_av.columns)), mdfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0, 1, len(mdfVT_av.columns)), mdfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mdfDT_av.columns)), mdfDT_av.loc[x], color='#E7298A', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (z-score)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_line_graph_DT_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()



#BM-specific spatially reproducible genes
xgs = ['ENSMUSG00000035557_Krt17_ProteinCoding__11','ENSMUSG00000057967_Fgf18_ProteinCoding__11',
       'ENSMUSG00000021242_Npc2_ProteinCoding__12','ENSMUSG00000074743_Thbd_ProteinCoding__2',
       'ENSMUSG00000026141_Col19a1_ProteinCoding__1','ENSMUSG00000040452_Cdh12_ProteinCoding__15',
       'ENSMUSG00000034177_Rnf43_ProteinCoding__11','ENSMUSG00000020891_Alox8_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#D95F02', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (transcripts)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_norm_graph_BM_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


#HP-specific spatially reproducible genes
xgs = ['ENSMUSG00000026586_Prrx1_ProteinCoding__1','ENSMUSG00000075602_Ly6a_ProteinCoding__15',
       'ENSMUSG00000002602_Axl_ProteinCoding__7','ENSMUSG00000025809_Itgb1_ProteinCoding__8',
       'ENSMUSG00000074344_Tmigd3_ProteinCoding__3','ENSMUSG00000000562_Adora3_ProteinCoding__3']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#1B9E77', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (transcripts)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_norm_graph_HP_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()

#VT-specific spatially reproducible genes
xgs = ['ENSMUSG00000005087_Cd44_ProteinCoding__2','ENSMUSG00000037279_Ovol2_ProteinCoding__2',
       'ENSMUSG00000059668_Krt4_ProteinCoding__15','ENSMUSG00000070047_Fat1_ProteinCoding__8',
       'ENSMUSG00000003032_Klf4_ProteinCoding__4','ENSMUSG00000008999_Bmp7_ProteinCoding__2',
       'ENSMUSG00000041415_Dicer1_ProteinCoding__12','ENSMUSG00000023041_Krt6b_ProteinCoding__15',
       'ENSMUSG00000031757_Mt4_ProteinCoding__8','ENSMUSG00000054252_Fgfr3_ProteinCoding__5',
       'ENSMUSG00000001497_Pax9_ProteinCoding__12','ENSMUSG00000058354_Krt6a_ProteinCoding__15',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000021732_Fgf10_ProteinCoding__13',
       'ENSMUSG00000030849_Fgfr2_ProteinCoding__7','ENSMUSG00000032125_Robo4_ProteinCoding__9',
       'ENSMUSG00000015468_Notch4_ProteinCoding__17','ENSMUSG00000027496_Aurka_ProteinCoding__2',
       'ENSMUSG00000020105_Lrig3_ProteinCoding__10','ENSMUSG00000000078_Klf6_ProteinCoding__13',
       'ENSMUSG00000070047_Fat1_ProteinCoding__8','ENSMUSG00000068196_Col8a1_ProteinCoding__16',
       'ENSMUSG00000013663_Pten_ProteinCoding__19','ENSMUSG00000029999_Tgfa_ProteinCoding__6',
       'ENSMUSG00000055333_Fat2_ProteinCoding__11','ENSMUSG00000027878_Notch2_ProteinCoding__3']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#f294c4')
    plt.plot(np.linspace(0,1,len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#7570B3', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (transcripts)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_norm_graph_VT_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()


#DT-specific spatially reproducible genes
xgs = ['ENSMUSG00000030029_Lrig1_ProteinCoding__6','ENSMUSG00000032179_Bmp5_ProteinCoding__9',
       'ENSMUSG00000020916_Krt36_ProteinCoding__11','ENSMUSG00000053797_Krt16_ProteinCoding__11',
       'ENSMUSG00000049382_Krt8_ProteinCoding__15','ENSMUSG00000042284_Itga1_ProteinCoding__13',
       'ENSMUSG00000062352_Itgb1bp1_ProteinCoding__12','ENSMUSG00000032179_Bmp5_ProteinCoding__9',
       'ENSMUSG00000000126_Wnt9a_ProteinCoding__11']

for x in xgs:
    plt.figure(figsize=(2, 2))
    plt.plot(np.linspace(0,1,len(mnfHP_av.columns)), mnfHP_av.loc[x], color='#74e7c5')
    plt.plot(np.linspace(0,1,len(mnfBM_av.columns)), mnfBM_av.loc[x], color='#fdac6e')
    plt.plot(np.linspace(0, 1, len(mnfVT_av.columns)), mnfVT_av.loc[x], color='#bab7d9')
    plt.plot(np.linspace(0,1,len(mnfDT_av.columns)), mnfDT_av.loc[x], color='#E7298A', linewidth=3)
    plt.xlabel('Relative sections (Epi -> LP)')
    plt.ylabel('Normalised expression (transcripts)')
    plt.title(x.rsplit('_')[1])
    plt.savefig(path + '/Unique_genes_norm_graph_DT_vs_All_' + x.rsplit('_')[1] + '.pdf', bbox_inches='tight')
plt.show()
