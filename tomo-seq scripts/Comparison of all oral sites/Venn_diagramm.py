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

#Read in correlated genes from all oral sites

#BM data
corrGenes_BM = pd.read_table(path + '/corr_genes_BM.tsv', sep = '\t', index_col=0)
#corrGenes_BM = pd.DataFrame(corrGenes_BM)

pvth = 0.05

gl_BM_dict = {k: corrGenes_BM[k][corrGenes_BM['adj-pv'] < pvth] for k in corrGenes_BM}
gl_BM = []
for k in gl_BM_dict:
    gl_BM  += list(gl_BM_dict[k].index)

from collections import Counter

gl_BM_cnt = Counter(gl_BM)
gl_BM = [g for g in gl_BM_cnt]

#HP data
corrGenes_HP = pd.read_table(path + '/corr_genes_HP.tsv', sep = '\t', index_col=0)
#corrGenes_HP = pd.DataFrame(corrGenes_HP)

pvth = 0.05

gl_HP_dict = {k: corrGenes_HP[k][corrGenes_HP['adj-pv'] < pvth] for k in corrGenes_HP}
gl_HP = []
for k in gl_HP_dict:
    gl_HP  += list(gl_HP_dict[k].index)

from collections import Counter

gl_HP_cnt = Counter(gl_HP)
gl_HP = [g for g in gl_HP_cnt]

#VT data
corrGenes_VT = pd.read_table(path + '/corr_genes_VT.tsv', sep = '\t', index_col=0)
#corrGenes_VT = pd.DataFrame(corrGenes_VT)

pvth = 0.05

gl_VT_dict = {k: corrGenes_VT[k][corrGenes_VT['adj-pv'] < pvth] for k in corrGenes_VT}
gl_VT = []
for k in gl_VT_dict:
    gl_VT  += list(gl_VT_dict[k].index)

from collections import Counter

gl_VT_cnt = Counter(gl_VT)
gl_VT = [g for g in gl_VT_cnt]

#DT data
corrGenes_DT = pd.read_table(path + '/corr_genes_DT.tsv', sep = '\t', index_col=0)
#corrGenes_DT = pd.DataFrame(corrGenes_DT)

pvth = 0.05

gl_DT_dict = {k: corrGenes_DT[k][corrGenes_DT['adj-pv'] < pvth] for k in corrGenes_DT}
gl_DT = []
for k in gl_DT_dict:
    gl_DT  += list(gl_DT_dict[k].index)

from collections import Counter

gl_DT_cnt = Counter(gl_DT)
gl_DT = [g for g in gl_DT_cnt]


gl_all = list(set(list(gl_BM) + list(gl_HP) + list(gl_DT) + list(gl_VT)))

import venn

labels = venn.get_labels([set(gl_BM),
                       set(gl_HP),
                     set(gl_DT),
                      set(gl_VT)],
                        fill=['number', 'percent'])

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

fig, ax = venn.venn4(labels, names = ['BM','HP','DT','VT'], colors = oralcolors)

plt.savefig(path + '/venn_genes_p_val_0.05_all.pdf', bbox_inches='tight')
plt.show()
