#!/usr/bin/env python3
import sys, os
import pandas as pd
import numpy as np
import glob
import sc_aautils as scaa
import Colors
import matplotlib
import matplotlib.pyplot as plt

# Default sans-serif font is HELVETICA
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
#Always use sans-serif fonts
matplotlib.rcParams['font.family'] = "sans-serif"


path= '/Users/kaikretzschmar/Dropbox/_TomoSeq/'


files = glob.glob('corr_genes_*tsv')

dfs = {f.rsplit('.tsv')[0].rsplit('_')[-1]: pd.read_table(f, sep = '\t', index_col = 0) for f in files}

for i, k in enumerate(dfs):
    dfs[k].columns = [col + '_' + k for col in dfs[k].columns]
    if i == 0:
        mdf = dfs[k].copy()
    else:
        mdf = mdf.merge(dfs[k], how = 'outer', left_index = True, right_index = True)

fig, ax = plt.subplots()
for k in dfs:
    ax.hist(mdf['adj-pv_'+k], bins = 100, label = k, alpha = 0.1)

smdf = mdf[['adj-pv_'+k for k in dfs]]
fsmdf = smdf.fillna(1)

Z, dg = scaa.hierarchicalClustering(fsmdf)

smdf = smdf.loc[fsmdf.index[dg['leaves']]]

fig, axs = plt.subplots(ncols = 2, figsize = (3*1.6, 2*3), width_ratios=[1, 5])
ax = axs[-1]
im = ax.imshow(smdf, aspect = 'auto', interpolation = 'nearest')
cbc = fig.colorbar(im, ax = ax, label = 'adj. P-value', shrink = 0.25)
ax.set_xticks([0,1,2,3])
ax.set_xticklabels([c.rsplit('_')[-1] for c in smdf.columns], rotation = 90)
ax.set_yticklabels('')
#fig, ax = plt.subplots()
ax = axs[0]
im = ax.imshow((smdf<0.05).astype(int), aspect = 'auto', interpolation = 'nearest', cmap = 'Reds')
#cbc = fig.colorbar(im, ax = ax, label = 'Reproducible')
ax.set_xticks([0,1,2,3])
ax.set_xticklabels([c.rsplit('_')[-1] for c in smdf.columns], rotation = 90)
#ax.axis('off')
ax.set_ylabel('genes')
fig.savefig('./heatmap_adj_pvals.pdf', bbox_inches = 'tight')

fig, axs = plt.subplots(nrows = 4, figsize = (3*1.6, 2*3))
for i, (ax, k) in enumerate(zip(axs, dfs)):
    col = 'adj-pv_'+k
    ax.hist(smdf[col], bins = 100, label = k, color = Colors.colors[i])
    ax.legend()
    ax.axvline(0.05, ls = '--', c = 'k')
    ax.set_ylabel('frequency')
ax.set_xlabel('adj. P-value')
fig.savefig('./histo_adj_pvals.pdf', bbox_inches = 'tight')