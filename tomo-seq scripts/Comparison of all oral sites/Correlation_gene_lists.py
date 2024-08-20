#!/usr/bin/env python3
#import sys, os
import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
import matplotlib
from itertools import combinations
from itertools import product
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt

# Default sans-serif font is HELVETICA
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
#Always use sans-serif fonts
matplotlib.rcParams['font.family'] = "sans-serif"


path= '/Users/kaikretzschmar/Dropbox/_TomoSeq/'

# Load the gene data for each sample
HP_df = pd.read_table(path + 'clusters_HP_reprGenes1.tsv', sep = '\t')
BM_df = pd.read_table(path + 'clusters_BM_reprGenes1.tsv', sep = '\t')
VT_df = pd.read_table(path + 'clusters_VT_reprGenes1.tsv', sep = '\t')
DT_df = pd.read_table(path + 'clusters_DT_reprGenes1.tsv', sep = '\t')


def group_genes_by_list(df):
    gene_sets = {}
    for list_number in df['hcl'].unique():
        gene_sets[list_number] = set(df[df['hcl'] == list_number]['Gene'])
    return gene_sets

# Process the data for each sample
hp_gene_sets = group_genes_by_list(HP_df)
bm_gene_sets = group_genes_by_list(BM_df)
vt_gene_sets = group_genes_by_list(VT_df)
dt_gene_sets = group_genes_by_list(DT_df)

# Combine all samples into a single dictionary
all_samples = {
    "HP": hp_gene_sets,
    "BM": bm_gene_sets,
    "VT": vt_gene_sets,
    "DT": dt_gene_sets
}

# Function to calculate Jaccard similarity
def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union != 0 else 0

# Create a list of all list names for labeling
all_list_names = [f'{sample} cl {list_num}' for sample, lists in all_samples.items() for list_num in lists]

# Initialize an empty similarity matrix
num_lists = len(all_list_names)
similarity_matrix = np.zeros((num_lists, num_lists))

# Fill the similarity matrix with Jaccard similarities
list_index_map = [(sample, list_num) for sample, gene_sets in all_samples.items() for list_num in gene_sets]

for i, (sample1, list1) in enumerate(list_index_map):
    for j, (sample2, list2) in enumerate(list_index_map):
        similarity_matrix[i, j] = jaccard_similarity(all_samples[sample1][list1], all_samples[sample2][list2])


np.fill_diagonal(similarity_matrix, np.nan)

# Convert the matrix to a DataFrame for easier labeling
similarity_df = pd.DataFrame(similarity_matrix, index=all_list_names, columns=all_list_names)

similarity_df.to_csv(path + '/Jaccard_correlation.tsv', sep = '\t')

# Plot the heatmap using seaborn
plt.figure(figsize=(10, 8))
sns.heatmap(similarity_df, annot=similarity_df.round(1), fmt='.1f', cmap="viridis", linewidths=0.5)
plt.title('Jaccard Similarity Heatmap of Gene Lists')
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.savefig(path + '/Jaccard_heatmap.pdf', bbox_inches = 'tight')
plt.show()


# Replace diagonal NaN values with 0 to avoid clustering issues
np.fill_diagonal(similarity_matrix, 0)

# Plot the clustermap using seaborn with the viridis color palette
sns.clustermap(similarity_df, annot=similarity_df.round(1), fmt='.1f', cmap="viridis", linewidths=0.5, mask=similarity_df.isnull(), figsize=(14, 12))
plt.title('Jaccard Similarity Clustermap of Gene Lists', y=1.05)
plt.show()