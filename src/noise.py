import numpy as np
import pandas as pd
from scipy.stats import binom

def add_noise(seq_df, n_noise, n_noise_clones, n_positions, noise_type='random', edges_df=None):
    variant_df = seq_df.copy()
    mutations = variant_df['Variant_Position'].unique()
    noise_mutations = np.random.choice(mutations, n_noise, replace=False)
    variant_df['Noise'] = variant_df['Variant_Position'].isin(noise_mutations) 
    max_clone = variant_df[['Parent', 'Child']].max().max()
    noisy_clones = [i + max_clone for i in range(1, n_noise_clones+1)]
    ccfs = {clone : {sample: np.random.randint(1, 100)/100 for sample in variant_df['Sample'].unique()} for clone in noisy_clones}
    noise_clusters_assign = dict(zip(variant_df['Variant_Position'].unique(), np.random.choice(noisy_clones, n_positions)))
    variant_df['Child'] = variant_df.apply(lambda row: row['Child'] if not row['Noise'] else noise_clusters_assign[row['Variant_Position']], axis=1)
    variant_df['CCF'] = variant_df.apply(lambda row: row['CCF'] if not row['Noise'] else ccfs[row['Child']][row['Sample']], axis=1)
    variant_df['Fract_Variant'] = variant_df.apply(lambda row: row['Fract_Variant'] if not row['Noise'] else np.maximum(row['Fract_Variant'], 1), axis=1)
    variant_df['True_VAF'] = np.minimum(1.0, (variant_df['Fract_Variant'] * variant_df['CCF'] * variant_df['Purity']) / variant_df['Total_Fract_CN'])
    
    def test(row):
        return row['Variant_Reads'] if not row['Noise'] else binom.rvs(row['Total_Reads'], p=row['True_VAF'])

    variant_df['Variant_Reads'] = variant_df.apply(test, axis=1)
    variant_df['Est_VAF'] = variant_df['Variant_Reads']/variant_df['Total_Reads']
    return variant_df, noisy_clones


def get_tree_edges(tree):

    edges = list(tree[['Parent', 'Child']].apply(tuple, axis=1))
    parents, children = np.array(list(zip(*edges)))
    leaves = children[~np.isin(children, parents)]

    descendants = {edge[-1]: get_descendants(edge, edges, leaves, [])[1:] for edge in edges}
    return descendants


def get_descendants(edge, edges, leaves, descendants=[]):
    if edge[-1] in leaves:
        descendants.append(edge[-1])
        return descendants
    else:
        descendants.append(edge[-1])
        edge_children = [i for i in edges if i[0] == edge[1]]
        for j in edge_children:
            descendants = get_descendants(j, edges, leaves, descendants)
        return descendants