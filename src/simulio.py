import pandas as pd
from graphviz import Digraph
import os



def tree_graphing_topology(nodes, edges):

    tree = Digraph(comment='Phylogenetic tree', node_attr={'shape':'box'})
    for j in range(len(nodes)):
        tree.node(str(nodes[j]), 'Cell ' + str(nodes[j]), style='filled')

    for j in range(len(edges)):
        tree.edge(str(edges[j][0]), str(edges[j][1]))

    tree = tree.unflatten(stagger=2)

    return tree 


def tree_graphing_profiles(node_profiles, edge_profiles, root, n_positions, max_positions, save_as_png=False):

    node_labels, node_profiles = zip(*node_profiles.items())
    node_profiles = [list(tuple(j) for j in i.values()) for i in node_profiles]
    edge_labels, edge_profiles = zip(*edge_profiles.items())
    edge_profiles = [''.join([str(i) + ': ' + str(j) +' \n' for i, j in k.items()]) for k in edge_profiles]

    tree = Digraph(comment='Phylogenetic tree', filename=str(save_as_png), format='png', node_attr={'shape':'box'})

    if n_positions <= max_positions:
        for j in range(len(node_labels)):
            A, B, M = list(zip(*node_profiles[j]))
            tree.node(str(node_labels[j]), 'Clone ' + str(node_labels[j]) + ':\n Allele A: {}\n Allele B: {}\n     SNVs: {}'.format(A, B, M), style='filled')

        for j in range(len(edge_labels)):
            if edge_labels[j][0] == root:
                tree.edge(str(edge_labels[j][0]), str(edge_labels[j][1]), color='darkmagenta')
            else:
                tree.edge(str(edge_labels[j][0]), str(edge_labels[j][1]))

    else:
        for j in range(len(node_labels)):
            tree.node(str(node_labels[j]), 'Clone ' + str(node_labels[j]), style='filled')

        for j in range(len(edge_labels)):
            if edge_labels[j][0] == root:
                tree.edge(str(edge_labels[j][0]), str(edge_labels[j][1]), label=edge_profiles[j], color='darkmagenta')
            else:
                tree.edge(str(edge_labels[j][0]), str(edge_labels[j][1]), label=edge_profiles[j])

    tree = tree.unflatten(stagger=2)

    if save_as_png != False:
        tree.render()

    return tree


def create_edge_snv_df(edges_snv_assignment):
    df = pd.DataFrame(edges_snv_assignment.items(), columns=['Edges', 'SNV'])
    df = df.set_index('Edges')['SNV'].apply(pd.Series).stack().reset_index(level=-1, drop=True).astype(int).reset_index()
    df['Parent'], df['Child'] = zip(*df.Edges)
    df = df[[0, 'Parent', 'Child']].sort_values(by=0)
    df.columns = ['Variant_Position','Parent', 'Child']
    return df
