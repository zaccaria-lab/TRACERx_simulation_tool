import pandas as pd
import numpy as np


from simulio import *



def split(leaves):
    leaves = np.array(leaves)
    n_leaves = len(leaves)
    n_ones = np.random.randint(1, n_leaves)
    n_zeros = n_leaves - n_ones
    mask = np.array(n_ones * [1] + n_zeros * [0], dtype=bool)
    np.random.shuffle(mask)
    left_leaves = leaves[mask].tolist()
    right_leaves = leaves[~mask].tolist()
    return left_leaves, right_leaves


def make_tree(leaves):
    parents = [max(leaves)] 
    tree_edges, parents, current_parent = make_tree_internal(leaves, parents)
    root = parents[-1] + 1
    initial_cc = parents[1]
    tree_edges.append((root, initial_cc))
    tree_nodes = leaves + parents[1:]
    tree_nodes.append(root)
    return tree_nodes, tree_edges, root, initial_cc


def make_tree_internal(leaves, parents): # ABI: this function is only working because you provide the first parent which is a number higher than all leaves. We will need to improve it.
    if len(leaves) == 1:
        return [], parents, leaves[0]     
    else:
        current_parent = max(parents) + 1 if len(parents) > 0 else 0 # ABI: is this line ever used?
        parents.append(current_parent)
        left_leaves, right_leaves = split(leaves)
        tree, parents, parent_left = make_tree_internal(left_leaves, parents)
        tree_right, parents, parent_right = make_tree_internal(right_leaves, parents)
        tree.extend(tree_right)
        tree.extend([(current_parent, parent_left), (current_parent, parent_right)])
        return tree, parents, current_parent


def prune_tree(nodes, edges, n_nodes):
    parents, children = np.array(list(zip(*edges)))
    leaves = np.array([i for i in children if i not in parents])

    if len(nodes) == n_nodes+1:
        return nodes, edges

    else:
        remove = np.random.choice(leaves, 1)[0]
        nodes.remove(remove)
        remove_idx = np.where(children == remove)[0][0]
        children = np.delete(children, remove_idx)
        parents = np.delete(parents, remove_idx)
        edges = list(zip(parents, children))
        nodes, edges = prune_tree(nodes, edges, n_nodes)
        return nodes, edges


def simulate_topology(n_nodes):

    nodes = list(np.arange(n_nodes)+1)
    n, e, root, initial_cc = make_tree(nodes)

    nodes, edges = prune_tree(n, e, n_nodes)

    return nodes, edges, root, initial_cc
