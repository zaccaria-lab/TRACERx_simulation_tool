import numpy as np


def tree_profiles(nodes, edges, leaves, root, initial_cc, n_positions, prop_truncal_snvs, minimum_branch_length, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, clonal_wgd, subclonal_wgd, n_subclonal_wgd=1, constant_multiplicity=True, prop_mutloss_positions=0, gain_events=None, loss_events=None, mutloss_truncal=False):
    positions = np.arange(n_positions)
    
    overall_snvs, snv_edge_labels, snv_edge_assign = snv_assignment(edges, positions, prop_truncal_snvs, minimum_branch_length)
    overall_cnas, subclone_nodes, cna_edge_labels, prop_mutation_loss, warnings = cna_assignment(positions, nodes, edges, leaves, root, initial_cc, prop_truncal_snvs, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, constant_multiplicity, snv_edge_assign, prop_mutloss_positions, gain_events, loss_events, mutloss_truncal)
    overall_wgd, wgd_edge_labels = wgd_assignment(positions, edges, clonal_wgd, subclonal_wgd, n_subclonal_wgd, initial_cc, subclone_nodes)
    l = (lambda e_dict, key : e_dict[key] if key in e_dict.keys() else {})
    profile_edge_labels = {i: l(cna_edge_labels, i) | l(snv_edge_labels, i) | l(wgd_edge_labels, i) for i in edges}

    profile_events = concat_all_events(overall_cnas, overall_snvs, overall_wgd)

    event_functions = generate_events()
    total_profiles = mutation_node_assignment(n_positions, profile_events, nodes, initial_cc, edges, event_functions)
    dcf_profiles = dcf_mutation_node_assignment(n_positions, profile_events, nodes, initial_cc, edges, event_functions)

    return total_profiles, dcf_profiles, snv_edge_assign, profile_edge_labels, profile_events, prop_mutation_loss, warnings 


def conditioned_choice(elements, minimum_distance, n_partitions, chosen_partition_indexes=[]):
    if n_partitions == 0:
        return chosen_partition_indexes
    else:
        available = [i for i in elements if len(chosen_partition_indexes) == 0 or min(abs(i-j) for j in chosen_partition_indexes) >= minimum_distance]
        if len(available) > 0:
            choice = np.random.choice(available, 1)[0]
            chosen_partition_indexes.append(choice)
            return conditioned_choice(available, minimum_distance, n_partitions-1, chosen_partition_indexes)
        else:
            return chosen_partition_indexes


def partition(edges, subclonal, minimum_distance):
    if minimum_distance == 1:
        assert len(subclonal) >= len(edges), 'You have more subclonal node branches than subclonal SNV events'
        partitions = np.sort(np.random.choice(np.arange(1, len(subclonal)-1), len(edges)-2, replace=False))
        subclonal_partitions = np.split(subclonal, partitions)
        assert all(len(i) >= 1 for i in subclonal_partitions), 'Empty partitions'
    else:
        partitions = np.sort(conditioned_choice(subclonal, minimum_distance, (len(edges)-2), []))
        subclonal_partitions = np.split(subclonal, partitions)
    assert np.array_equal(subclonal, np.concatenate(subclonal_partitions, axis=0))
    return subclonal_partitions


def snv_assignment(edges, mutations, prop_truncal, minimum_distance=1):
    subclonal = np.random.choice(mutations, int(len(mutations)*(1-prop_truncal)), replace=False)
    truncal = mutations[~np.isin(mutations, subclonal)]
    subclonal_partitions = partition(edges, subclonal, minimum_distance)
    if len(subclonal_partitions) < (len(edges)-1):
        assert minimum_distance > 1
        for i in range((len(edges)-1) - len(subclonal_partitions)):
            subclonal_partitions.append(np.array([]))
    subclonal_partitions.append(truncal)
    assert len(subclonal_partitions) == len(edges)
    edge_labels = dict(zip(edges, subclonal_partitions))
    overall_mutations = {i[1]: {p : 'snv' for p in j} for i, j in edge_labels.items()}
    snv_edge_labels = {i: {'SNVs': len(j)} for i, j in edge_labels.items()}
    return overall_mutations, snv_edge_labels, edge_labels


def create_event_dict(positions, n_samples, events=None):
    dict_idx = np.random.choice(positions, n_samples, replace=False)
    remaining = positions[~np.isin(positions, dict_idx)]
    if events is None:
        return dict_idx, remaining
    else:
        event_dict = {i : np.random.choice(events) for i in dict_idx}
        remaining_dict = {i : np.random.choice(events) for i in remaining}
        return event_dict, remaining_dict


def wgd_assignment(mutations, edges, clonal_wgd=False, subclonal_wgd=False, n_subclonal_wgd=1, initial_cc=None, subclone_nodes=None):

    wgd = dict()
    if clonal_wgd == True:
        assert initial_cc is not None, 'Please input truncal cancer cell (initial_cc)'
        wgd[initial_cc] = {i: 'wgd' for i in mutations}
    if subclonal_wgd == True:
        assert subclone_nodes is not None, 'Please input subclonal_nodes list (subclone_nodes)'
        wgd_subclones = np.random.choice(subclone_nodes, n_subclonal_wgd)
        for j in wgd_subclones:
            wgd[j] = {i: 'wgd' for i in mutations}
    wgd_edge_labels = {l: {'WGD': 1} for i, j in wgd.items() for l in edges if l[1] == i}

    return wgd, wgd_edge_labels


def concat_all_events(cna_events, snv_events, wgd_events):

    node_keys = (lambda dict, key: set(dict[key].keys()) if key in dict else set([]))
    node_values = (lambda dict, key, pos: [dict[key][pos]] if key in dict and pos in dict[key] else [])
    node_addition = (lambda cna, snv, wgd, node, pos: node_values(cna, node, pos) + node_values(snv, node, pos) + node_values(wgd, node, pos))
    node_union = (lambda cna, snv, wgd: set(cna.keys()).union(set(snv.keys()).union(set(wgd.keys()))))
    node_key_union = (lambda cna, snv, wgd, node: node_keys(cna, node).union(node_keys(snv, node).union(node_keys(wgd, node))))

    profile_events = {node : {pos : node_addition(cna_events, snv_events, wgd_events, node, pos) for pos in node_key_union(cna_events, snv_events, wgd_events, node)} for node in node_union(cna_events, snv_events, wgd_events)}

    for node, mutations in profile_events.items():
        for pos, events in mutations.items():
            np.random.shuffle(events)

    return profile_events


def apply_events(profile, profile_events, event_functions):
    for i in profile_events:
        profile = event_functions[i](profile)
    return profile


def mutation_node_assignment(n_positions, profile_events, nodes, initial_cc, edges, event_functions):
    profiles = {node : {pos: np.array([1, 1, 0]) for pos in range(n_positions)} for node in nodes}

    def assign_total_profiles(profiles, initial_cc, profile_events):
        profiles = set_profiles(profiles, initial_cc, profile_events[initial_cc])
        node_children = [c for p, c in edges if p == initial_cc]
        for node_child in node_children:
            profiles = assign_total_profiles(profiles, node_child, profile_events)
        return profiles

    def set_profiles(profiles, parent, node_events):
        for i in node_events.keys():
            profiles[parent][i] = apply_events(profiles[parent][i], node_events[i], event_functions)
        children = [c for p, c in edges if p == parent]
        for child in children:
            profiles = set_profiles(profiles, child, node_events)
        return profiles

    return assign_total_profiles(profiles, initial_cc, profile_events)


def dcf_apply_events(profile, profile_events, event_functions):
    for i in profile_events:
        if i in ['loss_A', 'cnloh_A']:
            i = 'dcf_' + i
        profile = event_functions[i](profile)
    return profile


def dcf_mutation_node_assignment(n_positions, profile_events, nodes, initial_cc, edges, event_functions):
    profiles = {node : {pos: np.array([1, 1, 0]) for pos in range(n_positions)} for node in nodes}

    def assign_total_profiles(profiles, initial_cc, profile_events):
        profiles = set_profiles(profiles, initial_cc, profile_events[initial_cc])
        node_children = [c for p, c in edges if p == initial_cc]
        for node_child in node_children:
            profiles = assign_total_profiles(profiles, node_child, profile_events)
        return profiles

    def set_profiles(profiles, parent, node_events):
        for i in node_events.keys():
            profiles[parent][i] = dcf_apply_events(profiles[parent][i], node_events[i], event_functions)
        children = [c for p, c in edges if p == parent]
        for child in children:
            profiles = set_profiles(profiles, child, node_events)
        return profiles

    return assign_total_profiles(profiles, initial_cc, profile_events)


def generate_events():

    cn_mut = (lambda x : 1 if x[-1] > 0 else 0)
    snv_cn = (lambda x: 1 if x[0] != 0 else 0 if x[-1] == 0 else -1)
    wgd = (lambda x: x * 2)
    gain_A = (lambda x : x + np.array([1, 0, cn_mut(x)]))
    gain_B = (lambda x : x + np.array([0, 1, 0]))
    gain_both = (lambda x : x + np.array([1, 1, cn_mut(x)]))
    snv = (lambda x: x + np.array([0, 0, snv_cn(x)]))
    loss_B = (lambda x: x - np.array([0, 1, 0]))
    cnloh_B = (lambda x: x + np.array([1, -1, cn_mut(x)]))
    cnloh_A = (lambda x: x + np.array([-1, 1, -cn_mut(x)]))
    loss_A = (lambda x: x - np.array([1, 0, cn_mut(x)]))
    dcf_loss_A = (lambda x: x - np.array([1, 0, 0]))
    dcf_cnloh_A = (lambda x: x + np.array([-1, 1, 0]))

    events = {'wgd': wgd, 'gain_A': gain_A, 'gain_B': gain_B, 'gain_both': gain_both, 'snv' : snv, 'loss_B' : loss_B, 'loss_A': loss_A, 'cnloh_B' : cnloh_B, 'cnloh_A': cnloh_A, 'dcf_loss_A' : dcf_loss_A, 'dcf_cnloh_A': dcf_cnloh_A}
    return events


def cna_assignment(positions, nodes, edges, leaves, root, initial_cc, prop_snv_truncal, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, constant_multiplicity=True, snv_edge_assign=None, prop_mutloss_positions=0, gain_events=None, loss_events=None, mutloss_truncal=False):

    if np.round((prop_pos_lost + prop_pos_gained),2) > 1:
        prop_pos_gained = 1 - prop_pos_lost

    assert np.round((prop_pos_lost + prop_pos_gained),2) <= 1, 'Prop pos gained + lost > 1, = {}'.format(np.round((prop_pos_lost + prop_pos_gained),2))

    if gain_events is None:
        gain_events = ['gain_A', 'gain_B', 'gain_both']
    if loss_events is None:
        loss_events = ['loss_B', 'cnloh_B']

    warnings = []

    subclone_nodes = np.array([i for i in nodes if i != root and i != initial_cc])

    n_truncal_gain, n_truncal_loss, n_subclonal_gain, n_subclonal_loss, n_mutloss_positions, warnings = get_npos_cnas(positions, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, prop_mutloss_positions, warnings)
    
    root_mutations = [pos for edge, pos in snv_edge_assign.items() if edge[0] == root][0]
    subclonal_mutations = positions[~np.isin(positions, root_mutations)]
    mutloss_cnas_assign, root_mutations, subclonal_mutations, prop_mutloss, warnings = get_mutloss(positions, subclone_nodes, root_mutations, subclonal_mutations, n_mutloss_positions, snv_edge_assign, edges, root, leaves, warnings, mutloss_truncal)
    overall_cnas, warnings = assign_cnas(positions, root_mutations, subclonal_mutations, subclone_nodes, snv_edge_assign, n_truncal_gain, n_truncal_loss, n_subclonal_gain, n_subclonal_loss, mutloss_cnas_assign, edges, leaves, root, initial_cc, gain_events, loss_events, warnings, constant_multiplicity)

    cna_edge_labels = {l: {'CNA gains': len([q for p, q in j.items() if q in gain_events]), 'CNA losses': len([q for p, q in j.items() if q in loss_events + ['loss_A', 'cnloh_A']])} for i, j in overall_cnas.items() for l in edges if l[1] == i}
    prop_mutation_loss = np.round(sum({l: len([q for p, q in j.items() if q in ['loss_A', 'cnloh_A']]) for i, j in overall_cnas.items() for l in edges if l[1] == i}.values())/len(positions),2)

    return overall_cnas, subclone_nodes, cna_edge_labels, prop_mutation_loss, warnings


def get_npos_cnas(positions, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, prop_mutloss_positions, warnings):

    possible_subclone_losses = prop_pos_lost*(1-prop_lost_truncal)

    if possible_subclone_losses < prop_mutloss_positions:
        prop_lost_truncal = np.max([0, np.round(1 - (prop_mutloss_positions/prop_pos_lost),2)])
        use_prop_mutloss_positions = prop_pos_lost*(1-prop_lost_truncal)
        warning = 'WARNING: Due to % mutation losses specified, the proportion of truncal losses is capped to {}'.format(prop_lost_truncal)
        print(warning)
        warnings.append(warning)
    else:
        use_prop_mutloss_positions = prop_mutloss_positions

    if prop_mutloss_positions > use_prop_mutloss_positions:
        warning = 'WARNING: Due to % losses, the proportion of mutations lost is capped to {}'.format(use_prop_mutloss_positions)
        print(warning)
        warnings.append(warning)

    n_mutloss_positions = round(len(positions)*use_prop_mutloss_positions)

    n_truncal_gain = round((len(positions)*prop_pos_gained*prop_gained_truncal))
    n_subclonal_gain = round(len(positions)*prop_pos_gained*(1-prop_gained_truncal))
    n_truncal_loss = round(len(positions)*prop_pos_lost*prop_lost_truncal)
    n_subclonal_loss = round(len(positions)*prop_pos_lost*(1-prop_lost_truncal)) - n_mutloss_positions

    if sum([n_truncal_gain, n_subclonal_gain, n_truncal_loss, n_subclonal_loss]) > len(positions):
        if n_truncal_loss > 0:
            n_truncal_loss = n_truncal_loss - 1
        elif n_truncal_gain > 0:
            n_truncal_gain = n_truncal_gain - 1
        else:
            n_subclonal_gain = n_subclonal_gain - 1

    return n_truncal_gain, n_truncal_loss, n_subclonal_gain, n_subclonal_loss, n_mutloss_positions, warnings


def get_mutloss(positions, subclone_nodes, root_mutations, subclonal_mutations, n_mutloss_positions, snv_edge_assign, edges, root, leaves, warnings, mutloss_truncal=False):
    if n_mutloss_positions > 0:
        leaf_positions = [i for edge, pos in snv_edge_assign.items() for i in pos if edge[1] in leaves]
        possible_loss_positions = positions[~np.isin(positions, leaf_positions)]
        possible_subclonal_loss_positions = possible_loss_positions[~np.isin(possible_loss_positions, root_mutations)]

        if mutloss_truncal == True:
            if (n_mutloss_positions <= len(root_mutations)):
                mutloss_positions = np.random.choice(root_mutations, n_mutloss_positions, replace=False)
            else:
                if len(possible_subclonal_loss_positions) >= n_mutloss_positions-len(root_mutations):
                    mutloss_positions = np.concatenate([root_mutations, np.random.choice(possible_subclonal_loss_positions, (n_mutloss_positions-len(root_mutations)), replace=False)])
                else:
                    warning = 'WARNING: Due to % truncal/subclonal snvs, the proportion of mutations lost is capped to {}'.format(np.round(len(possible_loss_positions)/len(positions),2))
                    print(warning)
                    warnings.append(warning)
                    mutloss_positions = possible_loss_positions           
        else:
            use_loss_positions = possible_loss_positions
            if len(use_loss_positions) >= n_mutloss_positions:
                mutloss_positions = np.random.choice(use_loss_positions, n_mutloss_positions, replace=False)
            else:
                warning = 'WARNING: Due to % truncal/subclonal snvs, the proportion of mutations lost is capped to {}'.format(np.round(len(use_loss_positions)/len(positions),2))
                print(warning)
                warnings.append(warning)
                mutloss_positions = use_loss_positions
        mutloss_cnas = {i: np.random.choice(['loss_A', 'cnloh_A'], 1)[0] for i in mutloss_positions}
        mutloss_nodes = only_mutloss_cnas(subclone_nodes, mutloss_positions, snv_edge_assign, edges, leaves, root)
        mutloss_cnas_assign = {i: dict(list(mutloss_cnas.items())[j] for j in list(np.where(np.array(mutloss_nodes) == i)[0])) for i in mutloss_nodes}
        root_mutations = root_mutations[~np.isin(root_mutations, mutloss_positions)]
        subclonal_mutations = subclonal_mutations[~np.isin(subclonal_mutations, mutloss_positions)]
        prop_mutloss = len(mutloss_positions)/len(positions)
        return mutloss_cnas_assign, root_mutations, subclonal_mutations, prop_mutloss, warnings
    else:
        return {}, root_mutations, subclonal_mutations, 0, warnings

    
def only_mutloss_cnas(subclone_nodes, mutloss_positions, snv_edge_assignment, edges, leaves, root):
    descendants = {edge: get_descendants(edge, edges, leaves, [])[1:] for edge in edges}
    snv_edge_assignment = {i:j for i, j in snv_edge_assignment.items()}
    mutloss_cna_edges = {pos: np.random.choice(subclone_nodes[np.isin(subclone_nodes, descendants[edge])]) for edge in snv_edge_assignment.keys() for pos in snv_edge_assignment[edge] if pos in mutloss_positions}
    mutloss_edges = [mutloss_cna_edges[pos] for pos in mutloss_positions if pos in mutloss_cna_edges.keys()]
    return mutloss_edges 


def assign_cnas(positions, root_mutations, subclonal_mutations, subclone_nodes, snv_edge_assign, n_truncal_gain, n_truncal_loss, n_subclonal_gain, n_subclonal_loss, mutloss_cnas_assign, edges, leaves, root, initial_cc, gain_events, loss_events, warnings, constant_multiplicity=True):
    if constant_multiplicity == True:
        if (n_subclonal_gain + n_subclonal_loss) > len(subclonal_mutations):
            warning = 'WARNING: The % subclonal CNAs > % subclonal SNVs and the constant_multiplicity is set to True. The proportion of subclonal CNAs is capped to {}'.format(np.round(len(subclonal_mutations)/len(positions),2))
            print(warning)
            warnings.append(warning)
            diff = len(subclonal_mutations) / (n_subclonal_gain + n_subclonal_loss)
            n_subclonal_gain = round(n_subclonal_gain * diff)
            n_subclonal_loss = np.minimum(round(n_subclonal_loss * diff), (len(subclonal_mutations) - n_subclonal_gain))
        subclonal_gains, subclonal_remaining = create_event_dict(subclonal_mutations, n_subclonal_gain, gain_events)
        subclonal_losses, remaining = create_event_dict(np.array(list(subclonal_remaining.keys())), n_subclonal_loss, loss_events)

        possible_truncal_positions = np.concatenate([root_mutations, list(remaining.keys())])
        possible_truncal_positions

        if (n_truncal_gain + n_truncal_loss) > len(possible_truncal_positions):
            warning = 'WARNING: Capping truncal CNAs!! SEE ERROR'
            print(warning)
            warnings.append(warning)
            diff = len(possible_truncal_positions) / (n_truncal_gain + n_truncal_loss)
            n_truncal_loss = round(n_truncal_loss * diff)
            n_truncal_gain = np.minimum(round(n_truncal_gain * diff), (len(possible_truncal_positions) - n_truncal_loss))

        truncal_gains, truncal_remaining = create_event_dict(possible_truncal_positions, n_truncal_gain, gain_events)
        truncal_losses, remaining = create_event_dict(np.array(list(truncal_remaining.keys())), n_truncal_loss, loss_events)

        truncal_cnas = truncal_gains | truncal_losses
        subclonal_cnas = subclonal_gains | subclonal_losses
        subclonal_cnas_pos = list(subclonal_cnas.keys())
        subclonal_nodes = prevent_subclonal_cnas(subclone_nodes, subclonal_cnas_pos, snv_edge_assign, edges, leaves, root)

    else:
        use_positions = np.concatenate([root_mutations, subclonal_mutations])

        if (n_truncal_gain + n_truncal_loss + n_subclonal_gain + n_subclonal_loss) <= len(use_positions):
            truncal_gains, remaining = create_event_dict(use_positions, n_truncal_gain, gain_events)
            truncal_losses, remaining = create_event_dict(np.array(list(remaining.keys())), n_truncal_loss, loss_events)
            subclonal_gains, remaining = create_event_dict(np.array(list(remaining.keys())), n_subclonal_gain, gain_events)
            subclonal_losses, remaining = create_event_dict(np.array(list(remaining.keys())), n_subclonal_loss, loss_events)
        else:
            print('ERROR! ABI FIX!!')

        truncal_cnas = truncal_gains | truncal_losses
        subclonal_cnas = subclonal_gains | subclonal_losses
        subclonal_cnas_pos = list(subclonal_cnas.keys())
        subclonal_nodes = subclonal_cnas_assign = np.random.choice(subclone_nodes, len(subclonal_cnas_pos))

    subclonal_cnas_assign = {int(i): dict(list(subclonal_cnas.items())[j] for j in list(np.where(np.array(subclonal_nodes) == i)[0])) for i in subclonal_nodes}
    node_present = (lambda nodedict, i: nodedict[i] if i in nodedict.keys() else {})
    cna_nodes = np.union1d(list(subclonal_cnas_assign.keys()), list(mutloss_cnas_assign.keys()))
    overall_cnas = {int(i): node_present(subclonal_cnas_assign, i) | node_present(mutloss_cnas_assign, i) for i in cna_nodes}
    overall_cnas[initial_cc] = truncal_cnas

    return overall_cnas, warnings


def prevent_subclonal_cnas(subclone_nodes, subclonal_cnas, snv_edge_assignment, edges, leaves, root):
    descendants = {edge: get_descendants(edge, edges, leaves, [])[1:] for edge in edges}
    snv_edge_assignment = {i:j for i, j in snv_edge_assignment.items() if i[0] != root}
    subclonal_cna_edges = {pos: np.random.choice(subclone_nodes[~np.isin(subclone_nodes, descendants[edge])]) for edge in snv_edge_assignment.keys() for pos in snv_edge_assignment[edge] if pos in subclonal_cnas}
    cna_edges = [subclonal_cna_edges[pos] for pos in subclonal_cnas if pos in subclonal_cna_edges.keys()]
    return cna_edges

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

def create_event_dict(positions, n_samples, events=None):
    dict_idx = np.random.choice(positions, n_samples, replace=False)
    remaining = positions[~np.isin(positions, dict_idx)]
    if events is None:
        return dict_idx, remaining
    else:
        event_dict = {i : np.random.choice(events) for i in dict_idx}
        remaining_dict = {i : np.random.choice(events) for i in remaining}
        return event_dict, remaining_dict