import pandas as pd
import numpy as np
import os
from IPython.display import display

from sequencing import *
from topology import *
from evolution import *
from simulio import *
from conversion import *
from noise import *

def tx_simulation(patient_name, 
                n_nodes,
                n_positions,
                n_samples,
                min_clones_per_sample,
                max_clones_per_sample,
                max_positions=25,
                prop_truncal_snvs=0.5,
                prop_pos_gained=0.5,
                prop_pos_lost=0.3,
                prop_gained_truncal=1,
                prop_lost_truncal=1,
                clonal_wgd=False,
                subclonal_wgd=False,
                n_subclonal_wgd=1,
                minimum_branch_length=1,
                constant_multiplicity=True, 
                prop_mutloss_positions=0,
                gain_events=None,
                loss_events=None,
                mutloss_truncal=False,
                all_leaf_observed=True,
                coverage=400,
                min_tumour_purity=0.7,
                max_tumour_purity=0.9,
                fract_noise=0,
                n_noise_clones=0,
                noise_type='random',
                rs=None,
                display_plots=False,
                output_directory=None,
                ccf_vaf_plots=False,
                save_output=True,
                display_output=False):

    if save_output == True:
        assert output_directory is not None, 'Please specify an output_directory'
        save_dir = os.path.join(output_directory, patient_name, 'sim')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)


        save_as_png = os.path.join(save_dir, patient_name + '_tree')
    else:
        save_as_png = False

    n_noise = round(n_positions * fract_noise)
    n_pos = n_positions + n_noise



    variant_df, clone_df, edges_df, edge_position_profiles, prop_mutation_loss, warnings = make_a_simulation(n_nodes, n_pos, n_samples, min_clones_per_sample, max_clones_per_sample,
                                                                    max_positions, prop_truncal_snvs, prop_pos_gained, prop_pos_lost, prop_gained_truncal,
                                                                    prop_lost_truncal, clonal_wgd, subclonal_wgd, n_subclonal_wgd, minimum_branch_length, constant_multiplicity, prop_mutloss_positions,
                                                                    gain_events, loss_events, mutloss_truncal, all_leaf_observed, coverage, min_tumour_purity, max_tumour_purity, rs, display_plots, 
                                                                    save_as_png=save_as_png)

    if fract_noise > 0:
        assert n_noise_clones > 0, 'If noisy clusters are added to the simulation, the n_noise_clones must be > 0'
        variant_df, noise_clones = add_noise(variant_df, n_noise, n_noise_clones, n_pos, noise_type, edges_df)
    else:
        noise_clones = None
        variant_df['Noise'] = False

    variant_df, tx_input_df = convert_tree_building_format(variant_df, patient_name) 
    edge_position_events_df = variant_df[['VAR_POS', 'CHR', 'POS']].drop_duplicates().merge(edge_position_profiles, on='VAR_POS')

    missing_mut = variant_df[['VAR_POS', 'TRUE_CCF']].groupby('VAR_POS')['TRUE_CCF'].sum()
    prop_missing_mut = 0 if len(missing_mut[missing_mut == 0])/n_positions < 0.01 else len(missing_mut[missing_mut == 0])/n_positions
    if prop_missing_mut > 0:
        warning = 'WARNING: Due to sampling of clones and mutation removal, {} of mutations are not present at all in samples'.format(prop_missing_mut)
        warnings.append(warning)

    patient_info_df = make_patient_info_df([patient_name, n_nodes, n_positions, n_samples, min_clones_per_sample,
                                    max_clones_per_sample, prop_truncal_snvs, prop_pos_gained,
                                    prop_pos_lost, prop_gained_truncal, prop_lost_truncal, clonal_wgd, subclonal_wgd, gain_events, loss_events, mutloss_truncal,
                                    subclonal_wgd*n_subclonal_wgd, constant_multiplicity, prop_mutation_loss, min_tumour_purity, max_tumour_purity, 
                                    fract_noise, noise_clones, noise_type, rs, warnings])                                                         

    if ccf_vaf_plots == True:
        if save_output == True:
            vaf_save = os.path.join(save_dir, patient_name + '_estvaf.png')
            ccf_save = os.path.join(save_dir, patient_name + '_ccfphylo.png')
        else:
            vaf_save = False
            ccf_save = False
        
        plot_variant_data(variant_df, 'VAF_OBS', save_as=vaf_save, display=display_plots)
        plot_variant_data(variant_df, 'CCF_OBS', save_as=ccf_save, display=display_plots)

    if save_output == True:
        variant_df.to_csv(os.path.join(save_dir, patient_name + '.tsv'), index=False, sep='\t')
        tx_input_df.to_csv(os.path.join(save_dir, patient_name + '_TX.tsv'), index=False, sep='\t')
        clone_df.to_csv(os.path.join(save_dir, patient_name + '_cloneprops.tsv'), index=False, sep='\t')
        edges_df.to_csv(os.path.join(save_dir, patient_name + '_treeedges.tsv'), index=False, sep='\t')
        patient_info_df.to_csv(os.path.join(save_dir, patient_name + '_info.tsv'), index=False, sep='\t')
        edge_position_events_df.to_csv(os.path.join(save_dir, patient_name + '_events.tsv'), index=False, sep='\t')

    if display_output == True:
        return patient_info_df, variant_df, clone_df, edges_df, tx_input_df, edge_position_events_df



def make_patient_info_df(data):

    info = ['patient_name', 'n_nodes', 'n_positions', 'n_samples', 'min_clones_per_sample',
               'max_clones_per_sample', 'prop_truncal_snvs', 'prop_pos_gained', 'prop_pos_lost',
               'prop_gained_truncal', 'prop_lost_truncal', 'clonal_wgd', 'subclonal_wgd', 'gain_events', 'loss_events', 'mutloss_truncal',
               'n_subclonal_wgd', 'constant_multiplicity', 'prop_mut_lost','min_tumour_purity', 'max_tumour_purity', 
               'fract_noise', 'noise_clones','noise_type', 'random_seed', 'simulation_warnings']

    patient_info = pd.DataFrame(list(zip(info, data)))
    return patient_info
    


def make_a_simulation(n_nodes,
                    n_positions,
                    n_samples,
                    min_clones_per_sample,
                    max_clones_per_sample,
                    max_positions=25,
                    prop_truncal_snvs = 0.6,
                    prop_pos_gained = 0.2,
                    prop_pos_lost = 0.3,
                    prop_gained_truncal = 1,
                    prop_lost_truncal = 1,
                    clonal_wgd=False,
                    subclonal_wgd=False,
                    n_subclonal_wgd=1,
                    minimum_branch_length=1,
                    constant_multiplicity=True, 
                    prop_mutloss_positions=0,
                    gain_events=None,
                    loss_events=None,
                    mutloss_truncal=False,
                    all_leaf_observed=True,
                    coverage=400,
                    min_tumour_purity=0.7,
                    max_tumour_purity=0.9,
                    rs=None,
                    display=False,
                    save_as_png=False):
    
    np.random.seed(rs)
    
    clone_profiles, dcf_profiles, edges_df, snv_edge_df, edge_position_profiles, diploid, leaves, prop_mutation_loss, warnings = simulation(n_nodes, n_positions, max_positions, prop_truncal_snvs, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, clonal_wgd, subclonal_wgd, n_subclonal_wgd, minimum_branch_length, constant_multiplicity, prop_mutloss_positions, gain_events, loss_events, mutloss_truncal, display, save_as_png)
    variant_df, clone_prop_df = make_sequencing_data(clone_profiles, dcf_profiles, n_samples, min_clones_per_sample, max_clones_per_sample, diploid, all_leaf_observed, leaves, coverage, min_tumour_purity, max_tumour_purity)
    variant_df = variant_df.merge(snv_edge_df, on='Variant_Position')
    
    return variant_df, clone_prop_df, edges_df, edge_position_profiles, prop_mutation_loss, warnings


def simulation(n_nodes, n_positions, max_positions, prop_truncal_snvs, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, clonal_wgd=False, subclonal_wgd=False, n_subclonal_wgd=1, minimum_branch_length=1, constant_multiplicity=True, prop_mutloss_positions=0, gain_events=None, loss_events=None, mutloss_truncal=False, display_tree=False, save_as_png=False):
    
    nodes, edges, root, initial_cc = simulate_topology(n_nodes)
    parents, children = np.array(list(zip(*edges)))
    leaves = children[~np.isin(children, parents)]

    total_profiles, dcf_profiles, snv_edge_assign, edge_label_profiles, profile_events, prop_mutation_loss, warnings = tree_profiles(nodes, edges, leaves, root, initial_cc, n_positions, prop_truncal_snvs, minimum_branch_length, prop_pos_gained, prop_pos_lost, prop_gained_truncal, prop_lost_truncal, clonal_wgd, subclonal_wgd, n_subclonal_wgd, constant_multiplicity, prop_mutloss_positions, gain_events, loss_events, mutloss_truncal)
    if display_tree == True:
        display(tree_graphing_profiles(total_profiles, edge_label_profiles, root, n_positions, max_positions, save_as_png))
    else:
        tree_graphing_profiles(total_profiles, edge_label_profiles, root, n_positions, max_positions, save_as_png) 

    edges_df = pd.DataFrame(edge_label_profiles).T.reset_index().fillna(0).astype(int)

    if clonal_wgd == False and (subclonal_wgd == False or n_subclonal_wgd == 0):
        edges_df['WGD'] = 0
    edges_df.columns = ['Parent', 'Child', 'CNA_gains', 'CNA_losses', 'SNVs', 'WGD']
    edges_df['Truncal'] = 'No'
    edges_df.loc[edges_df['Parent'] == root, 'Truncal'] = 'Yes'
    snv_edge_assign_df = create_edge_snv_df(snv_edge_assign)

    edge_events = {str([i for i in edges if i[1] == k][0]): v for k, v in profile_events.items()}
    edge_position_profiles = pd.DataFrame(edge_events).reset_index().rename(columns={'index': 'VAR_POS'})

    return total_profiles, dcf_profiles, edges_df, snv_edge_assign_df, edge_position_profiles, root, leaves, prop_mutation_loss, warnings 