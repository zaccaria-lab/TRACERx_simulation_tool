import pandas as pd
import numpy as np
import sys

from simulate import *
from conversion import *


def make_tx_simulations(n_sims_per_sample_group=50, sample_group_boundaries=(3,7), output_dir=None):

    data = pd.read_csv('../data/simulated_TRACERx_data.tsv', sep='\t')
    
    if output_dir is None:
        save_output = False
        display_output = True
    else:
        save_output = True
        display_output = False

    if sample_group_boundaries is None:
        groups = [data]
        node_groups = [(8, 30)]
    else:
        samples_low = sample_group_boundaries[0]
        samples_med = sample_group_boundaries[1]
        low_data = data[data.NSamples <= samples_low]  
        med_data = data[(data.NSamples > samples_low) & (data.NSamples <= samples_med)]
        high_data = data[data.NSamples > samples_med]
        groups = [low_data, med_data, high_data]
        node_groups = [(8, 16), (12, 24), (22, 30)]

    sample_parameters = {'NSamples': sample, 'NMutations': nmutations, '%MutTruncal' :muttruncal, 'Gain_Loss': cna, '%GainTruncal' : sample, '%LossTruncal': sample, 'Purity': purity, 'Coverage': coverage, 'WGD': wgd}
    seeds = np.random.choice(10 * n_sims_per_sample_group * len(groups), n_sims_per_sample_group * len(groups), replace=False).reshape(len(groups), -1)
    results = {}
    for i in range(len(groups)):
        sample_group_df = groups[i]
        use_seeds = seeds[i]
        nodes = np.random.choice(np.arange(node_groups[i][0],node_groups[i][1]), n_sims_per_sample_group)
        mut_loss = np.random.choice(np.arange(0.1, 0.26, 0.01), n_sims_per_sample_group)
        noise = np.random.choice(np.arange(0.05, 0.11, 0.01), n_sims_per_sample_group)
        mutloss_trunc = np.random.choice([True, False], n_sims_per_sample_group)
        for j in range(n_sims_per_sample_group):
            rs = use_seeds[j]
            patient_name = 'LTXSIM{}'.format(str((i*n_sims_per_sample_group)+j+1).zfill(3))
            print(patient_name, end=': ')
            sim_params = {feature: func(sample_group_df, feature) for feature, func in sample_parameters.items()}
            output = tx_simulation(patient_name, 
                        n_nodes = nodes[j],      
                        n_positions = sim_params['NMutations'],
                        n_samples = sim_params['NSamples'],
                        min_clones_per_sample = 3,
                        max_clones_per_sample = 8,
                        max_positions=20,
                        prop_truncal_snvs=sim_params['%MutTruncal'],
                        prop_pos_gained=sim_params['Gain_Loss'][0],
                        prop_pos_lost=sim_params['Gain_Loss'][1],    
                        prop_gained_truncal=sim_params['%GainTruncal'],
                        prop_lost_truncal=sim_params['%LossTruncal'],
                        clonal_wgd=sim_params['WGD']['clonal'],
                        subclonal_wgd=sim_params['WGD']['subclonal'],
                        n_subclonal_wgd=sim_params['WGD']['subclonal']*1,
                        minimum_branch_length=1,
                        constant_multiplicity=True, 
                        prop_mutloss_positions=mut_loss[j],                  
                        gain_events=None,
                        loss_events=None,
                        mutloss_truncal=mutloss_trunc[j],
                        all_leaf_observed=True,
                        coverage=sim_params['Coverage'],
                        min_tumour_purity=sim_params['Purity'][0],
                        max_tumour_purity=sim_params['Purity'][1],
                        fract_noise=noise[j],                               
                        n_noise_clones=1,                              
                        noise_type='random',
                        rs=rs,
                        display_plots=display_output,
                        output_directory=output_dir,
                        ccf_vaf_plots=False,
                        save_output=save_output,
                        display_output=True)
            if display_output == True:
                df_names = ['patient_info_df', 'variant_df', 'clone_df', 'edges_df', 'tx_input_df', 'edge_position_events_df']
                results[patient_name] = dict(zip(df_names, output))
    if display_output == True:
        return results
    else:
        return None


purity = (lambda data, feature: sorted(np.random.choice(list(data[data[feature] > 0.2][['Patient', 'Region', feature]].drop_duplicates()[feature]), 2)))
coverage = (lambda data, feature: np.random.choice(list(data[['Patient', 'Region', feature]].drop_duplicates()[feature])))
muttruncal = (lambda data, feature: np.random.choice(list(data[data[feature] < 0.90][['Patient', feature]].drop_duplicates()[feature])))
wgd_return = (lambda choice: {'clonal': True, 'subclonal': False} if choice == 'clonal' else {'clonal': False, 'subclonal': True} if choice == 'subclonal' else {'clonal': False, 'subclonal': False})
wgd = (lambda data, feature: wgd_return(np.random.choice(list(data[['Patient', feature]].drop_duplicates()[feature]))))
gainloss = (lambda patient_data: [patient_data['%PosGained'].values[0], patient_data['%PosLost'].values[0]])
cna = (lambda data, feature: gainloss(data[data.Patient == np.random.choice(list(data['Patient']))][['%PosGained', '%PosLost']].drop_duplicates()))
nmutations = (lambda data, feature: np.random.choice(list(data[data[feature] > 150][['Patient', feature]].drop_duplicates()[feature])))
sample = (lambda data, feature: np.random.choice(list(data[['Patient', feature]].drop_duplicates()[feature]))) 