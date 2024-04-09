import numpy as np
import pandas as pd

def convert_tree_building_format(variant_df, patient_name):
    df = variant_df.copy()
    nucleotides = ['A', 'C', 'T', 'G']
    n_pos = len(pd.unique(df['Variant_Position']))
    n_sample = len(pd.unique(df['Sample']))
    df['PATIENT'] = patient_name
    df['SAMPLE'] = df['Sample'].apply(lambda x: 'U_{}_SU_T1.R{}'.format(patient_name, x.split()[1]))
    df = assign_chr_pos_correct(df, n_pos)
    df['REF'], df['ALT'] = np.array([np.tile(np.random.choice(nucleotides, 2, replace=False), n_sample) for i in range(n_pos)]).reshape(-1,2).T
    df['REF_COUNT'] = df['Total_Reads'] - df['Variant_Reads']
    df['VAR_COUNT'] = df['Variant_Reads']
    df['DEPTH'] = df['Total_Reads']
    df['CLUSTER'] = df['Child']
    df['MUT_COPY'] = df['Fract_Variant']
    df['CCF_PHYLO'] = (df['Est_VAF']*df['Total_Fract_CN'])/(df['MUT_COPY']*df['Purity'])
    df['CCF_OBS'] = df['CCF_PHYLO']
    df['COPY_NUMBER_A'] = df['Fract_a']
    df['COPY_NUMBER_B'] = df['Fract_b']
    df['ACF'] = df['Purity']
    df['TOTAL_COPY'] = df['COPY_NUMBER_A'] + df['COPY_NUMBER_B']
    df['PLOIDY'] = df['SAMPLE'].map(df.groupby('SAMPLE')['TOTAL_COPY'].mean())
    df = df.fillna(0)

    df_tx = df[['PATIENT', 'SAMPLE', 'CHR', 'POS', 'REF', 'ALT', 'REF_COUNT', 'VAR_COUNT', 'DEPTH', 'CLUSTER', 'CCF_PHYLO', 'CCF_OBS', 'MUT_COPY', 'COPY_NUMBER_A', 'COPY_NUMBER_B', 'ACF', 'PLOIDY']]
    df_tx = df_tx.round(2)

    df_og = df[['PATIENT', 'SAMPLE', 'Variant_Position', 'CHR', 'POS', 'REF', 'ALT', 'REF_COUNT', 'VAR_COUNT', 'DEPTH', 'True_VAF', 'Est_VAF', 'CCF', 'DCF', 'CCF_OBS', 'MUT_COPY', 'COPY_NUMBER_A', 'COPY_NUMBER_B', 'TOTAL_COPY', 'CN_Prevalence', 'Purity', 'PLOIDY', 'CLUSTER', 'Noise']].copy()
    df_og.columns = ['PATIENT', 'SAMPLE', 'VAR_POS', 'CHR', 'POS', 'REF', 'ALT', 'REF_COUNT', 'VAR_COUNT', 'DEPTH', 'TRUE_VAF', 'VAF_OBS', 'TRUE_CCF', 'TRUE_DCF', 'CCF_OBS', 'MUT_CN', 'CN_A', 'CN_B', 'TOTAL_CN', 'CN_PREVALENCE', 'PURITY', 'PLOIDY', 'CLUSTER', 'NOISE']
    df_og = df_og.round(2)
    return df_og, df_tx

def assign_chr_pos_correct(df, n_positions):
    chr_lengths = np.array([247199719, 242751149, 199446827, 191263063, 180837866, 170896993, 158821424, 146274826, 140442298, 135374737, 134452384, 132289534, 114127980, 106360585, 100338915, 88822254, 78654742, 76117153, 63806651, 62435965, 46944323, 49528953])
    chr_lengths_dict = dict(zip(np.arange(1,23), chr_lengths))
    chr_prop = chr_lengths / sum(chr_lengths)
    chr = np.random.choice(np.arange(1,23), size=n_positions, p=chr_prop)
    pos = np.array([np.random.randint(1, chr_lengths[i-1]) for i in chr])
    chr_pos_df = pd.DataFrame([chr, pos], index=['CHR', 'POS']).T.reset_index().rename(columns={'index':'Variant_Position'})
    assert len(chr_pos_df.drop_duplicates()) == len(chr_pos_df)
    df = df.merge(chr_pos_df, on='Variant_Position', how='left')
    return df