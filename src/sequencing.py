import pandas as pd
import numpy as np
from scipy.stats import dirichlet, poisson, binom



def make_sequencing_data(clone_profiles, dcf_profiles, n_samples, min_clones_per_sample, max_clones_per_sample, diploid, all_leaf_observed=True, leaves=None, coverage=400, min_tumour_purity=0.1, max_tumour_purity=0.9):
    clones = list(clone_profiles.keys())
    sample_clone_proportions, sample_tumour_purities = all_samples(clones, leaves, n_samples, min_clones_per_sample, max_clones_per_sample, all_leaf_observed, diploid, min_tumour_purity, max_tumour_purity)
    fract_copy_num, vaf, fract_variant, fract_a, fract_b, cn_prevalence = calculate_vaf(clone_profiles, sample_clone_proportions, diploid, sample_tumour_purities)
    total_reads = get_total_reads(fract_copy_num, coverage)
    variant_reads = get_variant_reads(total_reads, vaf)
    ccf = calculate_ccf(sample_clone_proportions, sample_tumour_purities, clone_profiles)
    dcf = calculate_dcf(sample_clone_proportions, sample_tumour_purities, clone_profiles, dcf_profiles)
    variant_output_df = create_output_df(variant_reads, total_reads, vaf, ccf, dcf, fract_copy_num, fract_variant, fract_a, fract_b, sample_tumour_purities, cn_prevalence)
    clone_prop_df = create_clone_prop_df(sample_clone_proportions, diploid)
    return variant_output_df, clone_prop_df


def ensure_leaves_observed(n_samples, leaves):
    samples = np.arange(1,n_samples+1) 
    leaf_assignment = np.random.choice(samples, len(leaves))

    sample_leaf_assignment = {i : np.array([leaves[k] for k in range(len(leaves)) if leaf_assignment[k] == i]) for i in samples}
    return sample_leaf_assignment   


def simulate_sample_leaves(clones, assigned_clones, n_clones_per_sample, root, tumour_purity):
    n_clones_per_sample = np.maximum(n_clones_per_sample - len(assigned_clones), 0)
    if n_clones_per_sample > 0:
        sample_clones = np.random.choice(clones[~np.isin(clones,assigned_clones)], n_clones_per_sample, replace=False)
        sample_clones = np.concatenate([assigned_clones, sample_clones]).astype(int)
    else: 
        sample_clones = assigned_clones.astype(int)
    clone_proportions = dirichlet.rvs(np.ones(len(sample_clones)))[0]* tumour_purity
    sample_clone_proportions = dict(zip(sample_clones, clone_proportions))
    sample_clone_proportions[root] = (1 - tumour_purity)
    np.testing.assert_almost_equal(sum(sample_clone_proportions.values()),1), 'Sum of clone proportions is not 1'
    return sample_clone_proportions


def simulate_sample(clones, n_clones_per_sample, root, tumour_purity):
    sample_clones = np.random.choice(clones, n_clones_per_sample, replace=False)
    clone_proportions = dirichlet.rvs(np.ones(len(sample_clones)))[0] * tumour_purity
    sample_clone_proportions = dict(zip(sample_clones, clone_proportions))
    sample_clone_proportions[root] = (1 - tumour_purity)
    np.testing.assert_almost_equal(sum(sample_clone_proportions.values()),1), 'Sum of clone proportions is not 1'
    return sample_clone_proportions


def all_samples(clones, leaves, n_samples, min_clones_per_sample, max_clones_per_sample, all_leaf_observed, root, min_tumour_purity, max_tumour_purity):
    clones_use = clones.copy()
    clones_use.remove(root) 
    sample_tumour_purities = dict(zip(np.arange(1, n_samples+1), np.random.uniform(min_tumour_purity, max_tumour_purity, n_samples)))
    get_n_clones = np.random.choice(np.arange(min_clones_per_sample, max_clones_per_sample+1), n_samples)
    if all_leaf_observed == False:
        samples = {'Sample {}'.format(i): simulate_sample(clones_use, get_n_clones[i-1], root, sample_tumour_purities[i]) for i in range(1, n_samples+1)}
    else:
        assert leaves is not None, 'Please specify the leaf nodes'
        clones_use = np.array(clones_use)
        leaf_assignment = ensure_leaves_observed(n_samples, leaves)
        samples = {'Sample {}'.format(i): simulate_sample_leaves(clones_use, leaf_assignment[i], get_n_clones[i-1], root, sample_tumour_purities[i]) for i in range(1, n_samples+1)}
    return samples, sample_tumour_purities


def calculate_vaf(clone_profiles, sample_clone_prop, diploid, sample_tumour_purities):
    get_frac_copies = (lambda c, p, idx1, idx2 : np.array([sum(genotype[idx1: idx2])*p for pos, genotype in clone_profiles[c].items()]))
    sample_fract_copies = (lambda idx1, idx2 : {sample: sum(get_frac_copies(clone, prop, idx1, idx2) for clone, prop in proportions.items()) for sample, proportions in sample_clone_prop.items()})
    tumour_spec_fract_cn = (lambda idx1, idx2 : {sample: sum(get_frac_copies(clone, prop, idx1, idx2) for clone, prop in proportions.items() if clone != diploid) for sample, proportions in sample_clone_prop.items()})
    get_sample_purity_name = (lambda sample : int(sample.split(' ')[-1])) # ABI: This seems a bit hacky; can we save this as a dictionary instead?
    divide_by_purity = (lambda x: {sample : fract / sample_tumour_purities[get_sample_purity_name(sample)] for sample, fract in x.items()})
    get_normalisation = (lambda sample : sum(np.array([sample_clone_prop[sample][c] if genotype[-1] > 0 else 0 for pos, genotype in clone_profiles[c].items()]) for c in sample_clone_prop[sample]))
    normalised_division = (lambda fract, norm: np.divide(fract, norm, out=np.zeros_like(fract), where=(norm!=0)))
    normalise_mutations = (lambda x: {sample : normalised_division(fract, get_normalisation(sample)) for sample, fract in x.items()})
    
    get_cell_prevalence = (lambda c, p, idx1, idx2: np.array([p if (genotype[idx1] !=1) | (genotype[idx2] !=1) else 0 for pos, genotype in clone_profiles[c].items()]))
    sample_cell_prevalence = (lambda idx1, idx2 : {sample: sum(get_cell_prevalence(clone, prop, idx1, idx2) for clone, prop in proportions.items()) for sample, proportions in sample_clone_prop.items()})
    cn_prevalence = sample_cell_prevalence(0, 1)

    fract_copy_num = sample_fract_copies(0, 2)
    fract_variant_cn = tumour_spec_fract_cn(2, 3)
    assert_test = sample_fract_copies(2, 3)
    assert all(np.array_equal(fract_variant_cn[c],assert_test[c]) for c in fract_variant_cn)
    mutant_copy_number = normalise_mutations(fract_variant_cn)
    fract_a = tumour_spec_fract_cn(0, 1)
    fract_a = divide_by_purity(fract_a)
    fract_b = tumour_spec_fract_cn(1, 2)
    fract_b = divide_by_purity(fract_b)
    
    vaf = {sample: fract_variant_cn[sample]/fp for sample, fp in fract_copy_num.items()}

    return fract_copy_num, vaf, mutant_copy_number, fract_a, fract_b, cn_prevalence


def get_total_reads(fract_cn, coverage):
    fbar = (lambda sample_fp: sum(sample_fp)/len(sample_fp))
    lambdap = (lambda sample_fp, coverage: sample_fp * (coverage/fbar(sample_fp)))
    total_reads = {sample : poisson.rvs(lambdap(cn, coverage)) for sample, cn in fract_cn.items()}
    return total_reads


def get_variant_reads(total_reads, vaf):
    variant_reads = {sample_t : binom.rvs(n=sample_total_reads, p=sample_vaf) for sample_t, sample_total_reads in total_reads.items() for sample_v, sample_vaf in vaf.items() if sample_t == sample_v}
    return variant_reads


def calculate_ccf(clone_props, purities, clone_profiles):
    snvs = {clone: np.array([(profile[-1] > 0).astype(int) for pos, profile in c_profile.items()]) for clone, c_profile in clone_profiles.items()}

    get_sample_proportions = (lambda clones: np.sum(np.array([snvs[clone]*prop for clone, prop in clones.items()]), axis=0))
    get_sample_purity_name = (lambda sample: int(sample.split(' ')[-1]))
    ccf = {sample: get_sample_proportions(clones)/purities[get_sample_purity_name(sample)] for sample, clones in clone_props.items()}
    return ccf

def calculate_dcf(clone_props, purities, clone_profiles, dcf_profiles):
    ccf_muts = np.array([b[-1] for k, v in clone_profiles.items() for a, b in v.items()])
    dcf_muts = np.array([b[-1] for k, v in dcf_profiles.items() for a, b in v.items()])
    assert sum(ccf_muts > dcf_muts) == 0
    snvs = {clone: np.array([(profile[-1] > 0).astype(int) for pos, profile in c_profile.items()]) for clone, c_profile in dcf_profiles.items()}

    get_sample_proportions = (lambda clones: np.sum(np.array([snvs[clone]*prop for clone, prop in clones.items()]), axis=0))
    get_sample_purity_name = (lambda sample: int(sample.split(' ')[-1]))
    dcf = {sample: get_sample_proportions(clones)/purities[get_sample_purity_name(sample)] for sample, clones in clone_props.items()}
    return dcf

def make_df(data, col_name):
    df = pd.DataFrame(data)
    df = df.reset_index()
    final_df = pd.melt(df, id_vars='index')
    final_df.columns = ['Variant_Position', 'Sample', col_name]
    return final_df


def create_output_df(variant_reads, total_reads, vaf, ccf, dcf, fract_cn, fract_snv, fract_a, fract_b, purity, cn_prevalence):
    col_names = ['Variant_Reads', 'Total_Reads', 'True_VAF', 'CCF', 'DCF', 'Total_Fract_CN', 'Fract_Variant', 'Fract_a', 'Fract_b', 'CN_Prevalence']
    total_data = [variant_reads, total_reads, vaf, ccf, dcf, fract_cn, fract_snv, fract_a, fract_b, cn_prevalence]
    overall_data = [make_df(total_data[i], col_names[i]) for i in range(len(total_data))]
    output = overall_data[0]
    for i in range(1,len(total_data)):
        output = output.merge(overall_data[i], on=['Variant_Position', 'Sample'])
    purity_df = pd.DataFrame([purity]).T.reset_index()
    purity_df.columns = ['Sample', 'Purity']
    purity_df['Sample'] = purity_df['Sample'].apply(lambda x: 'Sample ' + str(x))
    output = output.merge(purity_df, on=['Sample'])
    output['Est_VAF'] = output['Variant_Reads']/ output['Total_Reads']
    return output


def create_clone_prop_df(clone_proportions, diploid):
    df = pd.DataFrame(clone_proportions)
    df = df.reset_index()
    final_df = pd.melt(df, id_vars='index')
    final_df.columns = ['Clone', 'Sample', 'Proportion']
    final_df = final_df[['Sample', 'Clone', 'Proportion']].dropna().sort_values(by=['Sample', 'Clone'])
    final_df['Is_Normal'] = 'No'
    final_df.loc[final_df['Clone'] == diploid, 'Is_Normal'] = 'Yes'
    return final_df
