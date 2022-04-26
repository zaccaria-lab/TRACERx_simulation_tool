# TRACERx evolutionary simulation tool


## Setup:
Creating simulations using the tool requires the following packages:

- Python = 3.9.5
- Numpy
- Pandas
- Jupyter 
- ipykernel
- ipython 
- Scipy 
- Python-graphviz 

These packages can be installed using anaconda into a new environment (simulate) and used to create simulations using the following commands:

```shell 
conda create -n simulate python=3.9.5 jupyter ipykernel pandas numpy ipython scipy python-graphviz
conda activate simulate
```

Then clone the current repo to run the tool using:

```shell
git clone https://github.com/zaccaria-lab/TRACERx_simulation_tool.git
```
## Running the tool:

Simulations can be created using the jupyter notebook 'make_tx_simulations.ipynb' within the scripts directory.
To create simulations, change the input parameter values to desired values and run all cells.


### Inputs:
The tool takes the following inputs -

(The default values are those used to create simulations used to benchmark the TRACERx tree reconstruction tool).

1. n_sims_per_sample_group - The number of simulations you would like to create. If the sample_group_boundaries parameter is not set to None, the number of simulations set will be created for each of the sample groups. (Default: 50). 

2. sample_group_boundaries - Enables the possibility to split the data into 3 groups (low, medium and high) in terms of the number of samples available for each patient. Setting this parameter to a tuple (a,b) will create simulations sampled from data of patients with samples:
- Low: <= a (# nodes sampled uniformly from 8-16)
- Medium: > a, <= b (# nodes sampled uniformly from 12-24)
- High: > b (# nodes sampled uniformly from 22-30)

  Alternatively, setting this value to None will use all patient data to create simulations (# nodes sampled uniformly from 8-30). (Default: (3,7))

3. output_dir - The output directory to save the simulation output. If this value is set to a file path, a separate directory will be created for each simulation with the prefix 'LTXSIM' within the output directory specified. Alternatively, if the output_directory is None, the simulation phylogenetic tree diagram will be printed and simulation output will be saved to a dictionary within the jupyter notebook where each simulation output can be accessed with 'LTXSIM{simulation number}'. (Default: None)


### Outputs:

Each simulations will be named 'LTXSIM{simulation number}'. In the case of separate sample groups, simulations will be numbered as follows: 
 - Low: (1, ...,  n_sims_per_sample_group)
 - Medium: (n_sims_per_sample_group + 1, ...,  2 * n_sims_per_sample_group)
 - High: ((2 * n_sims_per_sample_group) + 1, ..., 3 * n_sims_per_sample_group).


The tool provides the following outputs shown as: output file name if returned (output suffix if saved) -

1. patient_info_df (_info.tsv) - The parameters used to create the simulation.

2. variant_df (.tsv) - Simulated variant read count data for each mutation within the simulated phylogenetic tree for each simulated sample.

3. clone_df (_cloneprops.tsv) - Information regarding the clones present in each simulated sample, along with the clone proportions for each present clone.

4. edges_df (_treeedges.tsv) - Summary information for the number of events that are simulated on each phylogenetic tree branch (identified as Parent --> Child), including: SNVs, CNA gains, CNA losses, WGD and whether the branch is Truncal.

5. tx_input_df (_TX.tsv) - A version of the variant_df that has been converted to the correct format to input to the TRACERx phylogenetic tree reconstruction tool.

6. edge_position_events_df (_events.tsv) - A dataframe where each row is a simulated genomic position and each column is a phylogenetic tree branch, displaying specific events that occur on each branch at each genomic position. 

7. phylogenetic tree image (_tree.png) - An image version of the number of the phylogenetic tree showing a summary of the tree structure and number of events on each branch. This is printed directly rather than returned if output_dir parameter is None.

