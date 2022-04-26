# TRACERx evolutionary simulation tool:


## Running the tool:
Creating simulations using the tools requires the following packages:

1. Python >= 3.9
2. Numpy
3. Pandas
4. Jupyter 
5. ipykernel
6. ipython 
7. Scipy 
8. Python-graphviz 

These packages can be installed using anaconda into a new environment and used to create simulations using the following commands:

conda create -n simulate python=3.9 jupyter ipykernel pandas numpy ipython scipy python-graphviz seaborn matplotlib

conda activate simulate


## Inputs:
The tool takes the following inputs -

(The default values are those used to create simulations used to benchmark the TRACERx tree reconstruction tool).

1. n_sims_per_sample_group - The number of simulations you would like to create. If the sample_group_boundaries parameter is not set to None, the number of simulations set will be created for each of the sample groups. (Default: 50)

2. sample_group_boundaries - Enables the possibility to split the data into 3 groups (low, medium and high) in terms of the number of samples available for each patient. Setting this parameter to a tuple (a,b) will create simulations sampled from data of patients with samples:
- Low: <= a (No. nodes sampled uniformly from 8-16)
- Medium: > a, <= b (No. nodes sampled uniformly from 12-24)
- High: > b (No. nodes sampled uniformly from 22-30)
Alternatively, setting this value to None will use all patient data to create simulations (no. nodes sampled uniformly from 8-30). (Default: (3,7))

3. output_dir - The output directory to save the simulation output. If this value is not None, a separate directory will be created for each simulation with the prefix 'LTXSIM' within the output directory specified. Alternatively, if the output_directory is None, the simulation phylogenetic tree diagram will be printed and simulation output will be saved to a dictionary within the jupyter notebook where each simulation output can be accessed with 'LTXSIM<simulation number>'. (Default: None)


## Outputs:
The tool provides the following outputs -

1. patient_info_df - The parameters used to create the simulation.

2. variant_df - Simulated variant read count data for each mutation within the simulated phylogenetic tree for each simulated sample.

3. clone_df - Information regarding the clones present in each simulated sample, along with the clone proportions for each present clone.

4. edges_df - Summary information for the number of events that are simulated on each phylogenetic tree branch (identified as Parent --> Child), including: SNVs, CNA gains, CNA losses, WGD and whether the branch is Truncal.

5. tx_input_df - A version of the variant_df that has been converted to the correct format to input to the TRACERx phylogenetic tree reconstruction tool.

6. edge_position_events_df - A dataframe where each row is a simulated genomic position and each column is a phylogenetic tree branch, displaying specific events that occur on each branch at each genomic position. 

