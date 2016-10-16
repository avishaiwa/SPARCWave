####################################
Overview: 
=================
The SPARCWave package contains utilities for sparse and structured-sparse clustering of time series using  
time-freqeuncy signal representations obtained from the wavelets and scattering transforms, implementing the methods and experiments in [1]. Datasets for the experiments are also provided.
For advanced instructions, please refer to the README_ADVANCED file.
###################################


### Basic usage example:

Code for running the SPARCWave and SPARCWave Group-sparse clustering methods is contained in SPARCWave_run_real_data.py.

Here are two lines that run both methods (depending on the `group` argument), and evaluate the clustering result with
the Adjusted Rand Score if ground-truth labels are given:
```python
"""
AllDataMatrix: Dataset, numpy matrix
nclust: Number of clusters
s: Sparsity tuning parameter
niter: Number of iterations
group: Flase for SPARCWave, True for SPARCWave Group
groups_vector: A partitioning of time-frequency indices into groups
"""
best_sparse_kmeans_group,lgroup,w  = sparse_kmeans(AllDataMatrix = AllDataMatrix, nclust = nclust, s = s,niter = niter,
                                                        group = True,groups_vector =  groups_vector)
print adjusted_rand_score(labels,lgroup)                                                        
```

In the above, we assume a pre-defined tuning parameter `s`. To select it automatically, use a permutation-based method which can be parallelized with the `joblib` library:
```python
"""
machine_cores_to_use: How many cores to use in parallel
perm_num: How many permutations
s_list_group: A grid of s values
"""
machine_cores_to_use = 7 
perm_num = 5
s_list_group =  [1.25,1.5,2,3,4,5,5.5,6.5,7,8.5,10]

permuted_list = list()
    
for perm in range(perm_num):
    for i in range(AllDataMatrix_temp.shape[1]):
        np.random.shuffle(AllDataMatrix_temp[:,i])
    permuted_list.append(AllDataMatrix_temp)

results_group = Parallel(n_jobs = machine_cores_to_use)(delayed(get_gap_one_s_group)(i,AllDataMatrix = AllDataMatrix,permuted_list = permuted_list) for i in zip(s_list_group,[nclust]*len(s_list_group),[False]*len(s_list_group)))  
s = s_list_group[np.argmax(results_group)]
best_sparse_kmeans_group,lgroup,w = sparse_kmeans(AllDataMatrix = AllDataMatrix, nclust = nclust, s = s,niter = niter,
                                            group = True,groups_vector = groups_vector) print adjusted_rand_score(labels,lgroup)                                                        
print adjusted_rand_score(labels,lgroup)
```

### Data files:

**CSV files:**
data_growth_T32_Q2_new_ord.csv - Berkeley Growth dataset after scattering transform

labels_growth.csv - labels for Growth dataset

phoneme_scatter_T32_Q8_New_ord.csv - Phoneme dataset after scattering transform

labels_phoneme.csv - labels for Phoneme dataset

wheat_T32_Q8_new_ord.csv - Wheat dataset after scattering transform

labels_wheat.csv - labels for Wheat dataset

### Installation:

**Required Python libraries:** CVXPY, sklearn, joblib, numpy, pylab (and dependencies). 

### Acknowledgments:

SPARCWave was developed by Tom Hope, Avishai Wagner and Or Zuk, as part of work on the paper:

[1]  "Clustering Noisy Signals with Structured Sparsity Using Time-Frequency Representation", T. Hope, A. Wagner and O. Zuk, http://arxiv.org/abs/1510.05214 , 2015

Please cite the above paper if using the package.

For support, any questions or comments, please contact:

Tom Hope: tom.hope@huji.ac.il

Avishai Wagner: avishaiwa@gmail.com

Or Zuk: or.zuk@mail.huji.ac.il
