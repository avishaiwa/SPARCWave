# -*- coding: utf-8 -*-


import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from rpy2.robjects.numpy2ri import numpy2ri
from joblib import Parallel, delayed  
import sys
import rpy2.robjects as robjects
from SPACWave_functs import sparse_kmeans,get_gap_one_s_group

rkmeans  = robjects.r['kmeans']
package_src_path = '' ### ENTER THE PACKAGE FILE PATH HERE

robjects.r('source('+ "'" + package_src_path + 'R_functions.R' + "');")
R_sparse_kmeans = robjects.globalenv['R_sparse_kmeans']


if __name__ == '__main__':

    data_path = "" ### ENTERE THE DATA PATH HERE, E.G. "data_growth_T32_Q2_new_ord.csv"
    label_path = "" ### ENTERE THE LABEL PATH HERE, E.G. "labels_growth.csv"
    groups_vector = [2]*21  ### VECTOR DEFINING INDICES OF GROUPS, E.G. FOR SCATTERING
    nclust = 2 ###Number of clusters

    
    AllDataMatrix = np.loadtxt(open(data_path,"rb"),delimiter=",")
    labels = np.loadtxt(open(label_path,"rb"),delimiter=",")
    labels = labels.astype(int)

    machine_cores_to_use = 7 
    perm_num = 5; niter=5
    use_gap_statistic = False; s_list_group =  [1.25,1.5,2,3,4,5,5.5,6.5,7,8.5,10]
    s = 1.25
    print AllDataMatrix.shape
    print nclust
    use_group = True

    if use_group:
        if use_gap_statistic:
            "GROUP - GAP STATISTIC"
            sys.stdout.flush()

            print AllDataMatrix.shape
            AllDataMatrix_temp = np.copy(AllDataMatrix)
            permuted_list = list()
    
            for perm in range(perm_num):
                for i in range(AllDataMatrix_temp.shape[1]):
                    np.random.shuffle(AllDataMatrix_temp[:,i])
                permuted_list.append(AllDataMatrix_temp)
                
            results_group = Parallel(n_jobs = machine_cores_to_use)(delayed(get_gap_one_s_group)(i,AllDataMatrix = AllDataMatrix,permuted_list = permuted_list) for i in zip(s_list_group,[nclust]*len(s_list_group),[False]*len(s_list_group)))  
            s = s_list_group[np.argmax(results_group)]
            best_sparse_kmeans_group,lgroup,w = sparse_kmeans(AllDataMatrix = AllDataMatrix, nclust = nclust, s = s,niter = niter,
                                                        group = True,groups_vector = groups_vector) 
        else:
            print "GROUP - GIVEN S"
            sys.stdout.flush()

            best_sparse_kmeans_group,lgroup,w  = sparse_kmeans(AllDataMatrix = AllDataMatrix, nclust = nclust, s = s,niter = niter,
                                                        group = True,groups_vector =  groups_vector) 
        print adjusted_rand_score(labels,lgroup)
        
    else:
        if use_gap_statistic:
            print "GAP STATISTIC METHOD"
            sys.stdout.flush()

            results = R_sparse_kmeans(data = numpy2ri(AllDataMatrix),nclust = nclust,nperms  = perm_num, s = -1) 
            l =  np.array(results[0])
        
            print adjusted_rand_score(labels,l)                 
        else:
            print "S GIVEN AS INPUT"
            sys.stdout.flush()

            results = R_sparse_kmeans(data = numpy2ri(AllDataMatrix),nclust = nclust,nperms  = perm_num, s = s) 
            l =  np.array(results[0])
            w =  np.array(results[1])
            print adjusted_rand_score(labels,l) 
