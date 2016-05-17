# -*- coding: utf-8 -*-


import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from rpy2.robjects.numpy2ri import numpy2ri
from joblib import Parallel, delayed

import cvxpy as cvx
from sklearn.cluster import KMeans 
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.cluster import adjusted_rand_score
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import pylab as pl


rkmeans  = robjects.r['kmeans']
package_src_path = '' ### ENTER THE PACKAGE FILE PATH HERE
robjects.r('source('+ "'" + package_src_path + 'R_functions.R' + "');")

simulate_curves_pad = robjects.globalenv['simulate.curves.pad']
transform = robjects.globalenv['transform']
simulate_multi = robjects.globalenv['simulate.multidata']
R_sparse_kmeans = robjects.globalenv['R_sparse_kmeans']


from SPACWave_functs import sparse_kmeans,get_gap_one_s_group

if __name__ == '__main__':
    path_to_save_files = ''##
    ### UNIVARIATE/MULTIVARIATE simulations
    SIGNAL_TYPE = "univariate"
    #SIGNAL_TYPE = "multivariate"

    SIMULATION_TYPES = ["group","sparse"]
    SIMULATION_TYPES = ["sparse"]


    machine_cores_to_use = 7 
    nsim = 50; perm_num = 5; niter = 5
    s_list =  [3,4,5,5.5,6.5,7,8.5,10,12,15]

    if SIGNAL_TYPE == "univariate":
        nclust = 6
        N = 256
        pad = 128
        n_list = [2,3,4,5,8,10,12,15,20,25]
        sigma = 2.75
        multi=False
    if SIGNAL_TYPE == "multivariate":
        nclust = 5
        N = 128
        pad = 64
        n_list = [2,3,4,5,6,7,8,9,10]
        sigma = 1.75
        multi=True

    for n in n_list:
        ###CREATE GROUND-TRUTH LABELS
        c = range(nclust)
        labels = []
        for clst in c:
            labels = labels +[clst]*n

        sparse_group_res = np.zeros([nsim,1])
        sparse_res = np.zeros([nsim,1])
        
        for sim in range(nsim):
            print "****sim =" + str(sim) + " n=" + str(n)
            if SIGNAL_TYPE == "univariate":
            
                AllDataMatrix = simulate_curves_pad(N=N,sigma = sigma,n=n,pad=pad)
                AllDataMatrix = np.array(AllDataMatrix)
                AllDataMatrix = transform(numpy2ri(AllDataMatrix),"s8")
                AllDataMatrix = np.array(AllDataMatrix[0])
            if SIGNAL_TYPE == "multivariate":
            
                AllDataMatrix = simulate_multi(N=128,ncurves=n,sigma=sigma,pad=64)
                AllDataMatrix = np.array(AllDataMatrix[2])
            
            print AllDataMatrix.shape
            AllDataMatrix_temp = np.copy(AllDataMatrix)
            permuted_list = list()
    
            for perm in range(perm_num):
                for i in range(AllDataMatrix_temp.shape[1]):
                    np.random.shuffle(AllDataMatrix_temp[:,i])
                permuted_list.append(AllDataMatrix_temp)
            
            
            try:
                if "group" in SIMULATION_TYPES: 
                    results_group = Parallel(n_jobs = machine_cores_to_use)(delayed(get_gap_one_s_group)(i,AllDataMatrix = AllDataMatrix,permuted_list = permuted_list) for i in zip(s_list,[nclust]*len(s_list),[multi]*len(s_list)))  
     
                    best_sparse_kmeans_group,lgroup,wgroup  = sparse_kmeans(AllDataMatrix = AllDataMatrix, nclust = nclust, s = s_list[np.argmax(results_group)],
                                                                                                                   niter=niter,group = True, multi = multi) 
                    print(s_list[np.argmax(results_group)])
                    print(adjusted_rand_score(labels,lgroup))
    
                    sparse_group_res[sim,:] = adjusted_rand_score(labels,lgroup)
                    
                    path = path_to_save_files +"GROUP_n=" +str(n)+ "sigma=" +str(sigma) + "signal=" + SIGNAL_TYPE +".txt"
                    np.savetxt(path, sparse_group_res)
                if "sparse" in SIMULATION_TYPES: 

                    results = R_sparse_kmeans(data = numpy2ri(AllDataMatrix),nclust = nclust,nperms  = perm_num, s = -1) 
                    l =  np.array(results[0])
                    print adjusted_rand_score(labels,l)
                    sparse_res[sim,:] = adjusted_rand_score(labels,l)
                    
                    path = path_to_save_files +"sparse_n=" +str(n)+ "sigma=" +str(sigma) + "signal=" + SIGNAL_TYPE +".txt"
    
                    np.savetxt(path, sparse_res)                
            except:
                 continue

                