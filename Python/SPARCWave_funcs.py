# -*- coding: utf-8 -*-

import numpy as np
import cvxpy as cvx
from sklearn.cluster import KMeans    
from sklearn.metrics.pairwise import euclidean_distances

    
def sparse_kmeans(AllDataMatrix,nclust,s,niter,group=False,tree=False,multi=False,groups_vector = None):    
    w = [1/np.sqrt(AllDataMatrix.shape[1])]*AllDataMatrix.shape[1]
    wx = np.zeros((len(AllDataMatrix),AllDataMatrix.shape[1]))
    for j in range(AllDataMatrix.shape[1]):
      wx[:,j] = AllDataMatrix[:,j]*(np.sqrt(w)[j])
    alpha_group = s
    s_orig = s
    print nclust

    kmt = KMeans(n_clusters = nclust, init='k-means++', max_iter=200, n_init=100)
    kmt.fit(wx)
    kmlabels = np.array(kmt.labels_)
    for i in range(niter):
       
        aj_list = []
        for j in range(AllDataMatrix.shape[1]):
            dat_j = AllDataMatrix[:,j].reshape((len(AllDataMatrix),1))
            djall = euclidean_distances(dat_j, dat_j)
            sumd_all = np.sum(djall**2)/len(AllDataMatrix)
            nk_list = [];sumd_k_list = []
        
            for k in range(nclust):
                dat_j = AllDataMatrix[kmlabels==k,j]
                dat_j = dat_j.reshape((len(dat_j),1))
                if(len(dat_j)<1):
                    d = 0
                else:    
                    d = euclidean_distances(dat_j, dat_j)
                nk = len(dat_j)
                sumd_k = np.sum(d**2)
                nk_list.append(nk)
                sumd_k_list.append(sumd_k)
            
            nk_list = np.array(nk_list)
            sumd_k_list = np.array(sumd_k_list)
            #compute within-sum of squares over feature j
            nk_list[nk_list==0] = -1
            one_nk_list = 1./nk_list
            one_nk_list[np.sign(one_nk_list)== -1 ] = 0
            withinssj = np.sum(one_nk_list*sumd_k_list)
            
            aj = sumd_all - withinssj
            aj_list.append(aj)
        #2. get w
        a = np.array(aj_list)
        wvar = cvx.Variable(len(a))
        if tree:
            print "tree structure not supported"
            return -1
        else:    
           if group:
                obj = cvx.Minimize(sum(-1*a*wvar))
                if groups_vector is None:
                        lenseq = np.hstack((np.power(2, np.arange(np.log2(AllDataMatrix.shape[1])))[::-1],[1]))
                else:
                    lenseq = np.array(groups_vector)
                
                lenseq = lenseq.astype(int)
                    
                nlevels = len(lenseq)
                sqrtlenseq = np.sqrt(lenseq)
                indseq = np.cumsum(lenseq)
                t = cvx.Variable(nlevels)
                group0 = [sqrtlenseq[0]*cvx.norm(wvar[0:(indseq[0])],2)<=t[0]]
                group_constraints = group0
                for level in range(1,nlevels):
                    #print level
                    group_const = [sqrtlenseq[level]*cvx.norm(wvar[indseq[(level-1)]:(indseq[level])],2)<=t[level]]
                    group_constraints = group_constraints + group_const
 
                    constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0] +  group_constraints
                    constr = constr + [sum(t)<=alpha_group]
                if multi:
                    T = AllDataMatrix.shape[1]/3

                    t = cvx.Variable(T-1)
                    constr_list = []
                    for coeff in range((T-1)):
                        penalty = cvx.norm(wvar[coeff:(T*3):T],2)<=t[coeff]
                        constr_list.append(penalty)
                   
                    constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0] +constr_list
                    constr = constr + [sum(t)<=alpha_group]
           else:
            ####ORIGINAL SPARSE KMEANS PROBLEM
                print "ORIGINAL"
                obj = cvx.Minimize(sum(-1*a*wvar)) 
                constr = [cvx.square(cvx.norm2(wvar))<=1,cvx.norm1(wvar)<=s_orig, wvar>=0]  

        prob = cvx.Problem(obj, constr)
        try: 
            prob.solve()

        except:
            prob.solve(solver = cvx.SCS,verbose=False)#use_indirect=True
            
        w = wvar.value
        
        #3. update kmeans 
        wx = np.zeros((len(AllDataMatrix),AllDataMatrix.shape[1]))
        for j in range(AllDataMatrix.shape[1]):
            wj = np.sqrt(w[j][0,0])
            if np.isnan(wj):
               
                wj = 10**-30
            wx[:,j] = AllDataMatrix[:,j]*wj

        kmt = KMeans(n_clusters = nclust, init='k-means++', max_iter=200, n_init=100)
        kmt.fit(wx)
        kmlabels = np.array(kmt.labels_)
    
    return prob.value,kmlabels,w

def get_gap_one_s_group(s,AllDataMatrix,permuted_list):
    real_O_s,ll = sparse_kmeans(AllDataMatrix = AllDataMatrix,nclust = s[1], s = s[0],niter=5,group=True,multi=s[2]) 
    O_list_perm_s = list()
    i=0
    for AllDataMatrix_i in permuted_list:
        i+=1
        #print i    
        O,l,w = sparse_kmeans(AllDataMatrix = AllDataMatrix_i, nclust = s[1], s = s[0],niter = 5,group=True,multi=s[2])
        O_list_perm_s.append(-O)
    g = np.log(-real_O_s) - np.mean(np.log(O_list_perm_s))
    return g
    
def get_gap_one_s(s,AllDataMatrix,permuted_list):

    real_O_s,ll = sparse_kmeans(AllDataMatrix = AllDataMatrix,nclust = s[1], s = s[0],niter=5,group=False,multi=s[2]) 
    O_list_perm_s = list()
    i=0
    for AllDataMatrix_i in permuted_list:
        i+=1
        O,l,w = sparse_kmeans(AllDataMatrix = AllDataMatrix_i, nclust = s[1], s = s[0],niter=5,group=False,multi=s[2])
        O_list_perm_s.append(-O)
    g = np.log(-real_O_s) - np.mean(np.log(O_list_perm_s))
    return g

