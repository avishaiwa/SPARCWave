# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:22:54 2015

@author: tomhope
"""
import cvxpy as cvx

from sklearn.cluster import KMeans    
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.cluster import adjusted_rand_score
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as robjects
import pylab as pl
rkmeans  = robjects.r['kmeans']

#N=100
#d1 = np.random.multivariate_normal([0,0,0],np.eye(3),N)
#d2 = np.random.multivariate_normal([2,2,0],np.eye(3),N)
#AllDataMatrix = np.vstack((d1,d2))
AllDataMatrix = np.loadtxt(open("C://Users//tomhope//Dropbox//wavelet//trans_data_for_py.csv","rb"),delimiter=",")
labels = np.loadtxt(open("C://Users//tomhope//Dropbox//wavelet//labels_for_py.csv","rb"),delimiter=",")
w = [1/np.sqrt(AllDataMatrix.shape[1])]*AllDataMatrix.shape[1]

wx = np.zeros((len(AllDataMatrix),AllDataMatrix.shape[1]))
for j in range(AllDataMatrix.shape[1]):
  wx[:,j] = AllDataMatrix[:,j]*(np.sqrt(w)[j])
group = True
alpha_group = 25
s_orig = 25
nclust = 6
niter = 10
#kmt = KMeans(n_clusters=nclust, init='random', n_init=100,verbose=False, tol=0.0000000001)                    
#kmt.fit(wx)
#print kmt.labels_
#print adjusted_rand_score(labels, kmt.labels_)


kmt = rkmeans(x=numpy2ri(wx),centers=nclust,nstart=500)
kmlabels = np.array(kmt[0])
print adjusted_rand_score(labels, np.array(kmt[0]))
#overall iterations
for i in range(niter):
    print i
    #1.get bcssj

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
        #aj = totalss/n - wss/nk
        aj = sumd_all - withinssj
        aj_list.append(aj)
    #2. get w
    a = np.array(aj_list)
    lenseq = np.array([256,128,64,32,16,8,4,2,1,1])
    lenseq = np.array([256,128,64,32,16,8,8])
    nlevels = len(lenseq)

    sqrtlenseq = np.sqrt(lenseq)
    indseq = np.cumsum(lenseq)
    wvar = cvx.Variable(len(a))
    
    t = cvx.Variable(nlevels)
    ## Form objective.

    if group:
    ####GROUP SPARSE
        #obj = cvx.Minimize(sum(-1*a*wvar) + alpha_group*sum(t))
        obj = cvx.Minimize(sum(-1*a*wvar))

       
        group0 = [sqrtlenseq[0]*cvx.norm(wvar[0:(indseq[0]-1)],2)<=t[0]]
        group1 = [sqrtlenseq[1]*cvx.norm(wvar[indseq[0]:(indseq[1])],2)<=t[1]]
        group2 = [sqrtlenseq[2]*cvx.norm(wvar[indseq[1]:(indseq[2])],2)<=t[2]]
        group3 = [cvx.norm(wvar[indseq[2]:(indseq[3])],2)<=t[3]]
        group4 = [cvx.norm(wvar[indseq[3]:(indseq[4])],2)<=t[4]]
        group5 = [cvx.norm(wvar[indseq[4]:(indseq[5])],2)<=t[5]]
        group6 = [cvx.norm(wvar[indseq[5]:(indseq[6])],2)<=t[6]]
        

#        group0 = [cvx.norm(wvar[0:indseq[0]],2)<=t[0]]
#        group1 = [cvx.norm(wvar[indseq[0]:indseq[1]],2)<=t[1]]
        group0 = [sqrtlenseq[0]*cvx.norm(wvar[0:(indseq[0]-1)],2)<=t[0]]
        group1 = [sqrtlenseq[1]*cvx.norm(wvar[indseq[0]:(indseq[1])],2)<=t[1]]
        group2 = [sqrtlenseq[2]*cvx.norm(wvar[indseq[1]:(indseq[2])],2)<=t[2]]
        group3 = [sqrtlenseq[3]*cvx.norm(wvar[indseq[2]:(indseq[3])],2)<=t[3]]
        group4 = [sqrtlenseq[4]*cvx.norm(wvar[indseq[3]:(indseq[4])],2)<=t[4]]
        group5 = [sqrtlenseq[5]*cvx.norm(wvar[indseq[4]:(indseq[5])],2)<=t[5]]
        group6 = [sqrtlenseq[6]*cvx.norm(wvar[indseq[5]:(indseq[6])],2)<=t[6]]
#        
        #group7 = [cvx.norm(wvar[indseq[6]:(indseq[7])],2)<=t[7]]
#        group8 = [cvx.norm(wvar[indseq[7]:(indseq[8])],2)<=t[8]]
#        group9 = [cvx.norm(wvar[indseq[8]:(indseq[9])],2)<=t[9]]
        
    ###"correct" constraints
        #constr = [wvar>=0,sum(wvar)==1] + group0 + group1 + group2 + group3 + group4 + group5 + group6
    ##l2 constraints
        #constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0] + group0 + group1 + group2 + group3 + group4 + group5 + group6 + group7 + group8 + group9
        #constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0] + group0 + group1 + group2 + group3 + group4 + group5 + group6 + group7
        constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0] + group0 + group1 + group2 + group3 + group4 + group5 + group6
        constr = constr + [sum(t)<=alpha_group]#cvx.norm1(wvar)<=s_orig

####GROUP NORM AS IN LASSO
#        groupnormvec = [cvx.norm(wvar[0:(indseq[0]-1)],2),cvx.norm(wvar[indseq[0]:(indseq[1])],2),
#                    cvx.norm(wvar[indseq[1]:(indseq[2])],2),cvx.norm(wvar[indseq[2]:(indseq[3])],2),
#                    cvx.norm(wvar[indseq[3]:(indseq[4])],2),cvx.norm(wvar[indseq[4]:(indseq[5])],2),
#                    cvx.norm(wvar[indseq[5]:(indseq[6])],2)]
#        obj = cvx.Minimize(sum(-1*a*wvar) + alpha_group*sum(groupnormvec))
#        constr = [cvx.square(cvx.norm2(wvar))<=1,wvar>=0]  
    else:
    ####ORIGINAL SPARSE KMEANS PROBLEM
        #obj = cvx.Minimize(cvx.sum(cvx.mul_elemwise(-1*a,wvar))) 
        obj = cvx.Minimize(sum(-1*a*wvar)) 
        
        constr = [cvx.square(cvx.norm2(wvar))<=1,cvx.norm1(wvar)<=s_orig, wvar>=0]  
        #constr = [cvx.square(cvx.norm2(wvar))<=1, wvar>=0]  

    prob = cvx.Problem(obj, constr)
    #prob.solve()

    try: 
        prob.solve()
        print "default solver"
    except:
        
        print "SCS SOLVER"
        #prob.solve(solver =cvx.CVXOPT)
        prob.solve(solver = cvx.SCS,verbose=False)#use_indirect=True
        print prob.value
    w = wvar.value

    #3. update kmeans 
    wx = np.zeros((len(AllDataMatrix),AllDataMatrix.shape[1]))
    for j in range(AllDataMatrix.shape[1]):
        wj = np.sqrt(w[j][0,0])
        #wj = w[j][0,0]
        if np.isnan(wj):
            #print "bad"
            wj = 10**-20
#        else:
#            #print "yes"
#            #print wj
        wx[:,j] = AllDataMatrix[:,j]*wj

#    kmt = KMeans(n_clusters=nclust, init='random', n_init=100,verbose=False,tol=0.0000000001)                    
#    kmt.fit(wx)
    kmt = rkmeans(x=numpy2ri(wx),centers=nclust,nstart=100)
    kmlabels =  np.array(kmt[0])

print kmlabels
pl.plot(w)
from sklearn.metrics.cluster import adjusted_rand_score
print adjusted_rand_score(labels,kmlabels)