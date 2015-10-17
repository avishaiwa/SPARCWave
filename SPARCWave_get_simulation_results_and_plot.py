# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

n_list = [2,3,4,5,8,10,12,15,20,25]
n_list  = [27,30,33,36];sigma = 2.75

#### PLOT - YES/NO
MAKE_PLOT = False
#### GET UNIVARIATE/MULTIVARIATE RESULTS
SIGNAL_TYPE = "univariate"

#### ENTER PATH FOR RESULTS FILE CREATED BY SPARCWaveMAIN.py
RESULTS_PATH  = ""

sparse_avg_list =[]; group_avg_list = [] 
for n in n_list:
    print n
    f = RESULTS_PATH + 'sparse_n=' +str(n) + 'sigma=' +str(sigma) + 'signal=' + SIGNAL_TYPE + '.txt'
    res = np.loadtxt(open(f,"rb"),delimiter=",")
    avg = np.mean(res)
    sparse_avg_list.append(avg)
    
#    f = RESULTS_PATH + 'GROUP_n=' +str(n) + 'sigma=' +str(sigma) + 'signal=' + SIGNAL_TYPE + '.txt'
#    res2 = np.loadtxt(open(f,"rb"),delimiter=",")
#    avg2 = np.mean(res2)
#    group_avg_list.append(avg2)
    

print sparse_avg_list

if MAKE_PLOT:
    n = [2*6,3*6,4*6,5*6,8*6,10*6,12*6,15*6,20*6,25*6]
    n = [27*6,30*6,33*6,36*6]
    ax = plt.subplot(111)
    plt.plot(n,sparse_avg_list, marker ='o', ms=5, label="SPARCwave",color="green")
    #plt.plot(n,group_avg_list, marker ='o', ms=5, label="SPARCwave_GROUP",color="yellow")
    
    
    plt.xlabel('Number of signals', fontsize=12)
    plt.ylabel('Adjusted Rand Index', fontsize=12)
    plt.xlim(xmin=9,xmax=155)
    plt.ylim(ymin=0.2)
    
    plt.legend(loc="best",fancybox=True, framealpha=0.5,prop={'size':9.4})
    #plt.savefig("sim_plt.png", dpi = 500)
    #plt.savefig("sim_plt.eps", dpi = 500)
    plt.show()