library(clues);library(wmtsa);library(sparcl);library(flexclust);library(kernlab)
source('C:/Users/tomhope/Dropbox/wavelet/src/noise_NEW.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/transform_data.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/wave_sparse_Kmeans.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/HMMcode-Revised_Shima.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/similar_PCA_alg.R')
nclust=5
#n=15
#sigma=3
wavelet ="s8"
#data = simulate.curves(N=512,snr=0.2,n=n)
#sigma_list = seq(2.25,3.75,0.25)
#sigma_list = 3.75
#sigma_list = c(2,2.5,3)
sigma_list = c(1.25,1.75)
n_list = c(2,3,4,5,6,7,8,9,10)
#n_list = c(15,20)
# 
nsim=150
#sigma_list = 4
#n_list=45
for(sigma in sigma_list){
  
  for(ncurves in n_list){
  
    sim_res_kmeans = matrix(NA,nsim,5)
    sim_res_sparse_concat= matrix(NA,nsim,5)
    sim_res_sparse= matrix(NA,nsim,5)
    sim_res_top0 = matrix(NA,nsim,5)
    sim_res_top1 = matrix(NA,nsim,5)
    sim_res_top2 = matrix(NA,nsim,5)
    
    for (sim in 1:nsim){
      actual = c(rep(1,ncurves),rep(2,ncurves),rep(3,ncurves),
                 rep(4,ncurves),rep(5,ncurves))
      print(paste("*****SIM NUM ",sim)  )
      print(paste("*****sigma,n",paste(sigma,ncurves))  )
      
      #SIMULATE DATA
      d = simulate.multidata(N=128,ncurves=ncurves,sigma=sigma,pad=64)
      #d = simulate.multidata(N=256,ncurves=ncurves,sigma=sigma,pad=128)
      
      
      #trans_data = transform(data,wavelet)
      trans_data = d[[1]]
      l = d[[2]]
      concat = do.call(cbind,trans_data)
      lconcat = do.call(cbind,l)
      
      wbound = try(KMeansSparseCluster.permute(concat,K = nclust,nperms=5))
      if(class(wbound)!="KMeansSparseCluster.permute"){
        print("ERROR")
        next;
      }
      wbound = wbound$bestw
      km.out = KMeansSparseCluster(concat,K=nclust,wbound,maxiter=50,nstart=100)
      randnaive = zapsmall(adjustedRand(km.out[[1]]$Cs, as.integer(actual)))
      print(randnaive)
      sim_res_sparse_concat[sim,] = randnaive
      m = multi_sparse(trans_data,lconcat,actual=actual,5,5,s=3)
      randsparse = m[[1]]
      randk = m[[2]]
      randtop0 = m[[3]]
      randtop1 = m[[4]]
      randtop2 = m[[5]]
      print(randk)
      print(randsparse)
      print(randtop0)
      print(randtop1)
      print(randtop2)
      
      sim_res_kmeans[sim,] = randk
      sim_res_sparse[sim,] = randsparse
      sim_res_top0[sim,] = randtop0
      sim_res_top1[sim,] = randtop1
      sim_res_top2[sim,] = randtop2

    }
    
#     sim_res_sparse_concat = sim_res_sparse_concat[complete.cases(sim_res_sparse_concat),]
#     write.csv(sim_res_sparse_concat,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_concat.csv'))
#     write.csv(sim_res_kmeans,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_kmeans.csv'))
#     write.csv(sim_res_sparse,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_sparse.csv'))
#     write.csv(sim_res_top0,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_sparsetop0.csv'))
#     write.csv(sim_res_top1,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_sparsetop1.csv'))
#     write.csv(sim_res_top2,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'NEWmulti_sparsetop2.csv'))
#     
}
}

# n_list_hmm = c(2,3)
# #sigma_list = 4
# #n_list=45
# nsim_hmm = 150
# for(sigma in sigma_list){
#   
#   for(ncurves in n_list_hmm){
#     
#     sim_res_hmm = matrix(NA,nsim_hmm,5)
#     sim_res_PCA = matrix(NA,nsim_hmm,5)
#     
#     for (sim in 1:nsim_hmm){
#       actual = c(rep(1,ncurves),rep(2,ncurves),rep(3,ncurves),
#                  rep(4,ncurves),rep(5,ncurves))
#       print(paste("*****SIM NUM ",sim)  )
#       print(paste("*****sigma,n",paste(sigma,ncurves))  )
#       
#       #SIMULATE DATA
#       d = simulate.multidata(N=128,ncurves=ncurves,sigma=sigma,pad=64)
#       #d = simulate.multidata(N=256,ncurves=ncurves,sigma=sigma,pad=128)
#       
#       l = d[[2]]
#       print("run HMM")
#       h = try(hmm_multi(l,actual,states=2))
#       print(h)
#       if (class(h) != "numeric"){
#         next;
#       }else{
#         sim_res_hmm[sim,] = h
#         
#       }
#       print("run PCA")
#       p = pca_clust(l,actual)
#       print(p)
#       sim_res_PCA[sim,] = p
#       
#       
#     }
#     sim_res_hmm = sim_res_hmm[complete.cases(sim_res_hmm),]
#     
#     write.csv(sim_res_hmm,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'hmm.csv'))
#     write.csv(sim_res_PCA,paste0('../Dropbox/wavelet/sim_res/n',ncurves,'_sigma',sigma,'pca.csv'))
#     
#   }
# }
