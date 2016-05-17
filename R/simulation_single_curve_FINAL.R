library(clues);library(wmtsa);library(sparcl);library(flexclust)
source('C:/Users/tomhope/Dropbox/wavelet/src/noise_NEW.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/transform_data.R')
source('C:/Users/tomhope/Dropbox/wavelet/src/wave_sparse_Kmeans.R')
nclust=6
#n=15
#sigma=3
wavelet ="s8"
#data = simulate.curves(N=512,snr=0.2,n=n)

#sigma_list = seq(2.25,3.75,0.25)
#sigma_list = 3.75
sigma_list = (c(2.5,2.75))
n_list = (c(2,3,4,5,8,10,12,15,20,25))
n_list = c(10)
sigma_list = 2.5
#n_list = (c(2,3,4))
#n_list = 5
#n_list = c(15,20)
#n_list = c(25,30)
nsim = 1
#sigma_list = 2.75
#n_list=5
for(sigma in sigma_list){
  
for(n in n_list){

sim_res_kmeans = matrix(NA,nsim,5)
sim_res_sparse_orig = matrix(NA,nsim,5)

sim_res_sparse1 = matrix(NA,nsim,5)
sim_res_sparse2 = matrix(NA,nsim,5)
num_feats1 = matrix(NA,nsim,1)
num_feats2 = matrix(NA,nsim,1)

sim_res_top1 = matrix(NA,nsim,5)
sim_res_top2 = matrix(NA,nsim,5)

for (sim in 1:nsim){
  actual = c(rep(0,n),rep(1,n),rep(2,n),rep(3,n),
             rep(4,n),rep(5,n))
  
  print(paste("*****SIM NUM ",sim)  )
  print(paste("*****sigma,n",paste(sigma,n))  )
  
  #KMEANS
  data = simulate.curves.pad(N=256,sigma=sigma,n=n,pad=128)#600
  trans_data = transform(data,wavelet)
  dwt_coeffs_shrink = trans_data[[2]]
  #print("running kmeans")
  km = kmeans(data,centers = nclust,iter.max = 20,nstart = 100)
  #table(km$cluster,actual)
  randk = zapsmall(adjustedRand(km$cluster, as.integer(actual)))
  #print("RANDK")
  #print(randk)
  sim_res_kmeans[sim,] = randk

#   #SPARSE ORIGINAL
#   wbound = try(KMeansSparseCluster.permute(data,K = nclust,nperms=5))
#   if (class(wbound) != "KMeansSparseCluster.permute"){
#     next;
#   }
#   wbound_orig = wbound$bestw
#   sparse_orig =  try(KMeansSparseCluster(data,K=nclust,wbound=wbound_orig,maxiter=50,nstart=100))
#   if (class(sparse_orig) != "KMeansSparseCluster"){
#     next;
#   }
  #ws = sparse_orig[[1]]$ws
  
  #print("number of features");print(length(ws[ws>0]))
  
  #wb = sparse_orig[[1]]$wbound
  #randsparse_orig = zapsmall(adjustedRand(sparse_orig[[1]]$Cs, as.integer(actual)))
  #print("sparse ORIG");print(randsparse_orig)
  #sim_res_sparse_orig[sim,] = randsparse_orig

  #SPARSE
  #median(apply(trans_data[[1]],2,mad))

  #wbound = KMeansSparseCluster.permute(trans_data[[1]],K = nclust,nperms=5)
  #wbound = wbound$bestw
 
  f_list = c()
  grid = seq(1.1,4,0.5)
  sparse_list = list()
  index = 0
  for(s in grid){
    index=index+1
    #print(s)
    sparse = try(wave_sparse_kmeans(data = data, wavelet = wavelet,nclust=nclust,
                                    wbound = s))
    if (class(sparse) != "list"){
      next;
    }
    kmobj = sparse[[1]]
    ws = kmobj[[1]]$ws
    sparse_list[[index]] = sparse
    #print("number of features");print(length(ws[ws>0]))
    f_list = c(f_list,length(ws[ws>0]))
  }

  elbow_list = c()
  for(i in 1:(length(f_list)-1)){
    s = f_list[i+1]/f_list[i]
    elbow_list = c(elbow_list,s)
  }
  best = grid[which.max(elbow_list)]
  print(best)
  sparse = try(wave_sparse_kmeans(data = data, wavelet = wavelet,nclust=nclust,
                                  wbound = best))
  if(class(sparse) != "list"){
    next;
  }
  kmobj = sparse[[1]]; dwtcoeffs = sparse[[2]]
  ws = kmobj[[1]]$ws
  
  print("number of features");print(length(ws[ws>0]))
  num_feats1[sim,] = length(ws[ws>0])
  wb = kmobj[[1]]$wbound
  randsparse = zapsmall(adjustedRand(kmobj[[1]]$Cs, as.integer(actual)))
  #print("RANDS");print(randsparse)
  sim_res_sparse1[sim,] = randsparse
  sortedw = sort(ws,decreasing = TRUE)
 
  top_kfeat00 = which(sortedw>0)[length(which(sortedw>0))]
  topw2 = order(ws,decreasing = TRUE)[1:top_kfeat00]
  top_features2 = dwtcoeffs[,topw2]
  kmtop2 = kmeans(top_features2,centers = nclust,iter.max = 20, nstart = 100)

  randtop2 = zapsmall(adjustedRand(kmtop2$cluster, as.integer(actual)))
  #print("RANDTOP")
  print(randtop2)
  sim_res_top1[sim,] = randtop2
  
  #SPARSE
  sparse = try(wave_sparse_kmeans(data = data, wavelet = wavelet,nclust=nclust,
                                  wbound = NULL))
  if (class(sparse) != "list"){
    next;
  }
  kmobj = sparse[[1]]; dwtcoeffs = sparse[[2]]
  ws = kmobj[[1]]$ws
  #print("number of features");print(length(ws[ws>0]))
  num_feats2[sim,] = length(ws[ws>0])
  
  #table(kmobj[[1]]$Cs,actual)
  wb = kmobj[[1]]$wbound
  randsparse = zapsmall(adjustedRand(kmobj[[1]]$Cs, as.integer(actual)))
  print("RANDS");print(randsparse)
  sim_res_sparse2[sim,] = randsparse
  sortedw = sort(ws,decreasing = TRUE)
  
  top_kfeat00 = which(sortedw>0)[length(which(sortedw>0))]
  
  topw2 = order(ws,decreasing=TRUE)[1:top_kfeat00]
  top_features2 = dwtcoeffs[,topw2]
  kmtop2 = kmeans(top_features2,centers = nclust,iter.max = 20, nstart = 100)
  
  randtop2 = zapsmall(adjustedRand(kmtop2$cluster, as.integer(actual)))
  #print("RANDTOP")
  #print(randtop2)
  sim_res_top2[sim,] = randtop2
  #print("BOUND")
  #print(median(apply(data,2,sd)))
  #print(median(apply(trans_data[[1]],2,sd)))
  #print(min(apply(trans_data[[1]],2,sd)))
}

 write.csv(sim_res_kmeans,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'kmeans.csv'))
# write.csv(sim_res_sparse_orig,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'sparse_orig.csv'))
# write.csv(sim_res_sparse1,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'sparse_mad.csv'))
# write.csv(sim_res_sparse2,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'sparse_sd.csv'))
# write.csv(sim_res_top1,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'new_top1.csv'))
# write.csv(sim_res_top2,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'new_top2.csv'))
# write.csv(num_feats1,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'numfeats1.csv'))
# write.csv(num_feats2,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'numfeats2.csv'))
}

}
print(apply(sim_res_kmeans,2,mean)[3])
print(apply(sim_res_sparse2,2,mean)[3])
print(apply(sim_res_top1,2,mean)[3])
print(apply(sim_res_top2,2,mean)[3])



