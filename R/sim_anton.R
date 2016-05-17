library(clues);library(wmtsa);library(sparcl);library(flexclust); library(kernlab)

#source('noise_NEW.R')
#source('transform_data.R')
#source('kmeans.R')

nclust=6
#n=15
#sigma=3
wavelet ="s8"
#data = simulate.curves(N=512,snr=0.2,n=n)

#sigma_list = seq(2.25,3.75,0.25)
#sigma_list = 3.75
sigma_list = (c(2.5,2.75))
n_list = (c(2:5,10,12,15,20,25))
n_list=8
#n_list = c(15,20)
#n_list = c(25,30)

nsim = 150
#sigma_list = 2.5
#n_list=45
for(sigma in sigma_list){
  
  for(n in n_list){
    
    sim_res_abs_pam = matrix(NA,nsim,5)
    sim_res_rel_pam = matrix(NA,nsim,5)  
    sim_res_abs_km= matrix(NA,nsim,5)
           
    for (sim in 1:nsim){
      actual = c(rep(0,n),rep(1,n),rep(2,n),rep(3,n),
                 rep(4,n),rep(5,n))
      print(paste("*****SIM NUM ",sim)  )
      print(paste("*****sigma,n",paste(sigma,n))  )
      
      #KMEANS
      data = simulate.curves.pad(N=256,sigma=sigma,n=n,pad=128)
     
            
      mat_dwt <- t(apply(data,1, toDWT))                    # DWT over the rows
      mat_contr_abs <- t(apply(mat_dwt, 1, contrib))         
      mat_contr_rel <- t(apply(mat_dwt, 1, contrib, rel = TRUE))
      
        
      mat_pam_abs <- pam(mat_contr_abs, k = 6)
      mat_pam_rel <- pam(mat_contr_rel, k = 6)
      kms = kmeans(mat_contr_abs,centers = nclust,nstart = 100)
      
      randks = zapsmall(adjustedRand(kms$cluster, as.integer(actual)))
      randpam_abs = zapsmall(adjustedRand( mat_pam_abs$clustering, as.integer(actual)))
      randpam_rel = zapsmall(adjustedRand( mat_pam_rel$clustering, as.integer(actual)))
      
      print(randpam_abs); print(randpam_rel);       print(randks)
      sim_res_abs_pam[sim,] = randpam_abs; sim_res_rel_pam[sim,] = randpam_rel; sim_res_abs_km[sim,] = randks;
                     
    }
    
#     write.csv(sim_res_abs_pam,paste0('../sim_res/single/n',n,'_sigma',sigma,'abs_pam.csv'))
#     write.csv(sim_res_rel_pam,paste0('../sim_res/single/n',n,'_sigma',sigma,'res_pam.csv'))
#     write.csv(sim_res_abs_km,paste0('../sim_res/single/n',n,'_sigma',sigma,'abs_km.csv'))
    write.csv(sim_res_abs_pam,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'abs_pam.csv'))
    write.csv(sim_res_rel_pam,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'res_pam.csv'))
    write.csv(sim_res_abs_km,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'abs_km.csv'))
    
    
    
  }
}
