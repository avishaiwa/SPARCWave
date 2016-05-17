library(clues);library(wmtsa);library(sparcl);library(flexclust)
source('C:/Users/tomhope/Dropbox/wavelet/src/noise_NEW.R')
library(wavethresh)
for(f in list.files('C:/Users/tomhope/Dropbox/wavelet/src/Other Algorithms/curvclust/R/')){
  source(paste0('C:/Users/tomhope/Dropbox/wavelet/src/Other Algorithms/curvclust/R/',f))
  
}
nclust=6
sigma_list = (c(2.5,2.75))
n_list = (c(5,10,12,15,20,25))
n_list = (c(2,3,4))
n_list = 8
nsim = 150
for(sigma in sigma_list){
  
  for(n in n_list){
    sim_res_giacof = matrix(NA,nsim,5)
    
    for (sim in 1:nsim){
      actual = c(rep(0,n),rep(1,n),rep(2,n),rep(3,n),
                 rep(4,n),rep(5,n))
      print(paste("*****SIM NUM ",sim)  )
      print(paste("*****sigma,n",paste(sigma,n))  )
      
      data = simulate.curves.pad(N=256,sigma=sigma,n=n,pad=128)
      curvs = convert.curves.to.curvclust(data)
      CCD = new("CClustData",Y=curvs,filter.number=8)
      CCDred = getUnionCoef(CCD)
      # options setting
      CCO = new("CClustO")
      CCO["nbclust"] = nclust
      CCO["Gamma2.structure"] = "none"
      CCO["burn"] = 100
      CCR = getFCM(CCDred,CCO)
      CCRclust = apply(zapsmall(CCR@Tau),1,which.max)
      #table(CCRclust,actual)
      randCCR = zapsmall(adjustedRand(CCRclust, as.integer(actual)))
      print(randCCR)
      sim_res_giacof[sim,] = randCCR

    }
    write.csv(sim_res_giacof,paste0('../Dropbox/wavelet/sim_res/single/n',n,'_sigma',sigma,'giacof.csv'))
  }
}

