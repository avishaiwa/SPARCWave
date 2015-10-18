#http://www.inside-r.org/packages/cran/waveband/docs/test.data
#library(waveband)

#### SIMULATE
test.data<-function (type = "ppoly", n = 512, signal = 1, rsnr = 7, plotfn = FALSE)  {
  x <- seq(0, 1, length = n + 1)[1:n]
  if (type == "ppoly") {
    y <- rep(0, n)
    xsv <- (x <= 0.5)
    y[xsv] <- -16 * x[xsv]^3 + 12 * x[xsv]^2
    xsv <- (x > 0.5) & (x <= 0.75)
    y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 40 * x[xsv] + 28))/3 - 
      1.5
    xsv <- x > 0.75
    y[xsv] <- (x[xsv] * (16 * x[xsv]^2 - 32 * x[xsv] + 16))/3
  }
  else if (type == "blocks") {
    t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 
           0.76, 0.78, 0.81)
    h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
    y <- rep(0, n)
    for (i in seq(1, length(h))) {
      y <- y + (h[i] * (1 + sign(x - t[i])))/2
    }
  }
  else if (type == "bumps") {
    t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 
           0.76, 0.78, 0.81)
    h <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
    w <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 
           0.005, 0.008, 0.005)
    y <- rep(0, n)
    for (j in 1:length(t)) {
      y <- y + h[j]/(1 + abs((x - t[j])/w[j]))^4
    }
  }
  else if (type == "heavi") 
    y <- 4 * sin(4 * pi * x) - sign(x - 0.3) - sign(0.72 - 
                                                      x)
  else if (type == "doppler") {
    eps <- 0.05
    y <- sqrt(x * (1 - x)) * sin((2 * pi * (1 + eps))/(x + 
                                                         eps))
  }
  else {
    cat(c("test.data: unknown test function type", type, 
          "\n"))
    cat(c("Terminating\n"))
    return("NoType")
  }
  y <- y/sqrt(var(y)) * signal
  ynoise <- y + rnorm(n, 0, signal/rsnr)
  if (plotfn == TRUE) {
    if (type == "ppoly") 
      mlab = "Piecewise polynomial"
    if (type == "blocks") 
      mlab = "Blocks"
    if (type == "bumps") 
      mlab = "Bumps"
    if (type == "heavi") 
      mlab = "HeaviSine"
    if (type == "doppler") 
      mlab = "Doppler"
    plot(x, y, type = "l", lwd = 2, main = mlab, ylim = range(c(y, 
                                                                ynoise)),xlab="",ylab="",xaxt='n',yaxt='n')
    lines(x, ynoise, col = 2)
    lines(x, y)
  }
  return(list(x = x, y = y, ynoise = ynoise, type = type, rsnr = rsnr))
}

simulate.curves.pad<-function(N,sigma,n,pad){
  #N=512;snr=0.5
  #n=200
  #nclust=6
  st=0
  class0 = matrix(0,n,N+pad*2)
  class1 = matrix(0,n,N+pad*2);class2 =matrix(0,n,N+pad*2);class3 = matrix(0,n,N+pad*2);
  class4 = matrix(0,n,N+pad*2);class5 = matrix(0,n,N+pad*2);
  s=sigma
  for(i in 1:n){
    y0 = rnorm(N,0,s)
    shift = floor(runif(1)*st)
    
    y1 = test.data(type = "heavi", n = N, signal = 1, rsnr = Inf, plotfn = FALSE)
    y2 = test.data(type = "blocks", n = N, signal = 1, rsnr = Inf, plotfn = FALSE)
    y3 = test.data(type = "bumps", n = N, signal = 1, rsnr = Inf, plotfn = FALSE)
    y4 = test.data(type = "doppler", n = N, signal = 1, rsnr = Inf, plotfn = FALSE)
    y5 = test.data(type = "ppoly", n = N, signal = 1, rsnr = Inf, plotfn = FALSE)
    class0[i,] = c(rnorm(pad+shift,0,s), y0 + rnorm(N,0,s),rnorm(pad-shift,0,s))
    class1[i,] = c(rnorm(pad+shift,0,s),y1$ynoise + rnorm(N,0,s) ,rnorm(pad-shift,0,s))
    class2[i,] = c(rnorm(pad+shift,0,s),y2$ynoise + rnorm(N,0,s),rnorm(pad-shift,0,s))
    class3[i,] = c(rnorm(pad+shift,0,s),y3$ynoise + rnorm(N,0,s),rnorm(pad-shift,0,s))
    class4[i,] = c(rnorm(pad+shift,0,s),y4$ynoise + rnorm(N,0,s),rnorm(pad-shift,0,s))
    class5[i,] = c(rnorm(pad+shift,0,s),y5$ynoise + rnorm(N,0,s),rnorm(pad-shift,0,s))
    
  }
  
  data = rbind(class0,class1,class2,class3,class4,class5)
  return(data)
}

simulate.multiple.curves<-function(N,sigma,n,pad){
  
  A1 = simulate.curves.pad(N=N,sigma,n=n,pad=pad)
  A2 = simulate.curves.pad(N=N,sigma,n=n,pad=pad)
  A3 = simulate.curves.pad(N=N,sigma,n=n,pad=pad)
  A4 = simulate.curves.pad(N=N,sigma,n=n,pad=pad)
  A5 = simulate.curves.pad(N=N,sigma,n=n,pad=pad)
  
  l = list(); trans_list=list()
  l[[1]] = rbind(A1[1:n,],A2[1:n,],A3[(4*n+1):(5*n),],A4[(4*n+1):(5*n),],A5[(4*n+1):(5*n),])
  l[[2]] = rbind(A1[(n+1):(2*n),],A2[(n+1):(2*n),],A3[(n+1):(2*n),],A4[(n+1):(2*n),],A5[1:n,])
  l[[3]] = rbind(A1[(2*n+1):(3*n),],A2[(3*n+1):(4*n),],A3[(3*n+1):(4*n),],A4[(2*n+1):(3*n),],A5[(2*n+1):(3*n),])
  
  return(l)
}

simulate.multidata<-function(N,ncurves,sigma,wavelet = "s8",pad){
  
  trans_list = list()
  
  l = simulate.multiple.curves(N=N,sigma=sigma,n=ncurves,pad=pad)#20, 0.35 #snr=0.33 N=128
  for (g in 1:3){
    
    trans_data = transform(l[[g]],wavelet)
    trans_list[[g]] = trans_data[[1]]
  }
  concat = do.call(cbind,trans_list)
  
  return(list(trans_list,l,concat))
}

simulate.curves.list<-function(N,snr,n,g){
  curves_list = list()
  for(sim in 1:g){
    print(snr)
    curves_list[[sim]] = simulate.curves(N=N,sigma=snr,n=n)
  }
  return(curves_list)
}

multiple.curves.into.wavelets<-function(multiple_curves_list,wavelet){
  trans_list = list()
  trans_data = transform(multiple_curves_list[[1]],wavelet);  
  trans_list[[1]] = trans_data[[1]]
  trans_data = transform(multiple_curves_list[[2]],wavelet)
  trans_list[[2]] = trans_data[[1]]
  trans_data = transform(multiple_curves_list[[3]],wavelet)
  trans_list[[3]] = trans_data[[1]]
  return(trans_list)
}


###### SPARCL
library(sparcl)

R_sparse_kmeans = function(data,nclust,nperms,s=-1){
  
  if(s==-1){
    wbound = KMeansSparseCluster.permute(data,K = nclust,nperms=nperms,silent = TRUE)
    wbound = wbound$bestw
  }else{
    wbound=s
  }
  
  #print("sparse K-means")
  km.out = KMeansSparseCluster(data,K=nclust,wbounds=wbound,maxiter=20,nstart=100)
  predicted = km.out
  ret_list = list()
  ret_list[[1]] = km.out[[1]]$Cs
  ret_list[[2]] = km.out[[1]]$ws
  
  return(ret_list)  
}


################TRANSFORM

library(wmtsa)
soft.threshold<-function(z,a){
  sign(z)*ifelse(abs(z)>=a,(abs(z)-a),0)
}

transform <-function(data,wavelet){
  nrow = dim(data)[1]; ncol=dim(data)[2]
  N = dim(data)[2]
  dwt_coeffs = matrix(0,nrow,ncol)
  dwt_threshold_coeffs = matrix(0,nrow,ncol)
  for(i in 1:nrow(data)){
    x = data[i,]
    x = as.numeric(x)
    w = wavDWT(x,wavelet = wavelet)
    dwt_coeffs[i,] = unlist(w$data)
    sigma2 = (mad(w$data[[1]]))^2
    uni_thresh = sqrt(2*sigma2*log(N))
    for (s in 1:length(w$data)){
      w$data[[s]] = soft.threshold(w$data[[s]],uni_thresh)
    }
    dwt_threshold_coeffs[i,] = unlist(w$data)
  }
  trans_data = list(dwt_coeffs,dwt_threshold_coeffs)
  return(trans_data)
}

transform_list<-function(data,wavelet){
  trans_data = list()
  for(g in 1:length(data)){
    trans_data[[g]] = transform(data[[g]],wavelet)[[1]]
  }
  return(trans_data)
}

