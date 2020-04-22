
set.seed(1001)
data_simu_M1 = t(sapply(1:100000, function(k){
  Z  = rnorm(10)
  
  # generate FPC scores
  A1 = Z[1] + Z[2] - 4*Z[3] + Z[4] + Z[5] + 4*Z[7] + rnorm(1, 0, 2)
  A2 = 1.8*Z[1] + 0.9*Z[2] - 0.9*Z[4] - 1.8*Z[5] + rnorm(1, 0, 1)
  A3 = Z[1] - Z[2] -Z[4] + Z[5] + rnorm(1, 0, 0.5)
  A4 = 0.85*Z[3] - 0.85*Z[6] + 0.85*Z[7] + rnorm(1, 0, 0.25)
  A5 = 0.5*Z[1] + 0.5*Z[2] + 0.5*Z[3] + 0.5*Z[4] + 0.5*Z[5] + 0.5*Z[6] + rnorm(1, 0, 0.125)
  A6 = 0.3*Z[2] - 0.6*Z[2] + 0.6*Z[4] - 0.3*Z[5] + rnorm(1, 0, 0.0625)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
    0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)
  
  X = A1*sqrt(2)*sin(2*pi*tgrid) + A2*sqrt(2)*cos(2*pi*tgrid) + 
      A3*sqrt(2)*sin(4*pi*tgrid) + A4*sqrt(2)*cos(4*pi*tgrid) +
      A5*sqrt(2)*sin(6*pi*tgrid) + A6*sqrt(2)*cos(6*pi*tgrid)
  
  conf = Z[3:5]
  # response
  Y=1+X%*%beta/ntgrid + 0.2*Z[3] + 0.2*Z[4] + 0.2*Z[5] + rnorm(1, 0, 5)
  
  c(Y, X, conf)
}))


##
library(doParallel)

## Simplest Model
registerDoParallel(detectCores())
Model1_np = foreach(B = 1:1000, .combine = 'rbind') %dopar%{
  set.seed(B + 1000)
  data_raw = data_simu[sample.int(nrow(data_simu_M1), 100, replace = F), ]  # sample 200 rows from data_simu
  
  Y = data_raw[,1]
  X = data_raw[,2:52]
  conf = data_raw[,53:55]
  
  n=nrow(X)
  
  # get FPC scores
  demeanX=as.matrix(X-matrix(apply(X,2,mean),nrow=n, ncol=ncol(X), byrow=TRUE))
  
  eigenX=eigen(cov(demeanX))
  
  L=which.max(cumsum(eigenX$values)/sum(eigenX$values)>0.95) # truncate at FVE=95%
  
  eigenfunX=eigenX$vectors[,1:L]
  
  eigenfunX.rescale=eigenfunX*sqrt(ntgrid)
  
  treat=demeanX%*%eigenfunX.rescale/ntgrid
  
  ntreat=ncol(treat)
  
  nconf=ncol(conf)
  
  # PSW
  
  npfcbps=npFCBPS.Functional.1(treat, conf, rho = 0.01)
  
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz=cov.wt(xy,wt=c(npfcbps$w),cor=TRUE, center=FALSE)
  
  Covw=zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  Covu=cov(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # FLR
  
  dataset=data.frame(cbind(treat, conf, Y))
  
  colnames(dataset)[1+L+nconf]='Y'
  
  fitw=lm(Y~., data=dataset, weights=npfcbps$w) # with propensity score
  
  fitu=lm(Y~., data=dataset) # without weighting
  
  beta.hatw=eigenfunX.rescale%*%fitw$coefficients[2:(1+L)]
  
  beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+L)]
  
  c(Covw, Covu, beta.hatw, beta.hatu)
  
}

save(Model1_np, file = "Model1_np.rda")
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  