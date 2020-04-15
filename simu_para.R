
Q = 100
simu_para.1 = replicate(Q, list())
library(mvtnorm)

for(q in 1:Q)
{
  cat("q = ", q = q, "\n")
  
  Y=simudata.1[[q]]$Y
  X=simudata.1[[q]]$X
  conf=simudata.1[[q]]$conf
  
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
  
  fcbps=FCBPS(treat, conf)
  
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  
  Covw=zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  Covu=cov(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # FLR
  
  dataset=data.frame(cbind(treat, conf, Y))
  
  colnames(dataset)[1+L+nconf]='Y'
  
  fitw=lm(Y~., data=dataset, weights=fcbps$weights) # with propensity score
  
  fitu=lm(Y~., data=dataset) # without weighting
  
  beta.hatw=eigenfunX.rescale%*%fitw$coefficients[2:(1+L)]
  
  beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+L)]
  
  
  
  
  # fitting results
  
  simu_para.1[[q]]$Covw=Covw
  simu_para.1[[q]]$Covu=Covu
  
  simu_para.1[[q]]$beta.hatw=beta.hatw
  simu_para.1[[q]]$beta.hatu=beta.hatu
}

#---------------------------------
#Simu Setting 2
#---------------------------------
simu_para.2 = replicate(Q, list())
for(q in 1:Q)
{
  cat("q = ", q = q, "\n")
  
  Y=simudata.2[[q]]$Y
  X=simudata.2[[q]]$X
  conf=simudata.2[[q]]$conf
  
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
  
  fcbps=FCBPS(treat, conf)
  
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  
  Covw=zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  Covu=cov(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # FLR
  
  dataset=data.frame(cbind(treat, conf, Y))
  
  colnames(dataset)[1+L+nconf]='Y'
  
  fitw=lm(Y~., data=dataset, weights=fcbps$weights) # with propensity score
  
  fitu=lm(Y~., data=dataset) # without weighting
  
  beta.hatw=eigenfunX.rescale%*%fitw$coefficients[2:(1+L)]
  
  beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+L)]
  
  
  
  
  # fitting results
  
  simu_para.2[[q]]$Covw=Covw
  simu_para.2[[q]]$Covu=Covu
  
  simu_para.2[[q]]$beta.hatw=beta.hatw
  simu_para.2[[q]]$beta.hatu=beta.hatu
}

#---------------------------------
#Simu Setting 3
#---------------------------------

simu_para.3 = replicate(Q, list())
for(q in 1:Q)
{
  cat("q = ", q = q, "\n")
  
  Y=simudata.3[[q]]$Y
  X=simudata.3[[q]]$X
  conf=simudata.3[[q]]$conf
  
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
  
  fcbps=FCBPS(treat, conf)
  
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  
  Covw=zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  Covu=cov(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # FLR
  
  dataset=data.frame(cbind(treat, conf, Y))
  
  colnames(dataset)[1+L+nconf]='Y'
  
  fitw=lm(Y~., data=dataset, weights=fcbps$weights) # with propensity score
  
  fitu=lm(Y~., data=dataset) # without weighting
  
  beta.hatw=eigenfunX.rescale%*%fitw$coefficients[2:(1+L)]
  
  beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+L)]
  
  
  
  
  # fitting results
  
  simu_para.3[[q]]$Covw=Covw
  simu_para.3[[q]]$Covu=Covu
  
  simu_para.3[[q]]$beta.hatw=beta.hatw
  simu_para.3[[q]]$beta.hatu=beta.hatu
}



