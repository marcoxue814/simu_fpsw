###### A Small Simulation Study for Functional Treatment

#### data generation

library(MASS)
library(mvtnorm)

Q=100 # simulation runs


simudata2=replicate(Q, list())


n=100 # sample size

rk=6 # rank of the true covariance function

nconf=3 # number of confounders




tgrid=0:50/50

ntgrid=length(tgrid)

beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)

mu = rep(0, rk)
sigma = diag(1, rk)


amp=2 # amplifier

for(q in 1:Q)
{
  
  set.seed(n+q)
  
  zz=mvrnorm(n, mu, sigma)
  
  conf=cbind(zz[,1]+2*zz[,2], zz[,2]^2-zz[,3]^2, exp(zz[,3])-exp(1/2)) # confounders
  
  # FPC.stan=zz[,1:rk] # standardized FPC scores
  
  # X=outer(amp*FPC.stan[,1]/1,sqrt(2)*sin(2*pi*tgrid)) + outer(amp*FPC.stan[,2]/2,sqrt(2)*cos(2*pi*tgrid)) +
  #   outer(amp*FPC.stan[,3]/3,sqrt(2)*sin(4*pi*tgrid)) + outer(amp*FPC.stan[,4]/4,sqrt(2)*cos(4*pi*tgrid)) +
  #   outer(amp*FPC.stan[,5]/5,sqrt(2)*sin(6*pi*tgrid)) + outer(amp*FPC.stan[,1]/6,sqrt(2)*cos(6*pi*tgrid))
  
  X=outer(amp*zz[,1]/1,sqrt(2)*sin(2*pi*tgrid)) + outer(amp*zz[,2]/2,sqrt(2)*cos(2*pi*tgrid)) +
    outer(amp*zz[,3]/3,sqrt(2)*sin(4*pi*tgrid)) + outer(amp*zz[,4]/4,sqrt(2)*cos(4*pi*tgrid)) +
    outer(amp*zz[,5]/5,sqrt(2)*sin(6*pi*tgrid)) + outer(amp*zz[,1]/6,sqrt(2)*cos(6*pi*tgrid))
  
  set.seed(Q*n+q)
  
  Y=1+X%*%beta/ntgrid + conf%*%rep(0.2,nconf) + rnorm(n, 0, sqrt(0.5))
  
  
  
  simudata2[[q]]=list(Y=Y, X=X, conf=conf)
  
}
#### PSW


simufit2_para=replicate(Q, list())

for(q in 1:Q)
{
  cat("q = ", q = q, "\n")
  
  Y=simudata2[[q]]$Y
  X=simudata2[[q]]$X
  conf=simudata2[[q]]$conf
  
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
  
  t(fcbps$weights)
  
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  
  Covw=zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  #Covw
  
  #max(abs(c(zz$cov[1:ntreat,(ntreat+1):(ntreat+nconf)])))
  
  Covu=cov(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  #Covu
  
  # treat.star=scale(treat) # treatment is a set of FPC scores, which are uncorrelated.
  #
  # conf.star=conf%*%solve(chol(var(conf)))
  #
  # conf.star=scale(conf.star,center=TRUE, scale=TRUE)
  #
  # Sigma0.star=cov(cbind(treat.star,conf.star))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  #
  # Covw.star=t(treat.star)%*%diag(c(npfcbps$w))%*%conf.star
  #
  # xy.star=cbind(treat.star,conf.star)
  #
  # zz=cov.wt(xy.star,wt=c(npfcbps$w),cor=TRUE, center=FALSE)
  
  
  # FLR
  
  dataset=data.frame(cbind(treat, conf, Y))
  
  colnames(dataset)[1+L+nconf]='Y'
  
  fitw=lm(Y~., data=dataset, weights=fcbps$weights) # with propensity score
  
  summary(fitw)
  
  beta.hatw=eigenfunX.rescale%*%fitw$coefficients[2:(1+L)]
  
  # fitu=lm(Y~., data=dataset) # WITHOUT propensity score
  #
  # summary(fitu)
  #
  # beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+L)]
  # #
  # #
  # plot(tgrid, beta,ylim=c(-4,4), type='l', col=1)
  # points(tgrid, beta.hatw, type='l', col=2)
  # points(tgrid, beta.hatu, type='l', col=3)
  
  
  # fitting results
  
  simufit2_para[[q]]$Covw=Covw
  simufit2_para[[q]]$Covu=Covu
  
  simufit2_para[[q]]$beta.hatw=beta.hatw
  
}




#### summarize results

load('simuData.RData')

Q=100 # simulation runs

n=100 # sample size

tgrid=0:50/50 # time grid

ntgrid=length(tgrid)

beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)

L2.beta=mean(beta^2)

absCovw=rep(NA,Q)
#replicate(Q, list())

absCovu=rep(NA,Q)
#replicate(Q, list())

ISE=rep(NA,Q)


for(q in 1:Q)
{
  
  absCovw[q]=max(abs(simufit2_para[[q]]$Covw))
  
  absCovu[q]=max(abs(simufit2_para[[q]]$Covu))
  
  ISE[q]=mean((simufit2_para[[q]]$beta.hatw-beta)^2)
  
}

sapply(simufit2_para, function(x) nrow(x$Covw)) # L=4


grp.weight=factor(c(rep('Weighted',Q), rep('Unweighted', Q)))

absCov=c(absCovw, absCovu)

absCov=as.data.frame(absCov)

absCov.cb=cbind(absCov,grp.weight)

par(mfrow=c(1,2))
boxplot(absCov ~ grp.weight, data=absCov.cb, horizontal = TRUE, 
        main='Parametric Covariate Balancing')
boxplot(ISE/L2.beta, horizontal = TRUE, ylab='RISE',
        main='Causal Effect')

median(absCovw)

median(absCovu)
