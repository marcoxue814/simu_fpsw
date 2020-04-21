library(mvtnorm)


#------------------------------------------
# Setting 1 n = 100 
#------------------------------------------
Q = 1000
simudata1 = replicate(Q, list())

for(q in 1:Q){
  n = 100
  set.seed(n+q)
  
  Z = rmvnorm(n, rep(0, 10), diag(10))
  
  # generate FPC scores
  A1 = Z[,1] + Z[,2] - 4*Z[,3] + Z[,4] + Z[,5] + 4*Z[,7] + rnorm(n, 0, 2)
  A2 = 1.8*Z[,1] + 0.9*Z[,2] - 0.9*Z[,4] - 1.8*Z[,5] + rnorm(n, 0, 1)
  A3 = Z[,1] - Z[,2] -Z[,4] + Z[,5] + rnorm(n, 0, 0.5)
  A4 = 0.85*Z[,3] - 0.85*Z[,6] + 0.85*Z[,7] + rnorm(n, 0, 0.25)
  A5 = 0.5*Z[,1] + 0.5*Z[,2] + 0.5*Z[,3] + 0.5*Z[,4] + 0.5*Z[,5] + 0.5*Z[,6] + rnorm(n, 0, 0.125)
  A6 = 0.3*Z[,2] - 0.6*Z[,2] + 0.6*Z[,4] - 0.3*Z[,5] + rnorm(n, 0, 0.0625)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
    0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)
  
  # functional treatment
  X=outer(A1,sqrt(2)*sin(2*pi*tgrid)) + outer(A2,sqrt(2)*cos(2*pi*tgrid)) + 
    outer(A3,sqrt(2)*sin(4*pi*tgrid)) + outer(A4,sqrt(2)*cos(4*pi*tgrid)) + 
    outer(A5,sqrt(2)*sin(6*pi*tgrid)) + outer(A6,sqrt(2)*cos(6*pi*tgrid))
  
  # cofounders
  conf=Z[,3:5]
  nconf = ncol(conf)
  
  set.seed(Q*n + q)
  # response
  Y=1+X%*%beta/ntgrid+conf%*%rep(0.2,nconf) + rnorm(n, 0, 5)
  
  simudata1[[q]]=list(Y=Y, X=X, conf=conf)
}

FCBPS = function(treat, conf, iterations = 1000){
  
  treat = as.matrix(treat)
  
  conf = as.matrix(conf)
  
  n = nrow(treat)
  
  ntreat = ncol(treat)
  ncof = ncol(conf)
  treat.star = scale(treat)
  
  conf.star=conf%*%solve(chol(var(conf)))
  
  conf.star=scale(conf.star,center=TRUE, scale=TRUE)
  
  
  ##Run linear regression
  fit1 = lm(treat.star ~ -1 + conf.star)
  mcoef = as.matrix(coef(fit1))
  mcoef[is.na(mcoef)] = 0
  e.star = treat.star - conf.star%*%mcoef
  sigma = 1/n*t(e.star)%*%e.star
  params.curr = c(mcoef, sigma[upper.tri(sigma, diag = T)])
  
  probs.min = 1e-6
  stabilizers<-log(apply(treat.star,1, function(t) mean(pmin(pmax(dmvnorm(t, mean = rep(0,ntreat), sigma = diag(ntreat)), probs.min), 1-probs.min))))
  
  
  bal.func<-function(params.curr){
    
    
    beta.curr = matrix(params.curr[1:(ncof*ntreat)], nrow = ncof, byrow = FALSE)
    sigma.curr = diag(ntreat)
    sigma.curr[upper.tri(sigma.curr, diag = T)] = params.curr[(1+ncof*ntreat):length(params.curr)]
    sigma.curr[which(upper.tri(t(sigma.curr), diag = T),arr.ind = T)[,c(2,1)]] = params.curr[(1+ncof*ntreat):length(params.curr)]
    
    #probs.curr<-log(apply(treat.star, 1, function(t) dmvnorm(t, mean = conf.star%*%beta.curr, sigma = sigma.curr)))
    probs.curr<-rep(0, n)
    for(i in 1:n){
      probs.curr[i] <- dmvnorm(treat.star[i,], mean = (conf.star%*%beta.curr)[i,], sigma = sigma.curr, log = TRUE)
    }
    probs.curr<-pmin(log(1-probs.min),probs.curr)
    probs.curr<-pmax(log(probs.min),probs.curr)
    
    w.curr<-exp(stabilizers - probs.curr)
    
    ##Generate the vector of mean imbalance by weights.
    # w.curr.del<-1/n*t(wtXilde)%*%w.curr
    # w.curr.del<-as.matrix(w.curr.del)
    # w.curr<-as.matrix(w.curr)
    e.star = treat.star - conf.star%*%beta.curr
    g = NULL
    
    for(i in 1:n){
      g = c(g, c(w.curr[i]*as.matrix(treat.star[i,])%*%t(as.matrix(conf.star[i,]))))
    }
    
    g = matrix(g, nrow = 50, byrow = TRUE)
    gbar = apply(g, 2, mean)
    gbar = matrix(gbar, nrow = ntreat, byrow = TRUE)
    # gbar <- c(w.curr.del,
    #           1/n*t(sample.weights)%*%((Ttilde - Xtilde%*%beta.curr)^2/sigmasq - 1))
    
    ##Generate mean imbalance.
    loss1<-sum((1/n*t(as.matrix(e.star))%*%as.matrix(e.star) - sigma.curr)^2) + sum(gbar^2)
    out1<-list("loss"=loss1)
    out1
  }
  
  bal.loss<-function(x,...) bal.func(x,...)$loss
  
  
  gmm.init = params.curr
  
  opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS")
  params.opt = opt.bal$par
  beta.opt = matrix(params.opt[1:(ncof*ntreat)], nrow = ncof, byrow = FALSE)
  sigma.opt = diag(ntreat)
  sigma.opt[upper.tri(sigma.opt, diag = T)] = params.opt[(1+ncof*ntreat):length(params.opt)]
  sigma.opt[which(upper.tri(t(sigma.opt), diag = T),arr.ind = T)[,c(2,1)]] = params.opt[(1+ncof*ntreat):length(params.opt)]
  
  
  #probs.opt<-log(apply(treat.star, 1, function(t) dmvnorm(t, mean = conf.star%*%beta.opt, sigma = sigma.opt)))
  probs.opt<-rep(0,n)
  for(i in 1:n){
    probs.opt[i]<-dmvnorm(treat.star[i,], mean = (conf.star%*%beta.opt)[i,], sigma = sigma.opt, log = TRUE)
  }
  probs.opt<-pmin(log(1-probs.min),probs.opt)
  probs.opt<-pmax(log(probs.min),probs.opt)
  
  w.opt<-exp(stabilizers - probs.opt)
  w.opt.scaled = w.opt/sum(w.opt)
  value = opt.bal$value
  R = list()
  R$weights = w.opt.scaled 
  R$value = value
  class(R) = "FCBPS"
  
  return(R)
}

simu_para.1 = replicate(Q, list())
library(mvtnorm)

for(q in 1:Q)
{
  cat("q = ", q = q, "\n")
  
  Y=simudata1[[q]]$Y
  X=simudata1[[q]]$X
  conf=simudata1[[q]]$conf
  
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


save(simu_para.1, file = 'simu_para_1.rda')


