library(mvtnorm)
library(refund)
set.seed(1001)
data_simu_M3 = t(sapply(1:100000, function(k){
  Z  = rnorm(10)
  
  # generate FPC scores
  A1 = Z[1] + Z[2] - 4*Z[3] + Z[4] + Z[5] + 4*Z[7] + rnorm(1, 0, 2)
  A2 = 1.8*Z[1] + 0.9*Z[2] - 0.9*Z[4] - 1.8*Z[5] + rnorm(1, 0, 1)
  A3 = Z[1] - Z[2] -Z[4] + Z[5] + rnorm(1, 0, 0.5)
  A4 = 0.85*Z[3] - 0.85*Z[6] + 0.85*Z[7] + rnorm(1, 0, 0.25)
  A5 = 0.5*Z[1] + 0.5*Z[2] + 0.5*Z[3] + 0.5*Z[4] + 0.5*Z[5] + 0.5*Z[6] + rnorm(1, 0, 0.125)
  A6 = 0.3*Z[1] - 0.6*Z[2] + 0.6*Z[4] - 0.3*Z[5] + rnorm(1, 0, 0.0625)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
    0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)
  
  X = A1*sqrt(2)*sin(2*pi*tgrid) + A2*sqrt(2)*cos(2*pi*tgrid) + 
    A3*sqrt(2)*sin(4*pi*tgrid) + A4*sqrt(2)*cos(4*pi*tgrid) +
    A5*sqrt(2)*sin(6*pi*tgrid) + A6*sqrt(2)*cos(6*pi*tgrid)
  
  conf = c(A1+2*A2, A2-A3, A1+A3+A4)
  # response
  Y=1+X%*%beta/ntgrid + 0.2*conf[1] + 0.2*conf[2] + 0.2*conf[3] + rnorm(1, 0, 5)
  
  c(Y, X, conf)
}))


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

npFCBPS.Functional<-function(treat, conf, #method = 'exact', 
                             rho = 0.01, ...)
{
  
  treat=as.matrix(treat)
  
  conf=as.matrix(conf)
  
  n=nrow(treat)
  
  ## normalize treatment and confounder
  
  ntreat=ncol(treat) # dimension of treatment, i.e., number of FPC scores
  
  nconf=ncol(conf)# dimension of confounder
  
  treat.star=scale(treat) # treatment is a set of FPC scores, which are uncorrelated.
  
  
  conf.star=conf%*%solve(chol(var(conf)))
  
  conf.star=scale(conf.star,center=TRUE, scale=TRUE)
  
  #conf.star=as.matrix(conf)%*%t(chol(solve(cov(conf)))) # orthogonalize and standardize confounders
  
  K=ntreat + nconf + ntreat*nconf
  
  ## Empirical likelihood; 
  ## Equivalent form using Langrangian multiplier
  ## optim function
  # if(method=='exact')
  # {
  #   gfun.matrix=cbind(treat.star, conf.star,
  #                     t(sapply(1:n,function(i){outer(treat.star[i,],conf.star[i,])})))
  #   
  #   K=ncol(gfun.matrix)
  #   
  #   obj=function(gam, eps=1/n)
  #   {
  #     arg=n-gfun.matrix%*%gam
  #     
  #     idx=arg < eps
  #     
  #     ans=arg
  #     ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
  #     ans[!idx] = log(arg[!idx])
  #     
  #     -sum(ans)
  #   }
  #   
  #   gam.ini=rep(0,K)
  #   opt.gam=optim(par=gam.ini, fn=obj, method='BFGS')
  #   w=1/(n-gfun.matrix%*%opt.gam$par)
  #   sum(w)
  #   w/sum(w)
  # }
  # else
  # {
  
  Sigma0=cor(cbind(treat.star,conf.star))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # eta_prior_sd=rep(rho,ntreat*nconf)
  
  # K=ncol(gfun.matrix)
  
  eps=1/n
  
  hfun_given_theta=function(theta)
  {
    cbind(treat.star, conf.star,t(sapply(1:n,function(i)
    {
      outer(treat.star[i,],conf.star[i,])-theta*Sigma0
    })))
  }
  
  
  objgam_given_theta=function(gam, theta, eps)
  {
    arg=n-hfun_given_theta(theta)%*%gam
    
    idx=arg < eps
    
    ans=arg
    ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
    ans[!idx] = log(arg[!idx])
    
    -sum(ans)
  }
  
  getw_given_theta=function(theta, eps, sumw.tol=.001)
  {
    gam.ini=rep(0,K)
    gam.opt=optim(par=gam.ini, fn=objgam_given_theta, method='BFGS', theta=theta, eps=eps)
    w=1/(n-hfun_given_theta(theta)%*%gam.opt$par)
    sumw=sum(w) 
    w_rescaled=w/sumw
    
    if (abs(1-sumw)<=sumw.tol){log_el=-sum(log(w_rescaled))}
    if (abs(1-sumw)>=sumw.tol){log_el=-sum(log(w_rescaled))-10^4*(1+abs(1-sumw))}
    
    # log_el=-sum(log(w_rescaled))
    
    R=list()
    R$w=w_rescaled
    R$sumw=sumw
    R$log_el=log_el
    return(R)
  }
  
  objtheta = function(theta, eps=eps, sumw.tol=.001)
  { 
    
    el.inner=getw_given_theta(theta, eps=eps, sumw.tol=sumw.tol)
    
    c=1
    
    # el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*rho)
    
    el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*n^2*rho)
    
  }
  
  theta.opt=optimize(f=objtheta, interval=c(-1,1), eps=eps, sumw.tol=.001, maximum=TRUE)
  
  el.inner.opt=getw_given_theta(theta=theta.opt$maximum, eps=eps, sumw.tol=.001)
  
  w.opt=el.inner.opt$w
  sumw.opt=el.inner.opt$sumw
  inobj.opt=el.inner.opt$log_el
  # outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*rho)
  
  outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*n^2*rho)
  #objtheta(theta=theta.opt,eps=eps,sumw.tol=.001)
  
  
  R=list()
  R$w=w.opt
  R$sumw=sumw.opt
  R$inobj=inobj.opt
  R$outobj=outobj.opt
  
  class(R) = "npFCBPS"
  
  return(R)
  
}

npFCBPS.Functional.1<-function(treat, conf, #method = 'exact', 
                               rho = 0.01, ...)
{
  
  treat=as.matrix(treat)
  
  conf=as.matrix(conf)
  
  n=nrow(treat)
  
  ## normalize treatment and confounder
  
  ntreat=ncol(treat) # dimension of treatment, i.e., number of FPC scores
  
  nconf=ncol(conf)# dimension of confounder
  
  treat.star=scale(treat) # treatment is a set of FPC scores, which are uncorrelated.
  
  
  conf.star=conf%*%solve(chol(var(conf)))
  
  conf.star=scale(conf.star,center=TRUE, scale=TRUE)
  
  #conf.star=as.matrix(conf)%*%t(chol(solve(cov(conf)))) # orthogonalize and standardize confounders
  
  K=ntreat + nconf + 2*ntreat*nconf
  
  ## Empirical likelihood; 
  ## Equivalent form using Langrangian multiplier
  ## optim function
  # if(method=='exact')
  # {
  #   gfun.matrix=cbind(treat.star, conf.star,
  #                     t(sapply(1:n,function(i){outer(treat.star[i,],conf.star[i,])})))
  #   
  #   K=ncol(gfun.matrix)
  #   
  #   obj=function(gam, eps=1/n)
  #   {
  #     arg=n-gfun.matrix%*%gam
  #     
  #     idx=arg < eps
  #     
  #     ans=arg
  #     ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
  #     ans[!idx] = log(arg[!idx])
  #     
  #     -sum(ans)
  #   }
  #   
  #   gam.ini=rep(0,K)
  #   opt.gam=optim(par=gam.ini, fn=obj, method='BFGS')
  #   w=1/(n-gfun.matrix%*%opt.gam$par)
  #   sum(w)
  #   w/sum(w)
  # }
  # else
  # {
  
  Sigma0=cor(cbind(treat.star,conf.star))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # eta_prior_sd=rep(rho,ntreat*nconf)
  
  # K=ncol(gfun.matrix)
  
  eps=1/n
  
  hfun_given_theta=function(theta)
  {
    cbind(treat.star, conf.star, t(sapply(1:n,function(i)
    {
      outer(treat.star[i,],conf.star[i,])-theta*Sigma0
    })), t(sapply(1:n,function(i)
    {
      outer(treat.star[i,],conf.star[i,]^2)
    })))
  }
  
  
  
  objgam_given_theta=function(gam, theta, eps)
  {
    arg=n-hfun_given_theta(theta)%*%gam
    
    idx=arg < eps
    
    ans=arg
    ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
    ans[!idx] = log(arg[!idx])
    
    -sum(ans)
  }
  
  getw_given_theta=function(theta, eps, sumw.tol=.001)
  {
    gam.ini=rep(0,K)
    gam.opt=optim(par=gam.ini, fn=objgam_given_theta, method='BFGS', theta=theta, eps=eps)
    w=1/(n-hfun_given_theta(theta)%*%gam.opt$par)
    sumw=sum(w) 
    w_rescaled=w/sumw
    
    if (abs(1-sumw)<=sumw.tol){log_el=-sum(log(w_rescaled))}
    if (abs(1-sumw)>=sumw.tol){log_el=-sum(log(w_rescaled))-10^4*(1+abs(1-sumw))}
    
    # log_el=-sum(log(w_rescaled))
    
    R=list()
    R$w=w_rescaled
    R$sumw=sumw
    R$log_el=log_el
    return(R)
  }
  
  objtheta = function(theta, eps=eps, sumw.tol=.001)
  { 
    
    el.inner=getw_given_theta(theta, eps=eps, sumw.tol=sumw.tol)
    
    c=1
    
    # el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*rho)
    
    el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*n^2*rho)
  }
  
  theta.opt=optimize(f=objtheta, interval=c(-1,1), eps=eps, sumw.tol=.001, maximum=TRUE)
  
  el.inner.opt=getw_given_theta(theta=theta.opt$maximum, eps=eps, sumw.tol=.001)
  
  w.opt=el.inner.opt$w
  sumw.opt=el.inner.opt$sumw
  inobj.opt=el.inner.opt$log_el
  # outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*rho)
  outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*n^2*rho)
  #objtheta(theta=theta.opt,eps=eps,sumw.tol=.001)
  
  
  R=list()
  R$w=w.opt
  R$sumw=sumw.opt
  R$inobj=inobj.opt
  R$outobj=outobj.opt
  
  class(R) = "npFCBPS"
  
  return(R)
  
}

npFCBPS.Functional.2<-function(treat, conf, #method = 'exact', 
                               rho = 0.01, ...)
{
  
  treat=as.matrix(treat)
  
  conf=as.matrix(conf)
  
  n=nrow(treat)
  
  ## normalize treatment and confounder
  
  ntreat=ncol(treat) # dimension of treatment, i.e., number of FPC scores
  
  nconf=ncol(conf)# dimension of confounder
  
  treat.star=scale(treat) # treatment is a set of FPC scores, which are uncorrelated.
  
  
  conf.star=conf%*%solve(chol(var(conf)))
  
  conf.star=scale(conf.star,center=TRUE, scale=TRUE)
  
  #conf.star=as.matrix(conf)%*%t(chol(solve(cov(conf)))) # orthogonalize and standardize confounders
  
  K=ntreat + nconf + 3*ntreat*nconf
  
  ## Empirical likelihood; 
  ## Equivalent form using Langrangian multiplier
  ## optim function
  # if(method=='exact')
  # {
  #   gfun.matrix=cbind(treat.star, conf.star,
  #                     t(sapply(1:n,function(i){outer(treat.star[i,],conf.star[i,])})))
  #   
  #   K=ncol(gfun.matrix)
  #   
  #   obj=function(gam, eps=1/n)
  #   {
  #     arg=n-gfun.matrix%*%gam
  #     
  #     idx=arg < eps
  #     
  #     ans=arg
  #     ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
  #     ans[!idx] = log(arg[!idx])
  #     
  #     -sum(ans)
  #   }
  #   
  #   gam.ini=rep(0,K)
  #   opt.gam=optim(par=gam.ini, fn=obj, method='BFGS')
  #   w=1/(n-gfun.matrix%*%opt.gam$par)
  #   sum(w)
  #   w/sum(w)
  # }
  # else
  # {
  
  Sigma0=cor(cbind(treat.star,conf.star))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # eta_prior_sd=rep(rho,ntreat*nconf)
  
  # K=ncol(gfun.matrix)
  
  eps=1/n
  
  hfun_given_theta=function(theta)
  {
    cbind(treat.star, conf.star, t(sapply(1:n,function(i)
    {
      outer(treat.star[i,],conf.star[i,])-theta*Sigma0
    })), t(sapply(1:n,function(i)
    {
      outer(treat.star[i,], conf.star[i,]^2)
    })), t(sapply(1:n,function(i)
    {
      outer(treat.star[i,]^2,conf.star[i,])
    }))
    )
  }
  
  # if(condition = 0){
  #   hfun_given_theta=function(theta)
  #   {
  #     cbind(treat.star, conf.star, t(sapply(1:n,function(i)
  #     {
  #       outer(treat.star[i,],conf.star[i,])-theta*Sigma0
  #     })))
  #   }
  # }else if(condition = 1){
  #   hfun_given_theta=function(theta)
  #   {
  #     cbind(treat.star, conf.star, treat.star%*%t(treat.star), t(sapply(1:n,function(i)
  #     {
  #       outer(treat.star[i,],conf.star[i,])-theta*Sigma0
  #     })))
  #   }
  # }else if(condition = 2){
  #   hfun_given_theta=function(theta)
  #   {
  #     cbind(treat.star, conf.star, treat.star%*%t(treat.star), conf.star%*%t(conf.star), t(sapply(1:n,function(i)
  #     {
  #       outer(treat.star[i,],conf.star[i,])-theta*Sigma0
  #     })))
  #   }
  # }
  
  
  objgam_given_theta=function(gam, theta, eps)
  {
    arg=n-hfun_given_theta(theta)%*%gam
    
    idx=arg < eps
    
    ans=arg
    ans[idx] = log(eps) - 1.5 + 2 * arg[idx]/eps - 0.5 * (arg[idx]/eps)^2
    ans[!idx] = log(arg[!idx])
    
    -sum(ans)
  }
  
  getw_given_theta=function(theta, eps, sumw.tol=.001)
  {
    gam.ini=rep(0,K)
    gam.opt=optim(par=gam.ini, fn=objgam_given_theta, method='BFGS', theta=theta, eps=eps)
    w=1/(n-hfun_given_theta(theta)%*%gam.opt$par)
    sumw=sum(w) 
    w_rescaled=w/sumw
    
    if (abs(1-sumw)<=sumw.tol){log_el=-sum(log(w_rescaled))}
    if (abs(1-sumw)>=sumw.tol){log_el=-sum(log(w_rescaled))-10^4*(1+abs(1-sumw))}
    
    # log_el=-sum(log(w_rescaled))
    
    R=list()
    R$w=w_rescaled
    R$sumw=sumw
    R$log_el=log_el
    return(R)
  }
  
  objtheta = function(theta, eps=eps, sumw.tol=.001)
  { 
    
    el.inner=getw_given_theta(theta, eps=eps, sumw.tol=sumw.tol)
    
    c=1
    
    # el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*rho)
    
    el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*n^2*rho)
    
  }
  
  theta.opt=optimize(f=objtheta, interval=c(-1,1), eps=eps, sumw.tol=.001, maximum=TRUE)
  
  el.inner.opt=getw_given_theta(theta=theta.opt$maximum, eps=eps, sumw.tol=.001)
  
  w.opt=el.inner.opt$w
  sumw.opt=el.inner.opt$sumw
  inobj.opt=el.inner.opt$log_el
  # outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*rho)
  outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*n^2*rho)
  #objtheta(theta=theta.opt,eps=eps,sumw.tol=.001)
  
  
  R=list()
  R$w=w.opt
  R$sumw=sumw.opt
  R$inobj=inobj.opt
  R$outobj=outobj.opt
  
  class(R) = "npFCBPS"
  
  return(R)
  
}

##
library(doParallel)

## Simplest Model
registerDoParallel(detectCores())
Model3 = foreach(B = 1:1000, .combine = 'rbind') %dopar%{
  set.seed(B + 1000)
  data_raw = data_simu_M3[sample.int(nrow(data_simu_M3), 200, replace = F), ]  # sample 100 rows from data_simu
  
  Y = data_raw[,1]
  X = data_raw[,2:52]
  conf = data_raw[,53:55]
  
  n=nrow(X)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  
  # get FPC scores
  demeanX=as.matrix(X-matrix(apply(X,2,mean),nrow=n, ncol=ncol(X), byrow=TRUE))
  
  # eigenX=eigen(cov(demeanX))
  
  # L=which.max(cumsum(eigenX$values)/sum(eigenX$values)>0.95) # truncate at FVE=95%
  
  # eigenfunX=eigenX$vectors[,1:L]
  
  # eigenfunX.rescale=eigenfunX*sqrt(ntgrid)
  
  
  pcscores = fpca.sc(demeanX, pve = 0.95)$scores
  
  eigenfunX = fpca.sc(demeanX, pve = 0.95)$efunctions
  
  eigenfunX.rescale=eigenfunX*sqrt(ntgrid)
  
  eigenvalX = fpca.sc(demeanX, pve = 0.95)$evalues
  
  L = ncol(pcscores)
  
  treat=demeanX%*%eigenfunX.rescale/ntgrid
  
  ntreat=ncol(treat)
  
  nconf=ncol(conf)
  
  
  # PSW
  
  fcbps=FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.01)
  npfcbps.1 = npFCBPS.Functional.1(treat, conf, rho = 0.01)
  npfcbps.2 = npFCBPS.Functional.2(treat, cnof, rho = 0.01)
  
  # weighted correlation
  
  xy=cbind(treat,conf)
  
  zz.para=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  zz.np = cov.wt(xy, wt = c(npfcbps$w), cor=TRUE, center=FALSE)
  zz.np.1 = cov.wt(xy, wt = c(npfcbps.1$w), cor=TRUE, center=FALSE)
  zz.np.2 = cov.wt(xy, wt = c(npfcbps.2$w), cor=TRUE, center=FALSE)
  
  Corw.para=zz.para$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  Corw.np = zz.np$cor[1:ntreat,(ntreat+1):(ntreat+nconf)] 
  Corw.np.1 = zz.np.1$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  Corw.np.2 = zz.np.2$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  Coru=cor(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  
  # FLR
  
  
  fitw.para=lm(Y~treat, weights=fcbps$weights) # with propensity score of parametric method
  fitw.np=lm(Y~treat, weights = npfcbps$w)
  fitw.np.1=lm(Y~treat, weights = npfcbps.1$w)
  fitw.np.2=lm(Y~treat, weights = npfcbps.2$w)
  
  
  fitu=lm(Y~treat) # without weighting
  
  # treatment effect
  
  beta.hatw.para=eigenfunX.rescale%*%fitw.para$coefficients[2:(1+ntreat)]
  beta.hatw.np = eigenfunX.rescale%*%fitw.np$coefficients[2:(1+ntreat)]
  beta.hatw.np.1 = eigenfunX.rescale%*%fitw.np.1$coefficients[2:(1+ntreat)]
  beta.hatw.np.2 = eigenfunX.rescale%*%fitw.np.2$coefficients[2:(1+ntreat)]
  
  beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+ntreat)]
  
  
  # find the mean of the absolute value of the correlation between treat and conf
  absCorw.para = mean(abs(Corw.para))
  absCorw.np = mean(abs(Corw.np))
  absCorw.np.1 = mean(abs(Corw.np.1))
  absCorw.np.2 = mean(abs(Corw.np.2))
  
  absCoru = mean(abs(Coru))
  
  
  c(absCorw.para, absCorw.np, absCorw.np.1, absCorw.np.2, absCoru, beta.hatw.para, beta.hatw.np, beta.hatw.np.1, beta.hatw.np.2, beta.hatu)
  
}

save(Model3, file = "Model3.rda")