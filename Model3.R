library(mvtnorm)
library(refund)
set.seed(12358)
data_simu_M3 = t(sapply(1:100000, function(k){
  Z  = rnorm(6)
  
  # generate FPC scores
  A = 2*Z/(1:6)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  
  
  
  X = A[1]*sqrt(2)*sin(2*pi*tgrid) + A[2]*sqrt(2)*cos(2*pi*tgrid) +
    A[3]*sqrt(2)*sin(4*pi*tgrid) + A[4]*sqrt(2)*cos(4*pi*tgrid) +
    A[5]*sqrt(2)*sin(6*pi*tgrid) + A[6]*sqrt(2)*cos(6*pi*tgrid)
  
  C1 = 0.2*A[1] + rnorm(1, 0, 1)
  C2 = 0.1*A[2] + rnorm(1, 0, 1)
  C3 = 0.1*A[3]+ 0.1*A[4] + rnorm(1, 0, 1)
  
  conf = c(C1, C2, C3)
  # response
  Y = 1 + 2*A[1]+A[2]+0.5*A[3]+0.5*A[4] + C1 + C2 + C3 + rnorm(1, 0, 5)
  
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
    gbar = matrix(0, nrow = ntreat, ncol = ncof)
    
    for(i in 1:n){
      gbar = gbar + w.curr[i]*as.matrix(treat.star[i,])%*%t(as.matrix(conf.star[i,]))
    }
    
    # g = matrix(g, nrow = ntreat*ncof, byrow = TRUE)
    # gbar = apply(g, 2, mean)
    # gbar = matrix(gbar, nrow = ntreat, byrow = TRUE)
    # gbar <- c(w.curr.del,
    #           1/n*t(sample.weights)%*%((Ttilde - Xtilde%*%beta.curr)^2/sigmasq - 1))
    
    ##Generate mean imbalance.
    loss1<-sum((1/n*t(as.matrix(e.star))%*%as.matrix(e.star) - sigma.curr)^2) + sum((gbar/n)^2)
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
  J.opt = bal.loss(params.opt)
  R = list()
  R$weights = w.opt.scaled
  R$J = J.opt
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

FPCA = function(X, pve = 0.95){
  
  demeanX=as.matrix(X-matrix(apply(X,2,mean),nrow=nrow(X), ncol=ncol(X), byrow=TRUE))
  
  eigenX=eigen(cov(demeanX))
  
  L=which.max(cumsum(eigenX$values)/sum(eigenX$values)>0.95) # truncate at FVE=95%
  
  eigenfunX=eigenX$vectors[,1:L]
  
  eigenvalX=eigenX$values[1:L]/ncol(X)
  
  eigenfunX.rescale=eigenfunX*sqrt(ncol(X))
  
  scr = demeanX%*%eigenfunX.rescale/ncol(X)
  
  R = list()
  R$efn = eigenfunX.rescale
  R$eval = eigenvalX
  R$scr = scr
  return(R)
}

##
library(doParallel)

## Simplest Model
registerDoParallel(detectCores())
Model3 = foreach(B = 1:100, .combine = 'rbind') %dopar%{
  set.seed(B + 1000)
  data_raw = data_simu_M3[sample.int(nrow(data_simu_M3), 200, replace = F), ]  # sample 200 rows from data_simu
  
  Y = data_raw[,1]
  X = data_raw[,2:52]
  conf = data_raw[,53:55]
  
  n=nrow(X)
  
  tgrid=0:50/50
  ntgrid=length(tgrid)
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  
  ntreat=ncol(treat)
  
  nconf=ncol(conf)
  
  
  # PSW
  
  fcbps=FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.01)
  npfcbps.1 = npFCBPS.Functional.1(treat, conf, rho = 0.01)
  npfcbps.2 = npFCBPS.Functional.2(treat, conf, rho = 0.01)
  
  # # weighted correlation
  # 
  # xy=cbind(treat,conf)
  # 
  # zz.para=cov.wt(xy,wt=c(fcbps$weights),cor=TRUE, center=FALSE)
  # zz.np = cov.wt(xy, wt = c(npfcbps$w), cor=TRUE, center=FALSE)
  # zz.np.1 = cov.wt(xy, wt = c(npfcbps.1$w), cor=TRUE, center=FALSE)
  # zz.np.2 = cov.wt(xy, wt = c(npfcbps.2$w), cor=TRUE, center=FALSE)
  # 
  # Corw.para=zz.para$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  # Corw.np = zz.np$cor[1:ntreat,(ntreat+1):(ntreat+nconf)] 
  # Corw.np.1 = zz.np.1$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  # Corw.np.2 = zz.np.2$cor[1:ntreat,(ntreat+1):(ntreat+nconf)]
  # 
  # Coru=cor(cbind(treat,conf))[1:ntreat,(ntreat+1):(ntreat+nconf)]
  # 
  # # FLR
  # 
  # 
  # fitw.para=lm(Y~treat, weights=fcbps$weights) # with propensity score of parametric method
  # fitw.np=lm(Y~treat, weights = npfcbps$w)
  # fitw.np.1=lm(Y~treat, weights = npfcbps.1$w)
  # fitw.np.2=lm(Y~treat, weights = npfcbps.2$w)
  # 
  # 
  # fitu=lm(Y~treat) # without weighting
  # 
  # # treatment effect
  # 
  # beta.hatw.para=eigenfunX.rescale%*%fitw.para$coefficients[2:(1+ntreat)]
  # beta.hatw.np = eigenfunX.rescale%*%fitw.np$coefficients[2:(1+ntreat)]
  # beta.hatw.np.1 = eigenfunX.rescale%*%fitw.np.1$coefficients[2:(1+ntreat)]
  # beta.hatw.np.2 = eigenfunX.rescale%*%fitw.np.2$coefficients[2:(1+ntreat)]
  # 
  # beta.hatu=eigenfunX.rescale%*%fitu$coefficients[2:(1+ntreat)]
  # 
  # 
  # # find the mean of the absolute value of the correlation between treat and conf
  # absCorw.para = mean(abs(Corw.para))
  # absCorw.np = mean(abs(Corw.np))
  # absCorw.np.1 = mean(abs(Corw.np.1))
  # absCorw.np.2 = mean(abs(Corw.np.2))
  # 
  # absCoru = mean(abs(Coru))
  # 
  
  rbind(fcbps$weights, npfcbps$w, npfcbps.1$w, npfcbps.2)
  
}

save(Model3, file = "Model3.rda")