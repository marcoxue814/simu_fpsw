########## Nonparametric CBPS for Functional Treatment ###########
### INPUTS: 
### treat:    functional principal component (FPC) scores of the functional treatment, 
###               a n by L matrix, where n is the sample size, and L is the truncated # of FPCs
### conf:     confounder matrix, a n by D matrix, which may includes FPC scores of 
###               functional confounders.
### method:   'exact' (default) or 'penalized';



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
      
      el.inner$log_el - c*theta^2*sum(c(Sigma0)^2)/(2*rho)
      
    }
    
    theta.opt=optimize(f=objtheta, interval=c(-1,1), eps=eps, sumw.tol=.001, maximum=TRUE)
    
    el.inner.opt=getw_given_theta(theta=theta.opt$maximum, eps=eps, sumw.tol=.001)
    
    w.opt=el.inner.opt$w
    sumw.opt=el.inner.opt$sumw
    inobj.opt=el.inner.opt$log_el
    outobj.opt=inobj.opt - theta.opt$maximum^2*sum(c(Sigma0)^2)/(2*rho)
      #objtheta(theta=theta.opt,eps=eps,sumw.tol=.001)

    
    R=list()
    R$w=w.opt
    R$sumw=sumw.opt
    R$inobj=inobj.opt
    R$outobj=outobj.opt
    
    class(R) = "npFCBPS"
    
    return(R)
  
}
