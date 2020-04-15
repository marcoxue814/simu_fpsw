
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

# rm(list = ls())
# set.seed(123456)
# y = cbind(rnorm(50), rnorm(50), runif(50), rexp(50))
# x = cbind(2*rnorm(50), rnorm(50)^2, rnorm(50)*rexp(50))
# FCBPS(y,x)

