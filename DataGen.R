
library(MASS)

Q=100 # simulation runs


simudata.1=replicate(Q, list())
simudata.2=replicate(Q, list())
simudata.3=replicate(Q, list())

n=100 # sample size

rk=6 # rank of the true covariance function

nconf=3 # number of confounders

mu=rep(0, rk+nconf)

Var.FPC=diag(rep(1,rk))

Cov.FPCvsConf=matrix(0.3, nrow=rk, ncol=nconf)

Var.Conf=matrix(c(1,0.6,0.6,0.6,1,0.6,0.6,0.6,1), nrow=nconf, ncol=nconf)

# Cov.FPCvsConf=matrix(0.6, nrow=rk, ncol=nconf)
#
# Var.Conf=matrix(c(1,0.7,0.7,0.7,1,0.7,0.7,0.7,1), nrow=nconf, ncol=nconf)

Sigma=rbind(cbind(Var.FPC, Cov.FPCvsConf), cbind(t(Cov.FPCvsConf), Var.Conf))

eigen(Sigma)

tgrid=0:50/50

ntgrid=length(tgrid)

beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid)

#beta= 0.2*sqrt(2)*sin(2*pi*tgrid) -0.1* sqrt(2)*cos(2*pi*tgrid) +
#  0.2*sqrt(2)*sin(4*pi*tgrid) #+ 0.5*sqrt(2)*cos(4*pi*tgrid)


plot(tgrid,beta,type='l')
#----------------------------
# Simu Setting 1
#----------------------------
amp=2 # amplifier

for(q in 1:Q)
{
  
  set.seed(n+q)
  
  zz=mvrnorm(n, mu, Sigma)
  
  conf=zz[,(rk+1):(rk+nconf)] # confounders
  
  FPC.stan=zz[,1:rk] # standardized FPC scores
  
  X=outer(amp*FPC.stan[,1]/1,sqrt(2)*sin(2*pi*tgrid)) + outer(amp*FPC.stan[,2]/2,sqrt(2)*cos(2*pi*tgrid)) +
    outer(amp*FPC.stan[,3]/3,sqrt(2)*sin(4*pi*tgrid)) + outer(amp*FPC.stan[,4]/4,sqrt(2)*cos(4*pi*tgrid)) +
    outer(amp*FPC.stan[,5]/5,sqrt(2)*sin(6*pi*tgrid)) + outer(amp*FPC.stan[,1]/6,sqrt(2)*cos(6*pi*tgrid))
  
  
  set.seed(Q*n+q)
  
  Y=1+X%*%beta/ntgrid + conf%*%rep(0.2,nconf) + rnorm(n, 0, sqrt(0.5))
  
  
  
  simudata.1[[q]]=list(Y=Y, X=X, conf=conf)
  
}

#--------------------------
#Simu Setting 2
#--------------------------
for(q in 1:Q)
{
  set.seed(n+2*q)
  z = mvrnorm(n, rep(0,6), diag(6))
  
  conf = cbind(z[,1]+2*z[,2], z[,2]^2-z[,3]^2, exp(z[,3])-exp(1/2))
  
  X = outer(amp*z[,1]/1, sqrt(2)*sin(2*pi*tgrid)) + outer(amp*z[,2]/2, sqrt(2)*cos(2*pi*tgrid)) + 
      outer(amp*z[,3]/3, sqrt(2)*sin(4*pi*tgrid)) + outer(amp*z[,4]/4, sqrt(2)*cos(4*pi*tgrid)) +
      outer(amp*z[,5]/5, sqrt(2)*sin(6*pi*tgrid)) + outer(amp*z[,6]/6, sqrt(2)*cos(6*pi*tgrid))
    
  set.seed(Q*n + 2*q)
  Y=1+X%*%beta/ntgrid + conf%*%rep(0.2,nconf) + rnorm(n, 0, sqrt(0.5))
  
  simudata.2[[q]]=list(Y=Y, X=X, conf=conf)
  
}

#--------------------------
#Simu Setting 3
#--------------------------

for(q in 1:Q)
{
  set.seed(n+3*q)
  z = mvrnorm(n, rep(0,6), diag(6))
  
  conf = cbind(z[,1]+2*z[,2], z[,2]^2-z[,3]^2, exp(z[,3])-exp(1/2))
  
  X = outer(amp*z[,1]/1, sqrt(2)*sin(2*pi*tgrid)) + outer(amp*z[,2]/2, sqrt(2)*cos(2*pi*tgrid)) + 
      outer(amp*z[,3]/3, sqrt(2)*sin(4*pi*tgrid)) + outer(amp*z[,4]/4, sqrt(2)*cos(4*pi*tgrid)) +
      outer(amp*z[,5]/5, sqrt(2)*sin(6*pi*tgrid)) + outer(amp*z[,6]/6, sqrt(2)*cos(6*pi*tgrid))
  
  set.seed(Q*n + 3*q)
  Y=1+X%*%beta/ntgrid + z[,1:3]%*%rep(0.5, 3) + rnorm(n, 0, sqrt(0.5))
  
  simudata.3[[q]] = list(Y=Y, X=X, conf = conf)
}





