
#-----------------------------------
# Some functions needed
#-----------------------------------

# to get fpc scores
getFPC = function(X){
  N = nrow(X)
  A = matrix(NA, nrow = N, ncol = 6) #generate 6 FPC scores
  A[,1] = sqrt(16)*X[,1]
  A[,2] = sqrt(12)*X[,2]
  A[,3] = sqrt(8)*X[,3]
  A[,4] = sqrt(4)*X[,4]
  A[,5] = sqrt(1)*X[,5]
  A[,6] = sqrt(1/2)*X[,6]
  R = list()
  R$A = A
  return(R)
}

# to get functional variable
getX = function(FPC, tgrid){
  X = outer(FPC[,1], sqrt(2)*sin(2*pi*tgrid)) + outer(FPC[,2], sqrt(2)*cos(2*pi*tgrid)) +
    outer(FPC[,3], sqrt(2)*sin(4*pi*tgrid)) + outer(FPC[,4], sqrt(2)*cos(4*pi*tgrid)) +
    outer(FPC[,5], sqrt(2)*sin(6*pi*tgrid)) + outer(FPC[,6], sqrt(2)*cos(6*pi*tgrid))
  R = list()
  R$X = X
  return(R)
}

# to get covariates, nonlinear in fpc scores
getC_nonln = function(X){
  N = nrow(X)
  C = matrix(NA, nrow = N, ncol = 3) #generate 3 confounders
  C[,1] = (X[,1] + 0.5)^2 + rnorm(N, 0, 1)
  C[,2] = 0.2*X[,2] + rnorm(N, 0, 0.5)
  C[,3] = 0.2*X[,3] + rnorm(N, 0, 0.5)
  R = list()
  R$C = C
  return(R)
}

# to get covariates, linear in fpc scores
getC_ln = function(X){
  N = nrow(X)
  C = matrix(NA, nrow = N, ncol = 3) #generate 3 confounders
  C[,1] = X[,1] + rnorm(N, 0, 1)
  C[,2] = 0.2*X[,2] + rnorm(N, 0, 0.5)
  C[,3] = 0.2*X[,3] + rnorm(N, 0, 0.5)
  R = list()
  R$C = C
  return(R)
}

# to get outcome nonlinear in fpc scores 
getY_nonln = function(FPC, conf){
  #FPC: the six FPC scores
  #C: confounders
  N = nrow(C)
  Y = 1 + 2*FPC[,1] + 1*FPC[,2] + 0.5*FPC[,3] + 0.5*FPC[,4] + 2*conf[,1] + conf[,2]^2 +  5*rnorm(N)
  R = list()
  R$Y = Y
  return(R)
}

# to get outcome linear in fpc scores
getY_ln = function(FPC, conf){
  #FPC: the six FPC scores
  #C: confounders
  N = nrow(C)
  Y = 1 + 2*FPC[,1] + 1*FPC[,2] + 0.5*FPC[,3] + 0.5*FPC[,4] + 2*conf[,1] +  5*rnorm(N)
  R = list()
  R$Y = Y
  return(R)
}

# function to do fpca
FPCA = function(X, pve = 0.95){
  
  demeanX=as.matrix(X-matrix(apply(X,2,mean),nrow=nrow(X), ncol=ncol(X), byrow=TRUE))
  
  eigenX=eigen(cov(demeanX))
  
  L=which.max(cumsum(eigenX$values)/sum(eigenX$values)>pve) # truncate at FVE=95%
  
  eigenfunX=eigenX$vectors[,1:L]
  
  eigenvalX=eigenX$values[1:L]/ncol(X)
  
  eigenfunX.rescale=eigenfunX*sqrt(ncol(X))
  
  scr = demeanX%*%eigenfunX.rescale/ncol(X)
  
  R = list()
  R$perc = eigenX$values/sum(eigenX$values)
  R$efn = eigenfunX.rescale
  R$eval = eigenvalX
  R$scr = scr
  return(R)
}

#------------------------------------------
# Y|C is linear & X|C is linear pve = 0.95
#------------------------------------------
Q = 200
results_1 = replicate(Q, list())
simu_data_1 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu_data_1[[q]]$Y = Y
  simu_data_1[[q]]$X = X
  simu_data_1[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_1[[q]]$w.para = fcbps$weights
  results_1[[q]]$w.j = fcbps$J
  results_1[[q]]$w.np = npfcbps$w
  results_1[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_1, file = "simu_data_1.rda")
save(results_1, file = "results_1.rda")

#------------------------------------------
# Y|C is linear & X|C is linear pve = 0.99
#------------------------------------------
Q = 200
results_2 = replicate(Q, list())
simu_data_1 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu_data_1[[q]]$Y = Y
  simu_data_1[[q]]$X = X
  simu_data_1[[q]]$conf = conf
  
  pve = 0.99
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_2[[q]]$w.para = fcbps$weights
  results_2[[q]]$w.j = fcbps$J
  results_2[[q]]$w.np = npfcbps$w
  results_2[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_1, file = "simu_data_1.rda")
save(results_2, file = "results_2.rda")

#---------------------------------------------
# Y|C is linear & X|C is nonlinear pve = 0.95
#---------------------------------------------
Q = 200
results_3 = replicate(Q, list())
simu_data_2 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_nonln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu_data_2[[q]]$Y = Y
  simu_data_2[[q]]$X = X
  simu_data_2[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_3[[q]]$w.para = fcbps$weights
  results_3[[q]]$w.j = fcbps$J
  results_3[[q]]$w.np = npfcbps$w
  results_3[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_2, file = "simu_data_2.rda")
save(results_3, file = "results_3.rda")

#---------------------------------------------
# Y|C is linear & X|C is nonlinear pve = 0.99
#---------------------------------------------
Q = 200
results_4 = replicate(Q, list())
simu_data_2 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_nonln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu_data_2[[q]]$Y = Y
  simu_data_2[[q]]$X = X
  simu_data_2[[q]]$conf = conf
  
  pve = 0.99
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_4[[q]]$w.para = fcbps$weights
  results_4[[q]]$w.j = fcbps$J
  results_4[[q]]$w.np = npfcbps$w
  results_4[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_2, file = "simu_data_2.rda")
save(results_4, file = "results_4.rda")

#---------------------------------------------
# Y|C is nonlinear & X|C is linear pve = 0.95
#---------------------------------------------
Q = 200
results_5 = replicate(Q, list())
simu_data_3 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_nonln(FPC, conf)$Y
  simu_data_3[[q]]$Y = Y
  simu_data_3[[q]]$X = X
  simu_data_3[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_5[[q]]$w.para = fcbps$weights
  results_5[[q]]$w.j = fcbps$J
  results_5[[q]]$w.np = npfcbps$w
  results_5[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_3, file = "simu_data_3.rda")
save(results_5, file = "results_5.rda")

#---------------------------------------------
# Y|C is nonlinear & X|C is linear pve = 0.99
#---------------------------------------------
Q = 200
results_6 = replicate(Q, list())
simu_data_3 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_nonln(FPC, conf)$Y
  simu_data_3[[q]]$Y = Y
  simu_data_3[[q]]$X = X
  simu_data_3[[q]]$conf = conf
  
  pve = 0.99
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_6[[q]]$w.para = fcbps$weights
  results_6[[q]]$w.j = fcbps$J
  results_6[[q]]$w.np = npfcbps$w
  results_6[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_3, file = "simu_data_3.rda")
save(results_6, file = "results_6.rda")

#------------------------------------------------
# Y|C is nonlinear & X|C is nonlinear pve = 0.95
#------------------------------------------------
Q = 200
results_7 = replicate(Q, list())
simu_data_4 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_nonln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_nonln(FPC, conf)$Y
  simu_data_4[[q]]$Y = Y
  simu_data_4[[q]]$X = X
  simu_data_4[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_7[[q]]$w.para = fcbps$weights
  results_7[[q]]$w.j = fcbps$J
  results_7[[q]]$w.np = npfcbps$w
  results_7[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_4, file = "simu_data_4.rda")
save(results_7, file = "results_7.rda")

#------------------------------------------------
# Y|C is nonlinear & X|C is nonlinear pve = 0.99
#------------------------------------------------
Q = 200
results_8 = replicate(Q, list())
simu_data_4 = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for(q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 6
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_nonln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_nonln(FPC, conf)$Y
  simu_data_4[[q]]$Y = Y
  simu_data_4[[q]]$X = X
  simu_data_4[[q]]$conf = conf
  
  pve = 0.99
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_8[[q]]$w.para = fcbps$weights
  results_8[[q]]$w.j = fcbps$J
  results_8[[q]]$w.np = npfcbps$w
  results_8[[q]]$w.theta = npfcbps$theta
  
  setTxtProgressBar(pb, q)
}


save(simu_data_4, file = "simu_data_4.rda")
save(results_8, file = "results_8.rda")

#----------------------------------
# Simu Results Analysis
#----------------------------------

Q = 200
beta.hatw.para = replicate(Q, list())
beta.hatw.np = replicate(Q, list())
beta.hatu = replicate(Q, list())
# np.theta = rep(0, Q)
# para.J = rep(0, Q)
# L = rep(0,Q)
for(q in 1:Q){
  w.para = results_1[[q]]$w.para
  w.np = results_1[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  # L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  
  beta.hatw.para[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu[[q]] = efnX%*%fitu$coefficients[-1]
  
  # np.theta[q] = results_3[[q]]$w.theta
  # para.J[q] = results_3[[q]]$w.j
  
}

beta.hatw.para.m = NULL
for(q in 1:Q){
  beta.hatw.para.m = cbind(beta.hatw.para.m, beta.hatw.para[[q]])
}

beta.hatw.np.m = NULL
for(q in 1:Q){
  beta.hatw.np.m = cbind(beta.hatw.np.m, beta.hatw.np[[q]])
}

beta.hatu.m = NULL
for(q in 1:Q){
  beta.hatu.m = cbind(beta.hatu.m, beta.hatu[[q]])
}

beta.hatw.para.mean = apply(beta.hatw.para.m, 1, mean)
beta.hatw.np.mean = apply(beta.hatw.np.m, 1, mean)
beta.hatu.mean = apply(beta.hatu.m, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para = (beta.hatw.para.mean - beta)
bias.np = (beta.hatw.np.mean - beta)
bias.u = (beta.hatu.mean - beta)

ISE.para = apply(beta.hatw.para.m, 2, function(X) mean((X-beta)^2))
ISE.np = apply(beta.hatw.np.m, 2, function(X) mean((X-beta)^2))
ISE.u = apply(beta.hatu.m, 2, function(X) mean((X-beta)^2))



table1 = cbind(rbind(median(ISE.u),median(ISE.para), median(ISE.np)), 
               rbind(mean(ISE.u), mean(ISE.para), mean(ISE.np)),
               rbind(mean(bias.u^2), mean(bias.para^2), mean(bias.np^2)))
rownames(table1) = c("Unweighted", "para", "Np")
colnames(table1) = c("Median", "IMSE", "Ibias^2")

# boxplot.matrix(cbind(ISE.u, ISE.para,ISE.np),
#                names = c("Unweighted", "Para", "Np"),
#                ylab = "ISE", cex.main=2, cex.lab=1.5, cex.axis=1.5,
#                main = expression('PVE'[L]*' = 0.95,'*' PVE'[L^"*"]*' = 0.95'))


beta.hatw.para = replicate(Q, list())
beta.hatw.np = replicate(Q, list())
beta.hatu = replicate(Q, list())
# np.theta = rep(0, Q)
# para.J = rep(0, Q)
# L = rep(0,Q)
for(q in 1:Q){
  w.para = results_1[[q]]$w.para
  w.np = results_1[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.99)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  # L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  
  beta.hatw.para[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu[[q]] = efnX%*%fitu$coefficients[-1]
  
  # np.theta[q] = results_3[[q]]$w.theta
  # para.J[q] = results_3[[q]]$w.j
  
}

beta.hatw.para.m = NULL
for(q in 1:Q){
  beta.hatw.para.m = cbind(beta.hatw.para.m, beta.hatw.para[[q]])
}

beta.hatw.np.m = NULL
for(q in 1:Q){
  beta.hatw.np.m = cbind(beta.hatw.np.m, beta.hatw.np[[q]])
}

beta.hatu.m = NULL
for(q in 1:Q){
  beta.hatu.m = cbind(beta.hatu.m, beta.hatu[[q]])
}

beta.hatw.para.mean = apply(beta.hatw.para.m, 1, mean)
beta.hatw.np.mean = apply(beta.hatw.np.m, 1, mean)
beta.hatu.mean = apply(beta.hatu.m, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para = (beta.hatw.para.mean - beta)
bias.np = (beta.hatw.np.mean - beta)
bias.u = (beta.hatu.mean - beta)

ISE.para = apply(beta.hatw.para.m, 2, function(X) mean((X-beta)^2))
ISE.np = apply(beta.hatw.np.m, 2, function(X) mean((X-beta)^2))
ISE.u = apply(beta.hatu.m, 2, function(X) mean((X-beta)^2))

table2 = cbind(rbind(median(ISE.u),median(ISE.para), median(ISE.np)),
               rbind(mean(ISE.u), mean(ISE.para), mean(ISE.np)),
               rbind(mean(bias.u^2), mean(bias.para^2), mean(bias.np^2)))

rownames(table2) = c("Unweighted", "para", "Np")
colnames(table2) = c("Median", "IMSE", "Ibias^2")


# boxplot.matrix(cbind(ISE.u, ISE.para,ISE.np),
#                names = c("Unweighted", "Para", "Np"),
#                ylab = "ISE", cex.main=2, cex.lab=1.5, cex.axis=1.5,
#                main = expression('PVE'[L]*' = 0.95,'*' PVE'[L^"*"]*' = 0.99'))

beta.hatw.para = replicate(Q, list())
beta.hatw.np = replicate(Q, list())
beta.hatu = replicate(Q, list())
# np.theta = rep(0, Q)
# para.J = rep(0, Q)
# L = rep(0,Q)
for(q in 1:Q){
  w.para = results_2[[q]]$w.para
  w.np = results_2[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  # L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  
  beta.hatw.para[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu[[q]] = efnX%*%fitu$coefficients[-1]
  
  # np.theta[q] = results_3[[q]]$w.theta
  # para.J[q] = results_3[[q]]$w.j
  
}

beta.hatw.para.m = NULL
for(q in 1:Q){
  beta.hatw.para.m = cbind(beta.hatw.para.m, beta.hatw.para[[q]])
}

beta.hatw.np.m = NULL
for(q in 1:Q){
  beta.hatw.np.m = cbind(beta.hatw.np.m, beta.hatw.np[[q]])
}

beta.hatu.m = NULL
for(q in 1:Q){
  beta.hatu.m = cbind(beta.hatu.m, beta.hatu[[q]])
}

beta.hatw.para.mean = apply(beta.hatw.para.m, 1, mean)
beta.hatw.np.mean = apply(beta.hatw.np.m, 1, mean)
beta.hatu.mean = apply(beta.hatu.m, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para = (beta.hatw.para.mean - beta)
bias.np = (beta.hatw.np.mean - beta)
bias.u = (beta.hatu.mean - beta)

ISE.para = apply(beta.hatw.para.m, 2, function(X) mean((X-beta)^2))
ISE.np = apply(beta.hatw.np.m, 2, function(X) mean((X-beta)^2))
ISE.u = apply(beta.hatu.m, 2, function(X) mean((X-beta)^2))

table3 = cbind(rbind(median(ISE.u),median(ISE.para), median(ISE.np)),
               rbind(mean(ISE.u), mean(ISE.para), mean(ISE.np)),
               rbind(mean(bias.u^2), mean(bias.para^2), mean(bias.np^2)))

rownames(table3) = c("Unweighted", "para", "Np")
colnames(table3) = c("Median", "IMSE", "Ibias^2")

# boxplot.matrix(cbind(ISE.u, ISE.para,ISE.np),
#                names = c("Unweighted", "Para", "Np"),
#                ylab = "ISE", cex.main=2, cex.lab=1.5, cex.axis=1.5,
#                main = expression('PVE'[L]*' = 0.99,'*' PVE'[L^"*"]*' = 0.95'))

beta.hatw.para = replicate(Q, list())
beta.hatw.np = replicate(Q, list())
beta.hatu = replicate(Q, list())
# np.theta = rep(0, Q)
# para.J = rep(0, Q)
# L = rep(0,Q)
for(q in 1:Q){
  w.para = results_2[[q]]$w.para
  w.np = results_2[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.99)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  # L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  
  beta.hatw.para[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu[[q]] = efnX%*%fitu$coefficients[-1]
  
  # np.theta[q] = results_3[[q]]$w.theta
  # para.J[q] = results_3[[q]]$w.j
  
}

beta.hatw.para.m = NULL
for(q in 1:Q){
  beta.hatw.para.m = cbind(beta.hatw.para.m, beta.hatw.para[[q]])
}

beta.hatw.np.m = NULL
for(q in 1:Q){
  beta.hatw.np.m = cbind(beta.hatw.np.m, beta.hatw.np[[q]])
}

beta.hatu.m = NULL
for(q in 1:Q){
  beta.hatu.m = cbind(beta.hatu.m, beta.hatu[[q]])
}

beta.hatw.para.mean = apply(beta.hatw.para.m, 1, mean)
beta.hatw.np.mean = apply(beta.hatw.np.m, 1, mean)
beta.hatu.mean = apply(beta.hatu.m, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para = (beta.hatw.para.mean - beta)
bias.np = (beta.hatw.np.mean - beta)
bias.u = (beta.hatu.mean - beta)

ISE.para = apply(beta.hatw.para.m, 2, function(X) mean((X-beta)^2))
ISE.np = apply(beta.hatw.np.m, 2, function(X) mean((X-beta)^2))
ISE.u = apply(beta.hatu.m, 2, function(X) mean((X-beta)^2))

table4 = cbind(rbind(median(ISE.u),median(ISE.para), median(ISE.np)),
               rbind(mean(ISE.u), mean(ISE.para), mean(ISE.np)),
               rbind(mean(bias.u^2), mean(bias.para^2), mean(bias.np^2)))

rownames(table4) = c("Unweighted", "para", "Np")
colnames(table4) = c("Median", "IMSE", "Ibias^2")


# boxplot.matrix(cbind(ISE.u, ISE.para,ISE.np),
#                names = c("Unweighted", "Para", "Np"),
#                ylab = "ISE", cex.main=2, cex.lab=1.5, cex.axis=1.5,
#                main = expression('PVE'[L]*' = 0.99,'*' PVE'[L^"*"]*' = 0.99'))

TABLE1 = rbind(cbind(table1, table2), cbind(table3, table4))
# xtable(TABLE4, digits = 4)

fstat.para = matrix(0, nrow = Q, ncol = 6)
fstat.np = matrix(0, nrow = Q, ncol = 6)
fstat.u = matrix(0, nrow = Q, ncol = 6)
L = rep(0, Q)
for(q in 1:Q){
  w.para = results_2[[q]]$w.para
  w.np = results_2[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.99)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  bal.para = summary(lm(treat~conf, weights = w.para))
  bal.np = summary(lm(treat~conf, weights = w.np))
  bal.u = summary(lm(treat~conf))
  
  if(L[q] == 6){
    fstat.para[q,] = c(bal.para[[1]]$fstatistic[1], bal.para[[2]]$fstatistic[1], 
                       bal.para[[3]]$fstatistic[1], bal.para[[4]]$fstatistic[1],
                       bal.para[[5]]$fstatistic[1], bal.para[[6]]$fstatistic[1])
    fstat.np[q,] = c(bal.np[[1]]$fstatistic[1], bal.np[[2]]$fstatistic[1], 
                     bal.np[[3]]$fstatistic[1], bal.np[[4]]$fstatistic[1],
                     bal.np[[5]]$fstatistic[1], bal.np[[6]]$fstatistic[1])
    fstat.u[q,] = c(bal.u[[1]]$fstatistic[1], bal.u[[2]]$fstatistic[1], 
                    bal.u[[3]]$fstatistic[1], bal.u[[4]]$fstatistic[1],
                    bal.u[[5]]$fstatistic[1], bal.u[[6]]$fstatistic[1])
  }else{
    fstat.para[q,] = c(bal.para[[1]]$fstatistic[1], bal.para[[2]]$fstatistic[1], 
                       bal.para[[3]]$fstatistic[1], bal.para[[4]]$fstatistic[1],
                       bal.para[[5]]$fstatistic[1], NA)
    fstat.np[q,] = c(bal.np[[1]]$fstatistic[1], bal.np[[2]]$fstatistic[1], 
                     bal.np[[3]]$fstatistic[1], bal.np[[4]]$fstatistic[1],
                     bal.np[[5]]$fstatistic[1], NA)
    fstat.u[q,] = c(bal.u[[1]]$fstatistic[1], bal.u[[2]]$fstatistic[1], 
                    bal.u[[3]]$fstatistic[1], bal.u[[4]]$fstatistic[1],
                    bal.u[[5]]$fstatistic[1], NA)
  }
}


par(mfrow = c(3,2))
par(yaxt="n")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
nameset = c("Np", "Para","Unweighted")
boxplot(cbind( fstat.np[,1],fstat.para[,1], fstat.u[,1]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "First FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,2], fstat.para[,2],fstat.u[,2]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "Second FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,3], fstat.para[,3], fstat.u[,3]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "Third FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,4], fstat.para[,4], fstat.u[,4]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "Fourth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,5], fstat.para[,5], fstat.u[,5]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "Fifth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,6], fstat.para[,6], fstat.u[,6]), horizontal = TRUE,
        cex.main = 2, cex.lab = 1.5,
        main = "Sixth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

# table(L)



fstat.para = matrix(0, nrow = Q, ncol = 4)
fstat.np = matrix(0, nrow = Q, ncol = 4)
fstat.u = matrix(0, nrow = Q, ncol = 4)
L = rep(0, Q)
for(q in 1:Q){
  w.para = results_1[[q]]$w.para
  w.np = results_1[[q]]$w.np
  
  Y = simu_data_1[[q]]$Y
  X = simu_data_1[[q]]$X
  conf = simu_data_1[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  bal.para = summary(lm(treat~conf, weights = w.para))
  bal.np = summary(lm(treat~conf, weights = w.np))
  bal.u = summary(lm(treat~conf))
  
  fstat.para[q,] = c(bal.para[[1]]$fstatistic[1], bal.para[[2]]$fstatistic[1], 
                     bal.para[[3]]$fstatistic[1], bal.para[[4]]$fstatistic[1])
  fstat.np[q,] = c(bal.np[[1]]$fstatistic[1], bal.np[[2]]$fstatistic[1], 
                   bal.np[[3]]$fstatistic[1], bal.np[[4]]$fstatistic[1])
  fstat.u[q,] = c(bal.u[[1]]$fstatistic[1], bal.u[[2]]$fstatistic[1], 
                  bal.u[[3]]$fstatistic[1], bal.u[[4]]$fstatistic[1])
}


par(mfrow = c(2,2))
par(yaxt="n")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0)) 
boxplot(cbind(fstat.np[,1], fstat.para[,1], fstat.u[,1]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "First FPC", xlab = "F-statisitc")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)


boxplot(cbind(fstat.np[,2], fstat.para[,2], fstat.u[,2]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Second FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,3], fstat.para[,3], fstat.u[,3]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Third FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,4], fstat.para[,4], fstat.u[,4]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Fourth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)


#--------------------------------------------------------------------------------------  
# Repeat the code in Simu Results Analysis part will give the rest tables and figures
#-------------------------------------------------------------------------------------- 
