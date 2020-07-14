
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

