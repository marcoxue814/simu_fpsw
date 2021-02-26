#-----------------------------------
# Some functions needed
#-----------------------------------

# to get fpc scores
getFPC = function(X){
  N = nrow(X)
  A = matrix(NA, nrow = N, ncol = 8) #generate 6 FPC scores
  A[,1] = sqrt(1 ^ (-1.15))*X[,1]
  A[,2] = sqrt(2 ^ (-1.15))*X[,2]
  A[,3] = sqrt(3 ^ (-1.15))*X[,3]
  A[,4] = sqrt(4 ^ (-1.15))*X[,4]
  A[,5] = sqrt(5 ^ (-1.15))*X[,5]
  A[,6] = sqrt(6 ^ (-1.15))*X[,6]
  A[,7] = sqrt(7 ^ (-1.15))*X[,7]
  A[,8] = sqrt(8 ^ (-1.15))*X[,8]
  R = list()
  R$A = A
  return(R)
}

# to get functional variable
getX = function(FPC, tgrid){
  X = outer(FPC[,1], sqrt(2)*sin(2*pi*tgrid)) + outer(FPC[,2], sqrt(2)*cos(2*pi*tgrid)) +
    outer(FPC[,3], sqrt(2)*sin(4*pi*tgrid)) + outer(FPC[,4], sqrt(2)*cos(4*pi*tgrid)) +
    outer(FPC[,5], sqrt(2)*sin(6*pi*tgrid)) + outer(FPC[,6], sqrt(2)*cos(6*pi*tgrid)) + 
    outer(FPC[,7], sqrt(2)*sin(8*pi*tgrid)) + outer(FPC[,8], sqrt(2)*cos(8*pi*tgrid))
  R = list()
  R$X = X
  return(R)
}


# to get covariates, linear in fpc scores
getC_ln = function(X){
  N = nrow(X)
  C = matrix(NA, nrow = N, ncol = 3) #generate 3 confounders
  C[,1] = X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6] + X[,7] + X[,8] + rnorm(N, 0, 1)
  C[,2] = - X[,1] + X[,2] - X[,3] + X[,4] - X[,5] + X[,6] - X[,7] + X[,8] + rnorm(N, 0, 1)
  C[,3] = X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6] + X[,7] + rnorm(N, 0, 1)
  R = list()
  R$C = C
  return(R)
}

# to get outcome linear in fpc scores
getY_ln = function(FPC, conf){
  #FPC: the six FPC scores
  #C: confounders
  N = nrow(C)
  Y = 1 + 2*FPC[,1] + 1*FPC[,2] + 0.5*FPC[,3] + 0.5*FPC[,4] + 2*conf[,1] + conf[,2] + 0.5*conf[,3] +  3*rnorm(N)
  R = list()
  R$Y = Y
  return(R)
}


Q = 200
# results_1 = replicate(Q, list())
# simu_data_1 = replicate(Q, list())
# pb = txtProgressBar(min = 0, max = Q, style = 3)]

simu5_data = replicate(Q, list())
for (q in 1:Q){
  set.seed(123456+q)
  N = 200
  P = 8
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu5_data[[q]]$Y = Y
  simu5_data[[q]]$X = X
  simu5_data[[q]]$conf = conf
}

library(foreach)
library(parallel)

simu5_result = foreach(q = 1:Q) %dopar% {
  Y = simu5_data[[q]]$Y
  X = simu5_data[[q]]$X
  conf = simu5_data[[q]]$conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  list('w.para' = fcbps$weights, 'w.np' = npfcbps$w)
  
}


save(simu5_data, file = "simu5_data.rda")
save(simu5_result, file = "simu5_result.rda")


