#-----------------------------------
# Simulation Setting 5
#-----------------------------------


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

Q = 200
results_simu5 = replicate(Q, list())
simu5_data = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for (q in 1:Q){
  print(paste0("Working on iter=",q))
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
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_simu5[[q]]$w.para = fcbps$weights
  results_simu5[[q]]$w.np = npfcbps$w
}



save(simu5_data, file = "simu5_data.rda")
save(results_simu5, file = "results_simu5.rda")


#----------------------------------
# Simu Results Analysis
#----------------------------------

Q = 200
beta.hatw.para_simu5 = replicate(Q, list())
beta.hatw.np_simu5 = replicate(Q, list())
beta.hatu_simu5 = replicate(Q, list())
beta.hatu.cov_simu5 = replicate(Q, list())
# np.theta = rep(0, Q)
# para.J = rep(0, Q)
L = rep(0,Q)
for(q in 1:Q){
  w.para = results_simu5[[q]]$w.para
  w.np = results_simu5[[q]]$w.np
  
  Y = simu5_data[[q]]$Y
  X = simu5_data[[q]]$X
  conf = simu5_data[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  
  beta.hatw.para_simu5[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np_simu5[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu_simu5[[q]] = efnX%*%fitu$coefficients[-1]
  
  beta.hatu.cov_simu5[[q]] = efnX%*%fitu.cov$coefficients[2:(L[q]+1)]
  
  # np.theta[q] = results_3[[q]]$w.theta
  # para.J[q] = results_3[[q]]$w.j
  
}

beta.hatw.para.m_simu5 = NULL
for(q in 1:Q){
  beta.hatw.para.m_simu5 = cbind(beta.hatw.para.m_simu5, beta.hatw.para_simu5[[q]])
}

beta.hatw.np.m_simu5 = NULL
for(q in 1:Q){
  beta.hatw.np.m_simu5 = cbind(beta.hatw.np.m_simu5, beta.hatw.np_simu5[[q]])
}

beta.hatu.m_simu5 = NULL
for(q in 1:Q){
  beta.hatu.m_simu5 = cbind(beta.hatu.m_simu5, beta.hatu_simu5[[q]])
}

beta.hatu.cov.m_simu5 = NULL
for(q in 1:Q){
  beta.hatu.cov.m_simu5 = cbind(beta.hatu.cov.m_simu5, beta.hatu.cov_simu5[[q]])
}

beta.hatw.para.mean_simu5= apply(beta.hatw.para.m_simu5, 1, mean)
beta.hatw.np.mean_simu5 = apply(beta.hatw.np.m_simu5, 1, mean)
beta.hatu.mean_simu5 = apply(beta.hatu.m_simu5, 1, mean)
beta.hatu.cov.mean_simu5 = apply(beta.hatu.cov.m_simu5, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para_simu5 = (beta.hatw.para.mean_simu5 - beta)
bias.np_simu5 = (beta.hatw.np.mean_simu5 - beta)
bias.u_simu5 = (beta.hatu.mean_simu5 - beta)
bias.u.cov_simu5 = (beta.hatu.cov.mean_simu5 - beta)

ISE.para_simu5 = apply(beta.hatw.para.m_simu5, 2, function(X) mean((X-beta)^2))
ISE.np_simu5 = apply(beta.hatw.np.m_simu5, 2, function(X) mean((X-beta)^2))
ISE.u_simu5 = apply(beta.hatu.m_simu5, 2, function(X) mean((X-beta)^2))
ISE.u.cov_simu5 = apply(beta.hatu.cov.m_simu5, 2, function(X) mean((X-beta)^2))


table_simu5 = cbind(rbind(median(ISE.u.cov_simu5), median(ISE.u_simu5),median(ISE.para_simu5), median(ISE.np_simu5)), 
                    rbind(mean(ISE.u.cov_simu5), mean(ISE.u_simu5), mean(ISE.para_simu5), mean(ISE.np_simu5)),
                    rbind(mean(bias.u.cov_simu5^2),mean(bias.u_simu5^2), mean(bias.para_simu5^2), mean(bias.np_simu5^2)))
rownames(table_simu5) = c("dirAd", "Unweighted", "para", "Np")
colnames(table_simu5) = c("Median", "IMSE", "Ibias^2")

xtable(table_simu5, digits = 4)


#--------------------------
# balance performance
#--------------------------

fstat.para = matrix(NA, nrow = Q, ncol = 7)
fstat.np = matrix(NA, nrow = Q, ncol = 7)
fstat.u = matrix(NA, nrow = Q, ncol = 7)
L = rep(0, Q)
for(q in 1:Q){
  w.para = results_simu5[[q]]$w.para
  w.np = results_simu5[[q]]$w.np
  
  Y = simu5_data[[q]]$Y
  X = simu5_data[[q]]$X
  conf = simu5_data[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  bal.para = summary(lm(treat~conf, weights = w.para))
  bal.np = summary(lm(treat~conf, weights = w.np))
  bal.u = summary(lm(treat~conf))
  
  
  fstat.para[q,] = c(bal.para[[1]]$fstatistic[1], bal.para[[2]]$fstatistic[1],
                     bal.para[[3]]$fstatistic[1], bal.para[[4]]$fstatistic[1],
                     bal.para[[5]]$fstatistic[1], bal.para[[6]]$fstatistic[1],
                     bal.para[[7]]$fstatistic[1])
  fstat.np[q,] = c(bal.np[[1]]$fstatistic[1], bal.np[[2]]$fstatistic[1],
                   bal.np[[3]]$fstatistic[1], bal.np[[4]]$fstatistic[1],
                   bal.np[[5]]$fstatistic[1], bal.np[[6]]$fstatistic[1],
                   bal.np[[7]]$fstatistic[1])
  fstat.u[q,] = c(bal.u[[1]]$fstatistic[1], bal.u[[2]]$fstatistic[1],
                  bal.u[[3]]$fstatistic[1], bal.u[[4]]$fstatistic[1],
                  bal.u[[5]]$fstatistic[1], bal.u[[6]]$fstatistic[1],
                  bal.u[[7]]$fstatistic[1])
  
  
  
  
}


par(mfrow = c(4,2))
par(yaxt="n")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0))
nameset = c("Np", "Para","Unweighted")
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

boxplot(cbind(fstat.np[,5], fstat.para[,5], fstat.u[,5]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Fifth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,6], fstat.para[,6], fstat.u[,6]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Sixth FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset,
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)

boxplot(cbind(fstat.np[,7], fstat.para[,7], fstat.u[,7]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5,
        main = "Seventh FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset,
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)


#-----------------------------------
# Simulation Setting 6
#-----------------------------------


#-----------------------------------
# Some functions needed
#-----------------------------------

# to get fpc scores
getFPC = function(X){
  N = nrow(X)
  A = matrix(NA, nrow = N, ncol = 8) #generate 6 FPC scores
  A[,1] = sqrt(1 ^ (-3.5))*X[,1]
  A[,2] = sqrt(2 ^ (-3.5))*X[,2]
  A[,3] = sqrt(3 ^ (-3.5))*X[,3]
  A[,4] = sqrt(4 ^ (-3.5))*X[,4]
  A[,5] = sqrt(5 ^ (-3.5))*X[,5]
  A[,6] = sqrt(6 ^ (-3.5))*X[,6]
  A[,7] = sqrt(7 ^ (-3.5))*X[,7]
  A[,8] = sqrt(8 ^ (-3.5))*X[,8]
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

Q = 200
results_simu6 = replicate(Q, list())
simu6_data = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
for (q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 8
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu6_data[[q]]$Y = Y
  simu6_data[[q]]$X = X
  simu6_data[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_simu6[[q]]$w.para = fcbps$weights
  results_simu6[[q]]$w.np = npfcbps$w
  
  setTxtProgressBar(pb, q)
}

save(simu6_data, file = "simu6_data.rda")
save(results_simu6, file = "results_simu6.rda")



#----------------------------------
# Simu Results Analysis
#----------------------------------

Q = 200
beta.hatw.para_simu6 = replicate(Q, list())
beta.hatw.np_simu6 = replicate(Q, list())
beta.hatu_simu6 = replicate(Q, list())
beta.hatu.cov_simu6 = replicate(Q, list())
L = rep(0,Q)
for(q in 1:Q){
  w.para = results_simu6[[q]]$w.para
  w.np = results_simu6[[q]]$w.np
  
  Y = simu6_data[[q]]$Y
  X = simu6_data[[q]]$X
  conf = simu6_data[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  fitu.cov = lm(Y~treat + conf)
  
  beta.hatw.para_simu6[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np_simu6[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu_simu6[[q]] = efnX%*%fitu$coefficients[-1]
  beta.hatu.cov_simu6[[q]] = efnX%*%fitu.cov$coefficients[2:(L[q]+1)]
  
  
}

beta.hatw.para.m_simu6 = NULL
for(q in 1:Q){
  beta.hatw.para.m_simu6 = cbind(beta.hatw.para.m_simu6, beta.hatw.para_simu6[[q]])
}

beta.hatw.np.m_simu6 = NULL
for(q in 1:Q){
  beta.hatw.np.m_simu6 = cbind(beta.hatw.np.m_simu6, beta.hatw.np_simu6[[q]])
}

beta.hatu.m_simu6 = NULL
for(q in 1:Q){
  beta.hatu.m_simu6 = cbind(beta.hatu.m_simu6, beta.hatu_simu6[[q]])
}

beta.hatu.cov.m_simu6 = NULL
for(q in 1:Q){
  beta.hatu.cov.m_simu6 = cbind(beta.hatu.cov.m_simu6, beta.hatu.cov_simu6[[q]])
}

beta.hatw.para.mean_simu6= apply(beta.hatw.para.m_simu6, 1, mean)
beta.hatw.np.mean_simu6 = apply(beta.hatw.np.m_simu6, 1, mean)
beta.hatu.mean_simu6 = apply(beta.hatu.m_simu6, 1, mean)
beta.hatu.cov.mean_simu6 = apply(beta.hatu.cov.m_simu6, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid) +
  0.5*sqrt(2)*sin(4*pi*tgrid) + 0.5*sqrt(2)*cos(4*pi*tgrid) 

bias.para_simu6 = (beta.hatw.para.mean_simu6 - beta)
bias.np_simu6 = (beta.hatw.np.mean_simu6 - beta)
bias.u_simu6 = (beta.hatu.mean_simu6 - beta)
bias.u.cov_simu6 = (beta.hatu.cov.mean_simu6 - beta)

ISE.para_simu6 = apply(beta.hatw.para.m_simu6, 2, function(X) mean((X-beta)^2))
ISE.np_simu6 = apply(beta.hatw.np.m_simu6, 2, function(X) mean((X-beta)^2))
ISE.u_simu6 = apply(beta.hatu.m_simu6, 2, function(X) mean((X-beta)^2))
ISE.u.cov_simu6 = apply(beta.hatu.cov.m_simu6, 2, function(X) mean((X-beta)^2))

table_simu6 = cbind(rbind(median(ISE.u.cov_simu6), median(ISE.u_simu6),median(ISE.para_simu6), median(ISE.np_simu6)), 
                    rbind(mean(ISE.u.cov_simu6), mean(ISE.u_simu6), mean(ISE.para_simu6), mean(ISE.np_simu6)),
                    rbind(mean(bias.u.cov_simu6^2),mean(bias.u_simu6^2), mean(bias.para_simu6^2), mean(bias.np_simu6^2)))
rownames(table_simu6) = c("dirAd", "Unweighted", "para", "Np")
colnames(table_simu6) = c("Median", "IMSE", "Ibias^2")

xtable(table_simu6, digits = 4)

#--------------------------
# balance performance
#--------------------------

fstat.para = matrix(NA, nrow = Q, ncol = 2)
fstat.np = matrix(NA, nrow = Q, ncol = 2)
fstat.u = matrix(NA, nrow = Q, ncol = 2)
L = rep(0, Q)
for(q in 1:Q){
  w.para = results_simu6[[q]]$w.para
  w.np = results_simu6[[q]]$w.np
  
  Y = simu6_data[[q]]$Y
  X = simu6_data[[q]]$X
  conf = simu6_data[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  bal.para = summary(lm(treat~conf, weights = w.para))
  bal.np = summary(lm(treat~conf, weights = w.np))
  bal.u = summary(lm(treat~conf))
  
  
  fstat.para[q,] = c(bal.para[[1]]$fstatistic[1], bal.para[[2]]$fstatistic[1])
  fstat.np[q,] = c(bal.np[[1]]$fstatistic[1], bal.np[[2]]$fstatistic[1])
  fstat.u[q,] = c(bal.u[[1]]$fstatistic[1], bal.u[[2]]$fstatistic[1])
  
  
  
  
}


par(mfrow = c(1,2))
par(yaxt="n")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 4, 0, 0))
nameset = c("Np", "Para","Unweighted")
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


