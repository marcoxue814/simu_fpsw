#-----------------------------------
# Simulation Setting 5
#-----------------------------------


#-----------------------------------
# Some functions needed
#-----------------------------------

# to get fpc scores
getFPC = function(X){
  N = nrow(X)
  A = matrix(NA, nrow = N, ncol = 3) #generate 6 FPC scores
  A[,1] = sqrt(1 ^ (-3.5))*X[,1]
  A[,2] = sqrt(2 ^ (-3.5))*X[,2]
  A[,3] = sqrt(3 ^ (-3.5))*X[,3]
  R = list()
  R$A = A
  return(R)
}

# to get functional variable
getX = function(FPC, tgrid){
  X = outer(FPC[,1], sqrt(2)*sin(2*pi*tgrid)) + outer(FPC[,2], sqrt(2)*cos(2*pi*tgrid)) +
    outer(FPC[,3], sqrt(2)*sin(4*pi*tgrid))
  R = list()
  R$X = X
  return(R)
}


# to get covariates, linear in fpc scores
getC_ln = function(X){
  N = nrow(X)
  C = matrix(NA, nrow = N, ncol = 2) #generate 2 confounders
  C[,1] = X[,1]*ifelse(X[,1] > 0, 1, 0) #+ rnorm(N, 0, 1)
  C[,2] = X[,2]*ifelse(X[,2] < 0, 1, 0) #+ rnorm(N, 0, 1)
  # C[,3] = X[,3]* ifelse(X[,1] < 0, 1, 0) + rnorm(N, 0, 1)
  R = list()
  R$C = C
  return(R)
}

# to get outcome linear in fpc scores
getY_ln = function(FPC, conf){
  #FPC: the six FPC scores
  #C: confounders
  N = nrow(C)
  Y = 1 + 2*FPC[,1] + 1*FPC[,2] + 2*conf[,1] + conf[,2]  +  3*rnorm(N)
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
results_simu7 = replicate(Q, list())
simu7_data = replicate(Q, list())
pb = txtProgressBar(min = 0, max = Q, style = 3)
L = rep(0, Q)
for (q in 1:Q){
  print(paste0("Working on iter=",q))
  set.seed(123456+q)
  N = 200
  P = 3
  C = rmvnorm(N, mean = rep(0, P), sigma = diag(P))
  FPC =   getFPC(C)$A
  conf = getC_ln(C)$C
  
  tgrid = 0:50/50
  X = getX(FPC, tgrid = tgrid)$X
  Y = getY_ln(FPC, conf)$Y
  simu7_data[[q]]$Y = Y
  simu7_data[[q]]$X = X
  simu7_data[[q]]$conf = conf
  
  pve = 0.95
  treat = FPCA(X, pve = pve)$scr
  L[q] = ncol(treat)
  
  fcbps = FCBPS(treat, conf)
  npfcbps = npFCBPS.Functional(treat, conf, rho = 0.1/N)
  
  results_simu7[[q]]$w.para = fcbps$weights
  results_simu7[[q]]$w.np = npfcbps$w
  results_simu7[[q]]$sumw.para = fcbps$totalweights
  results_simu7[[q]]$sumw.np = npfcbps$sumw
  
  setTxtProgressBar(pb, q)
}

table(L)
save(simu7_data, file = "simu7_data.rda")
save(results_simu7, file = "results_simu7.rda")



#----------------------------------
# Simu Results Analysis
#----------------------------------

Q = 200
beta.hatw.para_simu7 = replicate(Q, list())
beta.hatw.np_simu7 = replicate(Q, list())
beta.hatu_simu7 = replicate(Q, list())
beta.hatu.cov_simu7 = replicate(Q, list())
L = rep(0,Q)
for(q in 1:Q){
  w.para = results_simu7[[q]]$w.para
  w.np = results_simu7[[q]]$w.np
  
  Y = simu7_data[[q]]$Y
  X = simu7_data[[q]]$X
  conf = simu7_data[[q]]$conf
  
  fpcaX = FPCA(X, pve = 0.95)
  treat = fpcaX$scr
  efnX = fpcaX$efn
  L[q] = ncol(treat)
  
  fitw.para = lm(Y~treat, weights = w.para)
  fitw.np = lm(Y~treat, weights = w.np)
  fitu = lm(Y~treat)
  fitu.cov = lm(Y~treat + conf)
  
  beta.hatw.para_simu7[[q]] = efnX%*%fitw.para$coefficients[-1]
  beta.hatw.np_simu7[[q]] = efnX%*%fitw.np$coefficients[-1]
  beta.hatu_simu7[[q]] = efnX%*%fitu$coefficients[-1]
  beta.hatu.cov_simu7[[q]] = efnX%*%fitu.cov$coefficients[2:(L[q]+1)]
  
  
}

beta.hatw.para.m_simu7 = NULL
for(q in 1:Q){
  beta.hatw.para.m_simu7 = cbind(beta.hatw.para.m_simu7, beta.hatw.para_simu7[[q]])
}

beta.hatw.np.m_simu7 = NULL
for(q in 1:Q){
  beta.hatw.np.m_simu7 = cbind(beta.hatw.np.m_simu7, beta.hatw.np_simu7[[q]])
}

beta.hatu.m_simu7 = NULL
for(q in 1:Q){
  beta.hatu.m_simu7 = cbind(beta.hatu.m_simu7, beta.hatu_simu7[[q]])
}

beta.hatu.cov.m_simu7 = NULL
for(q in 1:Q){
  beta.hatu.cov.m_simu7 = cbind(beta.hatu.cov.m_simu7, beta.hatu.cov_simu7[[q]])
}

beta.hatw.para.mean_simu7 = apply(beta.hatw.para.m_simu7, 1, mean)
beta.hatw.np.mean_simu7 = apply(beta.hatw.np.m_simu7, 1, mean)
beta.hatu.mean_simu7 = apply(beta.hatu.m_simu7, 1, mean)
beta.hatu.cov.mean_simu7 = apply(beta.hatu.cov.m_simu7, 1, mean)

tgrid=0:50/50
ntgrid=length(tgrid)
beta= 2*sqrt(2)*sin(2*pi*tgrid) + sqrt(2)*cos(2*pi*tgrid)

bias.para_simu7 = (beta.hatw.para.mean_simu7 - beta)
bias.np_simu7 = (beta.hatw.np.mean_simu7 - beta)
bias.u_simu7 = (beta.hatu.mean_simu7 - beta)
bias.u.cov_simu7 = (beta.hatu.cov.mean_simu7 - beta)

ISE.para_simu7 = apply(beta.hatw.para.m_simu7, 2, function(X) mean((X-beta)^2))
ISE.np_simu7 = apply(beta.hatw.np.m_simu7, 2, function(X) mean((X-beta)^2))
ISE.u_simu7 = apply(beta.hatu.m_simu7, 2, function(X) mean((X-beta)^2))
ISE.u.cov_simu7 = apply(beta.hatu.cov.m_simu7, 2, function(X) mean((X-beta)^2))

table_simu7 = cbind(rbind(median(ISE.u.cov_simu7), median(ISE.u_simu7),median(ISE.para_simu7), median(ISE.np_simu7)), 
                    rbind(mean(ISE.u.cov_simu7), mean(ISE.u_simu7), mean(ISE.para_simu7), mean(ISE.np_simu7)),
                    rbind(mean(bias.u.cov_simu7^2),mean(bias.u_simu7^2), mean(bias.para_simu7^2), mean(bias.np_simu7^2)))
rownames(table_simu7) = c("dirAd", "Unweighted", "para", "Np")
colnames(table_simu7) = c("Median", "IMSE", "Ibias^2")

xtable(table_simu7, digits = 4)

#--------------------------
# balance performance
#--------------------------

fstat.para = matrix(NA, nrow = Q, ncol = 2)
fstat.np = matrix(NA, nrow = Q, ncol = 2)
fstat.u = matrix(NA, nrow = Q, ncol = 2)
L = rep(0, Q)
for(q in 1:Q){
  w.para = results_simu7[[q]]$w.para
  w.np = results_simu7[[q]]$w.np
  
  Y = simu7_data[[q]]$Y
  X = simu7_data[[q]]$X
  conf = simu7_data[[q]]$conf
  
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
        cex.main=2, cex.lab=1.5, ylim = c(0, 500),
        main = "First FPC", xlab = "F-statisitc")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)


boxplot(cbind(fstat.np[,2], fstat.para[,2], fstat.u[,2]), horizontal = TRUE,
        cex.main=2, cex.lab=1.5, ylim = c(0, 500),
        main = "Second FPC", xlab = "F-statistic")
axis(2, at=seq(1,3, by=1), labels = FALSE)
text(y = 1:length(nameset), par("usr")[1], labels = nameset, 
     srt = 0, pos = 2, xpd = TRUE, cex = 1.5)









