# Simulation nxmxk n- genes, m- sample types, k replicate
setwd("D:\\Rochester books\\Research\\Sim_study\\")
library(xtable)
library(HTqPCR)
#############################################
########### Mean Imputation #################
#############################################

## dim of the data n, m, k;
n = 16
m = 6  
nsim <- 100

# LOAD FILES ONE BY ONE
load("res_16_6_6_100logit.RData"); k=6
# Mean Imputation to calculate sigmas and thetas from simulation results
#results for the mean of the original data
res.thetas <- array(NA, dim=c(dim(sim.output$res.table[[1]]$out.theta), 
                              length(sim.output$res.table)) )
dat.na.all<- array(NA, dim=c(dim(sim.output$res.table[[1]]$out.theta)[1], m*k, 
                             length(sim.output$res.table)) )
dimnames(dat.na.all)[c(1,2)] <- dimnames(res.thetas)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)
for (q in 1:length(sim.output$res.table)) { # by number of data sets
  thets <-NULL
  dat.na <- sim.output$res.table[[q]]$data.orig$data*(1-sim.output$res.table[[q]]$data.orig$nd)
  dat.na[which(dat.na<0.001)] <- NA
  dat.na.all[,,q] <- dat.na
  for (j in 1:dim(sim.output$res.table[[1]]$out.theta)[1]) { # be genes
    x <-NULL
    for (i in seq(1, m*k, by=k)) # by sample-type
    {x <- c(x, mean(dat.na[j,i:(i+k-1)], na.rm=TRUE))}
    thets <- rbind(thets, x) 
  }
  dat.na <- NULL
  res.thetas[,,q]  <- thets
}
# mean for genes with missing all sample-type was replaces by 40
res.thetas[which(is.na(res.thetas))]<- 40
# impute the means for NA in data set
# matchthe dim of mean.impute to dat.na (n gemes) by (m*k) by (nsim)
mean.impute<- array(NA, dim=c(dim(sim.output$res.table[[1]]$out.theta)[1], m*k, 
                              length(sim.output$res.table)) )
dimnames(mean.impute)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)
for (q in 1:length(sim.output$res.table)) { # by number of data sets
  x <- NULL
  for (i in 1:m) # by sample-type
  {
    repls <- matrix((rep(res.thetas[,i,q],k)), nrow = n, ncol=k)
    x <- cbind(x, repls)
  }
  mean.impute[,,q] <- x
}

# after mean.impute and dat.na.all are of the same dimentiones we can impute
na.ind <- which(is.na(dat.na.all), arr.ind = TRUE)
dat.na.all[na.ind] <- mean.impute[na.ind]

#results thetas
res.thetas <- array(NA, dim=c(n, m, 
                              length(sim.output$res.table)) )
dimnames(res.thetas)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)
for (q in 1:length(sim.output$res.table)) { # by number of data sets
  thets <-NULL
  for (j in 1:n) { # be genes
    x <-NULL
    for (i in seq(1, m*k, by=k)) # by sample-type
    {x <- c(x, mean(dat.na.all[j,i:(i+k-1),q], na.rm=TRUE))}
    thets <- rbind(thets, x) 
  }
  dat.na <- NULL
  res.thetas[,,q]  <- thets
}

# calculate mean shift to get variance properly
calc.mean.shift <- array(NA,dim=dim(dat.na.all))
dimnames(calc.mean.shift)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)
for (j in 1:nsim) { # by number of sim
  p=1
  shift <- mean.shift <- NULL
  for (i in 1:m) { # by sample type
    shift <- dat.na.all[,p:(p+k-1),j]-res.thetas[,i,j]
    mean.shift <- cbind(mean.shift, shift)
    p=p+k
  }
  calc.mean.shift[,,j] <- mean.shift
}
#results sigma2
res.sigma2 <- matrix(NA, ncol=n, 
                     nrow=length(sim.output$res.table))
colnames(res.sigma2) <- row.names(sim.output$true.table$true.thetas)
for (q in 1:length(sim.output$res.table)) {
  res.sigma2[q,]  <- apply(calc.mean.shift[,,q], 1, var)
}

###
# a function to get and format the results SI
out.table.SI <- function(n, m, k, res.thetas, res.sigma2, ind.big = NULL){
  # Bias
  Bias.sigma2 <- 1/k*(t(res.sigma2)) - 1/k*(tsigma2)
  bsigma2 <- apply(Bias.sigma2, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  names(Bias.sigma2) <- gsub("^.*?_","sigma 2 ", names(Bias.sigma2))
  bsigma2.mean <- apply(Bias.sigma2, 1, mean)
  quantile(bsigma2.mean, probs= c(0.25,0.5,0.75))
  
  # Thetas
  Bias.thetas <- array(NA, dim=dim(res.thetas) )
  dimnames(Bias.thetas)[c(1,2)] <- dimnames(res.thetas)[c(1,2)]
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-ind.big]}
  p=1
  for (q in (seq.loop)) {
    Bias.thetas[,,p]  <- res.thetas[,,q] - tthetas
    p=p+1
  }
  
  bthetas <- apply(Bias.thetas, c(1,2), mean)
  quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75))
  Bias.all=rbind(theta=quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75)), 
                 sigma2=quantile(bsigma2.mean, probs= c(0.25,0.5,0.75)))
  
  
  # Mean Square Bias
  MSE.sigma2 <- (bsigma2.mean)^2 + apply(Bias.sigma2, 1, sd)^2
  
  MSE.thetas <- (bthetas)^2 + apply(Bias.thetas, c(1,2), sd)^2
  
  MSE.bias.all=rbind(theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), 
                     sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  
  # combine all
  out <-cbind(Bias.all, MSE.bias.all)
  return(out)
}

tbeta   <- sim.output$true.table$true.beta
tsigma2 <- sim.output$true.table$true.sigma2
tthetas <- sim.output$true.table$true.thetas

resultsMean <- out.table.SI(n, m, k, res.thetas,res.sigma2)

#############################################
########### Single Imputation ###############
#############################################

load("res_SI_16_6_6_100logit.RData"); k=6

# SI to calculate sigmas and thetas from simmulation results
#results thetas
res.thetas<- array(NA, dim=c(dim(exprs(sim.output$res.table[[1]]$object))[1], m, 
                             length(sim.output$res.table)) )
dimnames(res.thetas)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)
for (q in 1:length(sim.output$res.table)) {
  thets <-NULL
  for (j in 1:dim(exprs(sim.output$res.table[[1]]$object))[1]) {
    x <-NULL
    for (i in seq(1, dim(exprs(sim.output$res.table[[q]]$object))[2], by=k)) 
    {x <- c(x, mean(exprs(sim.output$res.table[[q]]$object)[j,i:(i+k-1)]))}
    thets <- rbind(thets, x) 
  }
  res.thetas[,,q]  <- thets
}

calc.mean.shift <- array(NA,dim=c(dim(exprs(sim.output$res.table[[1]]$object)), nsim))
dimnames(calc.mean.shift)[c(1,2)] <- dimnames(exprs(sim.output$res.table[[1]]$object))
for (j in 1:nsim) {
  p=1
  shift <- mean.shift <- NULL
  for (i in 1:m) {
    shift <- exprs(sim.output$res.table[[j]]$object)[,p:(p+k-1)]-res.thetas[,i,j]
    mean.shift <- cbind(mean.shift, shift)
    p=p+k
  }
  calc.mean.shift[,,j] <- mean.shift
}
#results sigma2
res.sigma2 <- matrix(NA, ncol=dim(exprs(sim.output$res.table[[1]]$object))[1], 
                     nrow=length(sim.output$res.table))
colnames(res.sigma2) <- row.names(exprs(sim.output$res.table[[1]]$object))
for (q in 1:length(sim.output$res.table)) {
  res.sigma2[q,]  <- apply(calc.mean.shift[,,q],1,function(x) (m*k-1)/(m*k-m)*var(x))
}
###

tbeta   <- sim.output$true.table$true.beta
tsigma2 <- sim.output$true.table$true.sigma2
tthetas <- sim.output$true.table$true.thetas

# a function to get and format the results SI and mean imputation
out.table.SI.SEtheta <- function(n, m, k, res.thetas, res.sigma2, ind.big = NULL){
  # Bias
  # Sigma 2 m*k/(m*k-m)*
  Bias.sigma2 <- 1/k*(t(res.sigma2)) - 1/k*(tsigma2)
  bsigma2 <- apply(Bias.sigma2, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  names(Bias.sigma2) <- gsub("^.*?_","sigma 2 ", names(Bias.sigma2))
  bsigma2.mean <- apply(Bias.sigma2, 1, mean)
  quantile(bsigma2.mean, probs= c(0.25,0.5,0.75))
  
  # Thetas
  Bias.thetas <- array(NA, dim=dim(res.thetas) )
  dimnames(Bias.thetas)[c(1,2)] <- dimnames(res.thetas)[c(1,2)]
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-ind.big]}
  p=1
  for (q in (seq.loop)) {
    Bias.thetas[,,p]  <- res.thetas[,,q] - tthetas
    p=p+1
  }
  
  bthetas <- apply(Bias.thetas, c(1,2), mean)
  quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75))
  Bias.all=rbind(theta=quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75)), 
                 sigma2=quantile(bsigma2.mean, probs= c(0.25,0.5,0.75)))
  
  
  # Mean Square Bias
  MSE.sigma2 <- (bsigma2.mean)^2 + apply(Bias.sigma2, 1, sd)^2
  
  MSE.thetas <- (bthetas)^2 + apply(Bias.thetas, c(1,2), sd)^2
  
  #MSE.bias.all=c(MSE.beta, theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  MSE.bias.all=rbind(theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), 
                     sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  
  # combine all
  out <-cbind(Bias.all, MSE.bias.all)
  return(out)
}

resultsSI <- out.table.SI.SEtheta(n, m, k, res.thetas,res.sigma2)

#############################################
########### Multiple Imputation #############
#############################################

load("res_MI_16_6_6_100logit.RData"); k=6

### load true values
tbeta   <- sim.output$true.table$true.beta
tsigma2 <- sim.output$true.table$true.sigma2
tthetas <- sim.output$true.table$true.thetas

#results thetas
res.thetas<- array(NA, dim=c(dim(exprs(sim.output$res.table[[1]]$object[[1]]))[1], m, 
                             length(sim.output$res.table), 5))
dimnames(res.thetas)[c(1,2)] <- dimnames(sim.output$true.table$true.thetas)

for (z in 1:5) {
  for (q in 1:length(sim.output$res.table)) {
    thets <-NULL
    for (j in 1:dim(exprs(sim.output$res.table[[1]]$object[[z]]))[1]) {
      x <-NULL
      for (i in seq(1, dim(exprs(sim.output$res.table[[q]]$object[[z]]))[2], by=k)) 
      {x <- c(x, mean(exprs(sim.output$res.table[[q]]$object[[z]])[j,i:(i+k-1)]))}
      thets <- rbind(thets, x) 
    }
    res.thetas[,,q,z]  <- thets
  }
}

# MI to calculate sigmas and thetas from simmulation results

# combine MI results for thetas over 5 MI data sets by taking the average
mean.Q <- res.thetas.comb <- apply(res.thetas, c(1,2,3), mean) # dim = n m nsim

# calculate teh shift from mean.Q for each Qm
calc.Q.shift <- array(NA,dim=c(dim(res.thetas)))
dimnames(calc.Q.shift) <- dimnames(res.thetas)
for (i in 1:m) {
  for (j in 1:nsim) {
    #    p=1
    shift <- mean.shift <- NULL
    for (z in 1:5) {
      shift <- res.thetas[,i,j,z]-mean.Q[,i,j]
      mean.shift <- cbind(mean.shift, shift)
      #      p=p+k
    }
    calc.Q.shift[,i,j,] <- mean.shift
  }
}

# calculating vatiance of thetas combined from MI
Um.sigma <- t(matrix(rep(sim.output$res.table[[1]]$object[[6]]$sigma,each=6),ncol=16))
Um.stdev <- diag(sqrt(sim.output$res.table[[1]]$object[[6]]$cov.matrix))
# SE(Q.hat)
Um       <- Um.sigma*Um.stdev
# within imputation variability U
U.bar    <- (Um^2)[,1]

# calculating variance between MI, B
# be is of dimention 16x6x100 and 
# combining variance from MI, T
# in our case M=5
B <- T.var <- array(NA, dim=c(dim(res.thetas)[c(1,2)], 
                              length(sim.output$res.table)))
dimnames(B)[c(1)] <- dimnames(exprs(sim.output$res.table[[1]]$object[[z]]))[1]
for (i in 1:m) {
  for (q in 1:length(sim.output$res.table)) {
    # Sigma^2 (m*k-1)/(m*k-m)* constant is used
    # var() function gives s^2/(m*k-1), in this situation we don't need s^2/(m*k-m)
    B[,i,q]     <- apply(calc.Q.shift[,i,q,], 1, function(x) var(x))
    T.var[,i,q] <- U.bar + (1+1/5)*B[,i,q]
  }
}

# a function to get and format the results MI

out.table.MI.comb <- function(n, m, k, mean.Q, T.var, ind.big = NULL){
  # Bias
  # Sigmas
  Bias.sigma2 <- array(NA, dim=dim(T.var) )
  dimnames(Bias.sigma2)[c(1,2)] <- dimnames(mean.Q)[c(1,2)]
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-ind.big]}
  for (q in (seq.loop)) {
    Bias.sigma2[,,q]  <- T.var[,,q] - t(matrix(rep(tsigma2/k, each=m), nrow = m))
  }
  
  bsigma2 <- apply(Bias.sigma2, c(1,2), mean)
  SE.mean.sigma2 <- apply(Bias.sigma2, c(1,2), function(x) sd(x)/sqrt(length(x)))
  quantile(as.vector(bsigma2), prob=c(0.25, 0.5, 0.75))
  
  # Thetas
  Bias.thetas <- array(NA, dim=dim(mean.Q) )
  dimnames(Bias.thetas)[c(1,2)] <- dimnames(mean.Q)[c(1,2)]
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-ind.big]}
  for (q in (seq.loop)) {
    Bias.thetas[,,q]  <- mean.Q[,,q] - tthetas
  }
  
  bthetas <- apply(Bias.thetas, c(1,2), mean)
  SE.mean.thetas <- apply(Bias.thetas, c(1,2), function(x) sd(x)/sqrt(length(x)))
  quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75))
  Bias.all=rbind(theta=quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75)), 
                 sigma2=quantile(as.vector(bsigma2), probs= c(0.25,0.5,0.75)))
  
  # Mean Square Bias
  MSE.sigma2 <- (bsigma2)^2 + apply(Bias.sigma2, c(1,2), sd)^2
  
  MSE.thetas <- (bthetas)^2 + apply(Bias.thetas, c(1,2), sd)^2
  
  #MSE.bias.all=c(MSE.beta, theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  MSE.bias.all=rbind(theta =quantile(as.vector(MSE.thetas), probs=c(0.25, 0.5, 0.75)), 
                     sigma2=quantile(as.vector(MSE.sigma2), probs= c(0.25,0.5,0.75)))
  
  # combine all
  out <-cbind(Bias.all, MSE.bias.all)
  return(out)
}

# replaces res.thetas by mean.Q and res.sigma2 by T.var
resultsMI <- out.table.MI.comb(n, m, k, mean.Q, T.var)

#############################################
########### Direct Estimation ###############
#############################################

load("res_16_6_6_100logit.RData"); k=6

tbeta   <- sim.output$true.table$true.beta
tsigma2 <- sim.output$true.table$true.sigma2
tthetas <- sim.output$true.table$true.thetas

#results prop of non-detects
res.nd<- matrix(NA, ncol=1, 
                nrow=length(sim.output$res.table) )
colnames(res.nd) <- names(sim.output$res.table[[1]]$out.beta[c(1)])
for (q in 1:length(sim.output$res.table)) {
  res.nd[q,]  <- sim.output$res.table[[q]]$out.beta[c(1)]
}
summary(res.nd)

#results of beta
res.beta<- matrix(NA, ncol=length(sim.output$true.table$true.beta), 
                  nrow=length(sim.output$res.table) )
colnames(res.beta) <- names(sim.output$res.table[[1]]$out.beta[c(2,4)])
for (q in 1:length(sim.output$res.table)) {
  res.beta[q,]  <- sim.output$res.table[[q]]$out.beta[c(2,4)]
}
ind.big <- which(res.beta[,1]<3*tbeta[1])
ind.big
res.beta.update <- res.beta[-ind.big,]
#results sigma2
res.sigma2 <- matrix(NA, ncol=length(sim.output$res.table[[1]]$out.sigma2), 
                     nrow=length(sim.output$res.table) )
colnames(res.sigma2) <- names(sim.output$res.table[[1]]$out.sigma2)
for (q in 1:length(sim.output$res.table)) {
  res.sigma2[q,]  <- sim.output$res.table[[q]]$out.sigma2
}
res.sigma2.update <- res.sigma2[-ind.big,]
#results thetas
res.thetas<- array(NA, dim=c(dim(sim.output$res.table[[1]]$out.theta), 
                             length(sim.output$res.table)) )
dimnames(res.thetas)[c(1,2)] <- dimnames(sim.output$res.table[[1]]$out.theta)
for (q in 1:length(sim.output$res.table)) {
  res.thetas[,,q]  <- sim.output$res.table[[q]]$out.theta
}
res.thetas.update <- res.thetas[,,-ind.big]

out.table.DirEst.SE <- function(n, m, k, res.beta, res.thetas, res.sigma2, ind.big = NULL){
  # Bias
  Bias.beta <- t(res.beta) - tbeta
  bbeta <- apply(Bias.beta, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  bbeta.mean <- apply(Bias.beta, 1, mean)

  # Sigma 2
  Bias.sigma2 <- 1/k*t(m*k/(m*k-m)*res.sigma2) - 1/k*tsigma2
  bsigma2 <- apply(Bias.sigma2, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  names(Bias.sigma2) <- gsub("^.*?_","sigma 2 ", names(Bias.sigma2))
  bsigma2.mean <- apply(Bias.sigma2, 1, mean)
  quantile(bsigma2.mean, probs= c(0.25,0.5,0.75))
  
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-ind.big]}
  
  # Thetas
  Bias.thetas <- array(NA, dim=dim(res.thetas[,,seq.loop]) )
  dimnames(Bias.thetas)[c(1,2)] <- dimnames(sim.output$res.table[[1]]$out.theta)
  p=1
  for (q in (seq.loop)) {
    Bias.thetas[,,p]  <- sim.output$res.table[[q]]$out.theta - tthetas
    p=p+1
  }
  
  bthetas <- apply(Bias.thetas, c(1,2), mean)
  quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75))
  Bias.all=rbind(theta=quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75)), 
                 sigma2=quantile(bsigma2.mean, probs= c(0.25,0.5,0.75)))
  
  # Mean Square Bias
  MSE.beta <- (bbeta.mean)^2 + apply(Bias.beta, 1, sd)^2
  
  MSE.sigma2 <- (bsigma2.mean)^2 + apply(Bias.sigma2, 1, sd)^2
  
  MSE.thetas <- (bthetas)^2 + apply(Bias.thetas, c(1,2), sd)^2
  
  MSE.bias.all=rbind(theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), 
                     sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  
  # combine all
  out <-cbind(Bias.all, MSE.bias.all)
  return(out)
}

nsim=100

resultsDirEst <- out.table.DirEst.SE(n,m,k,res.beta, res.thetas, res.sigma2)
round(resultsMean, 3)
round(resultsSI, 3)
round(resultsMI, 3)
round(resultsDirEst, 3)
