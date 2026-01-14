# Simulation nxmxk n- genes, m- sample types, k replicate
setwd("D:\\Rochester books\\Research\\Sim_study\\")
library(HTqPCR)
#############################################
########### TABLE 3 supplement ##############
#############################################

#function that peforms all the calculations
get.results <- function(sim.output, n, m, k) {
  nsim=100
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
    Bias.sigma2 <- (t(res.sigma2)) - (tsigma2)
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
  
  results <- out.table.SI(n, m, k, res.thetas,res.sigma2)
  return(results)
}
n=16; m=6
k=4
load("res_16_6_4_100logit.RData"); 
t3logit16_6_4 <- get.results(sim.output, n, m, k)
k=6
load("res_16_6_6_100logit.RData"); 
t3logit16_6_6 <- get.results(sim.output, n, m, k)
k=10
load("res_16_6_10_100logit.RData"); 
t3logit16_6_10 <- get.results(sim.output, n, m, k)

#Table 3
round(t3logit16_6_4, 3)
round(t3logit16_6_6, 3)
round(t3logit16_6_10, 3)