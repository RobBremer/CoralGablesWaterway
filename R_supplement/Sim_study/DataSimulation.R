setwd("D:\\Rochester books\\Research\\Sim_study\\")
sessionInfo()

library(HTqPCR)
library(nondetects)
library(limma)
library(mvtnorm)
library(coda)
library(R2OpenBUGS)
library(mcgibbsit)
library(base)
library(tmvtnorm)
library(truncnorm)
library(xtable)
library(arm) 
library(doParallel)
sessionInfo()

# make a new qPCRset
qPCRdata<-function(DataName, ControlGene=NULL, Ct_nd, nsampt, nrepl) {
  # Ct_nd is a maximum Ct value
  # nsampt is a number of sample types (trt vs. normal)
  # nrepl is a number of replicates per sample type
  
  object<-DataName # data matrix
  ft <- rep("Target",nrow(object))
  # Coltrol gene
  if (!is.null(ControlGene)) { # COntrol Gene is a name of the gene
    ind<-which(rownames(object)==ControlGene)
    ft[ind] <- "Endogenous Control"
  }
  
  fc <- matrix("OK",nrow=nrow(object),ncol=ncol(object))
  fc[which(object>(Ct_nd-0.0001),arr.ind=TRUE)] <- "Undetermined"
  colnames(fc) <- colnames(object)
  rownames(fc) <- rownames(object)
  
  fl <- matrix("Passed",nrow=nrow(object),ncol=ncol(object))
  fl[which(object>(Ct_nd-0.0001),arr.ind=TRUE)] <- "Flagged"
  colnames(fl) <- colnames(object)
  rownames(fl) <- rownames(object)
  
  myData <- new("qPCRset", exprs=object, flag=fl)
  featureNames(myData) <- rownames(object)
  featureType(myData) <- ft
  featureCategory(myData) <- as.data.frame(fc)
  
  sType <- c(rep(LETTERS[1:nsampt], each=nrepl))
  #  PartID <- c(rep(ind, dim(object)[2]))
  #  tab <- data.frame(sampleName=colnames(object), sampleType=sType, dataPartID=PartID)
  tab <- data.frame(sampleName=colnames(object), sampleType=sType)
  phenoData(myData) <- AnnotatedDataFrame(data=tab)
  return(myData)
}

my.sim.data <- function(n, m, k, beta0, beta1, thetas1, true.sigma2, a, b) {
  
  # sample residuals 
  epsilon <- t(rmvnorm(n = m*k, mean = rep(0, n), 
                       sigma = diag(true.sigma2, n))) # use variance
  
  # gene 1 is a Control gene
  # small global shift a vector of dimention n by k*m
  #dj     <- runif(m*k, min=-0.7, max=0.7)
  # dj     <- thetas1[1,]
  #dj     <- mean(dj)-dj
  
  # create data matrix
  mymu   <- thetas1 + epsilon 
  
  # keep Control gene and add djto the rest of genes 
  #mydata <- rbind(Gene_1 = mymu[1,], sweep(mymu[-1,], 2, dj, FUN = "+") )
  mydata <- mymu
  
  # we can check is dj is calculated correctly by calculation dj1, 
  # it should be zero
  #dj1    <- mydata[1,]
  #dj1    <- dj1 - mean(dj1) # should be zero
  
  # now we need to calculate the pij and nd, dont need epsilon in mymu
  #mymu <- sweep(thetas1[-1,], 2, dj, FUN = "+")
  
  # make mymu into array gene by replicate by sample type
  mymu.arr <- array(mymu, c(n,k,m), 
                    dimnames = list(gene=rownames(mymu),
                                    repl=b, sampl=a))
  mymu.arr
  # thanspose the data
  #thetas.arr      <- aperm(thetas.arr, c(1,3,2))
  dimnames(mymu.arr)
  dim(mymu.arr)
  
  # calculate P(Z=0) and identify nondetects
  median.mu <- apply(mymu.arr, c(1,3), median)
  
  logit.pij <- beta0+beta1*mydata
  pij       <- exp(logit.pij)/(1+exp(logit.pij))
  # generate ND
  nd        <- matrix(rbinom(dim(mydata)[1]*dim(mydata)[2], 1, prob=pij), ncol=dim(mydata)[2])
  
  # transform nd back to matrix 
  dimnames(nd) <- dimnames(mymu)
  
  return(list(data=mydata, nd=nd))
}

my.res.beta <- function(object, tst1, data.orig) {
  object$fit$coefficients
  summary(object$fit)$cov.unscaled
  out <- list(out.beta   = c(out.miss.perc = sum(featureCategory(tst1)!="OK")/(sum(featureCategory(tst1)!="OK") + 
                                                                                 sum(featureCategory(tst1)=="OK") ),
                             beta0 = object$fit$coefficients[1], 
                             beta0var = summary(object$fit)$cov.unscaled[1,1],
                             beta1 = object$fit$coefficients[2], 
                             beta1var = summary(object$fit)$cov.unscaled[2,2],
                             covarb0b1 = summary(object$fit)$cov.unscaled[2,1]),
              out.theta  = object$theta,
              out.sigma2 = object$sigma2[,1],
              data.orig  = data.orig)
  return(out)
}

perform.sim.theta.bySample <- function(nsim, n, m, k, beta0, beta1, mean0, sd0, linkglm) {
  # generate the data
  mean.theta  <- rtruncnorm(m, a = 20, b = 40.5, mean = mean0, sd = sd0)
  sigma.theta <- runif(m, min=3, max=3) # variances
  # residuals of the epsilon, each gene has its own sd, genes are independent
  true.sigma2  <- runif(n, min=sqrt(0.06), max=sqrt(1.3)) # variances
  # generate m by n matrix of values
  thetas <- rtmvnorm(n, mean = mean.theta,
                     sigma = diag(sigma.theta, length(mean.theta)),
                     lower=rep(0, length = length(mean.theta)),
                     upper=rep( Inf, length = length(mean.theta)),
                     algorithm=c("rejection"))
  rownames(thetas) <- paste0("Gene_", c(1:dim(thetas)[1]))
  thetas1 <- t(apply((thetas), 1, function(x) rep(x, each=k)))
  # lables for sample type and replicate
  a <- paste0("samp_", c(1:m))
  b <- paste0("rep_", c(1:k))
  colnames(thetas1) <- apply(cbind(rep(a,each=k), rep(b, m)), 1, paste, collapse=":")
  true.table <- list(true.beta = c(beta0, beta1), 
                     true.thetas = (thetas),
                     true.sigma2 = true.sigma2)
  
  res.table <- list()
  #  for (h in 1:nsim) {
  foreach ( h =seq(1:nsim) ) %dopar% {
    set.seed(h)
    cat("\n ###### Working on data set number", h, "###### \n")
    repeat {
      mydata <- my.sim.data(n=n, m=m, k=k, 
                            beta0=beta0, beta1=beta1, 
                            thetas1 = thetas1, true.sigma2 = true.sigma2, a, b)
      cat("\n number of ND is",sum(mydata$nd), "\n") 
      cat("\n maximum number of NDs per gene is",max(apply(mydata$nd, 1, sum)), "\n") 
      if (sum(mydata$nd) > 0 & max(apply(mydata$nd, 1, sum)) < k*m){
        break
      } 
    }
    # mycontrol     <- c("Gene_1")  
    #  ControlGene <- mycontrol 
    mydata.nd    <- mydata$data * (1-mydata$nd)
    # replace by 40
    mydata.nd40  <- mydata.nd
    mydata.nd40[which(mydata$nd==1, arr.ind = TRUE)] <- 40
    # truncate at 40
    mydata.nd40[which(mydata$data>=40, arr.ind = TRUE)] <- 40
    mydata.nd40  <- rbind(mydata.nd40)
    
    DataName <- mydata.nd40
    Ct_nd = 40
    nsampt = m
    nrepl = k
    #create a qPCRset object
    tst1 <- qPCRdata(DataName=DataName, Ct_nd=40, nsampt = nsampt, nrepl=k)
    # impute nondetected values using qpcrImpute function from package "nondetects"
    tst1Impute<-qpcrImpute(tst1, groupVars=c("sampleType"), 
                            outform="Param", linkglm = linkglm)
    #save res of iter h in the res.table
    res.table[[ h ]] <- my.res.beta(object = tst1Impute, tst1 = tst1, data.orig= mydata)
    cat("\n ###### Program created data set number",  h, "###### \n")
    cat("\n", round(res.table[[h]]$out.beta[c(1,2,4)], 4) ,"\n")
  }
  sim.output <- list(true.table = true.table, res.table = res.table)
  #res.table
  save(sim.output, file = paste0("res_",n, "_", m, "_", k, "_", nsim, linkglm, ".RData"))
}
# the end of the perform.sim function

# function to do a SI
perform.sim.theta.bySample.SI <- function(nsim, n, m, k, beta0, beta1, mean0, sd0, linkglm) {
  # generate the data
  mean.theta  <- rtruncnorm(m, a = 20, b = 40.5, mean = mean0, sd = sd0)
  sigma.theta <- runif(m, min=3, max=3) # variances
  # residuals of the epsilon, each gene has its own sd, genes are independent
  true.sigma2  <- runif(n, min=sqrt(0.06), max=sqrt(1.3)) # variances
  # generate m by n matrix of values
  thetas <- rtmvnorm(n, mean = mean.theta,
                     sigma = diag(sigma.theta, length(mean.theta)),
                     lower=rep(0, length = length(mean.theta)),
                     upper=rep( Inf, length = length(mean.theta)),
                     algorithm=c("rejection"))
  rownames(thetas) <- paste0("Gene_", c(1:dim(thetas)[1]))
  thetas1 <- t(apply((thetas), 1, function(x) rep(x, each=k)))
  # lables for sample type and replicate
  a <- paste0("samp_", c(1:m))
  b <- paste0("rep_", c(1:k))
  colnames(thetas1) <- apply(cbind(rep(a,each=k), rep(b, m)), 1, paste, collapse=":")
  true.table <- list(true.beta = c(beta0, beta1), 
                     true.thetas = (thetas),
                     true.sigma2 = true.sigma2)
  
  res.table <- list()
  #  for (h in 1:nsim) {
  foreach ( h = seq(1:nsim), .export=ls(envir=globalenv())) %dopar% {
    set.seed(h)
    cat("\n ###### Working on data set number", h, "###### \n")
    repeat {
      mydata <- my.sim.data(n=n, m=m, k=k, 
                            beta0=beta0, beta1=beta1, 
                            thetas1 = thetas1, true.sigma2 = true.sigma2, a, b)
      cat("\n number of ND is",sum(mydata$nd), "\n") 
      cat("\n maximum number of NDs per gene is",max(apply(mydata$nd, 1, sum)), "\n") 
      if (sum(mydata$nd) > 0 & max(apply(mydata$nd, 1, sum)) < k*m){
        break
      } 
    }
    # mycontrol     <- c("Gene_1")  
    #  ControlGene <- mycontrol 
    mydata.nd    <- mydata$data * (1-mydata$nd)
    # replace by 40
    mydata.nd40  <- mydata.nd
    mydata.nd40[which(mydata$nd==1, arr.ind = TRUE)] <- 40
    # truncate at 40
    mydata.nd40[which(mydata$data>=40, arr.ind = TRUE)] <- 40
    mydata.nd40  <- rbind(mydata.nd40)
    
    DataName <- mydata.nd40
    Ct_nd = 40
    nsampt = m
    nrepl = k
    #create a qPCRset object
    tst1 <- qPCRdata(DataName=DataName, Ct_nd=40, nsampt = nsampt, nrepl=k)
    # impute nondetected values using qpcrImpute function from package "nondetects"
    tst1Impute<-qpcrImpute(tst1, groupVars=c("sampleType"), 
                            outform="Single", linkglm = linkglm)
    #save res of iter h in the res.table
    res.table[[ h ]] <- list(object = tst1Impute, tst1 = tst1, data.orig= mydata)
    cat("\n ###### Program created data set number",  h, "###### \n")
    #cat("\n", round(res.table[[h]]$out.beta[c(1,2,4)], 4) ,"\n")
  }
  sim.output <- list(true.table = true.table, res.table = res.table)
  #res.table
  save(sim.output, file = paste0("res_SI_",n, "_", m, "_", k, "_", nsim, linkglm, ".RData"))
}
# the end of the perform.sim function

# function to do a SI
perform.sim.theta.bySample.MI <- function(nsim, n, m, k, beta0, beta1, mean0, sd0, linkglm) {
  # generate the data
  mean.theta  <- rtruncnorm(m, a = 20, b = 40.5, mean = mean0, sd = sd0)
  sigma.theta <- runif(m, min=3, max=3) # variances
  # residuals of the epsilon, each gene has its own sd, genes are independent
  true.sigma2  <- runif(n, min=sqrt(0.06), max=sqrt(1.3)) # variances
  # generate m by n matrix of values
  thetas <- rtmvnorm(n, mean = mean.theta,
                     sigma = diag(sigma.theta, length(mean.theta)),
                     lower=rep(0, length = length(mean.theta)),
                     upper=rep( Inf, length = length(mean.theta)),
                     algorithm=c("rejection"))
  rownames(thetas) <- paste0("Gene_", c(1:dim(thetas)[1]))
  thetas1 <- t(apply((thetas), 1, function(x) rep(x, each=k)))
  # lables for sample type and replicate
  a <- paste0("samp_", c(1:m))
  b <- paste0("rep_", c(1:k))
  colnames(thetas1) <- apply(cbind(rep(a,each=k), rep(b, m)), 1, paste, collapse=":")
  true.table <- list(true.beta = c(beta0, beta1), 
                     true.thetas = (thetas),
                     true.sigma2 = true.sigma2)
  
  res.table <- list()
  #  for (h in 1:nsim) {
  foreach ( h = seq(1:nsim), .export=ls(envir=globalenv())) %dopar% {
    set.seed(h)
    cat("\n ###### Working on data set number", h, "###### \n")
    repeat {
      mydata <- my.sim.data(n=n, m=m, k=k, 
                            beta0=beta0, beta1=beta1, 
                            thetas1 = thetas1, true.sigma2 = true.sigma2, a, b)
      cat("\n number of ND is",sum(mydata$nd), "\n") 
      cat("\n maximum number of NDs per gene is",max(apply(mydata$nd, 1, sum)), "\n") 
      if (sum(mydata$nd) > 0 & max(apply(mydata$nd, 1, sum)) < k*m){
        break
      } 
    }
    # mycontrol     <- c("Gene_1")  
    #  ControlGene <- mycontrol 
    mydata.nd    <- mydata$data * (1-mydata$nd)
    # replace by 40
    mydata.nd40  <- mydata.nd
    mydata.nd40[which(mydata$nd==1, arr.ind = TRUE)] <- 40
    # truncate at 40
    mydata.nd40[which(mydata$data>=40, arr.ind = TRUE)] <- 40
    mydata.nd40  <- rbind(mydata.nd40)
    
    DataName <- mydata.nd40
    Ct_nd = 40
    nsampt = m
    nrepl = k
    #create a qPCRset object
    tst1 <- qPCRdata(DataName=DataName, Ct_nd=40, nsampt = nsampt, nrepl=k)
    # impute nondetected values using qpcrImpute function from package "nondetects"
    tst1Impute<-qpcrImpute(tst1, groupVars=c("sampleType"), 
                            outform="Multy", linkglm = linkglm)
    #save res of iter h in the res.table
    res.table[[ h ]] <- list(object = tst1Impute, tst1 = tst1, data.orig= mydata)
    cat("\n ###### Program created data set number",  h, "###### \n")
    #cat("\n", round(res.table[[h]]$out.beta[c(1,2,4)], 4) ,"\n")
  }
  sim.output <- list(true.table = true.table, res.table = res.table)
  #res.table
  save(sim.output, file = paste0("res_MI_",n, "_", m, "_", k, "_", nsim, linkglm, ".RData"))
}
# the end of the perform.sim.MI function

## dim of the data n, m, k; gene 1 is a Control gene
n = 16 # n=90
m = 6  
#k=3 # k= 4 # k=6 # k = 10
nsim <- 100 # nsim = 500
beta0 <- -35.7
beta1 <- 1

# set a random seed  #set.seed(291086)
### logit
#set.seed(10914)
#perform.sim.theta.bySample(nsim=100, n=16, m=6, k=3, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done

### probit
set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # done

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # done

### cloglog
set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # done

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # done

# set a random seed  #set.seed(291086)
### logit
# set.seed(10914)
# perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=3, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done

### probit
set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # done

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "probit")  # done

### cloglog
set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # done

set.seed(10914)
perform.sim.theta.bySample.SI(nsim=100, n=16, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "cloglog")  # done

# Multiple imput

### logit
set.seed(10914)
perform.sim.theta.bySample.MI(nsim=100, n=16, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

######################  90 genes  ####################
### logit

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=90, m=6, k=4, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # "logit" , "cloglog"  done    

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=90, m=6, k=6, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done

set.seed(10914)
perform.sim.theta.bySample(nsim=100, n=90, m=6, k=10, beta0, beta1, mean0 = 31, sd0 = 3.5, linkglm = "logit")  # done
