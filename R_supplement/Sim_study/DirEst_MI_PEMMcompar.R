# compar. with PEMM
#rm(list=ls())
########################## RUN FROM HERE    
# data and initials for make it work
#setwd("/axiom/home/vsherina/Downloads/Tanzy_new_model/ValJags")
library(HTqPCR)
library(nondetects)
library(limma)
library(mvtnorm)
library(coda)
library(rjags)
library(mcgibbsit)
library(base)
library(invgamma)
library(xtable)
library("mcmcse")
library(vioplot)
##########################################################################
# In this program 3 different data sets are used:
# read data, genes that have missing data in original geneset in nature2008 are removed ("real")
# data where values are masked at 30 ("dat")
# nondetects output ("MLE")
##########################################################################

###########################################################################
#load the data
data(nature2008)
groupVars=c("sampleType");    objectName <-"nature2008";    object <- get(objectName)

## calculate dj
## select normalization genes
dj <- formula <- batch <- NULL

# object is real data 
if(is.null(dj)){
  i.hk <- grep("control", featureType(object), ignore.case=TRUE)
  if(length(i.hk) == 0){
    message("No control genes found in featureType(object).
            Data will not be normalized.")
    dj <- rep(0, ncol(object))
  } else {
    if(length(i.hk) == 1){
      dj <- exprs(object)[i.hk,]
      dj <- dj-mean(dj)
    } else{
      dj <- colMeans(exprs(object)[i.hk,])
      dj <- dj-mean(dj)
    }
  }
}

## sample grouping variable(s)
if(is.null(formula)){
  if(is.null(groupVars)){
    nrep <- apply(pData(object), 1, paste, collapse=":")
  } else {
    if(length(groupVars)>1){
      nrep <- apply(pData(object)[,groupVars], 1, paste, collapse=":")
    } else nrep <- pData(object)[,groupVars]
  }
  if(class(nrep)!="character"){
    nrep <- droplevels(nrep)
  }
  
  if(length(levels(nrep))==1) formula=as.formula("~1") else formula=as.formula("~0+nrep")
  DesLM=model.matrix(formula,as.data.frame(nrep))
} else{
  formula <- as.formula(formula)
  DesLM=model.matrix(formula, data=pData(object))
}
#nrep
v          <- as.numeric(as.factor(nrep))
m             <- length(unique(v)) 

## function to remove genes with missing values and mask at cutoff
rem40.mask.cutoff <- function(cutoff1, cutoff2, object, objectName, grVar = v, id) {
  # remove all missing data and remove data >= 35 (as in McCall 2014)
  miss.ind <- unique(which(exprs(object) >= cutoff1, arr.ind = TRUE)[,1])
  vv       <- table(grVar)
  m        <- length(unique(grVar))
  ##### not %in% function
  #'%!in%' <- function(x,y)!('%in%'(x,y))
  #ind <- ind[ind %!in% miss.ind]
  ########
  ## creating a new object without 40 as missing and 34 as missing 
  Y_sub           <- object[-miss.ind,]
  es              <- exprs(object)[-miss.ind,]
  if (id == "dat") { es[es>=cutoff2] <- cutoff2 }
  ndArray         <- es
  ndArray[which( es >= cutoff2, arr.ind = T )] <- 1
  ndArray[which( es  < cutoff2, arr.ind = T )] <- 0
  apply(ndArray, 1, function(x) by(x, v, sum)) -> bla
  
  rem <- NULL
  for (i in 1:dim(bla)[2]) {
    for (j in 1:m) {
      if (bla[j,i] == vv[j]) {rem <- c(rem, i)}
    }
  }
  if (!is.null(rem)) {es <- es[-c(rem),]}
  
  colnames(es)   <- sampleNames(object)
  
  
  # create a new oject qPCRset for Yctl 
  ft <- rep("Target",nrow(es))
  ft[rownames(es)=="Becn1"] <- "Endogenous Control"
  
  fc <- matrix("OK",nrow=nrow(es),ncol=ncol(es))
  fc[which(es>=cutoff2,arr.ind=TRUE)] <- "Undetermined"
  colnames(fc) <- colnames(es)
  rownames(fc) <- rownames(es)
  
  fl <- matrix("Passed",nrow=nrow(es),ncol=ncol(es))
  fl[which(es>=cutoff2,arr.ind=TRUE)] <- "Flagged"
  colnames(fl) <- colnames(es)
  rownames(fl) <- rownames(es)
  
  object                  <- new("qPCRset", exprs=es, flag=fl)
  featureNames(object)    <- rownames(es)
  featureType(object)     <- ft
  featureCategory(object) <- as.data.frame(fc)
  phenoData(object)       <- phenoData(Y_sub)
  #object                  <- get(objectName)
  save(object, file = paste0(objectName,"_", cutoff1, cutoff2,"_", id, ".Rda"))
  
  return(object)
}
object <- get(objectName)
## create masked at 40 object ("real")
cutoff1 <- 40; cutoff2 <- 30 
#id = "real"
if(!file.exists(paste0(objectName, "_", cutoff1, cutoff2, "_real.Rda"))) {
  real <- rem40.mask.cutoff(cutoff1, cutoff2, object, objectName, v, id = "real")
} else {
  load(file=paste0(objectName,"_", cutoff1, cutoff2,"_real.Rda")) 
  real <- object}

object <- get(objectName)
## create masked at 30 object ("dat")
if(!file.exists(paste0(objectName, "_", cutoff1, cutoff2, "_dat.Rda"))) {
  dat <- rem40.mask.cutoff(cutoff1, cutoff2, object, objectName, v, id = "dat")
} else {
  load(file=paste0(objectName,"_", cutoff1, cutoff2, "_dat.Rda")) 
  dat <- object}

# percentage of missing 
sum(featureCategory(dat)!="OK")/(sum(featureCategory(dat)!="OK")+sum(featureCategory(dat)=="OK"))
sum(rowSums(featureCategory(dat)!="OK")>0)/dim(featureCategory(dat))[1]
# create "bayes"

# from nondetects package
# run qpcrImpute on YAMC vs p53/Ras mutation, only 2 groups comparison
if(!file.exists(paste0(objectName,cutoff2,"Params_All.Rda"))) {
  object <- dat
  MLE <- qpcrImpute(object, groupVars=c("sampleType"), outform="Param")
  save(MLE, file = paste0(objectName,cutoff2,"Params_All.Rda"))
} else {
  load(file=paste0(objectName,cutoff2,"Params_All.Rda")) 
}

## mean and variances for 2 groups comparison
thetas.MLE1                 <- MLE$theta[ ,"p53"]
thetas.MLE2                 <- MLE$theta[ ,"p53/Ras"]
thetas.MLE3                 <- MLE$theta[ ,"Ras"]
thetas.MLE4                 <- MLE$theta[ ,"YAMC"]

numer <- denom <- NULL  
for (k in 1:m) {denom <- c(denom, sum(v==k)-1) ;  numer <- c(numer, sum(v==k))}
#vv    <- numer 
numer <- sum(numer)
denom <- sum(denom)
my_const <- numer/denom
sigmas2.MLE           <- MLE$sigma2[ ,1]*my_const # constant to make var unbiased

# Data prep
## sample grouping variable(s)
get.nrep <- function(formula, groupVars) {
  if(is.null(formula)){
    if(is.null(groupVars)){
      nrep <- apply(pData(object), 1, paste, collapse=":")
    } else {
      if(length(groupVars)>1){
        nrep <- apply(pData(object)[,groupVars], 1, paste, collapse=":")
      } else nrep <- pData(object)[,groupVars]
    }
    if(class(nrep)!="character"){
      nrep <- droplevels(nrep)
    }
    if(length(levels(nrep))==1) formula=as.formula("~1") else formula=as.formula("~0+nrep")
    DesLM=model.matrix(formula,as.data.frame(nrep))
  } else{
    formula <- as.formula(formula)
    DesLM=model.matrix(formula, data=pData(object))
  }
  #print(formula)
  out <- list(DesLM=DesLM, nrep=nrep)
  return(out)
}

nrepAll <- get.nrep(NULL, groupVars = c("sampleType"))

## take out batch effect
object     <- dat
ind        <- grep("target", featureType(dat), ignore.case=TRUE)
Ct.nobatch <- exprs(object)
Y          <- sweep(Ct.nobatch[ind,], 2, dj)
tst        =  lmFit(Y, design=nrepAll$DesLM)

thetas.object.All <- tst$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigma2.object.All <- (tst$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

# multiple imputation for nature2008 
if(!file.exists(paste0(objectName,cutoff2,"Multy_All.Rda"))) {
  object <- dat
  MLE_MI <- qpcrImpute(object, groupVars=c("sampleType"), outform="Multy")
  save(MLE_MI, file = paste0(objectName,cutoff2,"Multy_All.Rda"))
} else {
  load(file=paste0(objectName,cutoff2,"Multy_All.Rda")) 
}

Ct.nobatch <- exprs(object)
Y          <- sweep(Ct.nobatch[ind,], 2, dj)

tstMI1 <- lmFit(sweep(exprs(MLE_MI[[1]])[ind,], 2, dj), design=nrepAll$DesLM)
tstMI2 <- lmFit(sweep(exprs(MLE_MI[[2]])[ind,], 2, dj), design=nrepAll$DesLM)
tstMI3 <- lmFit(sweep(exprs(MLE_MI[[3]])[ind,], 2, dj), design=nrepAll$DesLM)
tstMI4 <- lmFit(sweep(exprs(MLE_MI[[4]])[ind,], 2, dj), design=nrepAll$DesLM)
tstMI5 <- lmFit(sweep(exprs(MLE_MI[[5]])[ind,], 2, dj), design=nrepAll$DesLM)

thetas.MI1  <- tstMI1$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigmas2.MI1 <- (tstMI1$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

thetas.MI2  <- tstMI2$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigmas2.MI2 <- (tstMI2$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

thetas.MI3  <- tstMI3$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigmas2.MI3 <- (tstMI3$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

thetas.MI4  <- tstMI4$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigmas2.MI4 <- (tstMI4$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

thetas.MI5 <- tstMI5$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
sigmas2.MI5 <- (tstMI5$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

# observed data
object <- real
Ct.nobatch <- exprs(object)#-B
Y          <- sweep(Ct.nobatch[ind,], 2, dj)
tst        =  lmFit(Y, design=nrepAll$DesLM)
vv         <- table(v)
real.thetas.object.All <- tst$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
real.sigma2.object.All <- (tst$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #
real.All.n.object      <- vv

# theta bias
# real - dat
real_obs <- as.vector(round(real.thetas.object.All, 5) - round(thetas.object.All, 5))
# real - MLE
real_MLE <- as.vector(round(real.thetas.object.All, 5) - round(MLE$theta, 5))
# real - MLE_MI 1-5
# combine MI together
thetas.MI <- apply(array(c(thetas.MI1, thetas.MI2, thetas.MI3, thetas.MI4, thetas.MI5), dim=c(dim(thetas.MI1), 5)), c(1,2), mean)
dimnames(thetas.MI)[c(1,2)] <- dimnames(thetas.MI1)
real_MI  <- as.vector(round(real.thetas.object.All, 5) - round(thetas.MI, 5))
# combine data for boxplot without zeros
box.dat <- cbind(Observed = real_obs[real_obs!=0], DirEst = real_MLE[real_obs!=0], MI = real_MI[real_obs!=0])

### add PEMM
library(PEMM)
X = cutoff2-exprs(dat)
X[X==0] <- NA
## If data are MAR or MCAR, by only specifying phi=0,
## a penalized EM algorithm will be performed at default.
PEM.result  = PEMM_fun(X, phi=0)
## By specifying phi=0, lambda=0, K=0, an EM algorithm will be performed.
## Although when n is small, EM may not converge.
EM.result   = PEMM_fun(X, phi=0, lambda=0, K=0)
## Generate data with non-ignorable missingness -- observations with
## lower absolute values are more likely to be missing
PEMM.result = PEMM_fun(X, phi=1)
## Compare the mean estimates for data with MNAR from different methods
cor(colMeans(cutoff2-exprs(real)),PEMM.result$mu)
cor(colMeans(cutoff2-exprs(real)),PEM.result$mu)
cor(colMeans(cutoff2-exprs(real)),EM.result$mu)

real_PEMM <- as.vector(exprs(real) - round((cutoff2-PEMM.result$Xhat), 5 ))
real_PEM  <- as.vector(exprs(real) - (cutoff2-PEM.result$Xhat ) )
real_EM   <- as.vector(exprs(real) - (cutoff2-EM.result$Xhat  ) )

PEMM <- as.vector((cutoff2-PEMM.result$mu) )
PEM  <- as.vector((cutoff2-PEM.result$mu ) )
EM   <- as.vector((cutoff2-EM.result$mu  ) )

PEMMCt <- cutoff2-PEMM.result$Xhat[ind,]
Y            <- sweep(PEMMCt, 2, dj)
tst          =  lmFit(Y, design=nrepAll$DesLM)
vv           <- table(v)
PEMM.thetas  <- tst$coefficients  #aggregate(. ~ c(rep(c("YAMC","mp53/Ras"), each=k)), as.data.frame(t(exprs(cts40))), mean)
PEMM.sigmas2 <- (tst$sigma)^2 #(tst$stdev.unscaled*tst$sigma)^2 #

real_PEMM = as.vector(round(real.thetas.object.All, 5) - round(PEMM.thetas, 5))
# combine data for boxplot without zeros
box.comb <- cbind(Observed = real_obs[real_obs!=0], PEMM = real_PEMM[real_obs!=0], DirEst = real_MLE[real_obs!=0], MI = real_MI[real_obs!=0])
#boxplot(box.dat, main = "Difference in theta estimates from the known values")
#abline(h=0)
## points plot with mean
par(mfrow=c(1,2))
set.seed(2)
plot(x=jitter(rep(1:4, each=nrow(box.comb))), y=box.comb, pch=1, ylab = "", xlab = "", axes=FALSE, xlim = c(0.7, 4.3), ylim=c(-0.55, 1.2))
xticks <- seq(1, 4, 1); xlabls <- c("Truncated", "PEMM", "DirEst", "MI")
axis(1, at = xticks, labels = xlabls, col.axis="black", las=1, tck=-0.05)
axis(2)
box(lty = 1, col = 'black')
title(ylab = expression(paste({{hat(theta)}}[real] - {{hat(theta)}}[method]) ), line=2.5, cex.lab=1)
segments(x0=c(0.75, 1.75, 2.75, 3.75, 4.75), x1=c(1.25,2.25,3.25, 4.25, 5.25), y0=colMeans(box.comb), lwd=3)
abline(h=0, col="black", lwd=2, lty=2)
###### end of PEMM


# sigma bias 
# real - dat
real_obs_sigma <- as.vector(real.sigma2.object.All - sigma2.object.All)
summary(as.vector(real.sigma2.object.All - sigma2.object.All))
mean(as.vector(real.sigma2.object.All - sigma2.object.All))
# real - MLE
real_MLE_sigma <- as.vector(real.sigma2.object.All - sigmas2.MLE)
summary(as.vector(real.sigma2.object.All - sigmas2.MLE))
mean(as.vector(real.sigma2.object.All - sigmas2.MLE))
# real - MLE_MI 1-5
# combine MI together
sigmas2.MI <- apply(cbind(sigmas2.MI1, sigmas2.MI2, sigmas2.MI3, sigmas2.MI4, sigmas2.MI5), 1, mean)
dimnames(sigmas2.MI) <- dimnames(sigmas2.MI1)
real_MI_sigma  <- as.vector(real.sigma2.object.All - sigmas2.MI)
summary(as.vector(real_MI_sigma[real_MI_sigma!=0]))
mean(as.vector(real_MI_sigma[real_MI_sigma!=0]))
##PEMM sigma2
real_PEMM_sigma <-  as.vector(real.sigma2.object.All - PEMM.sigmas2)
#boxplot(real_MI_sigma[real_MI_sigma!=0])

# combine data for boxplot
box.dat <- cbind(Observed = real_obs_sigma[real_obs_sigma!=0], PEMM = real_PEMM_sigma[real_obs_sigma!=0], DirEst = real_MLE_sigma[real_obs_sigma!=0], MI = real_MI_sigma[real_obs_sigma!=0])
#boxplot(box.dat, main = "Difference in sigma2 estimates from the known values")
set.seed(2)
plot(x=jitter(rep(1:4, each=nrow(box.dat))), y=box.dat, pch=1, ylab = "", xlab = "", axes=FALSE, xlim = c(0.7, 4.3), ylim=c(-0.5, 1.1))
xticks <- seq(1, 4, 1); xlabls <- c("Truncated", "PEMM", "DirEst", "MI")
axis(1, at = xticks, labels = xlabls, col.axis="black", las=1, tck=-0.05)
axis(2)
box(lty = 1, col = 'black')
title(ylab = expression(paste( {{hat(sigma)}^2}[real] - {{hat(sigma)}^2}[method]) ), line=2.5, cex.lab=1)
segments(x0=c(0.75,1.75,2.75, 3.75), x1=c(1.25,2.25,3.25, 4.25), y0=colMeans(box.dat), lwd=3)
abline(h=0, col="black", lwd=2, lty=2)

