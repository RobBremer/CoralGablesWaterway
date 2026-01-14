library(HTqPCR)
library(nondetects)

#Dataset 2

GroupVar=c("sampleType"); objectName<-"nature2008"; data(nature2008)
setwd("D:/Rochester books/Research/Real_Data")

## calculate dj
object     <- get(objectName)
## select normalization genes
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
####
LoadFileName<-paste0(objectName,"_Single.RData")
load(LoadFileName)

Yctl       <- qPCRdata
table(pData(Yctl)$sampleType)

#function to normalize genes
qPCRnormalize<-function(object, dj=NULL){
  ## select normalization genes
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
  
  ## select target genes
  ind <- grep("target", featureType(object), ignore.case=TRUE)
  if(length(ind) == 0) stop("No target genes found in featureType(object).")
  Y <- sweep(exprs(object)[,], 2, dj)
}
######
## batch effects 
B <-NULL
######################################################################

load("Nature2008_Params.Rdata")
params<-qPCRdata
param_est<-params$sigma2

LoadFileName<-paste0(objectName,"_Single.RData")
load(LoadFileName)

object<-qPCRdata
if (!is.null(B)) {exprs(object) <- exprs(object)-B} #substract batch effect
#function to normalize genes
qPCRnormalize<-function(object, dj=NULL){
  ## select normalization genes
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
  
  ## select target genes
  ind <- grep("target", featureType(object), ignore.case=TRUE)
  if(length(ind) == 0) stop("No target genes found in featureType(object).")
  Y <- sweep(exprs(object)[,], 2, dj)
}
# apply the function
dj=NULL
genes<-qPCRnormalize(object, dj=dj)
exprs(object)<-genes

# normalize
dCt<-normalizeCtData(qPCRdata, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt <- object[-which(featureNames(dCt)=="Becn1"),]
## switch to negative dCt values
exprs(dCt) <- -exprs(dCt)

st                 =  paste(as.character(pData(dCt)$sampleType))
ust                <- unique(st)
tst2               <- lapply(ust, function(n){ exprs(dCt)[,st==n]})
names(tst2)        <- ust
sample.means       <- lapply(tst2, function(x) rowMeans(x))
no.trt.eff         <- lapply(ust, function(x) { tst2[[x]] - sample.means[[x]] })
names(no.trt.eff)  <- ust

by.genes <- matrix(unlist(no.trt.eff), nrow=nrow(no.trt.eff[[1]]))
# or
by.genes <- as.matrix(no.trt.eff[[1]])
for(k in 2:length(no.trt.eff)){
  by.genes <- cbind(by.genes, as.matrix(no.trt.eff[[k]]))
}

dim(by.genes)
par(mfrow=c(1,1))
variat  <- apply(by.genes, 1, function(x) var(as.vector(x))*(length(x)-1)/length(x))

# find the proportion of nondetects
featureCategory(object)[29,]

# for nature tol= 0.009
g.ind<-which(abs(sqrt(param_est[,1])-sqrt(variat))>0.009) # tolerance level
g.ind

# some of the genes without nondetects have lower var SI
genesNondet<-as.matrix(table(rownames(which(featureCategory(qPCRdata)!="OK", arr.ind = T))))
rownames(genesNondet)
sort(rownames(param_est)[g.ind])
sqrt(param_est[g.ind,1])-sqrt(variat[g.ind])
# more nondet - further away from the line 
round(sqrt(param_est[g.ind,1])-sqrt(variat[g.ind]),5)
var.dif.nature <- summary(round((param_est[g.ind,1])-(variat[g.ind]),3))
######################################################################

# Dataset 3
data(sagmb2011);    GroupVar=c("sampleType");             objectName<-"sagmb2011"
## calculate dj
object     <- get(objectName)
## select normalization genes
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
####
iCtrl      <- which(pData(sagmb2011)$sampleType=="Control")
filename   <- paste0("sagmb2011Coltrols_Single_batch_common_dj",".Rdata")
load(filename)
Yctl       <- qPCRdata
load("sagmb2011Pert.Rdata")
table(pData(Yctl)$sampleType)
table(pData(Ypert)$sampleType)
str(exprs(Ypert))

#function to normalize genes
qPCRnormalize<-function(object, dj=NULL){
  ## select normalization genes
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
  
  ## select target genes
  ind <- grep("target", featureType(object), ignore.case=TRUE)
  if(length(ind) == 0) stop("No target genes found in featureType(object).")
  Y <- sweep(exprs(object)[,], 2, dj)
}
######
## batch effects 
B <- qPCRnormalize(Yctl, dj=dj[iCtrl])
B <- B[,match(pData(Ypert)$batchID, pData(Yctl)$batchID)]
objectName <-"Ypert"
######################################################################

load("Ypert_Params_dj.RData")   
params<-qPCRdata
param_est<-params$sigma2

LoadFileName<-paste0(objectName,"_Single_dj.RData")  
load(LoadFileName)
object<-qPCRdata
if (!is.null(B)) {exprs(object) <- exprs(object)-B} #substract batch effect
#function to normalize genes
# apply the function
genes<-qPCRnormalize(object, dj=dj[-iCtrl])
exprs(object)<-genes

# for sagmb2011
# normalize
dCt<-normalizeCtData(qPCRdata, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt <- object[-which(featureNames(object)=="Becn1"),]
## switch to negative dCt values
exprs(dCt) <- -exprs(dCt)

st                 =  paste(as.character(pData(dCt)$sampleType))
ust                <- unique(st)
tst2               <- lapply(ust, function(n){ exprs(dCt)[,st==n]})
names(tst2)        <- ust
sample.means       <- lapply(tst2, function(x) rowMeans(x))
no.trt.eff         <- lapply(ust, function(x) { tst2[[x]] - sample.means[[x]] })
names(no.trt.eff)  <- ust

by.genes <- matrix(unlist(no.trt.eff), nrow=nrow(no.trt.eff[[1]]))
# or
by.genes <- as.matrix(no.trt.eff[[1]])
for(k in 2:length(no.trt.eff)){
  by.genes <- cbind(by.genes, as.matrix(no.trt.eff[[k]]))
}

dim(by.genes)

variat  <- apply(by.genes, 1, function(x) var(as.vector(x))*(length(x)-1)/length(x))
g.ind<-which(abs(sqrt(param_est[,1])-sqrt(variat))>0.009) # tolerance level
# for sagmb2011 
objectName = "sagmb2011"
g.ind

# some of the genes without nondetects have lower var SI
genesNondet<-as.matrix(table(rownames(which(featureCategory(qPCRdata)!="OK", arr.ind = T))))
rownames(genesNondet)
sort(rownames(param_est)[g.ind])
sqrt(param_est[g.ind,1])-sqrt(variat[g.ind])
# more nondet - further away from the line 
round(sqrt(param_est[g.ind,1])-sqrt(variat[g.ind]),5)
var.dif.sagmb <- summary(round((param_est[g.ind,1])-(variat[g.ind]),3))
var.dif.sagmb

#########################################################
GroupVar=c("sampleType","treatment"); objectName<-"oncogene2013"; data("oncogene2013")

## calculate dj
object     <- get(objectName)
load("oncogene2013_Params.RData")
params<-qPCRdata
param_est<-params$sigma2

## batch effects 
B <- NULL

LoadFileName<-paste0(objectName,"_Single.RData")
load(LoadFileName)

object<-qPCRdata
if (!is.null(B)) {exprs(object) <- exprs(object)-B} #substract batch effect
#function to normalize genes
qPCRnormalize<-function(object, dj=NULL){
  ## select normalization genes
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
  
  ## select target genes
  ind <- grep("target", featureType(object), ignore.case=TRUE)
  if(length(ind) == 0) stop("No target genes found in featureType(object).")
  Y <- sweep(exprs(object)[,], 2, dj)
}
# apply the function
dj=NULL
genes<-qPCRnormalize(object, dj=dj)
exprs(object)<-genes
dCt <- object[-which(featureNames(object)=="Becn1"),]
## switch to negative dCt values
exprs(dCt) <- -exprs(dCt)

#for oncogene2013
st                 =  paste(as.character(pData(dCt)$sampleType),
                            as.character(pData(dCt)$treatment), sep=":")
ust                <- unique(st)
tst2               <- lapply(ust, function(n){ exprs(dCt)[,st==n]})
names(tst2)        <- ust
sample.means       <- lapply(tst2, function(x) rowMeans(x))
no.trt.eff         <- lapply(ust, function(x) { tst2[[x]] - sample.means[[x]] })
names(no.trt.eff)  <- ust

by.genes <- matrix(unlist(no.trt.eff), nrow=nrow(no.trt.eff[[1]]))
dim(by.genes)

variat  <- apply(by.genes, 1, function(x) var(as.vector(x))*(length(x)-1)/length(x))

#########################################

# identefy genes below the line
# for oncogene 
tol=0.009
g.ind<-which(abs(sqrt(param_est[,1])-sqrt(variat))>0.009) # tolerance level
g.ind

# some of the genes without nondetects have lower var SI
genesNondet<-as.matrix(table(rownames(which(featureCategory(qPCRdata)!="OK", arr.ind = T))))
rownames(genesNondet)
sort(rownames(param_est)[g.ind])
sqrt(param_est[g.ind,1])-sqrt(variat[g.ind])
# more nondet - further away from the line 
round(sqrt(param_est[g.ind,1])-sqrt(variat[g.ind]),5)
var.dif.oncogene <- summary(round((param_est[g.ind,1])-(variat[g.ind]),3))

## Making a table of diff in var SI vs. MLE
## need to get values of  var.dif.nature and var.dif.sagmb from 2 other files first

t(cbind(var.dif.oncogene, var.dif.nature, var.dif.sagmb))
library(xtable)
xtable(data.frame(t(cbind(var.dif.oncogene, var.dif.nature, var.dif.sagmb)), row.names=c("Dataset 1", "Dataset 2", "Dataset 3")), digits=3)
######################################################################
t(cbind(var.dif.oncogene, var.dif.nature, var.dif.sagmb))
