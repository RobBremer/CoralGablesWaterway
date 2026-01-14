# set working dirrectory
setwd("D:\\Rochester books\\Research\\Real_Data")

# bioconductor packages
library(HTqPCR)
library(nondetects)
library(limma)
library(mvtnorm)
library(arm)
# load the data
# removing batch effect
# sagmb2011 SI for Control samples first then non-control samples
#############################################

data(sagmb2011);  GroupVar=c("sampleType"); objectName<-"sagmb2011"
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
Y_sub       <- sagmb2011[,iCtrl]
objectName <-"Yctl"
es <- exprs(Y_sub)
colnames(es) <- sampleNames(sagmb2011)[iCtrl]

# create a new oject qPCRset for Yctl 
ft <- rep("Target",nrow(es))
ft[rownames(es)=="Becn1"] <- "Endogenous Control"

fc <- matrix("OK",nrow=nrow(es),ncol=ncol(es))
fc[which(es>39.9999,arr.ind=TRUE)] <- "Undetermined"
colnames(fc) <- colnames(es)
rownames(fc) <- rownames(es)

fl <- matrix("Passed",nrow=nrow(es),ncol=ncol(es))
fl[which(es>39.9999,arr.ind=TRUE)] <- "Flagged"
colnames(fl) <- colnames(es)
rownames(fl) <- rownames(es)

Yctl <- new("qPCRset", exprs=es, flag=fl)
featureNames(Yctl) <- rownames(es)
featureType(Yctl) <- ft
featureCategory(Yctl) <- as.data.frame(fc)
phenoData(Yctl) <- phenoData(Y_sub)
object     <- get(objectName)
save(object, file = "sagmb2011Control.Rdata")

set.seed(0)
#single imputation for control genes
qPCRdata   <- qpcrImpute(Yctl, groupVars="sampleType", dj=dj[iCtrl], outform = "Single", batch =NULL)
filename   <- paste0("sagmb2011Coltrols_Single_batch_common_dj",".Rdata")
save(qPCRdata, file = filename)
load(filename)
Yctl       <- qPCRdata
Y_sub      <- sagmb2011[,-iCtrl]

es <- NULL
es <- exprs(Y_sub)
colnames(es) <- sampleNames(sagmb2011)[-iCtrl]

####
# create a new qPCRset object Ypert for non-control samples
ft <- rep("Target",nrow(es))
ft[rownames(es)=="Becn1"] <- "Endogenous Control"

fc <- matrix("OK",nrow=nrow(es),ncol=ncol(es))
fc[which(es>39.9999,arr.ind=TRUE)] <- "Undetermined"
colnames(fc) <- colnames(es)
rownames(fc) <- rownames(es)

fl <- matrix("Passed",nrow=nrow(es),ncol=ncol(es))
fl[which(es>39.9999,arr.ind=TRUE)] <- "Flagged"
colnames(fl) <- colnames(es)
rownames(fl) <- rownames(es)

Ypert <- new("qPCRset", exprs=es, flag=fl)
featureNames(Ypert) <- rownames(es)
featureType(Ypert) <- ft
featureCategory(Ypert) <- as.data.frame(fc)
phenoData(Ypert) <- phenoData(Y_sub)

save(Ypert, file = "sagmb2011Pert.Rdata")

# check and make sure that Yctl has Control samples only and 
# Ypert has no control samples
table(pData(Yctl)$sampleType)
table(pData(Ypert)$sampleType)
######
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
B          <- qPCRnormalize(Yctl, dj=dj[iCtrl])
B          <- B[,match(pData(Ypert)$batchID, pData(Yctl)$batchID)]
######################################################################
objectName <-"Ypert"
# creating data sets
# set seed
set.seed(0)
object   <- get(objectName)
#Parameters
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl], 
                       outform="Param", linkglm = "logit")
filename <- paste0(objectName,"_Params_dj",".Rdata")
save(qPCRdata, file = filename)
# set seed
set.seed(0)
#single
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, batch = B,  dj=dj[-iCtrl], 
                       outform="Single")
filename <- paste0(objectName,"_Single_dj",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy theta and noise
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                         outform="Multy", numsam=10, vary_fit=F, linkglm = "logit")
filename  <- paste0(objectName,"_theta_noise",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy all
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                        outform="Multy", numsam=10, linkglm = "logit")
filename  <- paste0(objectName,"_multy_10",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy noise
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                       outform="Multy", numsam=10, vary_model=F, vary_fit=F, linkglm = "logit")
filename  <- paste0(objectName,"_noise",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy theta 
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                       outform="Multy", numsam=10, vary_fit=F, add_noise=F, linkglm = "logit")
filename  <- paste0(objectName,"_theta",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy fit
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                       outform="Multy", numsam=10, vary_model=F, add_noise=F, linkglm = "logit")
filename  <- paste0(objectName,"_fit",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy fit theta
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, batch = B, dj=dj[-iCtrl],
                       outform="Multy", numsam=10, add_noise =F, linkglm = "logit")
filename  <- paste0(objectName,"_fit_theta",".Rdata")
save(qPCRdata, file = filename)
