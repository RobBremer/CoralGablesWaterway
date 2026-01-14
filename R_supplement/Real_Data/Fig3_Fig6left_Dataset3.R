library(HTqPCR)
library(nondetects)
### PLOTS
setwd("D:/Rochester books/Research/Real_Data")
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
#plot(variat,param_est[,1])
#lines(c(0,5), c(0,5))

#########################################
#########################################
#########################################

# MA like plot
g.ind<-which(abs(sqrt(param_est[,1])-sqrt(variat))>0.009) # tolerance level
# for sagmb2011 
objectName = "sagmb2011"
g.ind

#save plots as 500x500 png
par(mar=c(5.1, 4.1, 1, 2.1),cex=1,cex.axis=1,cex.main=1,cex.lab=1)
plot(((sqrt(param_est[g.ind,1])+sqrt(variat[g.ind]))/2), (sqrt(param_est[g.ind,1])-sqrt(variat[g.ind])) ,  
     ylim=c(0,max(sqrt(param_est[g.ind,1])-sqrt(variat[g.ind]))+.02), 
     xlim=c(min((sqrt(param_est[g.ind,1])+sqrt(variat[g.ind]))/2)-.01,max((sqrt(param_est[g.ind,1])+sqrt(variat[g.ind]))/2)+.1), 
     ylab="", xlab=expression(paste("(", hat(sigma), scriptscriptstyle(MLE)," + ", hat(sigma), scriptscriptstyle(SI), ") / 2")),
     cex.lab=1.5, pch=20,cex=1.5
)
title(ylab=expression(paste(hat(sigma), scriptscriptstyle(MLE)," - ", hat(sigma), scriptscriptstyle(SI))),mgp=c(2.5,1,0), cex.lab=1.5)

# some of the genes without nondetects have lower var SI
genesNondet<-as.matrix(table(rownames(which(featureCategory(qPCRdata)!="OK", arr.ind = T))))
rownames(genesNondet)
sort(rownames(param_est)[g.ind])
sqrt(param_est[g.ind,1])-sqrt(variat[g.ind])
# more nondet - further away from the line 
round(sqrt(param_est[g.ind,1])-sqrt(variat[g.ind]),5)
# Table 4, Dataset 3 
var.dif.sagmb <- summary(round((param_est[g.ind,1])-(variat[g.ind]),3))
var.dif.sagmb
which(round((param_est[g.ind,1])-(variat[g.ind]),5)==max(round((param_est[g.ind,1])-(variat[g.ind]),5)))
######################################################################


## PLOT 2 boxplots
#library(HTqPCR)
#library(nondetects)

data(sagmb2011);    GroupVar=c("sampleType");             objectName <-"Ypert"

LoadFileName<-paste0(objectName,"_Single_dj.RData")
load(LoadFileName)

cts <- qPCRdata

## Becn1 normalization
dCt <- normalizeCtData(cts, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt <- dCt[-which(featureNames(dCt)=="Becn1"),]

iD <- which(featureCategory(dCt)=="OK", arr.ind=T)
i40 <- which(featureCategory(dCt)=="Imputed", arr.ind=T)

## look at different conditions
if (length(GroupVar)>1) {
  conds <- apply(pData(dCt)[,GroupVar],1,paste,collapse=":")} else{
    conds <- pData(dCt)[,GroupVar]}

## switch to negative dCt values
exprs(dCt) <- -exprs(dCt)

#for plot 3 sagmb
condit <- "Wnt9a:1"
genename <- "Sema7a"

ys <- exprs(dCt)[genename,which(conds==condit)]

resids <- matrix(nrow=nrow(dCt), ncol=ncol(dCt))
for(i in 1:nrow(dCt)){
  for(j in 1:ncol(dCt)){
    ind <- which(conds==conds[j])
    resids[i,j] <- exprs(dCt)[i,j]-mean(exprs(dCt)[i,ind])
  }
}
#3 replaced by 40

iCtrl <- which(pData(sagmb2011)$sampleType=="Control")
cts40 <- sagmb2011[,-iCtrl]
dCt40 <- normalizeCtData(cts40, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt40 <- dCt40[-which(featureNames(dCt40)=="Becn1"),]

## switch to negative dCt values
exprs(dCt40) <- -exprs(dCt40)
#for plot 3
ys <- c(exprs(dCt40)[genename,which(conds==condit)], ys)

resids40 <- matrix(nrow=nrow(dCt40), ncol=ncol(dCt40))
for(i in 1:nrow(dCt40)){
  for(j in 1:ncol(dCt40)){
    ind <- which(conds==conds[j])
    resids40[i,j] <- exprs(dCt40)[i,j]-mean(exprs(dCt40)[i,ind])
  }
}
qPCRdata <- NULL

#4 Multiple imputation vary all
LoadFileName<-paste0(objectName,"_multy_10.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsMulty <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsMulty[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

#5 Multiple imputation, noise  
LoadFileName<-paste0(objectName,"_noise.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsNoise <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsNoise[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

#6  Multiple imputation, theta and noise  
LoadFileName<-paste0(objectName,"_theta_noise.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsThetaNoise <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsThetaNoise[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

#7  Multiple imputation, theta   
LoadFileName<-paste0(objectName,"_theta.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsTheta <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsTheta[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

#8  Multiple imputation, fit  
LoadFileName<-paste0(objectName,"_fit.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsFit <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsFit[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

#9  Multiple imputation, fir, theta 
LoadFileName<-paste0(objectName,"_fit_theta.RData")
load(LoadFileName)

ctsMulty <- qPCRdata[[2]]
dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys <- c(ys, exprs(dCtMulty)[genename,which(conds==condit)])

residsFitTheta <- matrix(nrow=nrow(dCtMulty), ncol=ncol(dCtMulty))
for(i in 1:nrow(dCtMulty)){
  for(j in 1:ncol(dCtMulty)){
    ind <- which(conds==conds[j])
    residsFitTheta[i,j] <- exprs(dCtMulty)[i,j]-mean(exprs(dCtMulty)[i,ind])
  }
}
qPCRdata <- NULL
dCtMulty <- NULL

# combine all in the list

# The following commented out part of the code contains the code for a plot with nondetects encoded as 40
# this plot is for your information only, it is not included in the publication

#boxesHDAC <- list("Observed"=resids[iD], "ND=40"=resids40[i40], "SI"=resids[i40], 
#                  "MI: noise"=residsNoise[i40], "MI: theta"=residsTheta[i40], 
#                  "MI: theta, noise"=residsThetaNoise[i40], "MI: fit"=residsFit[i40], 
#                  "MI: fit, theta"=residsFitTheta[i40], "MI: all"=residsMulty[i40])

#par(mar=c(7.1, 4.1, 1, 2.1),cex=1,cex.axis=1,cex.main=1,cex.lab=1)
#main="Resuduals boxplot, Sagmb2011"

#ylable <- expression(paste("-",Delta, Delta,"Ct residuals",sep=""))

#boxplot(boxesHDAC, xlab="", ylab=ylable, ylim=c(-12,12),xaxt="n", cex.lab=1.5)
#abline(h=0)
# x-axis lable 
#boxesHDAC <- c("Observed","ND=40", "SI", "MI: noise", expression(paste("MI: ",theta)), 
#               expression(paste("MI: ",theta, ", noise")), "MI: fit", 
#               expression(paste("MI: fit, ",theta)), "MI: all")

#axis(side=1, at=1:length(boxesHDAC), labels=boxesHDAC, las=2, padj=0.5, cex.axis=1.2)


### without ND=40
boxesHDAC <- list("Observed"=resids[iD], "SI"=resids[i40], 
                  "MI: noise"=residsNoise[i40], "MI: theta"=residsTheta[i40], 
                  "MI: theta, noise"=residsThetaNoise[i40], "MI: fit"=residsFit[i40], 
                  "MI: fit, theta"=residsFitTheta[i40], "MI: all"=residsMulty[i40])

par(mar=c(7.1, 4.1, 1, 2.1),cex=1,cex.axis=1,cex.main=1,cex.lab=1)

ylable <- expression(paste("-",Delta, Delta,"Ct residuals",sep=""))

boxplot(boxesHDAC, xlab="", ylab=ylable, ylim=c(-7,7),xaxt="n", cex.lab=1.5)
abline(h=0)
# x-axis lable 
boxesHDAC <- c("Observed", "SI", expression(paste("MI: ",epsilon)), expression(paste("MI: ",theta)), 
               expression(paste("MI: ",theta,",", epsilon)), "MI: fit", 
               expression(paste("MI: fit, ",theta)), "MI: all")

axis(side=1, at=1:length(boxesHDAC), labels=boxesHDAC, las=2, padj=0.5, cex.axis=1.5)

