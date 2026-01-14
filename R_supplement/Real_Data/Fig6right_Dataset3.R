# set working dirrectory
setwd("D:/Rochester books/Research/Real_Data")
library(nondetects)
library(HTqPCR)

# set seed
set.seed(0)

load("Ypert_Single_dj.RData")
GroupVar=c("sampleType")
data(sagmb2011)
cts       <- qPCRdata

## Becn1 normalization
dCt       <- normalizeCtData(cts, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt       <- dCt[-which(featureNames(dCt)=="Becn1"),]

iD        <- which(featureCategory(dCt)=="OK", arr.ind=T)
i40       <- which(featureCategory(dCt)=="Imputed", arr.ind=T)

## look at different conditions
if (length(GroupVar)>1) {
  conds <- apply(pData(dCt)[,GroupVar],1,paste,collapse=":")} else{
    conds <- pData(dCt)[,GroupVar]}

## switch to negative dCt values
exprs(dCt) <- -exprs(dCt)
#for plot 3
condit    <- "Wnt9a:1"
genename  <- "Sema7a"
ys        <- exprs(dCt)[genename,which(conds==condit)]
NDindic   <- as.matrix(featureCategory(dCt)[genename,which(conds==condit)], nrow=1)
xs        <- rep(2, length(ys))

#3 replaced by 40
iCtrl      <- which(pData(sagmb2011)$sampleType=="Control")
cts40     <- sagmb2011[,-iCtrl]
dCt40     <- normalizeCtData(cts40, norm = "deltaCt", deltaCt.genes = "Becn1")
dCt40     <- dCt40[-which(featureNames(dCt40)=="Becn1"),]

## switch to negative dCt values
exprs(dCt40) <- -exprs(dCt40)
#for plot 3
ys1       <- exprs(dCt40)[genename,which(conds==condit)]
ys        <- c(ys1, ys)
NDindic   <- c(as.matrix(featureCategory(dCt40)[genename,which(conds==condit)], nrow=1),NDindic)
xs        <- c(rep(1, length(ys1)), xs)
qPCRdata  <- NULL

#4 Multiple imputation vary all
load(file="Ypert_multy_10.RData")
#class(qPCRdata)
ctsMulty  <- qPCRdata[[1]]
dCtMulty  <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty  <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1       <- exprs(dCtMulty)[genename,which(conds==condit)]
ys        <- c(ys, ys1)
NDindic   <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(3, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

#5 Multiple imputation, noise  
load(file="Ypert_noise.RData")
#class(qPCRdata)
ctsMulty   <- qPCRdata[[1]]
dCtMulty   <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty   <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1        <- exprs(dCtMulty)[genename,which(conds==condit)]
ys         <- c(ys, ys1)
NDindic    <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(4, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

#6  Multiple imputation, theta and noise  
load(file="Ypert_theta_noise.RData")
#class(qPCRdata)
ctsMulty   <- qPCRdata[[1]]
dCtMulty   <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty   <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1        <- exprs(dCtMulty)[genename,which(conds==condit)]
ys         <- c(ys, ys1)
NDindic    <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(5, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

#7  Multiple imputation, theta and noise  
load(file="Ypert_theta.RData")
#class(qPCRdata)
ctsMulty   <- qPCRdata[[1]]
dCtMulty   <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty   <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1        <- exprs(dCtMulty)[genename,which(conds==condit)]
ys         <- c(ys, ys1)
NDindic    <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(6, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

#8  Multiple imputation, theta and noise  
load(file="Ypert_fit.RData")
#class(qPCRdata)
ctsMulty   <- qPCRdata[[1]]
dCtMulty   <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty   <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1        <- exprs(dCtMulty)[genename,which(conds==condit)]
ys         <- c(ys, ys1)
NDindic    <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(7, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

#9  Multiple imputation, theta and noise  
load(file="Ypert_fit_theta.RData")
#class(qPCRdata)
ctsMulty   <- qPCRdata[[1]]
dCtMulty   <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
dCtMulty   <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]

## switch to negative dCt values
exprs(dCtMulty) <- -exprs(dCtMulty)
#for plot 3
ys1        <- exprs(dCtMulty)[genename,which(conds==condit)]
ys         <- c(ys, ys1)
NDindic    <- c(NDindic, as.matrix(featureCategory(dCtMulty)[genename,which(conds==condit)], nrow=1))

for (k in 2:10) {
  ctsMulty <- qPCRdata[[k]]
  dCtMulty <- normalizeCtData(ctsMulty, norm = "deltaCt", deltaCt.genes = "Becn1")
  dCtMulty <- dCtMulty[-which(featureNames(dCtMulty)=="Becn1"),]
  ## switch to negative dCt values
  exprs(dCtMulty) <- -exprs(dCtMulty)
  indImp   <- which(conds==condit)[which(featureCategory(dCt)[genename,which(conds==condit)]=="Imputed")]
  newImp   <- exprs(dCtMulty)[genename,indImp]
  ys       <- c(ys, newImp)
  newInd   <- as.matrix(featureCategory(dCtMulty)[genename,indImp], nrow=1)
  NDindic  <- c(NDindic, newInd)
}
xs         <- c(xs, rep(8, length=length(ys)-length(xs)))
qPCRdata   <- NULL
dCtMulty   <- NULL

######################################################################
######################################################################
######################################################################
######################################################################

## PLOT 3 scatter plot
####################
indic      <- rep(NA,length=length(NDindic))
for (k in 1:length(NDindic))
{
  if (NDindic[k] == "OK") 
  {
    indic[k] = 21
  } else  {
    indic[k] = 8}
}

xdHDAC     <- xs[which(indic==8)]
ydHDAC     <- ys[which(indic==8)]
xndHDAC    <- xs[which(indic==21)]
yndHDAC    <- ys[which(indic==21)]

xrngHDAC   <- range(xs)
yrngHDAC   <- range(ys)

# plot add margines later
# x-axis lable 
boxesHDAC  <- c("ND=40",
                "SI",
                "MI: all", 
                expression(paste("MI: ",epsilon)),  
                expression(paste("MI: ",theta, "," , epsilon)),
                expression(paste("MI: ",theta)), 
                "MI: fit", 
                expression(paste("MI: fit, ",theta)))
par(mar=c(7.1, 4.1, 1, 2.1))
plot(x=jitter(xdHDAC,amount=0.1), y=ydHDAC, pch=8, xaxt="n", yaxt="n", 
     ylab=expression(paste("-",Delta, Delta,"Ct",sep="")), xlab="", main="", 
     xlim=xrngHDAC+c(-0.4,0.4), ylim=yrngHDAC, cex.lab=1.5)
points(x=jitter(xndHDAC,amount=0.1),y=yndHDAC,pch=21)
rng        <- xrngHDAC-c(-0.5,0.5)
abline(v=seq(from=rng[1],to=rng[2],by=1))
axis(side=2, at=c(-14,-12, -10, -8), cex.axis=1.5)
axis(side=1, at=1:(length(boxesHDAC)), labels=boxesHDAC, las=2, padj=0.5, cex.axis=1.5)
legend("bottomright",c("observed","non-detect"),pch=c(21,8), bg="white")
mtext("Perturbed Gene",side=1,line=6,cex=1.5)
