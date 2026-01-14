# set working dirrectory
setwd("D:\\Rochester books\\Research\\Real_Data")
library(HTqPCR)
library(nondetects)
library(limma)
library(mvtnorm)
library(arm)
#load the data
data(oncogene2013); GroupVar=c("sampleType","treatment"); objectName<-"oncogene2013"

# set seed
set.seed(0)
object<-get(objectName)
#Parameters
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Param", linkglm = "logit")
filename<-paste0(objectName,"_Params",".Rdata")
save(qPCRdata, file = filename)
# set seed
set.seed(0)
#single
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Single")
filename<-paste0(objectName,"_Single",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy all
qPCRdata  <- qpcrImpute(object, groupVars=GroupVar, 
                        outform="Multy", numsam=10, linkglm = "logit")
filename<-paste0(objectName,"_multy_10",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy noise
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Multy", numsam=10, vary_model=F, vary_fit=F, linkglm = "logit")
filename<-paste0(objectName,"_noise",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy theta and noise
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Multy", numsam=10, vary_fit=F, linkglm = "logit")
filename<-paste0(objectName,"_theta_noise",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy theta 
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Multy", numsam=10, vary_fit=F, add_noise=F, linkglm = "logit")
filename<-paste0(objectName,"_theta",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy fit
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Multy", numsam=10, vary_model=F, add_noise=F, linkglm = "logit")
filename<-paste0(objectName,"_fit",".Rdata")
save(qPCRdata, file = filename)

# set seed
set.seed(0)
#multy fit theta
qPCRdata <- qpcrImpute(object, groupVars=GroupVar, 
                       outform="Multy", numsam=10, add_noise =F, linkglm = "logit")
filename<-paste0(objectName,"_fit_theta",".Rdata")
save(qPCRdata, file = filename)
