# Simulation nxmxk n- genes, m- sample types, k replicate
setwd("D:\\Rochester books\\Research\\Sim_study\\")
library(HTqPCR)
#############################################
########### TABLE 2 main paper ##############
#############################################

#function that peforms all the calculations
get.results <- function(sim.output, n, m, k) {
nsim <- 100
#load true values fo parameters
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

out.table <- function(n, m, k, res.beta, res.thetas, res.sigma2, ind.big = NULL){
  if (is.null(ind.big)) {seq.loop <- seq(1,nsim)} else {seq.loop <- seq(1,nsim)[-c(ind.big)]}
  # Bias
  Bias.beta <- t(res.beta[seq.loop,]) - tbeta
  bbeta <- apply(Bias.beta, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  bbeta.mean <- apply(Bias.beta, 1, mean)
  
  # Sigma 2
  Bias.sigma2 <- t(m*k/(m*k-m)*res.sigma2) - tsigma2
  bsigma2 <- apply(Bias.sigma2, 1, function(x) quantile(x, c(0.25, 0.5, 0.75)))
  names(Bias.sigma2) <- gsub("^.*?_","sigma 2 ", names(Bias.sigma2))
  bsigma2.mean <- apply(Bias.sigma2, 1, mean)
  quantile(bsigma2.mean, probs= c(0.25,0.5,0.75))
  
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
  Bias.all=rbind(cbind(bbeta.mean,bbeta.mean,bbeta.mean), theta=quantile(as.vector(bthetas), prob=c(0.25, 0.5, 0.75)), 
                 sigma2=quantile(bsigma2.mean, probs= c(0.25,0.5,0.75)))
  # Mean Square Bias
  MSE.beta <- (bbeta.mean)^2 + apply(Bias.beta, 1, sd)^2
  
  MSE.sigma2 <- (bsigma2.mean)^2 + apply(Bias.sigma2, 1, sd)^2
  
  MSE.thetas <- (bthetas)^2 + apply(Bias.thetas, c(1,2), sd)^2
  
  MSE.bias.all=rbind(cbind(MSE.beta,MSE.beta,MSE.beta), 
                     theta=quantile(as.vector(MSE.thetas), prob=c(0.25, 0.5, 0.75)), 
                     sigma2=quantile(MSE.sigma2, probs= c(0.25,0.5,0.75)))
  
  # combine all
  out <-cbind(Bias.all, MSE.bias.all)
  return(out)
}

nsim=100
if (is.null(ind.big)|length(ind.big)==0) {
  results <- out.table(n,m,k,res.beta, res.thetas, res.sigma2)
} else {
  results <- out.table(n,m,k,res.beta, res.thetas, res.sigma2, ind.big)
}
return(results)
}

n=16; m=6; k=6
load("res_16_6_6_100logit.RData"); 
t2logit <- get.results(sim.output, n, m, k)
load("res_16_6_6_100probit.RData"); 
t2probit <- get.results(sim.output, n, m, k)
load("res_16_6_6_100cloglog.RData"); 
t2cloglog <- get.results(sim.output, n, m, k)
#Table 2
round(t2logit, 3)
round(t2probit, 3)
round(t2cloglog, 3)

#############################################
########### TABLE 3 main paper ##############
#############################################
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
n=90; m=6
k=4
load("res_90_6_4_100logit.RData"); 
t3logit90_6_4 <- get.results(sim.output, n, m, k)
k=6
load("res_90_6_6_100logit.RData"); 
t3logit90_6_6 <- get.results(sim.output, n, m, k)
k=10
load("res_90_6_10_100logit.RData"); 
t3logit90_6_10 <- get.results(sim.output, n, m, k)
#Table 3
round(t3logit16_6_4, 3)
round(t3logit16_6_6, 3)
round(t3logit16_6_10, 3)
round(t3logit90_6_4, 3)
round(t3logit90_6_6, 3)
round(t3logit90_6_10, 3)
