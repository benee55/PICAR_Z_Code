
rankSelect<-function(dimSeq,
                     obsMod,
                     obsCV,
                     GriffithOperatorEig,
                     XMatMod,
                     XMatCV,
                     typ){
#Inside Function

posInd<-which(obsMod!=0)
posInd_CV<-which(obsCV!=0)
dat<-list();dat_CV<-list()
dat[[1]]<-ifelse(obsMod==0,0,1)
dat[[2]]<-obsMod[posInd]
dat_CV[[1]]<-ifelse(obsCV==0,0,1)
dat_CV[[2]]<-obsCV[posInd_CV]

CVMSPE<-matrix(NA,nrow=2,ncol=length(dimSeq))
betaMat<-list()
betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=2,ncol=length(dimSeq))

for(process in 1:2){
  print(c("Occurrence","Prevalence")[process])
for(jk in 1:length(dimSeq)){
  if(jk%%10==0){print(jk)}
  # constant eigenvector
  zEig<-which.min(abs(GriffithOperatorEig$values))
  keepM<-c(1:dimSeq[jk],zEig)
  mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
  mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
  
  if(process==1){ # Occurrence
    X<-cbind(mBase,XMatMod)
    X_CV<-cbind(mBaseCV,XMatCV)
    startBeta<-optim(par= c(0,0),fn=dBern,dat=dat[[process]] , X= XMatMod, log=TRUE,control=list(fnscale=-1))$par
    coeff<-optim(par= c(rep(0,ncol(mBase)),startBeta),fn=dBern,dat=dat[[process]] , X= X, log=TRUE,control=list(fnscale=-1))$par
    betaMat[[process]][,jk]<-coeff[c(length(coeff)-1,length(coeff))]
    foo<-exp(X_CV%*%coeff)  
    predCV<-ifelse(foo/(1+foo)>0.5,1,0)
  }else{ # Prevalence
    X<-cbind(mBase,XMatMod[posInd,])
    X_CV<-cbind(mBaseCV,XMatCV[posInd_CV,])
    if(typ=="semi"){
      coeff<-as.numeric(solve(t(X)%*%X,(t(X)%*%log(dat[[process]]))))
      foo<-X_CV%*%coeff #XB
      tau2<-sum((log(dat[[process]])-foo)^2)/(length(dat[[process]])-2) #MSE 
      predCV<-exp(foo+(tau2/2))
      
    }else{
      startBeta<-optim(par= c(0,0),fn=dTpois,dat=dat[[process]] , X= XMatMod[posInd,], log=TRUE,control=list(fnscale=-1))$par
      coeff<-optim(par= c(rep(0,ncol(mBase)),startBeta),fn=dTpois,dat=dat[[process]] , X= X, log=TRUE,control=list(fnscale=-1))$par
      predCV<-as.numeric(exp(X_CV%*%coeff))
    }
  }
}
  betaMat[process,jk]<-coeff[c(length(coeff)-1,length(coeff))]
  CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
}
return(list(CVMSPE,betaMat))
}



plot(x=dimSeq[1:length(CVMSPE)],y=CVMSPE,typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[1,1:length(CVMSPE)],typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[2,1:length(CVMSPE)],typ="l")

dimSeq[which.min(CVMSPE)]
#8
betaMat[,which.min(CVMSPE)]

####################################################################################
#Bernouilli

rm(list=ls())

# Initialize
library(nimble);library(mvtnorm);library(fields)
setwd("~/work/PICAR_ZIP/code/run/") # Setup work directory
source(file = "source/sharedFunctions.R")
source(file = "source/Mixture_source.R")
load("../samples/simulatedExample/mix_n500_count_mat.RData")
load("../samples/simulatedExample/mesh/mesh_mix_n500_count_mat.RData")
#Select Rank for PICAR


dimSeq<-seq(2,100,by=2)

dat<-ifelse(ObsDat$obs.mod==0,0,1)
dat_CV<-ifelse(ObsDat$obs.cv==0,0,1)
CVMSPE<-vector("numeric")
betaMat<-matrix(NA,nrow=2,ncol=length(dimSeq))
for(jk in 1:length(dimSeq)){
  print(jk)
  zEig<-which.min(abs(GriffithOperatorEig$values))
  keepM<-c(1:dimSeq[jk],zEig)
  mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
  X<-cbind(mBase,ObsDat$XMat.mod)
  startBeta<-optim(par= c(0,0),fn=dBern,dat=dat , X= ObsDat$XMat.mod, log=TRUE,control=list(fnscale=-1))$par
  coeff<-optim(par= c(rep(0,ncol(mBase)),startBeta),fn=dBern,dat=dat , X= X, log=TRUE,control=list(fnscale=-1))$par
  
  betaMat[,jk]<-coeff[c(length(coeff)-1,length(coeff))]
  
  mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
  X_CV<-cbind(mBaseCV,ObsDat$XMat.cv)
  foo<-exp(X_CV%*%coeff)  
  predCV<-foo/(1+foo)
  predCV<-ifelse(predCV>0.5,1,0)
  CVMSPE[jk]<-mean((predCV-dat_CV)^2)
}

plot(x=dimSeq[1:length(CVMSPE)],y=CVMSPE,typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[1,1:length(CVMSPE)],typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[2,1:length(CVMSPE)],typ="l")

dimSeq[which.min(CVMSPE)] #4
betaMat[,which.min(CVMSPE)]

##########################################################################################
# Linear
rm(list=ls())

# Initialize
library(nimble);library(mvtnorm);library(fields)
setwd("~/work/PICAR_ZIP/code/run/") # Setup work directory
source(file = "source/sharedFunctions.R")
source(file = "source/Mixture_source.R")
load("../samples/simulatedExample/mix_n500_semi_mat.RData")
load("../samples/simulatedExample/mesh/mesh_mix_n500_semi_mat.RData")
#Select Rank for PICAR



dimSeq<-seq(2,100,by=2)
posInd<-which(ObsDat$obs.mod!=0)
posInd_CV<-which(ObsDat$obs.cv!=0)
dat<-ObsDat$obs.mod[posInd]
dat_CV<-ObsDat$obs.cv[posInd_CV]
CVMSPE<-vector("numeric")
betaMat<-matrix(NA,nrow=2,ncol=length(dimSeq))
for(jk in 1:length(dimSeq)){
  print(jk)
  zEig<-which.min(abs(GriffithOperatorEig$values))
  keepM<-c(1:dimSeq[jk],zEig)
  mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
  X<-cbind(mBase,ObsDat$XMat.mod[posInd,])
  coeff<-as.numeric(solve(t(X)%*%X,(t(X)%*%log(dat))))
  
  betaMat[,jk]<-coeff[c(length(coeff)-1,length(coeff))]
  
  mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
  X_CV<-cbind(mBaseCV,ObsDat$XMat.cv[posInd_CV,])
  predCV<-as.numeric(exp(X_CV%*%coeff))
  CVMSPE[jk]<-mean((predCV-dat_CV)^2)
}

plot(x=dimSeq[1:length(CVMSPE)],y=CVMSPE,typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[1,1:length(CVMSPE)],typ="l")
plot(x=dimSeq[1:length(CVMSPE)],y=betaMat[2,1:length(CVMSPE)],typ="l")

dimSeq[which.min(CVMSPE)]
betaMat[,which.min(CVMSPE)]


