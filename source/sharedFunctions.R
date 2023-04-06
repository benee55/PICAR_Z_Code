library(fields)
library(classInt)
library(nimble)
library(actuar)
#### 


plotRF<-function(dat,rangeDat,label,location,length.out=10,pch=16,cex=1){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1,alpha = 1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=pch,cex=cex,
       main=label)
}

################################
## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}


# Matern Cov Function + Acceptance Rate function
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}


# Summary 
summaryFunction<-function(mcmcDat,bmseThresh=0.01,time){
  
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                          apply(mcmcDat,2,hpd),
                          apply(mcmcDat,2,accRateFunc),
                          bmmat(mcmcDat)[,2],
                          abs(apply(mcmcDat,2,mean))*bmseThresh,
                          apply(mcmcDat,2,ess),
                          apply(mcmcDat,2,ess)/time)
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                                "Accept","BMSE",paste(bmseThresh,"x mean"),
                                "ESS","ESS/sec")
  return(summaryMat)
}

# NIMBLE FUNCTIONS

expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })


sqecov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-0.5*(dists[i,j]/phi)^2)
      }
    }
    
    return(result)
  })


matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    
    return(result)
  })


# Rank Selection 
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
  zEig<-which.min(abs(GriffithOperatorEig$values))
  
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    for(jk in 1:length(dimSeq)){
      if(jk%%10==0){print(jk)}
      # constant eigenvector
      
      keepM<-c(1:dimSeq[jk],zEig)
      if(process==1){ # Occurrence
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-cbind(mBase,XMatMod)
        X_CV<-cbind(mBaseCV,XMatCV)
        startBeta<-optim(par= c(0,0),fn=dBern,dat=dat[[process]] , X= XMatMod, log=TRUE,control=list(fnscale=-1))$par
        coeff<-optim(par= c(rep(0,ncol(mBase)),startBeta),fn=dBern,dat=dat[[process]] , X= X, log=TRUE,control=list(fnscale=-1))$par
        betaMat[[process]][,jk]<-coeff[c(length(coeff)-1,length(coeff))]
        foo<-exp(X_CV%*%coeff)  
        predCV<-ifelse(foo/(1+foo)>0.5,1,0)
      }else{ # Prevalence
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-cbind(mBase,XMatMod[posInd,])
        X_CV<-cbind(mBaseCV,XMatCV[posInd_CV,])
        if(typ=="semi"){
          coeff<-as.numeric(solve(t(X)%*%X,(t(X)%*%log(dat[[process]]))))
          bar<-foo<-as.numeric(X%*%coeff) #XB
          tau2<-sum((log(dat[[process]])-bar)^2)/(length(dat[[process]])-length(coeff)) #MSE 
          foo<-as.numeric(X_CV%*%coeff) #XB
          predCV<-exp(foo+(tau2/2))
          betaMat[[process]][,jk]<-coeff[c(length(coeff)-1,length(coeff))]
          
        }else{
          startBeta<-optim(par= c(0,0),fn=dTpois,dat=dat[[process]] , X= XMatMod[posInd,], log=TRUE,control=list(fnscale=-1))$par
          coeff<-optim(par= c(rep(0,ncol(mBase)),startBeta),fn=dTpois,dat=dat[[process]] , X= X, log=TRUE,control=list(fnscale=-1))$par
          predCV<-as.numeric(exp(X_CV%*%coeff))
          betaMat[[process]][,jk]<-coeff[c(length(coeff)-1,length(coeff))]
        }
      }
      CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
    }
    
    
  }
  return(list(CVMSPE,betaMat))
}



####################################################################
dTpois<-function(par,dat,X,log=FALSE){
  lambda<-as.numeric(exp(X%*%par))
  # llhd<-sum(dztpois(x=dat,lambda = lambda, log = log))
  llhd<-sum(dat*(X%*%par)-log(exp(lambda)-1)-log(factorial(dat)))
  if(log==TRUE){
    return(llhd)
  }else{
    exp(llhd)
  }
}

dBern<-function(par,dat,X,log=FALSE){
  foo<-as.numeric(exp(X%*%par))
  prob<-foo/(1+foo)
  llhd<-sum(dbinom(x=dat,size=1,prob=prob, log = log))
  if(log==TRUE){
    return(llhd)
  }else{
    exp(llhd)
  }
}




#####################
# Rank Selection 
glmRankSelect<-function(dimSeq,
                        obsMod,
                        obsCV,
                        GriffithOperatorEig,
                        XMatMod,
                        XMatCV,
                        typ,
                        tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]
  
  
  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]
  
  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  
  CVMSPE<-matrix(NA,nrow=2,ncol=length(dimSeq))
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=2,ncol=length(dimSeq))
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    for(jk in 1:length(dimSeq)){
      if(jk%%10==0){print(jk)}
      keepM<-c(1:dimSeq[jk])
      if(process==1){ # Occurrence
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV))
        
        lm1<-glm(dat[[process]]~0+X,family = "binomial")  
        coeffs<-lm1$coefficients
        estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
        betaMat[[process]][,jk]<-estMean
        foo<-exp(X_CV%*%coeffs)  
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
        CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
      }else{ # Prevalence
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X) 
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[c(length(coeffs)-1,length(coeffs))]
          CVMSPE[process,jk]<-mean((exp(predCV)-exp(dat_CV[[process]]))^2)
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
            CVMSPE[process,jk]<-mean(((predCV)-(dat_CV[[process]]))^2)
          }
          
        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")  
          coeffs<-lm1$coefficients
          estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
          betaMat[[process]][,jk]<-estMean
          predCV<-exp(X_CV%*%coeffs)  
          CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
        }
      }
      
    }
    
    
  }
  return(list(CVMSPE,betaMat))
}


glmRankSelect_Intercept<-function(dimSeq,
                                  obsMod,
                                  obsCV,
                                  GriffithOperatorEig,
                                  XMatMod,
                                  XMatCV,
                                  typ,
                                  tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]
  
  
  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]
  
  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  
  CVMSPE<-matrix(NA,nrow=2,ncol=length(dimSeq))
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=3,ncol=length(dimSeq))
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    for(jk in 1:length(dimSeq)){
      if(jk%%10==0){print(jk)}
      keepM<-c(1:dimSeq[jk])
      if(process==1){ # Occurrence
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,1,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,1,XMatCV))
        
        lm1<-glm(dat[[process]]~0+X,family = "binomial")  
        coeffs<-lm1$coefficients
        betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
        foo<-exp(X_CV%*%coeffs)  
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
        CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
      }else{ # Prevalence
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,1,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,1,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X) 
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          CVMSPE[process,jk]<-mean((exp(predCV)-exp(dat_CV[[process]]))^2)
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
            CVMSPE[process,jk]<-mean(((predCV)-(dat_CV[[process]]))^2)
          }
        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")  
          coeffs<-lm1$coefficients
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          predCV<-exp(X_CV%*%coeffs)  
          CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
        }
      }
    }
    
    
  }
  return(list(CVMSPE,betaMat))
}




#####################
# Rank Selection 
glmRankSelect_zeig<-function(dimSeq,
                        obsMod,
                        obsCV,
                        GriffithOperatorEig,
                        XMatMod,
                        XMatCV,
                        typ,
                        tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]
  
  
  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]
  
  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  
  CVMSPE<-matrix(NA,nrow=2,ncol=length(dimSeq))
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=2,ncol=length(dimSeq))
  zEig<-which.min(abs(GriffithOperatorEig$values))
  ########################################################
  # zEig<-which.min(abs(GriffithOperatorEig$values))
  ########################################################
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    for(jk in 1:length(dimSeq)){
      if(jk%%10==0){print(jk)}
      # constant eigenvector
      ########################################################
      # keepM<-c(1:dimSeq[jk],zEig)
      ########################################################
      keepM<-c(1:dimSeq[jk],zEig)
      if(process==1){ # Occurrence
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV))
        
        lm1<-glm(dat[[process]]~0+X,family = "binomial")  
        coeffs<-lm1$coefficients
        estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
        betaMat[[process]][,jk]<-estMean
        foo<-exp(X_CV%*%coeffs)  
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
      }else{ # Prevalence
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X) 
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[c(length(coeffs)-1,length(coeffs))]
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
          }
        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")  
          coeffs<-lm1$coefficients
          estMean<-coeffs[c(length(coeffs)-1,length(coeffs))]
          betaMat[[process]][,jk]<-estMean
          predCV<-exp(X_CV%*%coeffs)  
        }
      }
      CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
    }
    
    
  }
  return(list(CVMSPE,betaMat))
}


#####################
#####################
# Ice Sheet
glmRankSelect_ice<-function(dimSeqOcc,dimSeqPrev,
                            obsMod,
                            obsCV,
                            GriffithOperatorEig,
                            XMatMod,
                            XMatCV,
                            typ,
                            tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]
  
  
  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]
  
  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  maxdimLengt<-max(length(dimSeqOcc),length(dimSeqPrev))
  CVMSPE<-matrix(NA,nrow=2,ncol=maxdimLengt)
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=2,ncol=maxdimLengt)
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    if(process==1){ # Occurrence
      for(jk in 1:length(dimSeqOcc)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqOcc[jk])
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV))
        
        lm1<-glm(dat[[process]]~0+X,family = "binomial")  
        coeffs<-lm1$coefficients
        betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
        foo<-exp(X_CV%*%coeffs)  
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
        CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
      }
    }else{
      for(jk in 1:length(dimSeqPrev)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqPrev[jk])
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X) 
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          CVMSPE[process,jk]<-mean((exp(predCV)-exp(dat_CV[[process]]))^2)
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
            CVMSPE[process,jk]<-mean(((predCV)-(dat_CV[[process]]))^2)
          }
          
        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")  
          coeffs<-lm1$coefficients
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          predCV<-exp(X_CV%*%coeffs)  
          CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
        }
        
      }}
  }
  
  return(list(CVMSPE,betaMat))
}


glmRankSelect_Intercept_ice<-function(dimSeqOcc,
                                      dimSeqPrev,
                                      obsMod,
                                      obsCV,
                                      GriffithOperatorEig,
                                      XMatMod,
                                      XMatCV,
                                      typ,
                                      tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]


  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]

  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  maxdimLengt<-max(length(dimSeqOcc),length(dimSeqPrev))
  CVMSPE<-matrix(NA,nrow=2,ncol=maxdimLengt)
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=3,ncol=maxdimLengt)
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    if(process==1){ # Occurrence
      for(jk in 1:length(dimSeqOcc)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqOcc[jk])
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,1,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,1,XMatCV))

        lm1<-glm(dat[[process]]~0+X,family = "binomial")
        coeffs<-lm1$coefficients
        betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
        foo<-exp(X_CV%*%coeffs)
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
        CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
      }
    }else{
      for(jk in 1:length(dimSeqPrev)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqPrev[jk])
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,1,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,1,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X)
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          CVMSPE[process,jk]<-mean((exp(predCV)-exp(dat_CV[[process]]))^2)
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
            CVMSPE[process,jk]<-mean(((predCV)-(dat_CV[[process]]))^2)
          }

        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")
          coeffs<-lm1$coefficients
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          predCV<-exp(X_CV%*%coeffs)
          CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
        }

      }}
  }

  return(list(CVMSPE,betaMat))
}


glmRankSelect_bivalve<-function(dimSeqOcc,
                            dimSeqPrev,
                            obsMod,
                            obsCV,
                            GriffithOperatorEig,
                            XMatMod,
                            XMatCV,
                            typ,
                            tobit=FALSE){
  posInd<-which(obsMod!=0)
  posInd_CV<-which(obsCV!=0)
  dat<-list();dat_CV<-list()
  dat[[1]]<-ifelse(obsMod==0,0,1)
  dat[[2]]<-obsMod[posInd]
  
  
  dat_CV[[1]]<-ifelse(obsCV==0,0,1)
  dat_CV[[2]]<-obsCV[posInd_CV]
  
  if(typ=="semi"){
    if(tobit==FALSE){
      dat[[2]]<-log(dat[[2]])
      dat_CV[[2]]<-log(dat_CV[[2]])
    }
  }else{
    dat[[2]]<-dat[[2]]-1
    dat_CV[[2]]<-dat_CV[[2]]-1
  }
  maxdimLengt<-max(length(dimSeqOcc),length(dimSeqPrev))
  CVMSPE<-matrix(NA,nrow=2,ncol=maxdimLengt)
  betaMat<-list()
  betaMat[[1]]<-betaMat[[2]]<-matrix(NA,nrow=3,ncol=maxdimLengt)
  for(process in 1:2){
    print(c("Occurrence","Prevalence")[process])
    if(process==1){ # Occurrence
      for(jk in 1:length(dimSeqOcc)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqOcc[jk])
        mBase<-(AMat%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV))
        
        lm1<-glm(dat[[process]]~0+X,family = "binomial")  
        coeffs<-lm1$coefficients
        betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
        foo<-exp(X_CV%*%coeffs)  
        predCV<-ifelse((foo/(1+foo))>0.5,1,0)
        CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
      }
    }else{
      for(jk in 1:length(dimSeqPrev)){
        if(jk%%10==0){print(jk)}
        keepM<-c(1:dimSeqPrev[jk])
        mBase<-(AMat[posInd,]%*%GriffithOperatorEig$vectors[,keepM])
        mBaseCV<-(AMatCV[posInd_CV,]%*%GriffithOperatorEig$vectors[,keepM])
        X<-as.matrix(cbind(mBase,XMatMod[posInd,]))
        X_CV<-as.matrix(cbind(mBaseCV,XMatCV[posInd_CV,]))
        if(typ=="semi"){
          lm1<-lm(dat[[process]]~0+X) 
          coeffs<-lm1$coefficients
          predCV<-as.numeric(X_CV%*%coeffs) #XB
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          if(tobit==TRUE){
            predCV<-ifelse(predCV<0,0,predCV)
          }
          CVMSPE[process,jk]<-mean((exp(predCV)-exp(dat_CV[[process]]))^2)
        }else{
          lm1<-glm(dat[[process]]~0+X,family = "poisson")  
          coeffs<-lm1$coefficients
          betaMat[[process]][,jk]<-coeffs[(ncol(mBase)+1):ncol(X)]
          predCV<-exp(X_CV%*%coeffs)  
          CVMSPE[process,jk]<-mean((predCV-dat_CV[[process]])^2)
        }
        
      }}
  }
  
  return(list(CVMSPE,betaMat))
}