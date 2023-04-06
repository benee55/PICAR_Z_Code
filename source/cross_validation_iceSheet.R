# Cross Validation
#################################
# Generate Random Effects
##################################
predRE<-function(reVect,phi,sigma2,covFunc,fullDistMat,indList){
  # Occurrence RE
  covMat<-sigma2*covFunc(distMat = fullDistMat , phi=phi)
  reList<-covMat[indList[[2]],indList[[1]]]%*%solve(covMat[indList[[1]],indList[[1]]],reVect)
  return(reList)
}
########################################################################
# Crossvalidation

########################
ice_hurdle_CVMPSE_intercept<-function(mcmcDat,ObsDat,model,reList){
  
  n=length(ObsDat$obs.mod)
  n1=sum(ObsDat$obs.mod!=0)
  ##########################################################################################
  finalMSPE<-vector("numeric")
  predictionList<-list()
  
  ##########################################################################################
  #Occurrence
  ##########################################################################################
  # Random Effect
  cv.W<-reList[[1]]
  # Random Effect + Mean
  cv.beta<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[1]]))
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  if(model=="semi"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]][,-ncol(mcmcDat[[2]])]))
    mean.exp.cv.P<-exp(cv.W.P+cv.beta.P+(mcmcDat[[2]][,"tau2"]/2)) # E[X] for X~LN()
    pred.exp.cv.P<-apply(mean.exp.cv.P,2,mean)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  
  # RMSPE
  finalMSPE<-sqrt(mean((finalPrediction-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , trueObs)
  colnames(predictionList[[3]])<-c("pred","obs")
  
  return(list(finalMSPE=finalMSPE,predictionList))
}


########################
ice_hurdle_CVMPSE<-function(mcmcDat,ObsDat,model,reList){
  
  n=length(ObsDat$obs.mod)
  n1=sum(ObsDat$obs.mod!=0)
  ##########################################################################################
  finalMSPE<-vector("numeric")
  predictionList<-list()
  
  ##########################################################################################
  #Occurrence
  ##########################################################################################
  # Random Effect
  cv.W<-reList[[1]]
  # Random Effect + Mean
  cv.beta<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[1]]))
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  if(model=="semi"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]][,-ncol(mcmcDat[[2]])]))
    mean.exp.cv.P<-exp(cv.W.P+cv.beta.P+(mcmcDat[[2]][,"tau2"]/2)) # E[X] for X~LN()
    pred.exp.cv.P<-apply(mean.exp.cv.P,2,mean)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  
  # RMSPE
  finalMSPE<-sqrt(mean((finalPrediction-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction  , trueObs)
  colnames(predictionList[[3]])<-c("pred","obs")
  
  return(list(finalMSPE=finalMSPE,predictionList))
}



########################
ice_mix_CVMPSE_intercept<-function(mcmcDat,ObsDat,model,reList){
  
  n=length(ObsDat$obs.mod)
  n1=sum(ObsDat$obs.mod!=0)
  ##########################################################################################
  finalMSPE<-vector("numeric")
  predictionList<-list()
  
  ##########################################################################################
  #Occurrence
  ##########################################################################################
  # Random Effect
  cv.W<-reList[[1]]
  # Random Effect + Mean
  cv.beta<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[1]]))
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  if(model=="semi"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-apply(cv.beta.P+cv.W.P,2,function(x) {ifelse(x<0,0,x)})
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  
  # RMSPE
  finalMSPE<-sqrt(mean((finalPrediction-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , trueObs)
  colnames(predictionList[[3]])<-c("pred","obs")
  
  return(list(finalMSPE=finalMSPE,predictionList))
}


########################
ice_mix_CVMPSE<-function(mcmcDat,ObsDat,model,reList){
  
  n=length(ObsDat$obs.mod)
  ##########################################################################################
  finalMSPE<-vector("numeric")
  predictionList<-list()
  
  ##########################################################################################
  #Occurrence
  ##########################################################################################
  # Random Effect
  cv.W<-reList[[1]]
  # Random Effect + Mean
  cv.beta<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[1]]))
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  predictionList[[1]]<-cbind(cv.Results,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  if(model=="semi"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-apply(cv.beta.P+cv.W.P,2,function(x) {ifelse(x<0,0,x)})
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  
  # RMSPE
  finalMSPE<-sqrt(mean((finalPrediction-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction  , trueObs)
  colnames(predictionList[[3]])<-c("pred","obs")
  
  return(list(finalMSPE=finalMSPE,predictionList))
}
