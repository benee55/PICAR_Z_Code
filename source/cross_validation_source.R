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
val_zip_CVMPSE_intercept<-function(mcmcDat,ObsDat,model,reList){
  
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
  true.XB<-ObsDat$XMat.cv%*% ObsDat$parTruth[3:4]
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.truth<-true.cv<-exp(ObsDat$gpS.cv+true.XB)/(1+exp(ObsDat$gpS.cv+true.XB))
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.truth,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","truth","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
  true.XB.P<-ObsDat$XMat.cv%*% ObsDat$parTruth[7:8]
  # Mean FUnction Exp
  if(model=="zip"){
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.mean.exp.cv.P<-true.cv.P<-exp(ObsDat$gpW.cv+true.XB.P)
  }
  if(model=="semi"){
    exp.cv.P<-apply(cv.beta.P+cv.W.P,2,function(x) {ifelse(x<0,0,x)})
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.mean.exp.cv.P<-true.cv.P<- ifelse(ObsDat$gpW.cv+true.XB.P<0,0,ObsDat$gpW.cv+true.XB.P)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,true.mean.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","truth","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  finalPrediction.truth<-predictionList[[1]][,2]*predictionList[[2]][,2]
  
  # RMSPE
  finalMSPE[1]<-sqrt(mean((finalPrediction-trueObs)^2))
  finalMSPE[2]<-sqrt(mean((finalPrediction.truth-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , finalPrediction.truth , trueObs)
  colnames(predictionList[[3]])<-c("pred","truth","obs")
  
  # Precision: min(obs,pred)/pred
  precisionMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,1]) , 
                      sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,2]))
  # Recall: min(obs,pred)/obs
  recallMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,3]) , 
                   sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,3]))
  f1Mat<-2*(precisionMat*recallMat)/(precisionMat+recallMat)
  
  
  return(list(finalMSPE=finalMSPE,precisionMat=precisionMat,recallMat=recallMat,f1Mat=f1Mat,predictionList))
}


########################
val_zip_CVMPSE<-function(mcmcDat,ObsDat,model,reList){
  
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
  true.XB<-ObsDat$XMat.cv%*% ObsDat$parTruth[3:4]
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.truth<-true.cv<-exp(ObsDat$gpS.cv+true.XB)/(1+exp(ObsDat$gpS.cv+true.XB))
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.truth,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","truth","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean

  true.XB.P<-ObsDat$XMat.cv%*% ObsDat$parTruth[7:8]
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.mean.exp.cv.P<-true.cv.P<-exp(ObsDat$gpW.cv+true.XB.P)
  }
  if(model=="semi"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-apply(cv.beta.P+cv.W.P,2,function(x) {ifelse(x<0,0,x)})
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.mean.exp.cv.P<-true.cv.P<- ifelse(ObsDat$gpW.cv+true.XB.P<0,0,ObsDat$gpW.cv+true.XB.P)
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,true.mean.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","truth","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  finalPrediction.truth<-predictionList[[1]][,2]*predictionList[[2]][,2]
  
  # RMSPE
  finalMSPE[1]<-sqrt(mean((finalPrediction-trueObs)^2))
  finalMSPE[2]<-sqrt(mean((finalPrediction.truth-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , finalPrediction.truth , trueObs)
  colnames(predictionList[[3]])<-c("pred","truth","obs")
  
  # Precision: min(obs,pred)/pred
  precisionMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,1]) , 
                      sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,2]))
  # Recall: min(obs,pred)/obs
  recallMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,3]) , 
                   sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,3]))
  f1Mat<-2*(precisionMat*recallMat)/(precisionMat+recallMat)
  
  
  return(list(finalMSPE=finalMSPE,precisionMat=precisionMat,recallMat=recallMat,f1Mat=f1Mat,predictionList))
}


########################
val_hurdle_CVMPSE_intercept<-function(mcmcDat,ObsDat,model,reList){
  
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
  true.XB<-ObsDat$XMat.cv%*% ObsDat$parTruth[3:4]
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.truth<-true.cv<-exp(ObsDat$gpS.cv+true.XB)/(1+exp(ObsDat$gpS.cv+true.XB))
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.truth,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","truth","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean
  
  true.XB.P<-ObsDat$XMat.cv%*% ObsDat$parTruth[7:8]
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.cv.P<-exp(ObsDat$gpW.cv+true.XB.P)
    true.mean.exp.cv.P<-true.cv.P/(1-exp(-true.cv.P))
  }
  if(model=="semi"){
    cv.beta.P<-t(cbind(1,ObsDat$XMat.cv)%*%t(mcmcDat[[2]][,-ncol(mcmcDat[[2]])]))
    mean.exp.cv.P<-exp(cv.W.P+cv.beta.P+(mcmcDat[[2]][,"tau2"]/2)) # E[X] for X~LN()
    pred.exp.cv.P<-apply(mean.exp.cv.P,2,mean)
    tau2.p<-ObsDat$parTruth[9]
    true.mean.exp.cv.P<-exp(ObsDat$gpW.cv+true.XB.P+(tau2.p/2))# LogNormal
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,true.mean.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","truth","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  finalPrediction.truth<-predictionList[[1]][,2]*predictionList[[2]][,2]
  
  # RMSPE
  finalMSPE[1]<-sqrt(mean((finalPrediction-trueObs)^2))
  finalMSPE[2]<-sqrt(mean((finalPrediction.truth-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , finalPrediction.truth , trueObs)
  colnames(predictionList[[3]])<-c("pred","truth","obs")
  
  # Precision: min(obs,pred)/pred
  precisionMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,1]) , 
                      sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,2]))
  # Recall: min(obs,pred)/obs
  recallMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,3]) , 
                   sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,3]))
  f1Mat<-2*(precisionMat*recallMat)/(precisionMat+recallMat)
  
  
  return(list(finalMSPE=finalMSPE,precisionMat=precisionMat,recallMat=recallMat,f1Mat=f1Mat,predictionList))
}


########################
val_hurdle_CVMPSE<-function(mcmcDat,ObsDat,model,reList){
  
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
  true.XB<-ObsDat$XMat.cv%*% ObsDat$parTruth[3:4]
  # Link FUnction Logit
  exp.cv<-exp(cv.W+cv.beta)/(1+exp(cv.W+cv.beta))
  cv.Results<-apply(exp.cv,2,mean)
  cv.truth<-true.cv<-exp(ObsDat$gpS.cv+true.XB)/(1+exp(ObsDat$gpS.cv+true.XB))
  cv.Obs.bin<-ifelse(ObsDat$obs.cv==0,0,1) # Observations
  # occMSPE[1]<-mean((cv.Results-cv.Obs.bin)^2)
  # occMSPE[2]<-mean((cv.truth-cv.Obs.bin)^2)
  predictionList[[1]]<-cbind(cv.Results,cv.truth,cv.Obs.bin)
  colnames(predictionList[[1]])<-c("pred","truth","obs")
  ##########################################################################################
  #Prevalence
  ##########################################################################################
  # Random Effect
  cv.W.P<-reList[[2]]
  # Random Effect + Mean


  true.XB.P<-ObsDat$XMat.cv%*% ObsDat$parTruth[7:8]
  # Mean FUnction Exp
  if(model=="zip"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]]))
    exp.cv.P<-exp(cv.W.P+cv.beta.P)
    mean.exp.cv.P<-exp.cv.P/(1-exp(-exp.cv.P))
    pred.exp.cv.P<-apply(exp.cv.P,2,mean)
    true.cv.P<-exp(ObsDat$gpW.cv+true.XB.P)
    true.mean.exp.cv.P<-true.cv.P/(1-exp(-true.cv.P))
  }
  if(model=="semi"){
    cv.beta.P<-t((ObsDat$XMat.cv)%*%t(mcmcDat[[2]][,-ncol(mcmcDat[[2]])]))
    mean.exp.cv.P<-exp(cv.W.P+cv.beta.P+(mcmcDat[[2]][,"tau2"]/2)) # E[X] for X~LN()
    pred.exp.cv.P<-apply(mean.exp.cv.P,2,mean)
    tau2.p<-ObsDat$parTruth[9]
    true.mean.exp.cv.P<-exp(ObsDat$gpW.cv+true.XB.P+(tau2.p/2))# LogNormal
  }
  
  # CVMSPE
  trueObs<-ObsDat$obs.cv
  predictionList[[2]]<-cbind(pred.exp.cv.P,true.mean.exp.cv.P,trueObs)
  colnames(predictionList[[2]])<-c("pred","truth","obs")
  
  # Final Match
  finalPrediction<-predictionList[[1]][,1]*predictionList[[2]][,1]
  finalPrediction.truth<-predictionList[[1]][,2]*predictionList[[2]][,2]
  
  # RMSPE
  finalMSPE[1]<-sqrt(mean((finalPrediction-trueObs)^2))
  finalMSPE[2]<-sqrt(mean((finalPrediction.truth-trueObs)^2))
  predictionList[[3]]<-cbind(finalPrediction , finalPrediction.truth , trueObs)
  colnames(predictionList[[3]])<-c("pred","truth","obs")
  
  # Precision: min(obs,pred)/pred
  precisionMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,1]) , 
                      sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,2]))
  # Recall: min(obs,pred)/obs
  recallMat<-cbind(sum(apply(predictionList[[3]][,c(1,3)],1,min))/sum(predictionList[[3]][,3]) , 
                   sum(apply(predictionList[[3]][,c(2,3)],1,min))/sum(predictionList[[3]][,3]))
  f1Mat<-2*(precisionMat*recallMat)/(precisionMat+recallMat)
  
  
  return(list(finalMSPE=finalMSPE,precisionMat=precisionMat,recallMat=recallMat,f1Mat=f1Mat,predictionList))
}


