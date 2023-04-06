rm(list=ls())
# detach("package:nimble", unload = TRUE)

# Initialize
library(nimble);library(mvtnorm);library(fields)
setwd('/home/slee287/ZeroInflatedPICAR/code/realExample/bivalve/run') # Setup work directory
# setwd('~/Dropbox/PICAR_ZIP/code/realExample/bivalve/run/') # Setup work directory
source(file="../../../source/batchmeans.R")
source(file = "../../../source/sharedFunctions.R")
# Designate Observation Data
load("../samples/data.RData")
load("../samples/mesh.RData")
ObsDat<-list(obs.mod=zMod , grid.location.mod = gridLocationMod , XMat.mod= XMod,
             obs.cv=zCV , grid.location.cv = gridLocationCV , XMat.cv= XCV)
n<-length(ObsDat$obs.mod)
filename="CorPICAR_count_mix_run.RData"
CVfilename="CV_CorPICAR_count_mix_run.RData"
#Select Rank for PICAR
#Select Rank for PICAR
dimSeqOcc<-seq(2,100,by=1)
dimSeqPrev<-seq(2,100,by=1)
# ADD SEARCH

rsResults<-glmRankSelect_Intercept_bivalve(dimSeqOcc=dimSeqOcc,
                                           dimSeqPrev=dimSeqPrev,
                                           obsMod=ObsDat$obs.mod,
                                           obsCV=ObsDat$obs.cv,
                                           GriffithOperatorEig=GriffithOperatorEig,
                                           XMatMod=ObsDat$XMat.mod,
                                           XMatCV=ObsDat$XMat.cv,
                                           typ="count")


Z<-ObsDat$obs.mod
pBase_O<-dimSeqOcc[which.min(rsResults[[1]][1,])]; print(pBase_O)
pBase_P<-dimSeqPrev[which.min(rsResults[[1]][2,])]; print(pBase_P)
base_O<-GriffithOperatorEig$vectors[,c(1:pBase_O)]
base_P<-GriffithOperatorEig$vectors[,c(1:pBase_P)]
m_O<-ncol(base_O)
m_P<-ncol(base_P)
MBasis_O<-as.matrix(AMat%*%base_O)
MBasis_P<-as.matrix(AMat%*%base_P)
MQM_O<-diag(m_O)
MQM_P<-diag(m_P)

invMQM_O<-diag(m_O)
invMQM_P<-diag(m_P)
LA<-t(chol(MBasis_O%*%t(MBasis_O)+diag(n)*0.00001))
MBasis_P_full<-as.matrix(AMat%*%base_P)
LB<-t(chol(MBasis_P_full%*%t(MBasis_P_full)+diag(n)*0.00001))
LAtLB<-LA%*%t(LB)
matA<-solve(MBasis_P_full%*%t(MBasis_P_full)+diag(n)*0.00001, MBasis_P_full)
matB<-LAtLB%*%matA
baseCor<-t(MBasis_O)%*%solve((MBasis_O%*%t(MBasis_O)+diag(n)*0.00001),matB)


X_O<-cbind(1,ObsDat$XMat.mod)
X_P<-cbind(1,ObsDat$XMat.mod)
p_O<-ncol(X_O)
p_P<-ncol(X_P)
###############################################################
n<-length(Z)
keepInR <- list(pCov_O=100*diag(p_O),pCov_P=100*diag(p_P),
                MBasis_O = MBasis_O , MBasis_P = MBasis_P , 
                X_O=X_O, X_P=X_P)
consts   <- list(n=n, m_O=m_O , m_P=m_P , p_O=p_O , p_P=p_P)
data     <- list(Z=Z ,
                 MQM_O = MQM_O , MQM_P = MQM_P,
                 invMQM_O=invMQM_O,
                 baseCor=baseCor,
                 mn=rep(0,2*n) )
inits    <- list(beta_O=rnorm(p_O),tau_O=2,delta_O=rnorm(m_O,mean=0,sd=sqrt(2)),
                 beta_P=rnorm(p_P),tau_P=2,delta_P=rnorm(m_P,mean=0,sd=sqrt(2)),
                 rho=0.5)

# Build Model
source(file = "../../../source/Mixture_source.R")
Rmodel<-nimbleModel(code=Cor_PICAR_count_model_string, data = data,  
                    constants=consts, inits = inits)
cRmodel <- compileNimble(Rmodel)
RmodelConf <- configureMCMC(Rmodel,print = TRUE,useConjugacy = TRUE)
RmodelConf$printSamplers()

RmodelConf$addMonitors(c("beta_O", "tau_O","delta_O",
                         "beta_P", "tau_P","delta_P",
                         "rho"))
RModelMCMC <- buildMCMC(RmodelConf)
cRModelMCMC <- compileNimble(RModelMCMC)

# MCMC
niter=20000
pt<-proc.time()
cRModelMCMC$run(niter)
samples_extended <- as.matrix(cRModelMCMC$mvSamples)
ptFinal<-proc.time()-pt
ptFinal
save(ptFinal,samples_extended,file=filename)
for(i in 1:6){
  pt<-proc.time()
  cRModelMCMC$run(niter = niter, reset = FALSE)
  ptFinal<-ptFinal+(proc.time()-pt)
  samples_extended <- as.matrix(cRModelMCMC$mvSamples)
  save(ptFinal,samples_extended,file=filename)
}
rm(cRModelMCMC)
rm(RModelMCMC)
rm(RmodelConf)
rm(Rmodel)
rm(cRmodel)
save(ptFinal,samples_extended,base_O,base_P,file=filename)
######################################################################
# Plots
parInd<-c(grep(pattern="beta_O",x=colnames(samples_extended)),
          grep(pattern="beta_P",x=colnames(samples_extended)),
          grep(pattern="tau",x=colnames(samples_extended)),
          grep(pattern="rho",x=colnames(samples_extended)))
deltaInd<-c(grep(pattern="delta_O",x=colnames(samples_extended)),
            grep(pattern="delta_P",x=colnames(samples_extended)))

summaryMat<-round(summaryFunction(samples_extended[,parInd],
                                  time=ptFinal[3]),3)
print(summaryMat)
# thin<-round(seq(floor(nrow(samples_extended)*0.2),nrow(samples_extended),length.out = 2000))
# par(mfrow=c(2,4),mar=c(2,2,2,2))
# plot.ts(samples_extended[thin,parInd[1]]);abline(h=0,col='red')
# plot.ts(samples_extended[thin,parInd[2]]);abline(h=ObsDat$parTruth[3],col='red')
# plot.ts(samples_extended[thin,parInd[3]]);abline(h=ObsDat$parTruth[4],col='red')
# plot.ts(1/samples_extended[thin,parInd[7]])
# plot.ts(samples_extended[thin,parInd[4]]);abline(h=0,col='red')
# plot.ts(samples_extended[thin,parInd[5]]);abline(h=ObsDat$parTruth[7],col='red')
# plot.ts(samples_extended[thin,parInd[6]]);abline(h=ObsDat$parTruth[8],col='red')
# plot.ts(1/samples_extended[thin,parInd[8]]);
# plot.ts(samples_extended[thin,parInd[9]]);
##########################################################################################
# Cross Validation
source("../../../source/cross_validation_bivalve.R")
thin<-1:nrow(samples_extended)
CVSample<-samples_extended[thin,]
reMatList<-list()
dO_ind<-grep(pattern="delta_O", colnames(samples_extended))
dP_ind<-grep(pattern="delta_P", colnames(samples_extended))
reMatList[[1]]<-t(as.matrix(AMatCV%*%base_O)%*%t(samples_extended[thin,dO_ind]))
reMatList[[2]]<-t(as.matrix(AMatCV%*%base_P)%*%t(samples_extended[thin,dP_ind]))

mcmcDat<-list(CVSample[,c(parInd[1:4])],
              CVSample[,c(parInd[5:8])])

cvMSPE<-bivalve_CVMPSE_intercept(mcmcDat=mcmcDat,ObsDat=ObsDat,
                                     model="zip",reList=reMatList)
cvMSPE[[1]] #RMSPE
#Save CV and Summary
save(cvMSPE,summaryMat,file=CVfilename)

# PLot Random Effects
# burnin<-floor(nrow(samples_extended)*0.5)
# dO_ind<-grep(pattern="delta_O", colnames(samples_extended))
# dP_ind<-grep(pattern="delta_P", colnames(samples_extended))
# re<-apply(t(as.matrix(AMat%*%base_O)%*%t(samples_extended[-(1:burnin),dO_ind])),2,mean)
# par(mfrow=c(2,2))
# plotRF(dat=re,rangeDat=c(ObsDat$gpS.mod,re),label="O-Model",location=ObsDat$grid.location.mod,length.out=10)
# plotRF(dat=ObsDat$gpS.mod,rangeDat=c(ObsDat$gpS.mod,re),label="O-Truth",location=ObsDat$grid.location.mod,length.out=10)
# 
# re<-apply(t(as.matrix(AMat%*%base_P)%*%t(samples_extended[-(1:burnin),dP_ind])),2,mean)
# plotRF(dat=re,rangeDat=c(ObsDat$gpW.mod,re),label="P-Model",location=ObsDat$grid.location.mod,length.out=10)
# plotRF(dat=ObsDat$gpW.mod,rangeDat=c(ObsDat$gpW.mod,re),label="P-Truth",
#        location=ObsDat$grid.location.mod,length.out=10)

# cvMSPE_ZIP[[1]] #RMSPE
# [1] 1.873770 1.382939
# > cvMSPE_ZIP[[2]] #Precision
# [,1]      [,2]
# [1,] 0.4651761 0.6120646
# > cvMSPE_ZIP[[3]] #Recall
# [,1]      [,2]
# [1,] 0.4862519 0.5796417
# > cvMSPE_ZIP[[4]] #F1
# [,1]      [,2]
# [1,] 0.4754806 0.5954121