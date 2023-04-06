# rm(list=ls())
# detach("package:nimble", unload = TRUE)
# Initialize
library(nimble);library(mvtnorm);library(fields)
setwd('/home/slee287/ZeroInflatedPICAR/code/realExample/bivalve/run') # Setup work directory
# setwd('~/Dropbox/PICAR_ZIP/code/realExample/bivalve/run/') # Setup work directory
source(file="../../../source/batchmeans.R")
source(file = "../../../source/sharedFunctions.R")
load("../samples/BisquareBasisFunctions.RData")

# Designate Observation Data
load("../samples/data.RData")
ObsDat<-list(obs.mod=zMod , grid.location.mod = gridLocationMod , XMat.mod= XMod,
             obs.cv=zCV , grid.location.cv = gridLocationCV , XMat.cv= XCV)
n<-length(ObsDat$obs.mod)
filename="BS_count_hurdle_run.RData"
CVfilename="CV_BS_count_hurdle_run.RData"
posInd<-which(ObsDat$obs.mod!=0)
Z_O<-ifelse(ObsDat$obs.mod==0,0,1)
Z_P<-ObsDat$obs.mod[posInd]

##################
# Generate Basis Functions
##################
distA<-rdist(ObsDat$grid.location.mod,knotsA) 
distB<-rdist(ObsDat$grid.location.mod,knotsB)
distC<-rdist(ObsDat$grid.location.mod,knotsC)
basisA<-bisquareFunc(dist = distA , eps = epsA) # Basis Functions - Resolution A
basisB<-bisquareFunc(dist = distB , eps = epsB) # Basis Functions - Resolution B
basisC<-bisquareFunc(dist = distC , eps = epsC) # Basis Functions - Resolution C
rm(distA,distB,distC)
M<-cbind(basisA,basisB,basisC)
removeBasis<-which(apply(M,2,sum)==0)
M<-M[,-removeBasis]
rm(basisA,basisB,basisC)
distA_CV<-rdist(ObsDat$grid.location.cv,knotsA) 
distB_CV<-rdist(ObsDat$grid.location.cv,knotsB)
distC_CV<-rdist(ObsDat$grid.location.cv,knotsC)
basisA_CV<-bisquareFunc(dist = distA_CV , eps = epsA) # Basis Functions - Resolution A
basisB_CV<-bisquareFunc(dist = distB_CV , eps = epsB) # Basis Functions - Resolution B
basisC_CV<-bisquareFunc(dist = distC_CV , eps = epsC) # Basis Functions - Resolution C
rm(distA_CV,distB_CV,distC_CV)
M_CV<-cbind(basisA_CV,basisB_CV,basisC_CV)
removeBasis<-which(apply(M_CV,2,sum)==0)
M_CV<-M_CV[,-removeBasis]
rm(basisA_CV,basisB_CV,basisC_CV)
##################

MBasis_O<-M
MBasis_P<-M[posInd,]
m_O<-ncol(MBasis_O)
m_P<-ncol(MBasis_P)
n_O<-length(Z_O)
n_P<-length(Z_P)
p_O<-ncol(ObsDat$XMat.mod)
p_P<-ncol(ObsDat$XMat.mod)
MQM_O<-diag(m_O)
MQM_P<-diag(m_P)

###############################################################
n=length(ObsDat$obs.mod)
keepInR <- list(pCov_O=100*diag(p_O),pCov_P=100*diag(p_P),
                MBasis_O = MBasis_O , MBasis_P = MBasis_P , 
                X_O=ObsDat$XMat.mod, X_P=ObsDat$XMat.mod[posInd,])
consts   <- list(n_O=n_O, n_P=n_P, m_O=m_O , m_P=m_P , p_O=p_O , p_P=p_P)
data     <- list(Z_O=Z_O, Z_P= Z_P, 
                 MQM_O = MQM_O , MQM_P = MQM_P,
                 mn=rep(0,n) )
inits    <- list(beta_O=rnorm(3),tau_O=2,delta_O=rnorm(m_O,mean=0,sd=sqrt(2)),
                 beta_P=rnorm(3),tau_P=2,delta_P=rnorm(m_P,mean=0,sd=sqrt(2)))

# Build Model
source(file = "../../../source/HURDLE_source.R")
Rmodel<-nimbleModel(code=BS_count_model_string, data = data,  
                    constants=consts, inits = inits)
cRmodel <- compileNimble(Rmodel)
RmodelConf <- configureMCMC(Rmodel,print = TRUE,useConjugacy = TRUE)
RmodelConf$printSamplers()

RmodelConf$addMonitors(c("beta_O", "tau_O","delta_O",
                         "beta_P", "tau_P","delta_P"))
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
print("End MCMC")

######################################################################
# Plots
parInd<-c(grep(pattern="beta_O",x=colnames(samples_extended)),
          grep(pattern="beta_P",x=colnames(samples_extended)),
          grep(pattern="tau",x=colnames(samples_extended)))
deltaInd<-c(grep(pattern="delta_O",x=colnames(samples_extended)),
            grep(pattern="delta_P",x=colnames(samples_extended)))
print("Summary")
summaryMat<-round(summaryFunction(samples_extended[,parInd],
                                  time=ptFinal[3]),3)
print(summaryMat)
##########################################################################################
# Cross Validation
print("Cross Validation")
source("../../../source/cross_validation_bivalve.R")
thin<-1:nrow(samples_extended)
CVSample<-samples_extended[thin,]
reMatList<-list()
dO_ind<-grep(pattern="delta_O", colnames(samples_extended))
dP_ind<-grep(pattern="delta_P", colnames(samples_extended))
reMatList[[1]]<-t(M_CV%*%t(CVSample[,dO_ind]))
reMatList[[2]]<-t(M_CV%*%t(CVSample[,dP_ind]))

mcmcDat<-list(CVSample[,c(parInd[1:3])],
              CVSample[,c(parInd[4:6])])

cvMSPE<-bivalve_CVMPSE(mcmcDat=mcmcDat,ObsDat=ObsDat,
                       model="hurdle",reList=reMatList)
cvMSPE[[1]] #RMSPE
#Save CV and Summary
save(cvMSPE,summaryMat,file=CVfilename)
