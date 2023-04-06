# Truncated Poisson
dtpois <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if(x==0){
      if(log){
        return(-Inf)
      }else{
        return(0)  
      }
      
    }else{
      if (log){
        return(dpois(x, lambda, log = TRUE) - log(1-exp(-lambda)))
      }else{
        return(dpois(x, lambda, log = FALSE)/(1-exp(-lambda)))
      }
    }
    
  })

rtpois <- nimbleFunction(
  run = function(n = integer(), lambda = double()) {
    returnType(integer())
    return(rpois(1, lambda)+1)
  })

registerDistributions(list(
  dtpois = list(
    BUGSdist = "dtpois(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()')
  )))


dTruncPoissonVector <- nimbleFunction(
  run = function(x    = double(1),
                 lambda = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    
    
    returnType(double(0))
    logProb <- sum(dtpois(x, lambda = lambda, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)

dBernoulliVector <- nimbleFunction(
  run = function(x    = double(1),
                 prob = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    returnType(double(0))
    logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)

## Use R's lexical scoping to create R functions to be called
## via nimbleRcall that will find large objects like
## X and pCov in their enclosing environment.



makeRfuns <- function(keepInR) {
  X_O <- keepInR$X_O
  X_P <- keepInR$X_P
  M_O <- keepInR$MBasis_O
  M_P <- keepInR$MBasis_P
  pCov_O <- keepInR$pCov_O
  chol_pCov_O <- chol(pCov_O)
  pCov_P <- keepInR$pCov_P
  chol_pCov_P <- chol(pCov_P)
  
  RcalculateXB_O <- function(beta) {
    as.numeric(X_O %*% beta)
  }
  RcalculateXB_P <- function(beta) {
    as.numeric(X_P %*% beta)
  }
  
  RcalculateW_O <- function(delta) {
    as.numeric(M_O %*% delta)
  }
  RcalculateW_P <- function(delta) {
    as.numeric(M_P %*% delta)
  }
  
  dRdmnormB_O <- function(x, mean, log) {
    nimble::dmnorm_chol(x, mean, cholesky = chol_pCov_O, prec_param = FALSE, log = log)
  }
  
  dRdmnormB_P <- function(x, mean, log) {
    nimble::dmnorm_chol(x, mean, cholesky = chol_pCov_P, prec_param = FALSE, log = log)
  }
  list(RcalculateXB_O = RcalculateXB_O, 
       RcalculateW_O=RcalculateW_O, 
       RcalculateXB_P = RcalculateXB_P, 
       RcalculateW_P=RcalculateW_P, 
       dRdmnormB_O = dRdmnormB_O,
       dRdmnormB_P = dRdmnormB_P)
}

## Create the functions and put them in .GlobalEnv
Rfuns <- makeRfuns(keepInR=keepInR)
RcalculateXB_O_internal <- Rfuns$RcalculateXB_O ## I think these need to be in .GlobalEnv
RcalculateW_O_internal <- Rfuns$RcalculateW_O ## I think these need to be in .GlobalEnv
RcalculateXB_P_internal <- Rfuns$RcalculateXB_P ## I think these need to be in .GlobalEnv
RcalculateW_P_internal <- Rfuns$RcalculateW_P ## I think these need to be in .GlobalEnv
dRdmnormB_O_internal <- Rfuns$dRdmnormB_O
dRdmnormB_P_internal <- Rfuns$dRdmnormB_P

## Define the nimbleRcall nimbleFunctions to call the R functions.
RcalculateXB_O <- nimbleRcall(function(beta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateXB_O_internal')
RcalculateW_O <- nimbleRcall(function(delta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateW_O_internal')
RcalculateXB_P <- nimbleRcall(function(beta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateXB_P_internal')
RcalculateW_P <- nimbleRcall(function(delta = double(1)) {}, returnType = double(1), Rfun = 'RcalculateW_P_internal')
dRdmnormB_O <- nimbleRcall(function(x = double(1), mean = double(1), log = integer(0, default = 0)) {}, returnType = double(), Rfun = 'dRdmnormB_O_internal')
dRdmnormB_P <- nimbleRcall(function(x = double(1), mean = double(1), log = integer(0, default = 0)) {}, returnType = double(), Rfun = 'dRdmnormB_P_internal')

## Write the model using the nimbleRcalls.  Note that X and pCov do not appear in the model.


####################################################
# PICAR Count
PICAR_count_model_string <- nimbleCode({
  
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dtpois(lambda[j])
  }
  lambda[1:n_P] <- exp(W_P[1:n_P]+XB_P[1:n_P])
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])

  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  precMat_P[1:m_P,1:m_P]<-tau_P * MQM_P[1:m_P,1:m_P]
  # Process Model
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = mn[1:m_P], prec = precMat_P[1:m_P,1:m_P])
  # Parameter Model
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
  })


####################################################
# PICAR SemiContinuous
PICAR_semi_model_string <- nimbleCode({
  # Data Model
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dlnorm(meanlog=muLog[j],varlog=tau2)
  }
  
  muLog[1:n_P] <- W_P[1:n_P]+XB_P[1:n_P]
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])
  
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  precMat_P[1:m_P,1:m_P]<-tau_P * MQM_P[1:m_P,1:m_P]
  # Process Model
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = mn[1:m_P], prec = precMat_P[1:m_P,1:m_P])
  # Parameter Model
  tau2   ~   dinvgamma(shape=0.2,rate=0.2) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})

########################################################################################
# Gold Standard - Count
Rep_count_model_string <- nimbleCode({
  
  
  # Data Model
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dtpois(lambda[j])
  }
  lambda[1:n_P] <- exp(W_P[1:n_P]+XB_P[1:n_P])
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  
  # Process Model
  W_O[1:n_O]  <- chol_O[1:n_O,1:n_O]%*%delta_O[1:n_O]
  W_P[1:n_P]  <- chol_P[1:n_P,1:n_P]%*%delta_P[1:n_P]
  chol_O[1:n_O,1:n_O]  <- t(chol(covMat_O[1:n_O,1:n_O]))
  chol_P[1:n_P,1:n_P]  <- t(chol(covMat_P[1:n_P,1:n_P]))
  covMat_O[1:n_O,1:n_O]<- expcov(dists = dists_O[1:n_O,1:n_O],phi = phi_O)
  covMat_P[1:n_P,1:n_P]<- expcov(dists = dists_P[1:n_P,1:n_P],phi = phi_P)
  
  for(i in 1:n_O) {
    delta_O[i] ~ dnorm(mean=0,  var=sigma2_O)  
  }
  for(i in 1:n_P) {
    delta_P[i] ~ dnorm(mean=0,  var=sigma2_P)  
  }
  # Parameter Model
  sigma2_O   ~  dinvgamma(shape=2,rate=2) 
  phi_O   ~  dunif(0.01,1.3)
  sigma2_P   ~  dinvgamma(shape=2,rate=2) 
  phi_P   ~  dunif(0.01,1.3)
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})

########################################################################################
# Gold Standard - Semi Continuous
Rep_semi_model_string <- nimbleCode({
  # Data Model
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dlnorm(meanlog=muLog[j],varlog=tau2)
  }
  
  muLog[1:n_P] <- W_P[1:n_P]+XB_P[1:n_P]
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  
  # Process Model
  W_O[1:n_O]  <- chol_O[1:n_O,1:n_O]%*%delta_O[1:n_O]
  W_P[1:n_P]  <- chol_P[1:n_P,1:n_P]%*%delta_P[1:n_P]
  chol_O[1:n_O,1:n_O]  <- t(chol(covMat_O[1:n_O,1:n_O]))
  chol_P[1:n_P,1:n_P]  <- t(chol(covMat_P[1:n_P,1:n_P]))
  covMat_O[1:n_O,1:n_O]<- expcov(dists = dists_O[1:n_O,1:n_O],phi = phi_O)
  covMat_P[1:n_P,1:n_P]<- expcov(dists = dists_P[1:n_P,1:n_P],phi = phi_P)
  
  
  for(i in 1:n_O) {
    delta_O[i] ~ dnorm(mean=0,  var=sigma2_O)  
  }
  for(i in 1:n_P) {
    delta_P[i] ~ dnorm(mean=0,  var=sigma2_P)  
  }
  # Parameter Model
  tau2   ~  dinvgamma(shape=2,rate=2) 
  sigma2_O   ~  dinvgamma(shape=2,rate=2) 
  phi_O   ~  dunif(0.01,1.3)
  sigma2_P   ~  dinvgamma(shape=2,rate=2) 
  phi_P   ~  dunif(0.01,1.3)
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
  })


####################################################
# Bisquare basis 
BS_count_model_string <- nimbleCode({
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dtpois(lambda[j])
  }
  
  # Z_P[1:n_P] ~ dTruncPoissonVector(lambda[1:n_P])# Prevalence
  lambda[1:n_P] <- exp(W_P[1:n_P]+XB_P[1:n_P])
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])
  
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  precMat_P[1:m_P,1:m_P]<-tau_P * MQM_P[1:m_P,1:m_P]
  # Process Model
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = mn[1:m_P], prec = precMat_P[1:m_P,1:m_P])
  # Parameter Model
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})

####################################################
# Bisquare basis 
BS_semi_model_string <- nimbleCode({
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dlnorm(meanlog=muLog[j],varlog=tau2)
  }
  
  muLog[1:n_P] <- W_P[1:n_P]+XB_P[1:n_P]
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])
  
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  precMat_P[1:m_P,1:m_P]<-tau_P * MQM_P[1:m_P,1:m_P]
  # Process Model
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = mn[1:m_P], prec = precMat_P[1:m_P,1:m_P])
  # Parameter Model
  tau2   ~   dinvgamma(shape=0.2,rate=0.2) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})

####################################################
# PICAR Count - Cross-correlation
Cor_PICAR_count_model_string <- nimbleCode({
  
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dtpois(lambda[j])
  }
  lambda[1:n_P] <- exp(W_P[1:n_P]+XB_P[1:n_P])
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])
  
  # Process Model
  corMat[1:m_O,1:m_P]<-rho*(1/sqrt(tau_O))*(1/sqrt(tau_O))*baseCor[1:m_O,1:m_P]
  tcorMat[1:m_P,1:m_O]<-t(corMat[1:m_O,1:m_P])
  # Correlated O
  # condCov_O[1:m_O,1:m_O]<-(1/tau_O)*MQM_O[1:m_O,1:m_O]-(tau_O)*corMat[1:m_O,1:m_P]%*%(invMQM_P[1:m_P,1:m_P]%*%tcorMat[1:m_P,1:m_O])
  # condMean_O[1:m_O]<-corMat[1:m_O,1:m_P]%*%((tau_P)*invMQM_P[1:m_P,1:m_P]%*%delta_P[1:m_P])
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  # Correlated P
  condCov_P[1:m_P,1:m_P]<-(1/tau_P)*MQM_P[1:m_P,1:m_P]-(tau_O)*tcorMat[1:m_P,1:m_O]%*%(invMQM_O[1:m_O,1:m_O]%*%corMat[1:m_O,1:m_P])
  condMean_P[1:m_P]<-tcorMat[1:m_P,1:m_O]%*%((tau_O)*invMQM_O[1:m_O,1:m_O]%*%delta_O[1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = condMean_P[1:m_P], cov = condCov_P[1:m_P,1:m_P])
  # Parameter Model
  rho ~  dunif(-1,1) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})
# PICAR Semi - Cross-correlation
Cor_PICAR_semi_model_string <- nimbleCode({
  
  Z_O[1:n_O] ~ dBernoulliVector(prob[1:n_O]) # Occurrence
  prob[1:n_O] <- 1/(1+exp(-(XB_O[1:n_O]+W_O[1:n_O])))
  
  for(j in 1:n_P){
    Z_P[j] ~ dlnorm(meanlog=muLog[j],varlog=tau2)
  }
  muLog[1:n_P] <- W_P[1:n_P]+XB_P[1:n_P]
  # XB and Random Effects
  XB_O[1:n_O]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n_P]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n_O] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n_P] <- RcalculateW_P(delta_P[1:m_P])
  
  # Process Model
  corMat[1:m_O,1:m_P]<-rho*(1/sqrt(tau_O))*(1/sqrt(tau_O))*baseCor[1:m_O,1:m_P]
  tcorMat[1:m_P,1:m_O]<-t(corMat[1:m_O,1:m_P])
  # Correlated O
  # condCov_O[1:m_O,1:m_O]<-(1/tau_O)*MQM_O[1:m_O,1:m_O]-(tau_O)*corMat[1:m_O,1:m_P]%*%(invMQM_P[1:m_P,1:m_P]%*%tcorMat[1:m_P,1:m_O])
  # condMean_O[1:m_O]<-corMat[1:m_O,1:m_P]%*%((tau_P)*invMQM_P[1:m_P,1:m_P]%*%delta_P[1:m_P])
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  # Correlated P
  condCov_P[1:m_P,1:m_P]<-(1/tau_P)*MQM_P[1:m_P,1:m_P]-(tau_O)*tcorMat[1:m_P,1:m_O]%*%(invMQM_O[1:m_O,1:m_O]%*%corMat[1:m_O,1:m_P])
  condMean_P[1:m_P]<-tcorMat[1:m_P,1:m_O]%*%((tau_O)*invMQM_O[1:m_O,1:m_O]%*%delta_O[1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = condMean_P[1:m_P], cov = condCov_P[1:m_P,1:m_P])
  # Parameter Model
  tau2   ~   dinvgamma(shape=0.2,rate=0.2) 
  rho ~  dunif(-1,1) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})

