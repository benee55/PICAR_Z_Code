# Vectorized Approach
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


dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 nonzeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(nonzeroProb))
      ## or the probability if log = FALSE
      else return(nonzeroProb * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- (1-nonzeroProb) + nonzeroProb * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    else return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), nonzeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = nonzeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, nonzeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'nonzeroProb = double()')
  )))
# Zero-inflated Poisson

PICAR_count_model_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    Z[i] ~ dZIP(lambda = lambda[i] , nonzeroProb = prob[i])
    }
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  lambda[1:n] <- exp(W_P[1:n]+XB_P[1:n])
  
  # XB and Random Effects
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n] <- RcalculateW_P(delta_P[1:m_P])
  
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
# PICAR Count - Cross-correlation
Cor_PICAR_count_model_string <- nimbleCode({
  for(i in 1:n){
    Z[i] ~ dZIP(lambda = lambda[i] , nonzeroProb = prob[i])
  }
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  lambda[1:n] <- exp(W_P[1:n]+XB_P[1:n])
  
  # XB and Random Effects
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n] <- RcalculateW_P(delta_P[1:m_P])
  
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

####################################################
# Reparameterized Spatial Random Effects (Gold Standard)
Rep_count_model_string <- nimbleCode({
  # eigenListFunction <- eigenListFunctionGenerator()
  # Data Model
  for(i in 1:n){
    Z[i] ~ dZIP(lambda = lambda[i] , nonzeroProb = prob[i])
  }
  
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  lambda[1:n] <- exp(XB_P[1:n]+W_P[1:n])
  
  # XB and Random Effects
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  
  # Intermediate
  W_O[1:n] <- chol_O[1:n,1:n]%*%delta_O[1:n]
  W_P[1:n] <- chol_P[1:n,1:n]%*%delta_P[1:n]
  chol_O[1:n,1:n]  <- t(chol(covMat_O[1:n,1:n]))
  chol_P[1:n,1:n]  <- t(chol(covMat_P[1:n,1:n]))
  covMat_O[1:n,1:n]<- expcov(dists = dists[1:n,1:n],phi = phi_O)
  covMat_P[1:n,1:n]<- expcov(dists = dists[1:n,1:n],phi = phi_P)
  
  # Process Model
  for(i in 1:n) {
    delta_O[i] ~ dnorm(mean=0,  var=sigma2_O)  # manual entry of linear predictors
    delta_P[i] ~ dnorm(mean=0,  var=sigma2_P)  # manual entry of linear predictors
  }
  
  # Parameter Model
  sigma2_O   ~  dinvgamma(shape=2,rate=2) 
  phi_O   ~  dunif(0.01,1.3)
  sigma2_P   ~  dinvgamma(shape=2,rate=2) 
  phi_P   ~  dunif(0.01,1.3)
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
  
})

####################################################
####################################################
# Zero-inflated Tobit (ZIT) Model 
####################################################
####################################################

dTobit <- nimbleFunction(
  run = function(x = integer(), 
                 mu = double(), 
                 nonzeroProb = double(), 
                 tau2=double(),
                 log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(log(nonzeroProb) +  dnorm(x, mu, sd=sqrt(tau2),log = TRUE) )
      ## or the probability if log = FALSE
      else return(nonzeroProb * dnorm(x, mu, sd=sqrt(tau2),log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- (1-nonzeroProb) + nonzeroProb * pnorm(0, mean= mu, sd=sqrt(tau2),log = FALSE)
    if (log) return(log(totalProbZero))
    else return(totalProbZero)
  })

rTobit <- nimbleFunction(
  run = function(n = integer(), mu = double(), nonzeroProb = double() , tau2=double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = nonzeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rnorm(1, mu, sd=sqrt(tau2)))
  })

registerDistributions(list(
  dTobit = list(
    BUGSdist = "dTobit(mu, nonzeroProb, tau2)",
    discrete = FALSE,
    range = c(0, Inf),
    types = c('value = integer()', 'mu = double()', 'nonzeroProb = double()', 'tau2 = double()')
  )))


# Nimble Model
PICAR_semi_model_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    Z[i] ~ dTobit(mu = mu[i] , nonzeroProb = prob[i], tau2=tau2)
  }
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  mu[1:n] <- W_P[1:n]+XB_P[1:n]
  
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n] <- RcalculateW_P(delta_P[1:m_P])
  
  precMat_O[1:m_O,1:m_O]<-tau_O * MQM_O[1:m_O,1:m_O]
  precMat_P[1:m_P,1:m_P]<-tau_P * MQM_P[1:m_P,1:m_P]
  # Process Model
  delta_O[1:m_O] ~ dmnorm(mean = mn[1:m_O], prec = precMat_O[1:m_O,1:m_O])
  delta_P[1:m_P] ~ dmnorm(mean = mn[1:m_P], prec = precMat_P[1:m_P,1:m_P])
  
  # Parameter Model
  tau2   ~  dinvgamma(shape=2,rate=2) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})


Cor_PICAR_semi_model_string <- nimbleCode({
  for(i in 1:n){
    Z[i] ~ dTobit(mu = mu[i] , nonzeroProb = prob[i], tau2=tau2)
  }
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  mu[1:n] <- W_P[1:n]+XB_P[1:n]
  
  # XB and Random Effects
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  W_O[1:n] <- RcalculateW_O(delta_O[1:m_O])
  W_P[1:n] <- RcalculateW_P(delta_P[1:m_P])
  
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
  tau2   ~  dinvgamma(shape=2,rate=2) 
  tau_O   ~  dgamma(shape=0.5,rate=2000) 
  tau_P   ~  dgamma(shape=0.5,rate=2000) 
  beta_O[1:p_O] ~ dRdmnormB_O(mean = mn[1:p_O])
  beta_P[1:p_P] ~ dRdmnormB_P(mean = mn[1:p_P])
})


# Semi-continuous Model 

Rep_semi_model_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    Z[i] ~ dTobit(mu = mu[i] , nonzeroProb = prob[i], tau2=tau2)
  }
  prob[1:n] <- 1/(1+exp(-(XB_O[1:n]+W_O[1:n])))
  mu[1:n] <- W_P[1:n]+XB_P[1:n]
  
  # XB and Random Effects
  XB_O[1:n]<-RcalculateXB_O(beta_O[1:p_O])
  XB_P[1:n]<-RcalculateXB_P(beta_P[1:p_P])
  
  # Intermediate
  W_O[1:n] <- chol_O[1:n,1:n]%*%delta_O[1:n]
  W_P[1:n] <- chol_P[1:n,1:n]%*%delta_P[1:n]
  chol_O[1:n,1:n]  <- t(chol(covMat_O[1:n,1:n]))
  chol_P[1:n,1:n]  <- t(chol(covMat_P[1:n,1:n]))
  covMat_O[1:n,1:n]<- expcov(dists = dists[1:n,1:n],phi = phi_O)
  covMat_P[1:n,1:n]<- expcov(dists = dists[1:n,1:n],phi = phi_P)
  
  # Process Model
  for(i in 1:n) {
    delta_O[i] ~ dnorm(mean=0,  var=sigma2_O)  # manual entry of linear predictors
    delta_P[i] ~ dnorm(mean=0,  var=sigma2_P)  # manual entry of linear predictors
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