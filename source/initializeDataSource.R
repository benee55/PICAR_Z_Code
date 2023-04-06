

####################################################################################################

initializeData.hurdle<-function(n, # Sample Size for Model-fitting samples
                         sigma2.o, # Partial sill- Occurrence 
                         phi.o, # Range- Occurrence 
                         beta.o, # Regression Coefficent - Occurrence 
                         tau2.o=NULL, # Nugget Variance
                         sigma2.p, # Partial sill- Prevalence 
                         phi.p, # Range- Prevalence 
                         beta.p, # Regression Coefficent - Prevalence 
                         tau2.p, # Nugget Variance
                         rho, # Cross-correlation between Occurrence and Prevalence processes
                         seed, # Random Seed
                         cvPct=0.4, # Proportion of n for validation sample size
                         covFunc, # Put in actual covariance function (sqeCov, expCov, matCov)
                         xBounds=c(-1,1), # Bounds for covariates
                         model, # Count or semi-continuous: "count" or "semi"
                         type, # Hurdle or mixture mode: "hurdle" or "mixture"
                         covStructure, # Naming convention for Covariance function
                         uniformSpace=FALSE, # Sample from a uniform grid
                         maxGrid, # For uniform grid
                         pdf=FALSE # Pdf of data
){
  origN<-n
  n<-ceiling(n*(1+cvPct)) # Cross-Validation n 
  set.seed(seed)
  if(model=="count"){parTruth<-c(phi.o , sigma2.o , beta.o , tau2.o , 
                               phi.p , sigma2.p , beta.p , tau2.p, rho)}
  if(model=="semi"){parTruth<-c(phi.o , sigma2.o , beta.o , tau2.o , 
                                     phi.p , sigma2.p , beta.p , tau2.p, rho)}
  ######################################################################
  if(uniformSpace==TRUE){
    y<-x<-seq(0,maxGrid,length.out = floor(sqrt(n))) #0-1 Grid
    grid.location<-expand.grid(x,y) # Equally spaced points on grid
  }else{
    grid.location<-cbind(runif(n,0,maxGrid),runif(n,0,maxGrid))
  }
  # Covariates
  X<-cbind(runif(n,xBounds[1],xBounds[2]),
           runif(n,xBounds[1],xBounds[2])) 
  ########################################################################################
  # Process 1: Occurrence
  distS <- as.matrix(dist(grid.location))
  if(is.null(tau2.o)){
    SigmaS = sigma2.o*covFunc(distMat=distS,phi=phi.o)  
  }else{
    SigmaS = sigma2.o*covFunc(distMat=distS,phi=phi.o)  + tau2.o*diag(nrow(distS))
  }
  ########################################################################################
  # Process 2: Prevalence
  distW <- as.matrix(dist(grid.location))
  if(is.null(tau2.p)){
    SigmaW = sigma2.p*covFunc(distMat=distW,phi=phi.p)  
  }else{
    SigmaW = sigma2.p*covFunc(distMat=distW,phi=phi.p)
  }
  
  LS<-t(chol(SigmaS))
  LW<-t(chol(SigmaW))
  rhoLStLW<-rho*LS%*%t(LW)
  SigmaCombined<-rbind(cbind(SigmaS,rhoLStLW),
                       cbind(t(rhoLStLW),SigmaW))
  
  ########################################################################################
  gpAll<-t(chol(SigmaCombined))%*%rnorm(2*n)
  gpS<-gpAll[1:n,]
  meanS<-X%*%c(beta.o)
  Ax<-meanS+gpS
  pX<-exp(Ax)/(1+exp(Ax))
  Odat<-apply(pX,1,rbinom,n=1,size=1)
  NonZeroInd<-which(Odat!=0)
  ZeroInd<-which(Odat==0)
  ########################################################################################
  
  # Conditional Mean
  gpW.orig <-gpAll[(n+1):(2*n)]
  gpW <- gpW.orig[NonZeroInd]
  meanW<-X[NonZeroInd,]%*%c(beta.p)
  n1<-length(NonZeroInd)
  if(model=="count"){
    Bx<-meanW+gpW
    if(type=="hurdle"){Zdat<-rztpois(n = n1,lambda=exp(Bx))}
    if(type=="mixture"){Zdat<-rpois(n = n1,lambda=exp(Bx))}
    }
  if(model=="semi"){
    Bx<-meanW+gpW
    if(type=="hurdle"){Zdat<-exp(rnorm(n1,mean=Bx,sd=sqrt(tau2.p)))}
    if(type=="mixture"){Zdat<-rnorm(n1,mean = Bx,sd = sqrt(tau2.p)) ; Zdat<-ifelse(Zdat<0,0,Zdat)}
  }
  
  # ObsDat
  finalObs<-rep(0,n)
  for(h in 1:n1){
    finalObs[NonZeroInd[h]]<-Zdat[h]
  }
  
  #####################################################################
  keepIndex<-sample(1:n,origN)
  cvIndex<-(1:n)[!(1:n)%in%keepIndex]
  #####################################################################
  #Model Fitting Sample
  XMat.mod=X[keepIndex,]
  gpS.mod=gpS[keepIndex]
  gpW.mod=gpW.orig[keepIndex]
  obs.mod=finalObs[keepIndex]
  grid.location.mod=grid.location[keepIndex,]
  #####################################################################
  #Cross-Validation Sample
  XMat.cv=X[cvIndex,]
  gpS.cv=gpS[cvIndex]
  gpW.cv=gpW.orig[cvIndex]
  obs.cv=finalObs[cvIndex]
  grid.location.cv=grid.location[cvIndex,]
  ######################################################################################
  if(pdf){pdf(file=paste("Map_",model,"_",covStructure,"_",type,"_n",n,"_seed",seed,".pdf",sep=""))}
  
  visDat1<-list(gpS,meanS,Ax,Odat,gpW,meanW,Bx,Zdat)
  rangeY<-c(10,10,10,3,10,10,10,min(10,length(unique(Zdat))))
  mainTitle<-c("Occur: GP","Occur: Mean",
               "Occur: Mean + GP","Occur: Observation",
               "Prev: GP", "Prev: Mean",
               "Prev: Mean + GP","Prev: Observation")
  par(mfrow=c(2,4))
  for(i in 1:length(visDat1)){
    
    sim=visDat1[[i]]
    
    if(i==4){
      plot(x=grid.location[,1],y=grid.location[,2],pch=16,col=(Odat),cex=1.5,main=mainTitle[i],
           xlab="X", ylab="Y")
    }
    else if(i %in% 5:8){
      breaks <- seq(range(sim)[1],range(sim)[2],length.out=rangeY[i])
      pal <- tim.colors(length(breaks)-1)
      fb <- classIntervals(sim, n = length(pal), 
                           style = "fixed", fixedBreaks = breaks)
      col <- findColours(fb, pal)
      plot(x=grid.location[NonZeroInd,1],y=grid.location[NonZeroInd,2],col=col,pch=16,cex=1.5,
           main=mainTitle[i], xlab="X", ylab="Y")

    }
    
    else{breaks <- seq(range(sim)[1],range(sim)[2],length.out=rangeY[i])
    pal <- tim.colors(length(breaks)-1)
    fb <- classIntervals(sim, n = length(pal), 
                         style = "fixed", fixedBreaks = breaks)
    col <- findColours(fb, pal)
    plot(x=grid.location[,1],y=grid.location[,2],col=col,pch=16,cex=1.5,main=mainTitle[i],
        xlab="X", ylab="Y")
    }
    
  }
  
  dev.off()  
  
  return(list(parTruth=parTruth,
              seed=seed,
              XMat.mod=XMat.mod,
              gpS.mod=gpS.mod,
              gpW.mod=gpW.mod,
              obs.mod=obs.mod,
              grid.location.mod=grid.location.mod,
              XMat.cv=XMat.cv,
              gpS.cv=gpS.cv,
              gpW.cv=gpW.cv,
              obs.cv=obs.cv,
              grid.location.cv=grid.location.cv))
}

