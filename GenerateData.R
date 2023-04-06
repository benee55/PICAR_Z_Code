rm(list=ls())
library(nimble)
setwd("~/Dropbox/PICAR_ZIP/code/realExample/bivalve/") # Setup work directory
source(file="../../source/batchmeans.R")
source(file = "../../source/sharedFunctions.R")
dat<-as.matrix(read.csv("samples/datsc.csv")[,-c(1,5)])
coords.all<-dat[,c("x","y")]
range.coords<-apply(coords.all,2,range)
coords.all<-coords.all-matrix(abs(range.coords[1,]),nrow=nrow(coords.all),ncol=2,byrow = TRUE)
range.coords<-apply(coords.all,2,range)
gridLocation<-as.matrix(coords.all*matrix(1/abs(range.coords[2,]),nrow=nrow(coords.all),ncol=2,byrow = TRUE))
apply(gridLocation,2,range)

# Split Training and Test sample
set.seed(12345)
totN<-nrow(gridLocation)
modInd<-sample(1:totN,floor(totN*0.8))
cvInd<-(1:totN)[-modInd]

z<-dat[,"macoma"]
X<-dat[,c("mgs","silt","depth")]
zMod<-dat[modInd,"macoma"]
zCV<-dat[cvInd,"macoma"]
XMod<-X[modInd,]
XCV<-X[cvInd,]
gridLocationMod<-gridLocation[modInd,]
gridLocationCV<-gridLocation[cvInd,]

save(z,zMod,zCV,X,XMod,XCV,gridLocationMod,gridLocationCV,gridLocation,
     file="samples/data.RData")
par(mfrow=c(1,2))
posInd<-which(z!=0&z<30)
# Plot Locations with values and Zeros
plotRF(dat=z[posInd],rangeDat = z[posInd], label ="Obs",location = gridLocation[posInd,],length.out = 10,cex=0.5)
points(x=gridLocation[-posInd,1],y=gridLocation[-posInd,2],col="gray",pch=16,cex=0.2)

plot(x=gridLocation[posInd,1],y=gridLocation[posInd,2],col="red",pch=16,cex=0.3, main="Zeros vs. Non-Zeros")
points(x=gridLocation[-posInd,1],y=gridLocation[-posInd,2],col="black",pch=16,cex=0.3)


library(INLA)
  # Make Mesh for INLA
# rdomain <- inla.nonconvex.hull(gridLocationMod, # More wiggle room
#                                convex = -0.1,
#                                concave = -0.2,
#                                resolution = 100)
# mesh <- inla.mesh.2d(boundary = rdomain,
#                      max.edge = c(0.0175, 0.1), cutoff = 0.005,
#                      offset=c(-0.2, -0.01))

rdomain <- inla.nonconvex.hull(gridLocationMod,
                               convex = -0.05,
                               concave = 0.05,
                               resolution = 100)
plot(rdomain$loc,main="")
points(x=gridLocationMod[,1], y=gridLocationMod[,2],col="blue",pch=16,cex=0.3)
points(x=gridLocationCV[,1], y=gridLocationCV[,2],col="red",pch=18,cex=0.3)
mtext("Mesh",cex=2)

mesh <- inla.mesh.2d(loc=gridLocationMod,
                       boundary = rdomain,
                       # max.edge = c(0.0175, 0.1), cutoff = 0.005,
                     max.edge = c(0.0175, 0.3), cutoff = 0.005,
                     offset=c(-0.2, -0.02))
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(mesh,main="")
points(x=gridLocationMod[,1], y=gridLocationMod[,2],col="blue",pch=16,cex=0.3)
points(x=gridLocationCV[,1], y=gridLocationCV[,2],col="red",pch=18,cex=0.3)
mtext("Mesh",cex=2)


AMat <- inla.spde.make.A(mesh, loc=gridLocationMod) ; dim(AMat) 
AMatCV <- inla.spde.make.A(mesh, loc=gridLocationCV)

# prmesh2 <- inla.mesh.2d(boundary = rdomain,
#                         max.edge = c(0.45, 1), cutoff = 0.2)


  # Created Mesh
wp<-ncol(AMat)
DMat<-diag(apply(mesh$graph$vv,1,sum))
WeightMat<-mesh$graph$vv
PrecMat<-DMat-WeightMat  
Nnew<-nrow(WeightMat)
OrthSpace<-diag(Nnew)-(rep(1,Nnew)%*%t(rep(1,Nnew)))/Nnew
GriffithOperator<-OrthSpace%*%(WeightMat%*%OrthSpace)
GriffithOperatorEig<-eigen(GriffithOperator)
save(mesh,AMat,AMatCV,DMat,GriffithOperatorEig,PrecMat,
     file="samples/mesh.RData")


# par(mfrow=c(5,5), mar=c(2,2,2,2))
# for(k in 1:25){
#   plotRF(location = mesh$loc, dat = GriffithOperatorEig$vectors[,k],
#          rangeDat = GriffithOperatorEig$vectors[,k],label = paste("Eig",k),
#          length.out = 100)
#   }
