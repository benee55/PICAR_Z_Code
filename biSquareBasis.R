#Code to generate data
rm(list=ls())
################################################################
# Create Bisquare Basis Functions
################################################################
# Initialize
setwd("~/Dropbox/PICAR_ZIP/code/realExample/bivalve/") # Setup work directory
# setwd('/home/slee287/ZeroInflatedPICAR/code/')
source(file = "../../source/sharedFunctions.R")
load("samples/data.RData")

# Spatial Basis Functions
# Split domain into 4 partitions
# Resolution 1

pointA<-c(-0.025,0.65)# Topleft
pointC<-c(-0.025,-0.1) # bottom Left
pointB<-c(1.025,1.4) # Top Right
pointD<-c(1.025,0.65) # Bottom Right

# Begin Basis Functions
# 9 lines with 9 points
# Top Line: A to B
topCoords<-cbind(seq(pointA[1],pointB[1], length.out = 17),
                 seq(pointA[2],pointB[2], length.out = 17))
# Bottom Line: C to D
bottomCoords<-cbind(seq(pointC[1],pointD[1], length.out = 17),
                    seq(pointC[2],pointD[2], length.out = 17))
tempCoords<-list()
for(indexVal in 1:17){
  tempCoords[[indexVal]]<-cbind(seq(topCoords[indexVal,1],bottomCoords[indexVal,1], length.out = 17),
                                seq(topCoords[indexVal,2],bottomCoords[indexVal,2], length.out = 17))
}

# LEft Line: A to C
leftCoords<-cbind(seq(pointA[1],pointC[1], length.out = 17),
                  seq(pointA[2],pointC[2], length.out = 17))
# Right Line: B to D
rightCoords<-cbind(seq(pointB[1],pointD[1], length.out = 17),
                   seq(pointB[2],pointD[2], length.out = 17))
newCoords<-list()
for(indexVal in 1:17){
  newCoords[[indexVal]]<-cbind(seq(leftCoords[indexVal,1],rightCoords[indexVal,1], length.out = 17),
                               seq(leftCoords[indexVal,2],rightCoords[indexVal,2], length.out = 17))
}
# Resolution 1: 4 points
knotsA<-rbind(newCoords[[5]][c(5,13),],
              newCoords[[13]][c(5,13),])


# Resolution 2: 16 points
knotsB<-rbind(newCoords[[3]][c(3,7, 11, 15),],
              newCoords[[7]][c(3,7, 11, 15),],
              newCoords[[11]][c(3,7, 11, 15),],
              newCoords[[15]][c(3,7, 11, 15),])

# Resolution 3: 64 points
foo<-2*(1:8)
knotsC<-rbind(newCoords[[2]][foo,],
              newCoords[[4]][foo,],
              newCoords[[6]][foo,],
              newCoords[[8]][foo,],
              newCoords[[10]][foo,],
              newCoords[[12]][foo,],
              newCoords[[14]][foo,],
              newCoords[[16]][foo,])


# Full Knots
fullKnots<-rbind(knotsA,knotsB, knotsC)
# Visualization of Knots
par(mfrow=c(1,1))
plot(gridLocation, pch=16, col="black" , cex=0.25, 
     ylim=c(-0.1, 1.5), xlim=c(-0.1, 1.1))
lines(x=c(pointB[1],pointD[1]) , y=c(pointB[2],pointD[2]))
lines(x=c(pointA[1],pointB[1]) , y=c(pointA[2],pointB[2]))
lines(x=c(pointA[1],pointC[1]) , y=c(pointA[2],pointC[2]))
lines(x=c(pointD[1],pointC[1]) , y=c(pointD[2],pointC[2]))
points(knotsC,col="blue", pch=16, cex=1)
for(j in c(1,3,5,7,9,11,13,15,17)){
  lines(newCoords[[j]], lwd=2, col="blue")
  lines(tempCoords[[j]], lwd=2, col="blue")
}

points(knotsB,col="black", pch=16, cex=1.5)
for(j in c(1,5,9,13,17)){
  lines(newCoords[[j]], lwd=2, col="black")
  lines(tempCoords[[j]], lwd=2, col="black")
}
points(knotsA,col="red", pch=16, cex=2)
for(j in c(1,9,17)){
  lines(newCoords[[j]], lwd=2, col="red")
  lines(tempCoords[[j]], lwd=2, col="red")
}

# Basis Functions - Bandwidth
epsA<-min(as.numeric(dist(knotsA)))*1.5
epsB<-min(as.numeric(dist(knotsB)))*1.5
epsC<-min(as.numeric(dist(knotsC)))*1.5

# Basis Functions - Bisquare Functions
bisquareFunc<-function(dist, eps){ # Function to create Bisquare basis function
  (dist<eps)*((1-(dist/eps)^2)^2)
}


save( knotsA , knotsB,knotsC  , epsA , epsB ,epsC, bisquareFunc , 
      file="samples/BisquareBasisFunctions.RData")

