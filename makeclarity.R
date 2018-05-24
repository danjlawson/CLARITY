q("no")

library("devtools")
library(roxygen2)
setwd("/Users/madjl/code/Clarity")
#create("Clarity")
#create("ClaritySim")

setwd("Clarity")
document()
check()

setwd("..")

setwd("ClaritySim")
document()
check()


setwd("/Users/madjl/code/Clarity")
install("ClaritySim")
install("Clarity")

q("no")
library("Clarity")
library("ClaritySim")

##############################
## How the Clarity data are made

library("ClaritySim")
set.seed(1)
mysim<-simulateCoalescent(100,10) # Simlulate 100 objects in a 10 dimensional latent space
myroworder=order(apply(mysim$A,1,which.max))
mysim$A=mysim$A[myroworder,]
mysim$Y=mysim$Y[myroworder,myroworder]
mysim$Y0=mysim$Y0[myroworder,myroworder]

similarsim<-transformCoalescent(mysim)
alternatesim<-mixCoalescent(mysim)
dataraw=mysim$Y
datarep=similarsim$Y
datamix=alternatesim$Y

system("mkdir -p Clarity/data")
save(dataraw,file="Clarity/data/dataraw.RData")
save(datarep,file="Clarity/data/datarep.RData")
save(datamix,file="Clarity/data/datamix.RData")

##############################
## Basic Clarity usage
set.seed(1)
scanraw=Clarity_Scan(dataraw)
#scanraw=Clarity_Scan(dataraw,niter=20,objectivedistmin=1e-5,clist=scanraw)

scanrepfromraw=Clarity_Predict(datarep,scanraw)
scanmixfromraw=Clarity_Predict(datamix,clist=scanraw)

scanplot=Clarity_ComparisonPlot(scanraw[[1]],scanrepfromraw[[1]])
zlim=c(-1,1)
scanplot=Clarity_ComparisonPlot(scanraw[[19]],scanrepfromraw[[19]],cex.axis=0,zlimresidual=zlim)
scanplot=Clarity_ComparisonPlot(scanraw[[19]],scanmixfromraw[[19]],cex.axis=0,zlimresidual=zlim)


## From /Users/madjl/code/afs/testmm.R
