
library("devtools")
library(roxygen2)

## Document and check ClaritySim
setwd("~/code/Clarity")
setwd("ClaritySim")
document()
check()


## Document and check Clarity
setwd("~/code/Clarity")
setwd("Clarity")
document()
check()

## TO FIX c_Bootstrap example

## Install both packages
setwd("~/code/Clarity")
install("ClaritySim")
install("Clarity")

##############################
## How the Clarity data are made with ClaritySim

q("no")
library("ClaritySim")

set.seed(1)
n=100 ; k=10
mysim=simulateCoalescent(n,k,
                         sigma0=0.0001,
                         Amodel="uniform",
                         alpha=0,
                         minedge=0.1) # Simlulate 100 objects in a 10 dimensional latent space

similarsim<-transformCoalescent(mysim)
alternatesim<-mixCoalescent(mysim)
dataraw=mysim$Y
datarep=similarsim$Y
datamix=alternatesim$Y
datarawD=mysim$D
datarepD=similarsim$D
datamixD=alternatesim$D

system("mkdir -p Clarity/data")
save(dataraw,file="Clarity/data/dataraw.RData")
save(datarep,file="Clarity/data/datarep.RData")
save(datamix,file="Clarity/data/datamix.RData")
save(datarawD,file="Clarity/data/datarawD.RData")
save(datarepD,file="Clarity/data/datarepD.RData")
save(datamixD,file="Clarity/data/datamixD.RData")

