library("devtools")
library(roxygen2)

args=commandArgs(TRUE)

if("Clarity" %in% args){
    setwd("Clarity")
    d=document()
    print(d)
    c=check()
    print(c)
    setwd("..")
}
if("ClaritySim" %in% args){
    setwd("ClaritySim")
    d=document()
    print(d)
    c=check()
    print(c)
    setwd("..")
}
if("data" %in% args){
    library("ClaritySim")
    set.seed(1)
    n=100 ; k=10
    mysim=simulateCoalescent(n,k,sigma0=0.0001,Amodel="uniform",alpha=0,
                             minedge=0.1) # Simlulate 100 objects in a 10 dimensional latent space
    similarsim<-transformCoalescent(mysim)
    alternatesim<-mixCoalescent(mysim)
    dataraw=mysim$Y
    datarep=similarsim$Y
    datamix=alternatesim$Y
    
    system("mkdir -p Clarity/data")
    save(dataraw,file="Clarity/data/dataraw.RData")
    save(datarep,file="Clarity/data/datarep.RData")
    save(datamix,file="Clarity/data/datamix.RData")
}
if("install" %in% args){
    ## Make distributable packages
    system("rm Clarity.tar.gz ClaritySim.tar.gz")
    system("tar -czvf Clarity.tar.gz Clarity")
    system("tar -czvf ClaritySim.tar.gz ClaritySim")
    install.packages("ClaritySim.tar.gz",repos = NULL, type="source")
    install.packages("Clarity.tar.gz",repos = NULL, type="source")
}
