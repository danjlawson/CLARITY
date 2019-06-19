library("devtools")
library(roxygen2)

args=commandArgs(TRUE)

run_dont_test = FALSE
if("test" %in% args) {
    print("Running full test")
    run_dont_test =TRUE
}

if("Clarity" %in% args){
    setwd("Clarity")
    d=document()
    if(!any(is.null(d))){
        print(d)
    }else{
        print("No document problems returned")
    }
    c=check(run_dont_test = run_dont_test)
    print(c)
    setwd("..")
}
if("ClaritySim" %in% args){
    setwd("ClaritySim")
    d=document()
    if(!any(is.null(d))){
        print(d)
    }else{
        print("No document problems returned")
    }
    c=check(run_dont_test = run_dont_test)
    print(c)
    setwd("..")
}
if("data" %in% args){
    library("ClaritySim")
    set.seed(1)
    # Simlulate 100 objects in a 10 dimensional latent space
    mysim=simulateCoalescent(N=100,K=10,L=500,
                             alpha=0,
                             sigma0=0.0001,Amodel="uniform",
                             minedge=0.1)
    repsim<-transformCoalescent(mysim)
    mixsim<-mixCoalescent(mysim)
    datarawA=mysim$A
    datarawD=mysim$D
    dataraw=mysim$Y
    datarep=repsim$Y
    datamix=mixsim$Y
    
    system("mkdir -p Clarity/data")
    save(datarawA,file="Clarity/data/datarawA.RData",compress="xz")
    save(datarawD,file="Clarity/data/datarawD.RData",compress="xz")
    save(dataraw,file="Clarity/data/dataraw.RData",compress="xz")
    save(datarep,file="Clarity/data/datarep.RData",compress="xz")
    save(datamix,file="Clarity/data/datamix.RData",compress="xz")
}
if("install" %in% args){
    ## Make distributable packages
    system("rm Clarity.tar.gz ClaritySim.tar.gz")
    system("tar -czvf Clarity.tar.gz Clarity")
    system("tar -czvf ClaritySim.tar.gz ClaritySim")
    install.packages("ClaritySim.tar.gz",repos = NULL, type="source")
    install.packages("Clarity.tar.gz",repos = NULL, type="source")
}
