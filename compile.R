
args=commandArgs(TRUE)
if(length(args)==0 | args[1]=="-h" | args[1]=="--help"){
    cat("Usage: Rscript compile.R\n")
    stop("Must provide one of Clarity, ClaritySim, Install")
}

checkclarity=FALSE
checkclaritysim=FALSE
install=FALSE
if("clarity"%in%tolower(args))checkclarity=TRUE
if("claritysim"%in%tolower(args))checkclaritysim=TRUE
if("install"%in%tolower(args))install=TRUE
if("all"%in%tolower(args)){checkclarity=checkclaritysim=install=TRUE}

library("devtools")
library(roxygen2)
setwd("/Users/madjl/code/Clarity")
#create("Clarity")
#create("ClaritySim")

if(checkclarity){
    setwd("Clarity")
    document()
    check()
    setwd("..")
}

if(checkclaritysim){
    setwd("ClaritySim")
    document()
    check()
    setwd("..")
}


## Make distributable packages
if(install){
    system("rm Clarity.tar.gz ClaritySim.tar.gz")
    system("tar -czvf Clarity.tar.gz Clarity")
    system("tar -czvf ClaritySim.tar.gz ClaritySim")
    install.packages("ClaritySim.tar.gz",repos = NULL, type="source")
    install.packages("Clarity.tar.gz",repos = NULL, type="source")
}
