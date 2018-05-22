library("devtools")
library(roxygen2)
setwd("/Users/madjl/code/SimSim")
#create("SimSim")

setwd("SimSim")
document()
check()

setwd("..")
q()
install("SimSim")

library("SimSim")
