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

## Make distributable packages
system("rm Clarity.tar.gz ClaritySim.tar.gz")
system("tar -czvf Clarity.tar.gz Clarity")
system("tar -czvf ClaritySim.tar.gz ClaritySim")
install.packages("ClaritySim.tar.gz",repos = NULL, type="source")
install.packages("Clarity.tar.gz",repos = NULL, type="source")

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

library("Clarity")
library("ClaritySim")
set.seed(1)
     scanraw=Clarity_Scan(dataraw) # Generate an initial run
     scanraw=Clarity_Scan(dataraw,clist=scanraw) # Run it for longer
     # Apply it to a new dataset with the same structure
     scanrepfromraw=Clarity_Predict(datarep,clist=scanraw) 
     # Apply it to a new dataset with slightly different structure
     scanmixfromraw=Clarity_Predict(datamix,clist=scanraw)

## Plot the case where the residuals shouldn't matter
Clarity_ComparisonPlot(scanraw$scan[[19]],scanrepfromraw$scan[[19]],
           name2="Predicted (no structural change)",zlimresidual=c(-1,1))


Clarity_ComparisonPlot(scanraw$scan[[19]],scanrepfromraw$scan[[19]],
                       name2="Predicted (no structural change)",zlimresidual=c(-1,1),
                       order=NA)

## Plot the case where the residuals should matter
Clarity_ComparisonPlot(scanraw$scan[[19]],scanmixfromraw$scan[[19]],
           name2="Predicted (one structural change)",zlimresidual=c(-1,1))

set.seed(1)
scanraw=Clarity_Scan(dataraw)
#scanraw=Clarity_Scan(dataraw,niter=20,objectivedistmin=1e-5,clist=scanraw)

scanrepfromraw=Clarity_Predict(datarep,scanraw)
scanmixfromraw=Clarity_Predict(datamix,clist=scanraw)

scanplot=Clarity_ComparisonPlot(scanraw[[1]],scanrepfromraw[[1]])
zlim=c(-1,1)
scanplot=Clarity_ComparisonPlot(scanraw[[19]],scanrepfromraw[[19]],cex.axis=0,zlimresidual=zlim)
scanplot=Clarity_ComparisonPlot(scanraw[[19]],scanmixfromraw[[19]],cex.axis=0,zlimresidual=zlim)


##################
## How to handle a second matrix
## Y = l^2 A X A^T + (1-l)^2 B W B^T
## where l (lambda) is a scalar multiplying "rows" of A and B so that they together sum to 1.

c_objectivefunction=function(A,X,Y,lambda=1,B=NULL,W=NULL){
    Yhat=0
    if(lambda>0) Yhat = Yhat + lambda^2 * A %*% X %*% t(A)
    if(lambda<1) Yhat = Yhat + (1-lambda)^2 * B %*% W %*% t(B) 
    sum((Y - Yhat)^2)
}

## Learning X
## Set Y' = (1/l^2) Y - ((1-l)^2/l^2) B W B^T
## Then learning X proceeds as before.
## Learning W is completely symmetric to learning X; just swap arguments.
c_updateX=function(A,X,Y,lambda=1,B=NULL,W=NULL,fixed=NULL,norm=TRUE) {
    if(lambda==1){
        Yprime=Y
    }else{
        Yprime=(1/lambda^2) * Y - ((1-lambda)^2/lambda^2) * B %*% W %*% t(B)
    }
    Xsol=try(c_updateX_solve(A,Yprime),silent=TRUE)
    if(class(Xsol)=="try-error"){
        tAA=(t(A) %*% A) 
        num = t(A) %*% (Yprime) %*% A
        denom = tAA %*% X %*% tAA
        Xsol=X * (num / denom)
    }
    Xsol
}

## Learning A
## Set Y' = (1/l^2) Y - ((1-l)^2/l^2) B W B^T
## Then learning A proceeds as before.
## Learning B is completely symmetric to learning X; just swap arguments.
c_updateA<-function(A,X,Y,lambda=1,B=NULL,W=NULL,fixed=NULL,norm=TRUE) {
    ## This is probably the most efficient way to implement this, based on the choices of ordering available
    if(lambda==1){
        Yprime=Y
    }else{
        Yprime=(1/lambda^2) * Y - ((1-lambda)^2/lambda^2) * B %*% W %*% t(B)
    }
    tAA=(t(A) %*% A) 
    num = t(Yprime) %*% A %*% X + Yprime %*% A %*% t(X)
    denom = A %*% X %*% tAA %*% t(X) + A %*% t(X) %*% tAA %*% X
    ret=A * (num / denom)
    if(norm)ret=ret/rowSums(ret)
    ret
}

## Learning lambda
## f = tr {(Y - l^2 (AXA^T) - (1-l)^2 B W B^T )(...)^T }
##   = a0 - l^2 a1 - (1-l^2) a2 + l ^4 b1 + l^2(1-l^2) b2 + (1-l)^4 b3
##   = [a0-a2+b3] + l [ 2 a2 - 4 b3 ] + l^2 [-a1 -a2 + b2 + 6 b3 ] + l^3 [ -2 b2 - 4 b3 ] + l^4 [b1 +b2 + b3]
## with:
## a0 = tr{ Y Y^T }
## a1 = tr{ Y A [X + X^T] A^T }
## a2 = tr{ Y B [W + W^T] B^T }
## b1 = tr{ A X A^T A X^T A^T }
## b2 = tr{ A X A^T A X^T A^T + B W B^T B W^T B^T }
## b3 = tr{ B W B^T B W^T B^T }
## So we are seeking places where the derivative of f wrt lambda = 0, i.e. maxima and minima.
## dfdl = [ 2a2 - 4b3] + 2 l [-a1 -a2 +b2 + 6b3 ] + 3 l^2 [-2 b2 - 4 b3 ] + 4 l^3 [ b1 + b2 + b3 ] = 0
## This can be solved with the R function polyroot(c(C0,C1,C2,C3))
## where:
## C0 = [ 2a2 - 4b3 ]
## C1 = 2 [-a1 -aa2 +b2 + 6b3 ]
## C2 = 3 [-2 b2 - 4 b3 ]
## C3 = 4 [ b1 + b2 + b3 ]
## We are interested in the lambda that minimises f on the real line, (0,1)
## Or if no solutions lie here, argmin(f(0),f(1))
## (probably wise to insist on a bound epsilon<lambda<1-epsilon for epsilon small.)

## NB would be very easy to do this numerically no?
## By precomputing Y, AXA^T, BWB^T, all crosses and all traces.
## Then
tr<-function(x)sum(diag(x))

c_dolambda_numericalscan=function(A,X,Y,lambda=1,B=NULL,W=NULL,fixed=NULL,norm=TRUE) {
    trYY=tr(Y %*% t(Y))
    trYA=tr(Y %*% A %*% t(X) %*% t(A))
    tr
}

##
A=scanraw[[19]]$A
Y=scanraw[[19]]$Y
X=scanraw[[19]]$X
lambda=0.99
Yprime=(1/lambda^2) * Y - ((1-lambda)^2/lambda^2) * A %*% X %*% t(A)
c_test=Clarity_fixedK(Yprime,2)
B=c_test$A
W=c_test$X

c_objectivefunction(A,X,Y)
c_objectivefunction(A,X,Y,lambda,B,W)

Clarity_ComparisonPlot(scanraw[[19]],scanmixfromraw[[19]],
           name2="Predicted (one structural change)",zlimresidual=c(-1,1))

## Some of this is from examples written into:
## /Users/madjl/code/afs/testmm.R
