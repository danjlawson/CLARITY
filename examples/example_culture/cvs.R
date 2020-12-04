## Comparing Cultural Values to World Bank economic data
##
## Author: Daniel Lawson (dan.lawson@bristol.ac.uk)
## Date: 2020
## Licence: GPLv3 


library("Clarity")
### Read the CVs data
cvs=read.csv("Cultural_Values.csv.gz",row.names=1)
## Get the Y for culture
cvsdist=c_dist(cvs)

### Read the WB data
wb=as.matrix(read.csv("World_Bank_Unnormalised.csv.gz",row.names=1))
feature_sanitise=function(wb){
    ## standardize, and remove outliers
    tmean=apply(wb,2,mean)
    wb=t(t(wb)-tmean)
    tsd=apply(wb,2,sd)
    wb=t(t(wb)/tsd)
    wb=wb[,tsd<100]
    wb[wb< -10] = -10
    wb[wb> 10] = 10
    wb
}
wb=feature_sanitise(wb)

## Get the Y for world bank economics
wbdist=c_dist(wb)

#######################
## Save the results
matrixAsGraph=function(smatrix){
    resmat=matrix(0,nrow=dim(smatrix)[1]*(dim(smatrix)[1])/1/2,ncol=3)
    ii=1
    for(i in 1:(dim(smatrix)[1]-1))
        for(j in (i+1):dim(smatrix)[1]){
            resmat[ii,]=c(i,j,smatrix[i,j])
            ii=ii+1
        }
    resmat
}

system("mkdir -p cultureeconomics")
write.table(matrixAsGraph(wbdist),
            file="cultureeconomics/economics_distance.csv",sep=",",
            row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrixAsGraph(cvsdist),
            file="cultureeconomics/culture_distance.csv",sep=",",
            row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(wb,
            file="cultureeconomics/economics_features.csv",sep=",",
            row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(cvs,
            file="cultureeconomics/culture_features.csv",sep=",",
            row.names=FALSE,col.names=TRUE,quote=FALSE)
labeldf=data.frame(index=1:dim(cvs)[1],
                   country=rownames(cvs))
write.table(labeldf,"cultureeconomics/labels.csv",sep=",",
            row.names=FALSE,col.names=TRUE,quote=FALSE)

######################
## Running clarity
wbscan = Clarity_Scan(wbdist,105)
cvsscan = Clarity_Scan(cvsdist,105)

pred_wb_from_cvs = Clarity_Predict(wbdist,cvsscan)
pred_cvs_from_wb = Clarity_Predict(cvsdist,wbscan)

## Pvalues for the cultural values data
wb_from_cvs=Clarity_Compare(cvsscan,D=cvs,Dnew=wb,nbs=200,H0="structure") 
cvs_from_wb=Clarity_Compare(wbscan, D=wb,Dnew=cvs,nbs=200,H0="structure") 

## Pvalues for the cultural values data
scale_wb_from_cvs=Clarity_Compare(cvsscan,D=cvs,Dnew=wb,nbs=200,H0="scale") 
scale_cvs_from_wb=Clarity_Compare(wbscan, D=wb,Dnew=cvs,nbs=200,H0="scale") 

######################
## Pvalues for the residuals
myk=15
R15_wb_from_cvs=Clarity_Compare(cvsscan,k=myk,D=cvs,Dnew=wb,nbs=200,H0="scale") 
R15_cvs_from_wb=Clarity_Compare(wbscan,k=myk,D=wb,Dnew=cvs,nbs=200,H0="scale") 

######################
## plotting ##########
## Some color palettes
colwr=colorRampPalette(c("white","yellow","orange","red","darkred","black"))(100)
colbwg=colorRampPalette(c("darkblue","blue","white","green","darkgreen"))(100)
 
## A sensible ordering
tdist=dist(cvsdist)
tdist=dist(Clarity_Persistence(cvs_from_wb$pred))
tdist=dist(Clarity_Extract(cvs_from_wb$pred$scan[[15]],summary=I))
tdend=(as.dendrogram(hclust(tdist)))
tdend=reorder(tdend,colMeans(as.matrix(tdist)))
to=labels(tdend)

######################
pdf("Fig5_CvsVsWb.pdf",height=10,width=10)
cex.axis.X=0.7
cex.axis.Y=0.5
cex.axis=1.5
XPthresh=0.01
zmin=0;zmax=24000
scalefun=sqrt
scalemar=c(5,4,5,1)
layout(matrix(1:6,nrow=3,ncol=2,byrow=T),
       heights=c(0.9,0.6,0.9),widths=c(1,0.15))
###### Similarity
tmp=cvsdist[to,to] #pred_cvs_from_wb$scan[[myk]]$Y[to,to]
#rownames(tmp)=colnames(tmp)=rep("",109)
zmax=4#max(abs(range(tmp)))
zr=seq(0,zmax,length.out=100)
Clarity_Chart(tmp,mar=c(5,1,5,1),las.axis=2,
        scalefun=I,text=FALSE,zlim=range(zr),col=(colwr),
        cex.axis.X=cex.axis.X,cex.axis.Y=0)
abline(a=0,b=1,col="white")
mtext("a) Dissimilarity Y(Culture)",line=1,adj=0,cex=cex.axis)
## Similarity scale
par(mar=scalemar)
image(1,zr,matrix((zr),nrow=1),col=(colwr),axes=FALSE,xlab="",ylab="")
axis(2,las=1,cex.axis=cex.axis)
mtext("Distant",3,cex=0.75)
mtext("Close",1,cex=0.75)
####### Persistences
tmp=plot(scale_cvs_from_wb,
         order=to,
#    zlim=scalefun(c(zmin,zmax)),
    mar=c(5,1,5,1),scalefun=function(x){x^(1/4)},
    text=FALSE,#cols=(colwr),
    cex.axis.X=cex.axis.X,cex.axis.Y=0,las=2)
c_legend(100,115,size=c(2,4),cex.text=1,laboffset=1,gap=1)
#title(xlab="Model dimension k",line=2.5,cex.lab=cex.axis)
mtext("b) Persistence P(Culture from Economics)" ,line=1,adj=0,cex=cex.axis)
abline(h=myk,lwd=2,lty=2)
## Scale
zrl=seq(zmin,zmax,length.out=6)
zr=seq(zmin,zmax,length.out=100)
par(mar=scalemar)
image(1,zr,matrix(scalefun(zr),nrow=1),col=(colwr),axes=FALSE,xlab="",ylab="")
axis(2,las=1,cex.axis=cex.axis,at=zrl,labels=round(zrl))
title(main="")
#### residuals
rplot=plot(R15_cvs_from_wb,summary=I,order=to,
           zlim="symmetric",
           mar=c(5,1,5,1),las.axis=2,
        scalefun=I,text=FALSE,col=colbwg,
        cex.axis.X=cex.axis.X,cex.axis.Y=0)
abline(a=0,b=1,col="grey")
mtext(paste0("c) Residuals R(Culture from Economics,k=",myk,")"),line=1,adj=0,cex=cex.axis)
c_legend(100,115,size=c(2,2),cex.text=1,col=c("#0000FFFF"),laboffset=1)
###### Residuals scale
zmax=max(abs(rplot$R0))
zr=seq(-zmax,zmax,length.out=100)
par(mar=scalemar)
image(1,zr,matrix((zr),nrow=1),col=(colbwg),axes=FALSE,xlab="",ylab="")
axis(2,las=1,cex.axis=cex.axis)
                                        #mtext(expression(Y-hat(Y)),3,0)
mtext("Further\nthan\nexpected",3,cex=0.75)
mtext("Closer\nthan\nexpected",1,2.5,cex=0.75)
dev.off()
