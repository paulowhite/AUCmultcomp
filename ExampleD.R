rm(list=ls())

setwd("~/set/your/own/path")

## {{{ Technical Parameters
## important parameters for precision of p-value computation
## (passed in ... in pmvnorm(), see documentation about the algorithms of pmvnorm for details)
myabseps <- 1e-05
mymaxpts <- 2500000
## }}}

## {{{ source useful R functions
source("Rfunctions/PlotClosedTest.R")
source("Rfunctions/SingleStep.R")
## }}}

## {{{ load packages
library(riskRegression)
library(survival)
library(multcomp)
library(ggsci)
## }}}

## {{{ load data
d <- read.table(file="paquid.csv",
                header=TRUE,
                sep=";")
head(d)
## }}}


## {{{ Compute estimated AUC and var-cov
t1 <- Sys.time()   
Res <- Score(list("DSST"=-d$DSST,
                  "MMSE"=-d$MMSE,
                  "IST"=-d$IST,
                  "BVRT"=-d$BVRT
                  ),
             formula=Hist(time,status)~1,
             data=d,
             keep="vcov",
             null.model = FALSE,
             conf.int=TRUE,
             metrics=c("auc"),
             times=c(5,10),
             plots="ROC")
t2 <- Sys.time()
compTimeROC <- difftime(t2,t1)
## print main result
## compTimeROC
Res
## }}}


## {{{ Extract AUCs and vcov estimates
AUCs.D <- Res$AUC$score$AUC
vcovAUCs.D <- Res$AUC$vcov     # varcov matrix
names(AUCs.D) <- attributes(vcovAUCs.D)$dimnames[[1]]
## simpify the names (to improve the plot below)
rownames(vcovAUCs.D) <- colnames(vcovAUCs.D) <- names(AUCs.D)
NewOrderNames.D <- c("model=DSST, times=5", 
                     "model=MMSE, times=5", 
                     "model=IST, times=5", 
                     "model=BVRT, times=5",
                     "model=DSST, times=10", 
                     "model=MMSE, times=10", 
                     "model=IST, times=10", 
                     "model=BVRT, times=10"
                     )
AUCs.D <- AUCs.D[NewOrderNames.D]
vcovAUCs.D <- vcovAUCs.D[NewOrderNames.D,NewOrderNames.D]
## }}}



## {{{ Define the constrast matrix
AUCs.B <- AUCs.D[1:4]
names(AUCs.B) <- gsub(", times=5", "", gsub("model=", "", names(AUCs.B)))
Cmat.B <- multcomp:::contrMat(AUCs.B, type = "Tukey")
Cmat.D <- rbind(cbind(Cmat.B,matrix(0,ncol=4,nrow=6)),cbind(matrix(0,ncol=4,nrow=6),Cmat.B))
colnames(Cmat.D) <- gsub("times", "t", gsub("model=", "", names(AUCs.D)))
rownames(Cmat.D) <- c(paste0(rownames(Cmat.B), ", t=5"),paste0(rownames(Cmat.B), ", t=10"))
Cmat.D
## }}}


## {{{ Computation for Example D
## This needs 4.5 hours to run
## StartClock.D <- Sys.time()
## ResPlotfunc.D <- PlotClosedTest(theest=AUCs.D,
                                ## thevcov=vcovAUCs.D,
                                ## C=Cmat.D,
                                ## show= c("closed","single-step",
                                        ## "bonferroni","bonf-Holm",
                                        ## "westfall",
                                        ## "shaffer"),
                                ## CrossSize=5,
                                ## CircleSize=11,
                                ## computePvalues=TRUE,
                                ## xlim=c(-0.3,2.5),
                                ## ylim=c(-0.75,3.5))
## StopClock.D <- Sys.time()
## compTimeResPlot.D <- difftime(StopClock.D,StartClock.D)
## }}}

## {{{ A comment
## Note: computing the adjusted p-values using the function PlotClosedTest()
## is possible, but not recommended. This is because: (1) the computation is using this function is
## very slow (unlike the alternative code below) and (2) the motivation of using this function
## is unclear as a plot including >4000 intersection hypotheses is expected to be useless.
## }}}


## {{{ How to obtain the p-values using mulcomp
## Define some functions for compatibility with the use of multcomp
coef.AUCinference <- function(x,...){
    x$coef
}
vcov.AUCinference <- function(x,...){
    x$vcov
}
x.D <- list(coef=AUCs.D, vcov=vcovAUCs.D)
class(x.D) <- "AUCinference"
## --- use mulcomp package  ----
myglht.D <- glht(x.D,linfct=Cmat.D)
AllpMulcomp.D <- rbind(
    "Unadjusted"= summary(myglht.D,test=adjusted(type = "none"))$test$pvalues,
    "Closed max-test (constaints)" = summary(myglht.D,test=adjusted(type = "Westfall", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Shaffer" = summary(myglht.D,test=adjusted(type = "Shaffer"))$test$pvalues,
    "Closed max-test (free)"=summary(myglht.D,test=adjusted(type = "free", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni-Holm"=summary(myglht.D,test=adjusted(type = "holm"))$test$pvalues,
    "Single-step max-test"=summary(myglht.D,test=adjusted(type = "single-step", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni"=summary(myglht.D,test=adjusted(type = "bonferroni"))$test$pvalues)
round(AllpMulcomp.D,3)
## }}}


## {{{ Technicalities and parameters to create the plot of p-values
maxCardJ <- multcomp:::maxsets(Cmat.D[order(AllpMulcomp.D[1,]),]) 
matPval <- t(AllpMulcomp.D[,order(AllpMulcomp.D[1,])])
yatnames <- c(0,0.01,
              0.05,
              0.10,
              0.20,
              0.5,1)
transfo <- function(x){log(1+40*x)}
matPval0 <- matPval
matPval <- transfo(matPval)
yat <- transfo(yatnames)
mycols <- c("black","red","ForestGreen","orange","darkblue","purple","grey40")
mylwd <- 1.5
mypch <- 19
ylimlag <- 2
myzoomcoef <- 0.7
## }}}

## {{{ plot of p-values Example D
## pdf("Example-D-2.pdf",width=11.5*myzoomcoef,height=13*myzoomcoef)
par(mai=c(3.7, 2.65, 0.2, 0.2))
## --
plot(1:nrow(matPval),
     matPval[,1],
     type="b",axes=FALSE,ylim=c(0,transfo(1+ylimlag)),
     pch=mypch,
     xlim=c(1.25,nrow(matPval)),
     xlab="",ylab="p-value",lwd=mylwd,
     col=mycols[1])
## --
segments(x0=0,x1=nrow(matPval)*0.98,
         y0=yat[which(yatnames==0.05)],y1=yat[which(yatnames==0.05)],
         col="grey60",lty=2,lwd=2)
## --
axis(1,at=1:nrow(matPval),rownames(matPval),las=2)
axis(2,at=yat,las=2,labels=paste0(100*yatnames,"%"))
## --
for(i in 2:ncol(matPval)){
    lines(1:nrow(matPval),
          matPval[,i],col=mycols[i],lwd=mylwd,type="b",
          pch=mypch)
}
## --
text(x=1:nrow(matPval),y=transfo(1+ylimlag),
     do.call("c",lapply(maxCardJ,length)))
mtext(at=transfo(1+ylimlag),expression(paste("#",S[k]," :")),line=5,side=2,las=2,cex=1.2)
mtext(at=transfo(1+ylimlag/2),expression(paste(max("(#J)",J %subset% S[k])," :")),line=5,side=2,las=2,cex=1.2)
text(x=1:nrow(matPval),y=transfo(1+ylimlag/2),
     do.call("c",lapply(maxCardJ,function(x)max(do.call("c",lapply(x,length))))))
## --
for(j in 1:ncol(matPval)){
    mtext(at=c(-0,1:nrow(matPval)+0.3),
          side=1,
          text=c(colnames(matPval0)[j],
                 format(round(matPval0[1:2,j]*100,digits=2),nsmall=2),
                 format(round(matPval0[3:6,j]*100,digits=1),nsmall=1),
                 format(round(matPval0[7:nrow(matPval0),j]*100,digits=0),nsmall=0)
                 ),
          line=10+j-1,adj=1,col=mycols[j])
}
## dev.off()
## }}}
