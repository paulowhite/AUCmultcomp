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
             times=5,
             plots="ROC")
t2 <- Sys.time()
compTimeROC <- difftime(t2,t1)
## print main result
## compTimeROC
Res
## }}}

## {{{ Extract AUCs and vcov estimates
AUCs.A <- Res$AUC$score$AUC
vcovAUCs.A <- Res$AUC$vcov     # varcov matrix
names(AUCs.A) <- attributes(vcovAUCs.A)$dimnames[[1]]
## simpify the names (to improve the plot below)
names(AUCs.A) <- gsub(", times=5", "", gsub("model=", "", names(AUCs.A)))
rownames(vcovAUCs.A) <- colnames(vcovAUCs.A) <- names(AUCs.A)
## }}}



## {{{ Define the constrast matrix
Cmat.A <- multcomp:::contrMat(AUCs.A, type = "Dunnett",
                              base=which(names(AUCs.A)=="MMSE"))
Cmat.A
## }}}


## {{{ Figure of Example A
myzoomcoef <- 0.7
## pdf("Example-A.pdf",width=11.5*myzoomcoef,height=7.5*myzoomcoef)
par(mai=c(0.1,0.1,0.1,0.1))
ResPlotfunc.A <- PlotClosedTest(theest=AUCs.A,
                                thevcov=vcovAUCs.A,
                                C=Cmat.A,
                                CrossSize=5,
                                CircleSize=11,
                                computePvalues=TRUE,
                                mypalette=c("black","red","darkblue"),
                                xlim=c(-0.35,2.5),
                                ylim=c(-0.75,3.5))
## dev.off()
## }}}



## {{{ How to obtain the p-values using mulcomp
## Define some functions for compatibility with the use of multcomp
coef.AUCinference <- function(x,...){
    x$coef
}
vcov.AUCinference <- function(x,...){
    x$vcov
}
x.A <- list(coef=AUCs.A, vcov=vcovAUCs.A)
class(x.A) <- "AUCinference"
## --- use mulcomp package  ----
myglht.A <- glht(x.A,linfct=Cmat.A)
AllpMulcomp.A <- rbind(
    "Unadjusted"= summary(myglht.A,test=adjusted(type = "none"))$test$pvalues,
    "Closed max-test (constaints)" = summary(myglht.A,test=adjusted(type = "Westfall", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Shaffer" = summary(myglht.A,test=adjusted(type = "Shaffer"))$test$pvalues,
    "Closed max-test (free)"=summary(myglht.A,test=adjusted(type = "free", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni-Holm"=summary(myglht.A,test=adjusted(type = "holm"))$test$pvalues,
    "Single-step max-test"=summary(myglht.A,test=adjusted(type = "single-step", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni"=summary(myglht.A,test=adjusted(type = "bonferroni"))$test$pvalues)
round(AllpMulcomp.A,3)
## }}}



## {{{ Correlations between the test statistics
Vmat.A <- Cmat.A%*%vcovAUCs.A%*%t(Cmat.A)
A.A <- diag(sqrt(diag(Vmat.A)))
T.A <- solve(A.A)%*%Cmat.A%*%AUCs.A
Corr.A <- solve(A.A)%*%Vmat.A%*%solve(A.A)
rownames(Corr.A) <- colnames(Corr.A) <- rownames(Cmat.A)
round(Corr.A,2)
## }}}



## {{{ ROC curves Figure
## pdf("ROC-curves-5-years.pdf",width=6,height=6)
par(mai=c(0.9,0.9,0.1,0.1))
mycols <- c("blue",
            "red",
            "ForestGreen",
            "black")
plotROC(Res, time=5,
        legend=FALSE,
        col=mycols,xlab="False Positive Rate (FPR)",
        ylab="True Positive Rate (TPR)",
        )
AUCToplot <- riskRegression:::getLegendData(object = Res,
                                            times = 5, auc.in.legend = TRUE,
                                            brier.in.legend = FALSE,
                                            drop.null.model = TRUE)
for(i in 1:nrow(AUCToplot)){
    thy <- 0.250-(i/nrow(AUCToplot))*0.250
    text(x=0.5,y=thy,labels=AUCToplot[i,1], pos=4)
    text(x=0.5+0.12,y=thy,labels=AUCToplot[i,2], pos=4)
    segments(x0=0.5,x1=0.5-0.12/2,y0=thy,y1=thy,lwd=2,col=mycols[i])
}
dev.off()
## }}}
