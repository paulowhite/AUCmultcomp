rm(list=ls())

setwd("~/set/your/own/path")


## {{{ Technical Parameters
## important parameters for precision of p-value computation
## (passed in ... in pmvnorm(), see documentation about the algorithms of pmvnorm for details)
myabseps <- 1e-05
mymaxpts <- 2500000
## }}}

## {{{ source useful R functions
source("PlotClosedTest.R")
source("SingleStep.R")
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
AUCs.B <- Res$AUC$score$AUC
vcovAUCs.B <- Res$AUC$vcov     # varcov matrix
names(AUCs.B) <- attributes(vcovAUCs.B)$dimnames[[1]]
## simpify the names (to improve the plot below)
names(AUCs.B) <- gsub(", times=5", "", gsub("model=", "", names(AUCs.B)))
rownames(vcovAUCs.B) <- colnames(vcovAUCs.B) <- names(AUCs.B)
## }}}



## {{{ Define the constrast matrix
Cmat.B <- multcomp:::contrMat(AUCs.B, type = "Tukey")
Cmat.B
## }}}


## {{{ Plot Example B
myzoomcoef <- 0.85
## pdf("Example-B.pdf",width=18*myzoomcoef,height=10*myzoomcoef)
par(mai=c(0.1,0.1,0.1,0.1))
ResPlotfunc.B <- PlotClosedTest(theest=AUCs.B,
                                thevcov=vcovAUCs.B,
                                C=Cmat.B,
                                CrossSize=5,
                                CircleSize=11,
                                computePvalues=TRUE,
                                mypalette=c("black","red","ForestGreen","orange","darkblue","purple"),
                                xlim=c(-0.89,19))
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
x.B <- list(coef=AUCs.B, vcov=vcovAUCs.B)
class(x.B) <- "AUCinference"
## --- use mulcomp package  ----
myglht.B <- glht(x.B,linfct=Cmat.B)
AllpMulcomp.B <- rbind(
    "Unadjusted"= summary(myglht.B,test=adjusted(type = "none"))$test$pvalues,
    "Closed max-t test (constaints)" = summary(myglht.B,test=adjusted(type = "Westfall", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Shaffer" = summary(myglht.B,test=adjusted(type = "Shaffer"))$test$pvalues,
    "Closed max-t test (free)"=summary(myglht.B,test=adjusted(type = "free", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni-Holm"=summary(myglht.B,test=adjusted(type = "holm"))$test$pvalues,
    "Single-step max-t test"=summary(myglht.B,test=adjusted(type = "single-step", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni"=summary(myglht.B,test=adjusted(type = "bonferroni"))$test$pvalues)
round(AllpMulcomp.B,3)
## }}}



## {{{ Correlations between the test statistics
Vmat.B <- Cmat.B%*%vcovAUCs.B%*%t(Cmat.B)
A.B <- diag(sqrt(diag(Vmat.B)))
T.B <- solve(A.B)%*%Cmat.B%*%AUCs.B
Corr.B <- solve(A.B)%*%Vmat.B%*%solve(A.B)
round(Corr.B[order(-abs(T.B)),order(-abs(T.B))],2)
## }}}


