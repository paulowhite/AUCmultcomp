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
                  "MMSE"=-d$MMSE
                  ),
             formula=Hist(time,status)~1,
             data=d,
             keep="vcov",
             null.model = FALSE,
             conf.int=TRUE,
             metrics=c("auc"),
             times=c(3,5,10),
             plots="ROC")
t2 <- Sys.time()
compTimeROC <- difftime(t2,t1)
## print main result
## compTimeROC
Res
## }}}



## {{{ Extract AUCs and vcov estimates
AUCs.C <- Res$AUC$score$AUC
vcovAUCs.C <- Res$AUC$vcov     # varcov matrix
names(AUCs.C) <- attributes(vcovAUCs.C)$dimnames[[1]]
## simpify the names (to improve the plot below)
names(AUCs.C) <- gsub("times", "t", gsub("model=", "", names(AUCs.C)))
rownames(vcovAUCs.C) <- colnames(vcovAUCs.C) <- names(AUCs.C)
## }}}

## {{{ Define the constrast matrix
Cmat.C <- rbind(c(1,0,0,-1,0,0),
                c(0,1,0,0,-1,0),
                c(0,0,1,0,0,-1))
colnames(Cmat.C) <- names(AUCs.C)
rownames(Cmat.C) <- paste0("MMSE-DSST, ", paste0("t=",c(3,5,10)))
## }}}



## {{{ Plot Example C
myzoomcoef <- 0.7
## pdf("Example-C.pdf",width=11.5*myzoomcoef,height=7.5*myzoomcoef)
layout(matrix(1,1,1))
par(mai=c(0.1,0.1,0.1,0.1))
ResPlotfunc.C <- PlotClosedTest(theest=AUCs.C,
                                thevcov=vcovAUCs.C,
                                C=Cmat.C,
                                CrossSize=5,
                                CircleSize=11,
                                computePvalues=TRUE,
                                mypalette=c("deeppink3","black","skyblue2"),
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
x.C <- list(coef=AUCs.C, vcov=vcovAUCs.C)
class(x.C) <- "AUCinference"
## --- use mulcomp package  ----
myglht.C <- glht(x.C,linfct=Cmat.C)
AllpMulcomp.C <- rbind(
    "Unadjusted"= summary(myglht.C,test=adjusted(type = "none"))$test$pvalues,
    "Closed max-test (constaints)" = summary(myglht.C,test=adjusted(type = "Westfall", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Shaffer" = summary(myglht.C,test=adjusted(type = "Shaffer"))$test$pvalues,
    "Closed max-test (free)"=summary(myglht.C,test=adjusted(type = "free", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni-Holm"=summary(myglht.C,test=adjusted(type = "holm"))$test$pvalues,
    "Single-step max-test"=summary(myglht.C,test=adjusted(type = "single-step", abseps=myabseps, maxpts=mymaxpts))$test$pvalues,
    "Bonferroni"=summary(myglht.C,test=adjusted(type = "bonferroni"))$test$pvalues)
round(AllpMulcomp.C,4)
## }}}


## {{{ Estimated AUCs and se (%)
cbind(Res$AUC$score[,1:2],round(Res$AUC$score[,3:4]*100,1))
## }}}


## {{{ Correlations between the test statistics
Vmat.C <- Cmat.C%*%vcovAUCs.C%*%t(Cmat.C)
A.C <- diag(sqrt(diag(Vmat.C)))
T.C <- solve(A.C)%*%Cmat.C%*%AUCs.C
Corr.C <- solve(A.C)%*%Vmat.C%*%solve(A.C)
rownames(Corr.C) <- colnames(Corr.C) <- rownames(Cmat.C)
round(Corr.C,2)
## }}}



