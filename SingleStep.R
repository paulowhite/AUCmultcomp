### SingleStep.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 12 2018 (14:00) 
## Version: 
## Last-Updated: Mar 13 2018 (11:48) 
##           By: Paul Blanche
##     Update #: 26
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

SingleStep <- function(theest,thevcov,myabseps=1e-08, mymaxpts=25000){
    delta <- theest
    deltavcov <- thevcov
    # compute s.e. 
    sedelta <- sqrt(diag(deltavcov))
    # Correlation matrix
    diagInv <- diag(sedelta^(-1),ncol=length(sedelta))
    corrdelta <-  diagInv%*% deltavcov %*% t(diagInv)
    # compute individual hypothesis test statistic
    zdelta <- delta/sedelta
    # to compute single step adjusted individual p-values (max test)
    # main function
    func1 <- function(a, nstat, MatCorr) {
        p <- 1 - pmvnorm(lower = rep(-abs(a), nstat),
                         upper = rep(abs(a),nstat),
                         mean = rep(0, nstat),
                         sigma = MatCorr,
                         abseps = myabseps,
                         maxpts = mymaxpts)
        result <- p[1]
        return(as.numeric(result))
    }
    # function to loop over individual hypotheses
    func2 <- function(MatCorr, stat) {
        nstat <- length(stat)
        pval <- rep(NA, nstat)
        for (i in 1:nstat) {
            pval[i] <- func1(stat[i], nstat = nstat, MatCorr = MatCorr)
        }
        return(pval)
    }
    # compute single step adjusted individual p-values (max test)
    SingleStep <- func2(MatCorr=corrdelta, stat=zdelta)
    SingleStep   
}

##----------------------------------------------------------------------
### SingleStep.R ends here
