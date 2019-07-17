### PlotClosedTest.R --- 
##----------------------------------------------------------------------
## Author: Paul Blanche
## Created: Mar 12 2018 (13:45) 
## Version: 
## Last-Updated: Jul 15 2019 (13:49) 
##           By: Paul Blanche
##     Update #: 955
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


PlotClosedTest <- function(theest,thevcov,C,digits=3,
                           abseps = 1e-05,
                           maxpts = 2500000,
                           show=c("closed-free",
                                  "single-step",
                                  "bonferroni",
                                  "bonf-Holm",
                                  "shaffer",
                                  "closed-constraints"
                                  ),
                           Trace=TRUE,
                           CrossSize=5,
                           CircleSize=7,
                           computePvalues=TRUE,
                           mypch=1,                                             # how to display the "only" hypotheses for which it is sufficient to compute a p-value
                           mypalette=rev(pal_igv("default", alpha = 0.95)(10)), # from library("ggsci")   # color to display the subsets (i.e. one color per test)
                           xlim,
                           ylim
                           ){
    require(ggsci)
    ## {{{ To nicely display p-values
    specdec <- function(x, k) {
        if(!is.na(x)){
            xr <- round(x, k)
            if(xr==0){
                if(k>1){
                    mystring <- paste0(paste(rep("0",k-1),collapse=""),1)
                }else{
                    mystring <- "1"
                }
                paste0("<.",mystring)
            }else{
                ## format(xr, nsmall=k, scientific=FALSE)
                if(xr!=1){
                    gsub("0\\.","\\.", format(xr, nsmall=k, scientific=FALSE)) # trick to not print "0.123" but ".123" instead
                }else{
                    "1"
                }
            }
        }else{
            ""
        }
    }
    ## }}}

    ## {{{ Preliminaries
    theK <- nrow(C) # number of individual hypotheses
    ## compute differences
    delta <- C%*%theest
    ## compute varcov of differences
    deltavcov <- C%*%thevcov%*%t(C)
    ## compute s.e. of differences
    sedelta <- sqrt(diag(deltavcov))
    ## }}}

    ## {{{ REORDER null hypotheses by increasing p-value
    ## unadjusted p-values
    praw <- 2*(1-pnorm(abs(as.vector(delta))/sedelta))
    ## look how to order
    TheNewOrder <- order(praw,decreasing=FALSE)
    ## reorder everything
    praw <- praw[TheNewOrder]
    delta <- delta[TheNewOrder]
    deltavcov <- deltavcov[TheNewOrder,TheNewOrder]
    sedelta <- sedelta[TheNewOrder]
    C <- C[TheNewOrder,,drop=FALSE]
    ## }}}

    ## {{{ preliminaries
    ## just to facilitate names of object without
    ## loosing original labels to display on the plot
    notNulRownames <- !is.null(rownames(C)) 
    if(notNulRownames){
        LabelsHypo <- rownames(C)
        rownames(C) <- as.character(1:nrow(C))
    }    
    ## just to index the null hypotheses
    allhypo <- as.character(1:nrow(C))     
    ## }}}

    
    ## {{{ Compute all p-values

    ## This part of the code does not call multcomp
    mytest <- function(hypo){
        indices <- which(allhypo %in% hypo)
        p <- min(SingleStep(theest=delta[indices,drop=FALSE],
                            thevcov=deltavcov[indices,indices,drop=FALSE],
                            myabseps=abseps,
                            mymaxpts=maxpts))
        p
    }
    ## compute p-values for each hypothesis in the closure
    ## (assuming free combination)
    ## combn(allhypo,x) to display all subset of allhypo of x elements
    if(computePvalues){
        Allp <- lapply(theK:1,function(x) combn(allhypo,x,mytest))
    }else{
        Allp <- lapply(theK:1,function(x) combn(allhypo,x,FUN=function(x){NA}))
    }
    ## compute names of intersection hypothesis, corresponding to the pvalues
    AllIntersecHypo <- lapply(theK:1,function(x) combn(allhypo,x))
    Allnames <- lapply(1:length(AllIntersecHypo),function(y){apply(AllIntersecHypo[[y]],2,function(x){paste("H",paste(x,collapse="."),sep=".")})})    
    AllFancyNames <- lapply(1:length(AllIntersecHypo),function(y){apply(AllIntersecHypo[[y]],2,
                                                                        function(x){
                                                                            if(length(x)!=theK){
                                                                                myJ <- paste(x,collapse=",")
                                                                            }else{
                                                                                myJ <- paste0("1,...,",theK)
                                                                            }
                                                                            if(length(x)>1){
                                                                                eval(bquote(expression(H[0]^group("{",.(myJ),"}"))))
                                                                            }else{
                                                                                eval(bquote(expression(H[0]^group("",.(myJ),""))))
                                                                            }
                                                                        })})
    thep <- do.call("c",Allp)
    names(thep) <- do.call("c",Allnames)
    ## }}}

    ## {{{ single step p-values
    sspval <- SingleStep(theest=delta,
                         thevcov=deltavcov,
                         myabseps=abseps,
                         mymaxpts=maxpts)
    ## }}}

    ## {{{ bonferroni, bonferroni-holm
    # bonferroni
    pbonf <- pmin(praw*theK,1)   
    # bonferroni-holm
    pbonf
    pBH <- pmin(praw*(theK:1),1)
    # monotonicity enforcement
    for(i in 2:theK){
        pBH[i] <- min(max(pBH[i],pBH[i-1]),1)
    }
    pBH <- pBH[names(praw)]
    # similar to :
    ## p.adjust(praw, method = "holm")
    ## }}}
    

    ## {{{ technical details for plot: set colors of each intersection hypothesis
    # we also use that to compute the "free" combination p-values
    then <- do.call("c",lapply(theK:1,function(x) choose(theK,x)))
    nmax <- max(then)
    ## browser()
    Alltocolor <- lapply(1:length(AllIntersecHypo),function(y){
        apply(AllIntersecHypo[[y]],2,function(x){
            ## min(x) # wrong because > min("12","2") returns [1] "12"
            as.character(min(as.integer(as.character(x))))
        })})
    ## }}}

    ## {{{ pvalues using closed testing assuming no restriction between hypotheses (free)
    padjclosed <- rep(NA,theK)
    names(padjclosed) <- allhypo
    for(thehypo in allhypo){
        padjclosed[thehypo] <- max(thep[do.call("c",Alltocolor)==thehypo]) # Here is what is needed to be changed for Westfall (1997) method            
    }
    ## monotonicity enforcement
    ## (besides, we need this because otherwise the above maximun is not computed over a subset large enough)
    for(i in 2:length(padjclosed)){
        padjclosed[i] <- max(padjclosed[1:i])
    }
    ## }}}

    ## {{{ To display restriction, i.e. that closure is actually smaller that assumind no restriction and compute p-values as Westfall 1997
    pWestfall <- rep(NA,theK)
    names(pWestfall) <- allhypo
    ## Here we call multcomp to use the algorithm of Westfall (1997)
    ms <- multcomp:::maxsets(C)
    ListWestfallSeq <- vector("list",theK)
    MinListSeq <- vector("list",theK)
    names(ListWestfallSeq) <- allhypo
    names(MinListSeq) <- allhypo
    ShaffMult <- rep(0,theK)
    names(ShaffMult) <- allhypo   
    for(i in 1:theK){  # Loop over individual hypothesis
        thehypoi <- allhypo[i]
        msi <- ms[[i]]
        ## if(i!=1){ListWestfallSeq[[thehypoi]] <- ListWestfallSeq[[allhypo[i-1]]]}
        if(Trace) print(paste0("To reject H.",thehypoi,", we need to reject all those"))        
        for(j in 1:length(msi)){ # Loop over all the hypotheses needed to reject to reject this individual hypothesis
            nijinter <- length(msi[[j]]) # number of intercections 
            allmaxinterij <- AllIntersecHypo[[theK+1-nijinter]] # list of all intersection of that size            
            ShaffMult[thehypoi] <- max(ShaffMult[thehypoi],nijinter) # save multiplier for shaffer correction
                                        # which column to fing the corresponding hypothesis
            whichcolij <- apply(allmaxinterij,2,function(x){
                setequal(x,(1:theK)[msi[[j]]])
            }) # be careful, does not contained our name, but the index of our names!            
            Hstepij <- paste("H",paste(allmaxinterij[,whichcolij],collapse="."),sep=".")
            if(Trace) print(paste0("below ",Hstepij, " which are inlcuded in H.",thehypoi,",  that is"))
            MinListSeq[[thehypoi]] <- c(MinListSeq[[thehypoi]],Hstepij)
            ## add all the intercetion of hypo, for which one of the hypo is the null hypo which is assessed            
            if(nijinter>1){
                ## if(i==1){}
                ## Belowijcandidates <- unlist(Allnames[(theK+2-nijinter):theK])
                Belowijcandidates <- unlist(lapply((length(AllIntersecHypo):1)[1:(nijinter-1)],
                                                   function(y){apply(AllIntersecHypo[[y]],2,function(x){paste("H",paste(x,collapse="."),sep=".")})}))

                tokeepij1 <- unlist(lapply((length(AllIntersecHypo):1)[1:(nijinter-1)],
                                           function(y){apply(AllIntersecHypo[[y]],2,function(x){all(x %in% as.character((1:theK)[msi[[j]]]))})}))
                tokeepij2 <- unlist(lapply((length(AllIntersecHypo):1)[1:(nijinter-1)],
                                           function(y){apply(AllIntersecHypo[[y]],2,function(x){i %in% x})}))

                ToaddAtij <- Belowijcandidates[tokeepij1 & tokeepij2]                
                if(j==1){
                    ListWestfallSeq[[thehypoi]] <- unique(c(Hstepij,
                                                            ToaddAtij))
                }else{
                    ListWestfallSeq[[thehypoi]] <- unique(c(ListWestfallSeq[[thehypoi]],
                                                            Hstepij,
                                                            ToaddAtij))
                }
            }else{
                ## ListWestfallSeq[[thehypoi]] <- c(ListWestfallSeq[[thehypoi]],Hstepij) # ,paste("H",thehypoi,sep="."))
                if(j==1){
                    ListWestfallSeq[[thehypoi]] <- c(Hstepij)
                }else{
                    ListWestfallSeq[[thehypoi]] <- unique(c(ListWestfallSeq[[thehypoi]],Hstepij))
                }                               
            }
            ## if(j==length(msi)){
            ## ListWestfallSeq[[thehypoi]] <- unique(ListWestfallSeq[[thehypoi]])
            ## }
            if(Trace) print(ListWestfallSeq[[thehypoi]])
        }
        pWestfall[thehypoi] <- max(thep[ListWestfallSeq[[thehypoi]]])
    }
    ## monotonicity enforcement for westfall (closed testing with restriction)
    for(i in 2:length(praw)){
        pWestfall[i] <- max(pWestfall[1:i])
    }
   
    
    # List of Hypotheses not needed to be rejected because of constraint
    # i.e. because  "free combination" condition is not satistisfied.
    unlisallnames <- unlist(Allnames)
    ## 
    AllHpoWestfall <- unique(do.call("c",ListWestfallSeq))

    ## 
    NamesToCross <- unlisallnames[which(!(unlisallnames %in% AllHpoWestfall))]
    
    # list of name of hypo to test sequentially (maxset of Westfall 1997)
    NamesToCircle  <- unique(do.call("c",MinListSeq))
    
    
    ## }}}

    ## {{{ Shaffer p-value
    pShaffer <- pmin(praw*ShaffMult[1:theK],1)
    ## monotonicity enforcement
    for(i in 2:theK){
        pShaffer[i] <- min(max(pShaffer[i],pShaffer[i-1]),1)
    }
    pShaffer <- pShaffer[names(praw)]
    ## }}}

    
    ## {{{ plot

    ## {{{ initialize
    if(missing(xlim)){ xlim <- c(-1,nmax)}
    if(missing(ylim)){ ylim <- c(-0.75,theK+0.25)}    
    plot(NULL,
         xlim=xlim,
         ## ylim=c(-0.5,theK+1),
         ylim=ylim,
         axes=FALSE,
         ylab="",
         xlab="")
    ## }}}
    ## {{{ Loop over all null hypotheses
    for(i in 1:theK){              # loop from 1 to number of hyppotheses (theK)
        thetext <- Allnames[[i]]   # Names of intersection null hypotheses
        theFancytext <- AllFancyNames[[i]]   # Name of null hypotheses
        thepval <- Allp[[i]]       # p-values (unadjusted) of the null hypotheses
        ni <- length(thetext)      # How many null hypotheses
        allthex <- seq(from=(1/(then[i] + 1)), # compute x-axis values where to plot each null intercection hypothesis
                       to=(then[i]/(then[i] + 1)),
                       length.out=then[i])*nmax       
        if(i==theK){
            thehyp <- AllIntersecHypo[[i]]    
        }
        for(j in 1:ni){         # For each number i of hypotheses we loop over all the intersection of j hypotheses  
            thex <- allthex[j]  # where to draw: x-axis value
            mycolij <- mypalette[as.numeric(Alltocolor[[i]][j])] # which color to assign to the null intercection hypothesis
            text(eval(parse(text=theFancytext[j])),y=(theK-i+1)+0.1,x=thex,col=mycolij,cex=1.2) # plot name of the null intercection hypothesis
            ## should we cross?
            if(CrossSize>0 & thetext[j]%in% NamesToCross){
                points(y=(theK-i+1)+0.1,x=thex,cex=CrossSize, pch=4, lwd=3)
            }
            ## should we circle?
            if(CircleSize>0 & thetext[j]%in% NamesToCircle){
                points(y=(theK-i+1)+0.1,x=thex,cex=CircleSize, pch=mypch, lwd=3, col=mycolij)
            }
            text(specdec(thepval[j],digits),y=(theK-i+1),x=thex,pos=1,col=mycolij) # print p-value on the plot
            ## Now we print the adjusted p-values for each individual null hypothesis using the different methods
            if(i==theK){
                if("closed-constraints" %in% show ){
                    text(specdec(pWestfall[j],digits),y= 0.4, x=thex, pos=1,col=mycolij)
                }
                if("shaffer" %in% show ){
                    text(specdec(pShaffer[j],digits),y= 0.2, x=thex, pos=1,col=mycolij)
                }                 
                if("closed-free" %in% show ){
                    text(specdec(padjclosed[j],digits),y= 0, x=thex, pos=1,col=mycolij)
                }
                if("bonf-Holm" %in% show ){
                    text(specdec(pBH[j],digits),y= -0.2, x=thex, pos=1,col=mycolij)
                }        
                if("single-step" %in% show ){
                    text(specdec(sspval[j],digits), y=-0.4,x=thex,pos=1,col=mycolij)
                }
                if("bonferroni" %in% show ){
                    text(specdec(pbonf[j],digits), y= -0.6, x=thex, pos=1,col=mycolij)
                }
            }
        }                
    }
    ## add unadjusted p-values (in case we did not displays all p-val of the closure, otherwise already done)
    if(!computePvalues){
        for(i in 1:theK){
            text(specdec(praw[i],digits),
                 y=1,
                 x=allthex[i],
                 pos=1,
                 col=mypalette[Alltocolor[[theK]][i]]) 
        }
    }
    ## }}}
    ## {{{ Add legend
    text("Unadjusted: ",y= 1,x=0,pos=1)
    if("closed-constraints" %in% show ){text("Closed max-t test (constraints): ",y= 0.4,x=0,pos=1)}
    if("shaffer" %in% show ){text("Shaffer: ",y= 0.2,x=0,pos=1)}   
    if("closed-free" %in% show ){text("Closed max-t test (free): ",y= 0,x=0,pos=1)}
    if("bonf-Holm" %in% show ){text("Bonferroni-Holm: ",y= -0.2,x=0,pos=1)}
    if("single-step" %in% show ){text("Single-step max-t test: ",y= -0.4, x=0,pos=1)}
    if("bonferroni" %in% show ){text("Bonferroni: ",y= -0.6,x=0,pos=1)}
    ## add number of intersection, i.e. cardinality of J
    for(i in 1:(theK-1)){
        text(paste0("(#J=",i+1,")"),y= (i+1)+0.1,x=0,pos=2)
    }                            
    ## To add labels corresponding to the index of each null hypothesis
    if(notNulRownames){
        for(j in 1:nrow(C))
            text(paste0("(",LabelsHypo[j],")"), y= 0.8, x=allthex[j], pos=1)
    }
    ## }}}
    ## }}}
    
    
    ## {{{ Prepare some output  
    AllpvaluesOut <- cbind(praw,
                           pWestfall,
                           pShaffer,
                           padjclosed,
                           pBH,
                           sspval,
                           pbonf                           
                           )
    colnames(AllpvaluesOut) <- c("Unadjusted",
                                 "Closed max-t test (constraints)",
                                 "Shaffer",
                                 "Closed max-t test (free)",
                                 "Bonferroni-Holm",
                                 "Single-step max-t test",
                                 "Bonferroni"                                 
                                 )
    if(notNulRownames){
        rownames(AllpvaluesOut) <- LabelsHypo
    }
    ## }}}
    
    ## {{{ Uniform comparison of all methods
    if(computePvalues){
        CheckUnif <- matrix(NA,ncol=dim(AllpvaluesOut)[2],
                            nrow=dim(AllpvaluesOut)[2])
        colnames(CheckUnif) <- colnames(AllpvaluesOut)
        rownames(CheckUnif) <- colnames(AllpvaluesOut)   
        for(i in 1:ncol(CheckUnif)){
            for(j in 1:ncol(CheckUnif)){
                rownames(CheckUnif)[i]
                colnames(CheckUnif)[j]
                AreAllSmaller <- all(max(AllpvaluesOut[,i]-AllpvaluesOut[,j])<=abseps)
                AreAllBigger <-  all(min(AllpvaluesOut[,i]-AllpvaluesOut[,j])>=-abseps)
                if(AreAllSmaller){
                    CheckUnif[i,j] <- "<="
                }
                if(AreAllBigger){
                    CheckUnif[i,j] <- ">="
                }
                if( !AreAllSmaller & !AreAllBigger){
                    CheckUnif[i,j] <- "><"
                }
                if( all(max(abs(AllpvaluesOut[,i]-AllpvaluesOut[,j]))<=abseps) ){
                    CheckUnif[i,j] <- "="
                }
            }
        }
        diag(CheckUnif) <- NA
    }else{
        CheckUnif <- NA
    }
    ## }}}

    ## 
    ## {{{ Set of maximal sets that we
    ## look at with the stepwise (step-down) procedure of  Westfall (1997)
    ## in the format as the Table 2 of Westfall (1997)
    AllMSwithNames <- lapply(ms,function(x){lapply(x,function(y){
        rownames(C)[y]})})
    msnames <- paste0("{",sapply(lapply(AllMSwithNames,function(x){lapply(x, function(y){paste0("{",paste0(y, collapse=","),"}")})}),paste,collapse=", "),"}")
    ## }}}

    ## {{{ 
    IntHypoPvalues <- Allp
    for(i in 1:length(IntHypoPvalues)){
        names(IntHypoPvalues[[i]]) <- Allnames[[i]]
    }
    ## }}}
    
    ## {{{ output
    out <- list(pvalues=AllpvaluesOut,
                ListWestfallSeq=ListWestfallSeq,
                msraw=ms, 
                ms=msnames,
                AllHpoWestfall=AllHpoWestfall,
                MinListSeq=MinListSeq,
                UnifComp=CheckUnif,
                IntHypoPvalues=IntHypoPvalues)
    ## }}}
    out
}

##----------------------------------------------------------------------
### PlotClosedTest.R ends here

