source('scripts/tiltedCDF.R')
dnorminv<-function(y) sqrt(-2*log(sqrt(2*pi)*y))
NormalRiskDiff <- function(t,this.mean,this.sd,this.pi){
    risk <- pnorm(t-1,mean=this.mean,sd=this.sd,lower.tail=FALSE)
    prot <- pnorm(t,mean=this.mean,sd=this.sd,lower.tail=FALSE)
    risk-prot-this.pi
}

pPoisConv <- function(t,lambda,norm.sd,alphal=1,risk.allele=FALSE){
    my.range <- seq(qpois(0,lambda),qpois(1-1e-8,lambda))
    gen.dist <- dpois(my.range,lambda=lambda)
    prevs <- pnorm(t,alphal*(my.range+risk.allele),norm.sd,lower.tail=FALSE)
    gen.dist%*%prevs
}
dPoisConv <- function(t,lambda,norm.sd,alphal=1,risk.allele=FALSE){
    my.range <- seq(qpois(0,lambda),qpois(1-1e-8,lambda))
    gen.dist <- dpois(my.range,lambda=lambda)
    dens <- dnorm(t,alphal*(my.range+risk.allele),norm.sd)
    gen.dist%*%dens
}
PoissonRiskDiff <- function(t,mean,env.sd,gen.range,this.pi){
    prot <- pPoisConv(t,mean,env.sd,gen.range,risk.allele=FALSE)
    risk <- pPoisConv(t,mean,env.sd,gen.range,risk.allele=TRUE)
    risk - prot - this.pi

}
PoissonRiskDiffThr <- function(t,mean,env.sd,gen.range,this.pi,bigGamma=1,rcv=FALSE,verboseOutput=FALSE){
    if(rcv)
        recover()
    prot <- pPoisConv(t,mean,env.sd,gen.range,risk.allele=FALSE)
    risk <- pPoisConv(t,mean,env.sd,gen.range,risk.allele=TRUE)
    riskDiff <- risk-prot
    if(verboseOutput)
        return(c(riskDiff,this.pi)*bigGamma)
    else
        return((riskDiff - this.pi)*bigGamma)
}
getAsym <- function(y){
    ifelse(y>37,1,(exp(y)-1)/(exp(y)+1))
}
SolveNormal <- function(myTheta,myGamma,fitCost,Ne,L,h2,verbose.output=FALSE) {
    ## recover()

    my.pi <- myGamma/(4*Ne*fitCost)
    ## my.mean <- ifelse(
    ##     myGamma<100,
    ##     Re(WrightPhenoCGFprime(
    ##         s=0,
    ##         myTheta=myTheta,
    ##         myGamma=myGamma,
    ##         Ne=Ne,
    ##         L=L,
    ##         ploidy=2,
    ##         h2=h2
    ##     )),
    ##     myTheta*L/myGamma
    ## )
    b <- getAsym(myGamma)
    my.var <- myTheta*L*b/(myGamma*h2)
    ## WrightPhenoCGFdoubleprime(
    ##     s=0,
    ##     myTheta=myTheta,
    ##     myGamma=myGamma,
    ##     Ne=Ne,
    ##     L=L,
    ##     ploidy=2,
    ##     h2=h2
    ## )
    my.seq <- seq(0,0+20*sqrt(my.var),length.out=1000)
    blah <- NormalRiskDiff(my.seq,0,sqrt(my.var),my.pi)
    lower <- my.seq[max(which(blah>0))]
    upper <- my.seq[min(which(blah<0))]

    this.t <- uniroot(
        f=function(t)  NormalRiskDiff(t,0,sqrt(my.var),my.pi),
        lower=lower,
        upper=upper
    )$root


    my.dens <- dnorm(this.t,mean=0,sd=sqrt(my.var))
    my.prev <- pnorm(this.t,mean=0,sd=sqrt(my.var),lower.tail=FALSE)
    if(verbose.output) {
        list(
            cent.t=this.t,
            std.t=this.t/sqrt(my.var),
            ph.sd=sqrt(my.var),
            prev=my.prev,
            dens=my.dens
        )
    } else {
        return(my.prev)
    }

}


SolvePoisson <- function(myTheta,myGamma,fitCost,Ne,L,h2,allow.small.effect=FALSE,verbose.output=FALSE){

    if(allow.small.effect){
        if(myGamma<100) warning('Scaled selection coefficient may be too small')
    } else {
        stopifnot(myGamma>100)
    }
    my.pi <- myGamma/(4*Ne*fitCost)
    my.mean <- myTheta*L/myGamma
    upper <- qpois(1-1e-7,my.mean)
    lower <- qpois(1e-7,my.mean)
    env.sd <- sqrt(my.mean*(1-h2)/h2)
    my.range <- max((lower-1),0):(upper+1)

    this.t <- uniroot(
        f=function(t)  PoissonRiskDiff(t,my.mean,env.sd,my.range,my.pi),
        lower=my.mean,
        upper=2*upper-my.mean
    )$root
    my.dens <- dPoisConv(this.t,my.mean,env.sd,my.range)
    my.prev <- pPoisConv(this.t,my.mean,env.sd,my.range)
    if(verbose.output){
        list(
            prev=as.numeric(my.prev),
            t=this.t,
            mean=my.mean,
            env.sd=env.sd,
            range=my.range,
            prevs=my.prev,
            dens=my.dens
        )
    } else {
        as.numeric(my.prev)
    }


}

SolvePoissonGivenThr <- function(myTheta,thr,fitCost,Ne,L,h2,allow.small.effect=FALSE,verbose.output=FALSE){


    bigGamma <- 2*Ne*fitCost
    max.pi <- 1
    min.pi <- 1e-8
    thetaU <- myTheta*L


    getPoisRiskDiff <- function(THISPI) {
        my.mean = thetaU/(THISPI*bigGamma)
        env.sd = sqrt(my.mean*(1-h2)/h2)
        my.range <- seq(qpois(1e-8,my.mean),qpois(1-1e-8,my.mean),by=1)
        PoissonRiskDiffThr(thr,my.mean,env.sd,my.range,THISPI,bigGamma)
    }

    my.seq <- exp(seq(log(min.pi),log(max.pi),length.out=1000))
    my.diffs <- sapply(my.seq,getPoisRiskDiff)
    idx <- max(which(diff(sign(my.diffs))!=0))
    lower <- my.seq[idx]
    upper <- my.seq[idx+1]

    soln <- uniroot(
        f=function(THISPI) {
            my.mean = thetaU/(THISPI*bigGamma)
            env.sd = sqrt(my.mean*(1-h2)/h2)
            my.range <- seq(qpois(1e-8,my.mean),qpois(1-1e-8,my.mean),by=1)
            PoissonRiskDiffThr(thr,my.mean,env.sd,my.range,THISPI,bigGamma)
        },
        lower=lower,
        upper=upper,
        tol=.Machine$double.eps^0.5
    )


    this.pi <- soln$root
    my.mean <- thetaU/(this.pi*bigGamma)
    my.range <- seq(qpois(1e-8,my.mean),qpois(1-1e-8,my.mean),by=1)
    env.sd = sqrt(my.mean*(1-h2)/h2)
    ## PoissonRiskDiffThr(thr,my.mean,env.sd,my.range,this.pi,bigGamma,rcv=FALSE)
    my.sel <- this.pi*fitCost
    my.gamma <- my.sel* 4*Ne #bigGamma


    my.dens <- dPoisConv(thr,my.mean,env.sd,my.range)
    my.prev <- pPoisConv(thr,my.mean,env.sd,my.range)
    if(verbose.output){
        list(
            prev=as.numeric(my.prev),
            thr=thr,
            mean=my.mean,
            env.sd=env.sd,
            range=my.range,
            dens=my.dens,
            my.sel=my.sel,
            my.gamma=my.gamma,
            this.pi=this.pi,
            cost=fitCost,
            target=L,
            theta=myTheta,
            h2=h2
        )
    } else {
        as.numeric(my.prev)
    }


}



SolveTiltedSingleEffectForLGivenR <- function(myTheta,myGamma,fitCost,Ne,prev,h2,alpha=1,ploidy=2,env.mean=0,env.shape='normal') {
    ## work in progress
    ## recover()

    weights <- GetWeights(theta=myTheta,gamma=myGamma,Ne=Ne)

    my.Ls <- c(10^seq(3,8,by=1),10^seq(8.25,14,by=1/4))
    my.obj <- numeric()
    i <- 1
    go <- TRUE
    while(i==1|go){
        my.obj[i] <- prevObjective(
            gamma=myGamma,
            theta=myTheta,
            alpha=alpha,
            Ne=Ne,
            L=my.Ls[i],
            h2=h2,
            fitCost=fitCost,
            prev=prev,
            ploidy=ploidy,
            env.mean=env.mean,
            renew.weights=FALSE,
            env.shape=env.shape
        )
        if(i>1){
            go <- sign(my.obj[i])==sign(my.obj[i-1])
        }
        i <- i+1
    }
    my.ind <- length(my.obj)

    solnL <- uniroot(
        f=function(ELL){
            prevObjective(
                gamma=myGamma,
                theta=myTheta,
                alpha=alpha,
                Ne=Ne,
                L=ELL,
                h2=h2,
                fitCost=fitCost,
                prev=prev,
                ploidy=ploidy,
                env.mean=env.mean,
                renew.weights=FALSE,
                env.shape=env.shape
            )
        },
        interval=my.Ls[c(my.ind-1,my.ind)]#c(1e3,5e11)
    )

    return(
        list(
            eta=myEta,
            beta=myBeta,
            L=solnL$root
        )
    )

}


##################################################################################################
## main function to solve for threshold (and by proxy, prevalence), given selection coefficient ##
##################################################################################################
SolveTiltedSingleEffectForThr <- function(myTheta,myGamma,fitCost,Ne,L,h2,alpha=1,ploidy=2) {

    my.pi <- myGamma/(4*Ne*fitCost)

    my.mean <- WrightPhenoCGFprime(0,myTheta,myGamma,alpha,Ne,L,ploidy,h2)
    my.sd <- sqrt(WrightPhenoCGFdoubleprime(0,myTheta,myGamma,alpha,Ne,L,ploidy,h2))
    out.length <- 5
    thr.seq <- seq(my.mean,my.mean+20*my.sd,length.out=out.length)
    not.solved <- TRUE
    j <- 0
    while(not.solved){

        eta.seq <- sapply(
            thr.seq,
            function(THR) nleqslv(
                              x=0,
                              function(ETA) WrightPhenoCGFprime(ETA,myTheta,myGamma,alpha,Ne,L,ploidy,h2)-THR
                          )$x
        )
        beta.seq <- mapply(
            function(THR,ETA) nleqslv(
                                  x=0,
                                  function(BETA) WrightPhenoCGFprime(ETA-BETA,myTheta,myGamma,alpha,Ne,L,ploidy,h2)-(THR-alpha)
                              )$x,
            THR=thr.seq,
            ETA=eta.seq
        )
        risk.obj <- mapply(
            function(ETA,BETA)
                riskObjective(ETA,BETA,my.pi,myTheta,myGamma,alpha,Ne,L,ploidy,h2),
            ETA=eta.seq,
            BETA=beta.seq
        )
        idx <- which(sign(risk.obj)==1& c(diff(sign(risk.obj))!=0,FALSE))
        if(length(idx)==0) {
            out.length <- 2*out.length
            thr.seq <- seq(my.mean,tail(thr.seq,1),length.out=out.length)
            next
        }
        flip <- idx:(idx+1)

        sol.att <- nleqslv(
            x = c(mean(eta.seq[flip]),mean(beta.seq[flip])),
            fn=TwoObjectiveGammaFixed,
            my.pi=my.pi,
            myTheta=myTheta,
            myGamma=myGamma,
            L=L,
            Ne=Ne,
            h2=h2,
            alpha=alpha,
            ploidy=ploidy
        )

        ## sol.att <- nleqslv(
        ##     x = c(mean(eta.seq[flip]),mean(beta.seq[flip]),mean(thr.seq[flip])),
        ##     fn=MultiObjectiveGammaFixed
        ## )

        if(sol.att$termcd==1){
            not.solved <- FALSE
            eta <- sol.att$x[1]
            beta <- sol.att$x[2]
            ## thr <- sol.att$x[3]
        } else if (j>5){
            return(NULL)

        } else {
            out.length <- 2*out.length
            thr.seq <- seq(thr.seq[flip[1]],thr.seq[flip[2]],length.out=out.length)
            j <- j+1
        }

    }

    thr <- nleqslv(
        x=mean(thr.seq[flip]),
        fn=function(X) WrightGenoCGFprime(eta,myTheta,myGamma,alpha,Ne,L,ploidy)-X
    )$x

    return(
        list(
            eta=eta,
            beta=beta,
            thr=thr,
            L=L
        )
    )

}
######################################################################################################
## main function to solve for gamma (and by proxy, prevalence and risk effect), given the threshold ##
######################################################################################################
SolveTiltedSingleEffectForGamma <- function(myTheta,thr,fitCost,Ne,L,h2,alpha=1,ploidy=2,env.mean=0,env.var=NULL,eta.bounds=c(1e-8,6)) {

    ## recover()

    gamma.lb <- nleqslv(
        x=log(2*alpha*L/thr-1)/2,
        fn=function(GAMMA) WrightGenoCGFprime(s=0,myTheta,GAMMA,alpha,Ne,L,ploidy)-thr
    )$x


    if(gamma.lb>10){
        gamma.ub <- 1480
    } else {
        out <- nleqslv(
            x=gamma.lb,
            fn=function(GAMMA){
                gen.mean <- WrightGenoCGFprime(
                    s=0,
                    myTheta=myTheta,
                    myGamma=GAMMA,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=2
                )
                phen.sd <- WrightPhenoCGFdoubleprime(
                    s=0,
                    myTheta=myTheta,
                    myGamma=GAMMA,
                    alpha=alpha,
                    Ne=Ne,
                    L=L,
                    ploidy=2,
                    h2=h2
                )^(1/2)
                10-(thr-gen.mean)/phen.sd
            }
        )
        gamma.ub <- out$x
    }



    out.length <- 10
    gamma.seq <- seq(gamma.lb,gamma.ub,length.out=out.length)
    not.solved <- TRUE

    j=1
    while(not.solved){

        ## search for starting point ##
        eta.seq <- sapply(
            gamma.seq,
            function(GAMMA) nleqslv(
                                x=0,
                                function(ETA) WrightPhenoCGFprime(ETA,myTheta,GAMMA,alpha,Ne,L,ploidy,h2)-thr
                            )$x
        )
        beta.seq <- mapply(
            function(GAMMA,ETA) nleqslv(
                                    x=0,
                                    function(BETA) WrightPhenoCGFprime(ETA-BETA,myTheta,GAMMA,alpha,Ne,L,ploidy,h2)-(thr-alpha)
                                )$x,
            GAMMA=gamma.seq,
            ETA=eta.seq
        )
        risk.obj <- mapply(
            function(ETA,BETA,GAMMA)
                riskObjective(ETA,BETA,GAMMA/(4*Ne*fitCost),myTheta,GAMMA,alpha,Ne,L,ploidy,h2),
            ETA=eta.seq,
            BETA=beta.seq,
            GAMMA=gamma.seq
        )
        risk.obj


        ## check whether there is a good starting point in interval
        ## if not, expand interval appropriately, and try againi
        if(all(risk.obj<0)){
            gamma.seq <- seq(gamma.lb,gamma.seq[2],length.out=out.length)
            next
        } else if(tail(risk.obj,1)>0){
            gamma.ub <- 2*gamma.ub
            gamma.seq <- seq(gamma.lb,gamma.ub,length.out=out.length)
            next
        }


        idx <- which(sign(risk.obj)==1& c(diff(sign(risk.obj))!=0,FALSE))
        flip <- idx:(idx+1)


        ## solve ##
        start <- c(mean(eta.seq[flip]),mean(beta.seq[flip]),mean(gamma.seq[flip]))
        ##start <- c(min(eta.seq[flip]),min(beta.seq[flip]),min(gamma.seq[flip]))
        sol.att <- nleqslv(
            x = start,
            fn=ThreeObjectiveThrFixed,
            thr=thr,
            myTheta=myTheta,
            alpha=alpha,
            Ne=Ne,
            L=L,
            ploidy=ploidy,
            h2=h2,
            fitCost=fitCost
        )



        if(sol.att$termcd==1){
            not.solved <- FALSE
            eta <- sol.att$x[1]
            beta <- sol.att$x[2]
            myGamma <- sol.att$x[3]
        }  else {
            out.length <- out.length*2
            gamma.seq <- seq(gamma.seq[flip[1]],gamma.seq[flip[2]],length.out=out.length)
        }
        if(length(risk.obj)>200)
            return(NULL)

    }

    return(
        list(
            eta=eta,
            beta=beta,
            gamma=myGamma
        )
    )

}
