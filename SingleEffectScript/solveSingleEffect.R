 #source('scripts/tiltedCDF.R')
NormalRiskDiff <- function(t,this.mean,this.sd,this.pi){
    risk <- pnorm(t-1,mean=this.mean,sd=this.sd,lower.tail=FALSE)
    prot <- pnorm(t,mean=this.mean,sd=this.sd,lower.tail=FALSE)
    risk-prot-this.pi
}

pPoisConv <- function(t,lambda,env.sd,my.range,risk.allele=FALSE){
    gen.dist <- dpois(my.range,lambda=lambda)
    prevs <- pnorm(t,my.range+risk.allele,env.sd,lower.tail=FALSE)
    gen.dist%*%prevs
}
dPoisConv <- function(t,lambda,env.sd,my.range,risk.allele=FALSE){
    gen.dist <- dpois(my.range,lambda=lambda)
    dens <- dnorm(t,my.range+risk.allele,env.sd)
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
SolveNormal <- function(myTheta,myGamma,fitCost,Ne,L,h2,verbose.output=FALSE) {
    ## recover()
    my.pi <- myGamma/(4*Ne*fitCost)
    my.mean <- ifelse(
        myGamma<100,
        Re(WrightPhenoCGFprime(
            s=0,
            myTheta=myTheta,
            myGamma=myGamma,
            Ne=Ne,
            L=L,
            ploidy=2,
            h2=h2
        )),
        myTheta*L/myGamma
    )
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
    my.seq <- seq(my.mean,my.mean+20*sqrt(my.var),length.out=1000)
    blah <- NormalRiskDiff(my.seq,my.mean,sqrt(my.var),my.pi)
    lower <- my.seq[max(which(blah>0))]
    upper <- my.seq[min(which(blah<0))]

    this.t <- uniroot(
        f=function(t)  NormalRiskDiff(t,my.mean,sqrt(my.var),my.pi),
        lower=lower,
        upper=upper
    )$root
    

    my.dens <- dnorm(this.t,mean=my.mean,sd=sqrt(my.var))
    my.prev <- pnorm(this.t,mean=my.mean,sd=sqrt(my.var),lower.tail=FALSE)
    if(verbose.output) {
        list(
            raw.t=this.t,
            cent.t=this.t-my.mean,
            std.t=(this.t-my.mean)/sqrt(my.var),
            mean=my.mean,
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
    my.gamma <- my.sel*4*Ne

    
    my.dens <- dPoisConv(thr,my.mean,env.sd,my.range)    
    my.prev <- pPoisConv(thr,my.mean,env.sd,my.range)
    if(verbose.output){
        list(
            prev=as.numeric(my.prev),
            thr=thr,
            mean=my.mean,
            env.sd=env.sd,
            range=my.range,
            prevs=my.prev,
            dens=my.dens,
            my.sel=my.sel,
            my.gamma=my.gamma,
            this.pi=this.pi
        )
    } else {
        as.numeric(my.prev)
    }


}


args = commandArgs(trailingOnly=TRUE)
liaMutRate = as.numeric(args[1])
THR = as.numeric(args[2])
fitCost = as.numeric(args[3])
popSize = as.numeric(args[4])
liabilitySize = as.numeric(args[5])
target_h2 = as.numeric(args[6])
theta = 4 * popSize * liaMutRate

if (THR <1000){
   temp = SolvePoissonGivenThr(theta,THR, fitCost,popSize, liabilitySize, target_h2,verbose.output=T)
   poissonGamma = temp$my.gamma
   if (poissonGamma <50){
      b = (exp(poissonGamma)-1)/(exp(poissonGamma)+1)
   }else {b=1}
   envSD = sqrt(2.0 * liabilitySize * liaMutRate/temp$my.sel * b * (1-target_h2)/target_h2)
}else{
   rho= THR/(2*liabilitySize)
   gamma = log((1-rho)/rho)
   envSD = sqrt(2.0 * liabilitySize * liaMutRate/(gamma/(4*popSize)) * (exp(gamma)-1)/(exp(gamma)+1) * (1-target_h2)/target_h2)
}
print(envSD)