source('scripts/solveSingleEffect.R')


solveTwoEffect <- function(bs,bt,Ne,Ls,Ll,as,al,r2n=NULL,Ve=NULL,u,C,LL.soln=FALSE,var.ratio=NULL,equalize.observed.vars=FALSE){

    ## bs:  asymmetry for small effects
    ## Ne:  pop size
    ## Ls:  number of small effect loci
    ## Ll:  number of large effect loci
    ## as:  small effect size
    ## al:  large effect size
    ## r2n: ratio of Ve to Vas (small effect variance)
    ## u:   mutation rate
    ## T:   threshold position
    ## C:   cost
    ## h2l: ratio of large effect variance to total genetic variance (not currently used)

    ## bt <- bs

    if(recover.flag) recover()

    if(LL.soln)
        stopifnot(!is.null(var.ratio))

    if((!is.null(r2n))&(!is.null(Ve))){
        stop('either r2n or Ve must be NULL')
    }


    ## solve for ft
    ys <- log((1+bs)/(1-bs))
    ft <- ys*(4*Ne*C)^-1  #originally ys*(2*Ne*C)^-1

    ## get small effect additive genetic variance
    Vas <- 8*Ne*u*Ls*as^2*bs/ys  ## originally 4*Ne*u*Ls*as^2*bs/ys
    if(is.null(Ve))
        Ve <- Vas*(1-r2n)/r2n
    norm.sd <- sqrt(Vas+Ve)

    ## initial guess for deltal
    init.deltal <- 1##al*ft



    ## compute selection coefficient
    init.sl <- init.deltal*C
    init.yl <- 4*Ne*init.sl  #originally 2*Ne*init.sl

    if(init.yl<3)
        warning('initial scaled selection coefficient is too small')


    ## mean number of large effect alleles per individual
    init.mean.nl <- 2*Ll*u/init.sl

    ## solve for tstart conditional on initial guess of deltal
    init.soln <- uniroot(
        f=function(THR){
            ftild <- dPoisConv(
                t=THR,
                lambda=init.mean.nl,
                norm.sd=norm.sd,
                alphal=al
            )
            ft-ftild
        },
        interval=c(0,12*norm.sd)
    )
    init.tstar <- init.soln$root


    init.Ll <- Ll



    ## get 2d solution
    soln <- nleqslv(
        x=c(0,init.tstar),
        fn=function(X){
            deltal <- 1/(1+exp(X[1]))
            tstar <- X[2]

            my.s <- deltal*C
            my.y <- 4*Ne*my.s ## originally 2*Ne*my.s
            mean.nl <- 2*Ll*u/my.s
            my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
            prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
            risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)
            deltal.tild <- risk - prot
            ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
            Val <- 2*al*mean.nl
            diff.one <- deltal-deltal.tild
            diff.two <- ft-ft.tild
            ##diff.three <- Val-Vas
            ##diff.three <- 4*al^2*mean.nl - Vas*h2l/(1-h2l)
            c(diff.one,diff.two)
        },
        control=list(scalex=c(1,1/init.tstar),maxit=400)
    )
    trans.deltal <- soln$x[1]
    deltal <- 1/(1+exp(trans.deltal))
    tstar <- soln$x[2]


    if(LL.soln){

        ## added the normal jitters because giving the solution as the initial condition led
        ## it to somehow find other pathological solutions for some reason. Confused..
        init.trans.deltal <- trans.deltal + rnorm(1,0,abs(trans.deltal)/1000)
        new.init.tstar <- soln$x[2] + rnorm(1,0,tstar/1000)
        init.Ll <- Ll + rnorm(1,0,Ll/1000)
        init.trans.bs <- log((1-bs)/bs)



        ## get 3d solution
        soln <- nleqslv(
            x=c(init.trans.deltal,new.init.tstar,init.Ll,init.trans.bs),
            fn=function(X){
                if(FALSE){
                    deltal <- init.trans.deltal
                    tstar <- new.init.tstar
                    Ll <- init.Ll
                }
                ## recover()

                ## inputs
                deltal <- 1/(1+exp(X[1]))
                tstar <- X[2]
                Ll <- X[3]
                bs <- 1/(1+exp(X[4]))

                ## small effect stuff
                ys <- log((1+bs)/(1-bs))
                Vas <- 8*Ne*u*Ls*as^2*bs/ys # originall 4*Ne*u*Ls*as^2*bs/ys
                if(is.null(Ve))
                    Ve <- Vas*(1-r2n)/r2n
                norm.sd <- sqrt(Vas + Ve)

                ## large effect stuff
                my.s <- deltal*C
                my.y <- 4*Ne*my.s #originally 2*Ne*my.s 
                mean.nl <- 2*Ll*u/my.s
                my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
                prot <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
                risk <- pPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=TRUE)

                ## global stuff
                L <- Ll+Ls
                meana <- (as*Ls + al*Ll)/L
                maxg <- 2*meana*L



                ## differences
                deltal.tild <- risk - prot
                diff.one <- deltal-deltal.tild

                ft.tild <- dPoisConv(tstar,mean.nl,norm.sd,alphal=al,risk.allele=FALSE)
                diff.two <- ft-ft.tild

                if(equalize.observed.vars){
                    Vol <- 2*deltal^2*mean.nl
                    Vos <- Vas*ft^2
                    diff.three <- Vol-var.ratio*Vos
                } else {
                    Val <- 2*al^2*mean.nl
                    diff.three <- Val-var.ratio*Vas
                }

                bt.tild <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)
                diff.four <- bt - bt.tild
                c(diff.one,diff.two,diff.three,diff.four)
            },
            control=list(scalex=c(1,1/new.init.tstar,1/init.Ll,1),maxit=800,allowSingular=TRUE)
        )
        trans.deltal <- soln$x[1]
        deltal <- 1/(1+exp(trans.deltal))
        tstar <- soln$x[2]
        Ll <- soln$x[3]
        bs <- 1/(1+exp(soln$x[4]))

        ys <- log((1+bs)/(1-bs))
        ft <- ys*(4*Ne*C)^-1 #originally ys*(2*Ne*C)^-1

        ## get small effect additive genetic variance
        Vas <- 8*Ne*u*Ls*as^2*bs/ys  #originally 4*Ne*u*Ls*as^2*bs/ys
        if(is.null(Ve))
            Ve <- Vas*(1-r2n)/r2n
        norm.sd <- sqrt(Vas+Ve)
    }



    deltas <- as*ft
    mean.nl <- 2*Ll*u/(deltal*C)
    Val <- 2*al^2*mean.nl
    Vt <- Vas + Ve + Val

    ## risk scale variances
    Vos <- Vas*ft^2
    Vol <- 2*deltal^2*mean.nl
    Vrg <- Vos + Vol

    ## bulk stuff
    prev <- as.numeric(pPoisConv(tstar,mean.nl,sqrt(Vas+Ve),alphal=al))
    my.range <- seq(qpois(1e-8,mean.nl),qpois(1-1e-8,mean.nl))
    ldens <- dpois(my.range,mean.nl) ## density on large effect liability
    pdil <- pnorm(tstar, al*my.range,norm.sd,lower.tail=FALSE) ## prev among inds with i large effect alleles
    pidl <- ldens*pdil/prev ## density on number of large effect alleles conditional on having disease


    ## does mut-sel balance actually hold??
    mutl <- 2*Ll*al*u
    muts <- 2*Ls*as*u*bs
    mutall <- muts+mutl

    sell <- mean.nl*al*deltal*C
    sels <- Vas*ft*C
    selall <- sell+sels
    stopifnot(abs((mutall-selall)/(mutall+selall))<1e-8)



    tol <- 1e-5
    min.gl <- uniroot(
        function(X)
            tol-(1-pPoisConv(X,mean.nl,norm.sd,alphal=al)),
        interval=c(-10*tstar,10*tstar)
    )$root
    max.gl <- uniroot(
        function(X)
            tol-pPoisConv(X,mean.nl,norm.sd,alphal=al),
        interval=c(-10*tstar,10*tstar)
    )$root
    seq.li <- seq(min.gl,max.gl,length.out=1000)
    li.dense <- sapply(seq.li,function(G) dPoisConv(G,mean.nl,norm.sd,alphal=al))






    ## heritabilies on liability scale
    h2s <- Vas/(Val+Vas/r2n)
    h2l <- Val/(Val+Vas/r2n)
    h2 <- (Vas+Val)/(Val+Vas/r2n)



    ## liability scale h2 estimation via normal approx
    phit <- dnorm(qnorm(1-prev))
    risk.var <- (prev*(1-prev))

    ## observed scale
    h2os <- Vos/risk.var
    h2ol <- Vol/risk.var
    h2o <- Vrg/risk.var
    prs <- Vos/(Vos+Vol)

    ## naive estimates on liability scale
    h2s.est <- h2os*risk.var/Vt*phit^-2
    h2l.est <- h2ol*risk.var/Vt*phit^-2
    h2all.est <- h2o*risk.var/Vt*phit^-2

    implied.al.norm <- deltal/phit
    implied.al.true <- deltal/ft


    ## other stuff
    L <- Ll+Ls
    meana <- (as*Ls + al*Ll)/L
    maxg <- 2*meana*L
    bt <- (2*Ls*bs*as + 2*Ll*1*al)/(maxg)

    return(
        list(
            ft=ft,
            phit=phit,
            bs=bs,
            as=as,
            al=al,
            meana=meana,
            maxg=maxg,
            bt=bt,
            my.range=my.range,
            deltas=deltas,
            deltal=deltal,
            ss=deltas*C,
            sl=deltal*C,
            ys=4*Ne*deltas*C,  #originally 2*Ne*deltas*C
            yl=4*Ne*deltal*C,  #originally 2*Ne*deltal*C
            Ne=Ne,
            u=u,
            C=C,
            ors=(prev+deltas)/prev*(1-prev)/(1-prev-deltas),
            orl=(prev+deltal)/prev*(1-prev)/(1-prev-deltal),
            tstar=tstar,
            mean.nl=mean.nl,
            prev=prev,
            Vas=Vas,
            Ve=Ve,
            Val=Val,
            Vos=Vos,
            Vol=Vol,
            h2s=h2s,
            h2l=h2l,
            h2=h2,
            h2os=h2os,
            h2ol=h2ol,
            h2o=h2o,
            h2s.est=h2s.est, ## small effect h2 you would infer by naive application of ltm methods
            h2l.est=h2l.est, ## large effect h2 you would infer by naive application of ltm methods
            h2all.est=h2all.est, ## total h2 you would infer by naive application of ltm methods
            pgs.est=h2s.est/h2all.est,
            implied.al.norm=implied.al.norm,
            implied.al.true=implied.al.true,
            pgs=Vas/(Val+Vas),
            pgl=Val/(Val+Vas),
            Vos=Vos,
            Vol=Vol,
            pos=Vos/(Vos+Vol),
            pol=Vol/(Vos+Vol),
            ldens=ldens,
            pdil=pdil,
            pidl=pidl,
            Ls=Ls,
            Ll=Ll,
            li.dense=li.dense
        )
    )
}



