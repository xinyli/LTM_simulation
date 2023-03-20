#source('scripts/solveTwoEffect.R')

recover.flag <- FALSE

## params
Ls <- 1e6
Ne <- 5e3
u <- 1e-7
as <- 1
C <- 0.2
r2n <- 1/2
my.bt <- 0.3
my.als <- exp(seq(log(10),log(500),length.out=200))


last.Ll <- numeric()
last.Ll[1] <- 8000
last.bs <- numeric()
last.bs[1] <- 0.232

my.cols <- c('as','Ls','bs','rhos','al','Ll','bl','rhol','bt','rhot','Ne','u','C','Ve','h2s','h2l','prev')
output <- data.frame(matrix(ncol = length(my.cols), nrow = 0))
colnames(output) <- my.cols

solns <- list()
for( i in seq_along(my.als)){

    solns[[i]] <- solveTwoEffect(
        bs=last.bs[i],
        bt=my.bt,
        Ne=Ne,
        Ls=Ls,
        Ll=last.Ll[i],
        as=as,
        al=my.als[i],
        r2n=r2n,
        Ve=NULL,
        u=u,
        C=C,
        LL.soln=TRUE,
        var.ratio=1,
        equalize.observed.vars=TRUE
    )

    output[i,'as'] <- solns[[i]]$as  ## small effect size
    output[i,'Ls'] <- solns[[i]]$Ls  ## numer of small effect loci
    output[i,'bs'] <- solns[[i]]$bs  ## b for small effect loci
    output[i,'rhos'] <- 1/2-output[i,'bs']/2  ## rho for small effect loci
    output[i,'al'] <- solns[[i]]$al  ## large effect size
    output[i,'Ll'] <- solns[[i]]$Ll  ## numer of large effect loci
    output[i,'bl'] <- solns[[i]]$bs  ## b for large effect loci
    output[i,'rhol'] <- 0            ## rho for large effect loci
    output[i,'bt'] <- solns[[i]]$bt  ## total b
    output[i,'rhot'] <-  1/2 - solns[[i]]$bt/2  ## total rho

    output[i,'Ne'] <- solns[[i]]$Ne  ## population size
    output[i,'u'] <- solns[[i]]$u    ## mutation rate
    output[i,'C'] <- solns[[i]]$C    ## cost of disease
    output[i,'Ve'] <- solns[[i]]$Ve  ## environmental variance
    output[i,'h2s'] <- solns[[i]]$h2s  ## small effect h2 on liability scale
    output[i,'h2l'] <- solns[[i]]$h2l  ## large effect h2 on liability scale
    output[i,'prev'] <- solns[[i]]$prev ## prevalence


    last.Ll[i+1] <- output[i,'Ll']
    last.bs[i+1] <- output[i,'bs']

}
#save(output,file='twoEffectTestRun.Robj')
write.table(output, "params_test_twoeffect_Ls1e6_new.txt", col.names=T, row.names=F, quote=F)
#a=read.table("params_test_twoeffect_Ls1e6.txt")
