library(ggplot2)

## set results directory  
setwd("~/Documents/PhD/Project/Rotation3_Berg/JeremySimulation/prev_results/set4_nucleotide/") 
## load the parameter table for simulation or input the table 
load("../../param.table.Rdata")  
  
## the input is the parameter table used for simulation 
## rho is the all the rho range simulated 
merge_into_paramtable <- function(params_table, rho){
  for (r in rho)
  {
    index = which(params_table$rho== as.numeric(r))
    e = round(params_table$env.sd[index], 3) ## find the environmental variance in the table 
    temp_prefix = paste("PopSize5000_LiaSize100000_rho",r,"_envSD", e, "_all", sep="")
    temp_h2 = read.table(paste(temp_prefix,".h2", sep="")) 
    temp_prev = read.table(paste(temp_prefix,".prev", sep="")) 
    params_table$h2[index] = mean(temp_h2$V1)
    params_table$h2_sd[index] = sd(temp_h2$V1)/ sqrt(dim(temp_h2)[1])
    params_table$prev_emp[index] = mean(temp_prev$V1)
    params_table$prev_sd[index] = sd(temp_prev$V1)/ sqrt(dim(temp_prev)[1])
  }
  results = as.data.frame(params_table)
  return(results)
}

rho = c("0.0001","0.00025","0.0005","0.0025","0.005","1e-05","2.5e-05","5e-05")
results = merge_into_paramtable(new.table, rho)
ggplot() + geom_line(data= results, aes(x=rho, y =prev, col="predicted")) + 
  geom_line(data=results, aes(x=rho,y=prev_emp, col="empirical")) + geom_point(data=results, aes(x=rho,y=prev_emp)) +
  scale_x_continuous(trans='log10')
  
ggplot() + geom_line(data= results, aes(x=rho, y =h2, col="h2"))  + geom_point(data=results, aes(x=rho,y=h2)) +
  scale_x_continuous(trans='log10')
