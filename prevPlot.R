library(ggplot2)

## set results directory  
setwd("~/Documents/PhD/Project/Rotation3_Berg/JeremySimulation") 
## load the parameter table for simulation 
load("param.table.Rdata")  
  
rho = c("0.0001","0.00025","0.0005","0.0025","0.005","1e-05","2.5e-05","5e-05")
 
for (r in rho)
{
  index = which(new.table$rho== as.numeric(r))
  e = round(new.table$env.sd[index], 3) ## find the environmental variance in the table 
  temp_prefix = paste("prev_results/set4_nucleotide/PopSize5000_LiaSize100000_rho",r,"_envSD", e, "_all", sep="")
  temp_h2 = read.table(paste(temp_prefix,".h2", sep="")) 
  temp_prev = read.table(paste(temp_prefix,".prev", sep="")) 
  new.table$h2[index] = mean(temp_h2$V1)
  new.table$h2_sd[index] = sd(temp_h2$V1)/ sqrt(dim(temp_h2)[1])
  new.table$prev_emp[index] = mean(temp_prev$V1)
  new.table$prev_sd[index] = sd(temp_prev$V1)/ sqrt(dim(temp_prev)[1])
}
results = as.data.frame(new.table)

ggplot() + geom_line(data= results, aes(x=rho, y =prev, col="predicted")) + 
  geom_line(data=results, aes(x=rho,y=prev_emp, col="empirical")) + geom_point(data=results, aes(x=rho,y=prev_emp)) +
  scale_x_continuous(trans='log10')
  
ggplot() + geom_line(data= results, aes(x=rho, y =h2, col="h2"))  + geom_point(data=results, aes(x=rho,y=h2)) +
  scale_x_continuous(trans='log10')
