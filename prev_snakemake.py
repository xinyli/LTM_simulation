import numpy as np 
import pandas as pd

params_table = pd.read_csv("../params_test_extended.txt", delim_whitespace=True)
#params_table = params_table.iloc[8:10,].reset_index(drop=True)

## global parameter ( doesn't change) 
mu=1e-6
cyc = 100
fitCost = 0.5
sampleInt = 50

## simulation variable 
rep = list(np.arange(0,11))
rhos = np.array(params_table["rho"])
liaSizes = np.array((params_table["target.size"]).astype(int))
N = np.array(params_table["Ne"].astype(int))
#envSD = np.round(np.array(params_table["env.sd"]),3)

rule all:
  input: 
    expand(expand("PopSize{N}_LiaSize{liaSizes}_rho{rhos}_rep{{rep}}.prev",zip, N=N, liaSizes=liaSizes, rhos=rhos), rep=rep),
    #expand("PopSize{N}_LiaSize{liaSizes}_rho{rhos}_all.h2", zip, N=N, liaSizes=liaSizes, rhos=rhos)
    

rule slim_simulate_withsegregating:
  input:
    slim_script="LTM_prev_nucleotide.slim"
  params:
    mu=mu,
    fitCost=fitCost,
    cyc=cyc, 
    sampleInt = sampleInt, 
    time="50:00:00",
    partition="jnovembre",
    mem="5Gb"
  output:
    "PopSize{N}_LiaSize{liaSizes}_rho{rhos}_rep{rep}.prev"
  shell:
    """thr=`awk 'BEGIN {{print 1e5*2*{wildcards.rhos}}}'`;
    env=`Rscript --vanilla solveSingleEffect.R 1e-6 ${{thr}} 0.9 5000 1e5 0.5 | awk '{{print $2}}'`;
    set +u; slim  -d mu={params.mu} -d rho_input={wildcards.rhos} -d p={wildcards.N} -d liaSize={wildcards.liaSizes} -d f={params.fitCost}  -d e=${{env}} -d cyc={params.cyc} -d sampleInt={params.sampleInt} -d rep={wildcards.rep} {input.slim_script} > PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}.temp; set -u;""" 
    #set +u; rm PopSize{wildcards.N}_LiaSize{wildcards.liaSizes}_rho{wildcards.rhos}_rep{wildcards.rep}.temp; set -u"""


rule result_combined: 
   input: 
     h2=expand("PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_rep{rep}.h2", rep=rep),
     prev=expand("PopSize{{N}}_LiaSize{{liaSizes}}_rho{{rhos}}_rep{rep}.prev", rep=rep)
   output:
     h2="PopSize{N}_LiaSize{liaSizes}_rho{rhos}_all.h2",
     prev="PopSize{N}_LiaSize{liaSizes}_rho{rhos}_all.prev"
   shell: 
     """cat {input.h2} >> {output.h2}; cat {input.prev} >> {output.prev}""" 
 
