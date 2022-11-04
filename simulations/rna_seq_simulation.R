library(doParallel)
registerDoParallel(cores=3)
library(doRNG)
library(tidyverse)
library(EValueWeightsPaper)


alpha <- 0.1
nreps <- 100

data("meanDispPairs")
meanDispPairs_filtered <- filter(meanDispPairs, mean >= 1)


rnaseq_sim_combs <- expand.grid(K = 10000,
                                m = 40,
                                es = seq(0.3, to=0.9, by=0.1),
                                gamma = 0.5,
                                prop_signal_microarray = c(0, 0.5, 1.0),
                                alpha=alpha,
                                prior_dof = 3.64,
                                prior_s2 = 0.0144,
                                seed=1:nreps)


set.seed(10)

r <- foreach(i=1:nrow(rnaseq_sim_combs)) %dorng% {
  print(i);
  sim_df <- rnaseq_microarray_sim(
    K=rnaseq_sim_combs$K[i],
    m=rnaseq_sim_combs$m[i],
    mean_disp_pairs=meanDispPairs_filtered,
    es = rnaseq_sim_combs$es[i],
    gamma = rnaseq_sim_combs$gamma[i],
    prop_signal_microarray = rnaseq_sim_combs$prop_signal_microarray[i],
    prior_dof =  rnaseq_sim_combs$prior_dof[i],
    prior_s2 =  rnaseq_sim_combs$prior_s2[i]
  )
  mtp_eval <- apply_multiple_testing_methods(sim_df, alpha=rnaseq_sim_combs$alpha[i])
  cbind(mtp_eval, rnaseq_sim_combs[i,], row.names = NULL)
}




saveRDS(r, file="rnaseq_microarray_simulation_results.Rds")




