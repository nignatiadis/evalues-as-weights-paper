library(doParallel)
registerDoParallel(cores=4)
library(doRNG)
library(tidyverse)
library(EValueWeightsPaper)


alpha <- 0.1
nreps <- 2000

ttest_sim_combs <- expand.grid(m = 20000,
                               pi0 = 0.95,
                               effect_size = seq(0.5, to=1.5, by=0.2),
                               var_halfrange = c(0, 0.1, 0.3, 0.5),
                               evalue_ncp=10,
                               L=6,
                               n_samples=10,
                               alpha=alpha,
                               seed=1:nreps)


set.seed(10)

r <- foreach(i=1:nrow(ttest_sim_combs)) %dorng% {
  print(i);
  sim_df <- single_sample_ttest_sim(
    m = ttest_sim_combs$m[i],
    pi0 = ttest_sim_combs$pi0[i],
    var_halfrange = ttest_sim_combs$var_halfrange[i],
    effect_size = ttest_sim_combs$effect_size[i],
    L = ttest_sim_combs$L[i],
    evalue_ncp = ttest_sim_combs$evalue_ncp[i],
    n_samples = ttest_sim_combs$n_samples[i])

  mtp_eval <- apply_multiple_testing_methods(sim_df, alpha=ttest_sim_combs$alpha[i])
  cbind(mtp_eval, ttest_sim_combs[i,])
}

saveRDS(r, file="single_sample_ttest_simulation_results_with_storey.Rds")

