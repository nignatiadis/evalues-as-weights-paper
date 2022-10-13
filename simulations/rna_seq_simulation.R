library(doParallel)
registerDoParallel(cores=3)
library(doRNG)
library(tidyverse)
library(EValueWeightsPaper)


alpha <- 0.1
nreps <- 100

data("meanDispPairs.RData")
meanDispPairs_filtered <- filter(meanDispPairs, mean >= 1)# top_n(meanDispPairs, 10000, mean)

tmp <- rnaseq_microarray_sim(20000, 40, meanDispPairs_filtered, es=0.4, gamma=0.5, prop_signal_microarray = 1.0)
tmp_res <- apply_multiple_testing_methods(tmp, alpha)

sum(p.adjust(tmp$pvalue, method="BH") <= 0.1)
sum(p.adjust(tmp$filter_pvalue, method="BH") <= 0.1)
#sim_df <- ttest_sim(20000, 0.95, 2, evalue_ncp=10, n_samples=10)
#tmp <- apply_apply_multiple_testing_methods(sim_df, alpha)


rnaseq_sim_combs <- expand.grid(K = 20000,
                                m = 40,
                                es = seq(0.3, to=0.9, by=0.1),
                                gamma = 0.5,
                                prop_signal_microarray = c(0, 0.5, 1.0),
                                alpha=alpha,
                                seed=1:nreps)


set.seed(10)

r <- foreach(i=1:nrow(rnaseq_sim_combs)) %dorng% {
  print(i);
  sim_df <- rnaseq_microarray_sim(
    rnaseq_sim_combs$K[i],
    rnaseq_sim_combs$m[i],
    meanDispPairs_filtered,
    es = rnaseq_sim_combs$es[i],
    gamma = rnaseq_sim_combs$gamma[i],
    prop_signal_microarray = rnaseq_sim_combs$prop_signal_microarray[i]
  )
  mtp_eval <- apply_multiple_testing_methods(sim_df, alpha=rnaseq_sim_combs$alpha[i])
  mutate(cbind(mtp_eval, rnaseq_sim_combs[i,], row.names = NULL), alpha=rnaseq_sim_combs$alpha[i])
}




saveRDS(r, file="rnaseq_microarray_simulation_results.Rds")




