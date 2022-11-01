library(testthat)

ntest1 <- 500000
set.seed(1)
sim_6_df <- ttest_sim(ntest1, 0.9, 3.0, evalue_ncp=1.5, n_samples=6)
ev_6_fun <- attr(sim_6_df, "evalue_fun")


approx_evalue <- integrate(function(x) {ev_6_fun(x)*dchisq(x, df=5)},0.000001, 1000)

expect_equal(approx_evalue$value, 1)

expect_lt(abs(mean(sim_6_df[sim_6_df$H==0,]$evalue) - 1), 6*0.5/sqrt(500000))
expect_equal(sum(sim_6_df$H), 0.1*ntest1)

all_rejected <- rep(1, ntest1)
all_rejected_eval <- fdp_eval(sim_6_df$H, all_rejected)

expect_equal(all_rejected_eval$rjs, ntest1)
expect_equal(all_rejected_eval$pow, 1)
expect_equal(all_rejected_eval$FDP, 0.9)

none_rejected <- rep(0, ntest1)
none_rejected_eval <- fdp_eval(sim_6_df$H, none_rejected)
expect_equal(none_rejected_eval$rjs, 0)
expect_equal(none_rejected_eval$pow, 0)
expect_equal(none_rejected_eval$FDP, 0)

perfect_eval <- fdp_eval(sim_6_df$H, sim_6_df$H)
expect_equal(perfect_eval$rjs, 0.1*ntest1)
expect_equal(perfect_eval$pow, 1)
expect_equal(perfect_eval$FDP, 0)

first_half_rejected_eval <- fdp_eval(sim_6_df$H, rep(c(1,0), each=ntest1/2))
expect_equal(first_half_rejected_eval$rjs, ntest1/2)
expect_equal(first_half_rejected_eval$pow,  sum(sim_6_df$H[1:(ntest1/2)])/sum(sim_6_df$H))
expect_equal(first_half_rejected_eval$FDP,  sum((1-sim_6_df$H[1:(ntest1/2)]))/(ntest1/2))




alpha <- 0.1
sim_df <- ttest_sim(20000, 0.95, 2, evalue_ncp=10, n_samples=10)

expect_equal(sum(1-sim_df$H), 0.95*20000)

Storey <- FALSE
# Check unweighted BH
unweighted_res <- tau_weighted_bh(sim_df$pvalue, 1, Storey=Storey) <= alpha
expect_equal(unweighted_res, p.adjust(sim_df$pvalue,  method="BH") <= alpha)

# Check SIM BH
sim_res <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey)
sim_res_detail <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey, detail=TRUE)
expect_equal(sim_res_detail$total_rjs, max(sim_res_detail$rjs))

sim_res_pvals2unity <- sim_mtp(sim_df$pvalue, rep(1, 20000), alpha, Storey=Storey)
sim_res_pvals2unity_detail <- sim_mtp(sim_df$pvalue, rep(1, 20000), alpha, Storey=Storey, detail=TRUE)
expect_equal(sim_res_pvals2unity_detail$theta_opt,  0)
expect_equal(sim_res_pvals2unity_detail$p,  sim_df$pvalue)
expect_equal(sim_res_pvals2unity, p.adjust(sim_df$pvalue, method="BH") <= alpha)

sim_res_pvals1unity <- sim_mtp(rep(1, 20000), sim_df$pvalue, alpha, Storey=Storey)
sim_res_pvals1unity_detail <- sim_mtp(rep(1, 20000), sim_df$pvalue, alpha, Storey=Storey, detail=TRUE)
expect_equal(sim_res_pvals1unity_detail$theta_opt,  pi/2)
expect_equal(sim_res_pvals1unity_detail$p,  sim_df$pvalue)
expect_equal(sim_res_pvals1unity, p.adjust(sim_df$pvalue, method="BH") <= alpha)


# Check EP BH
ep_res <- ep_BH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
expect_equal(ep_res,  p.adjust(sim_df$pvalue / sim_df$evalue,  method="BH") <= alpha)

# Check WBH
normalized_res <- normalized_wBH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
ws_normalized <- sim_df$evalue / mean(sim_df$evalue)
expect_equal(normalized_res,  p.adjust(sim_df$pvalue / ws_normalized,  method="BH") <= alpha)

# Check IHW
library(IHW)
IHW_res <- ihw_nmeth_wrapper(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
IHW_res_detail <- ihw(sim_df$pvalue, sim_df$evalue, alpha, lambdas=Inf)
expect_equal(IHW_res, adj_pvalues(IHW_res_detail) <= alpha)

# Check unweighted Storey
unweighted_res_storey <- tau_weighted_bh(sim_df$pvalue, 1, Storey=TRUE) <= alpha
pi0 <- (1+sum(sim_df$pvalue > 0.5))*2/20000
expect_equal(unweighted_res_storey, p.adjust(sim_df$pvalue, method="BH") <=   alpha/pi0   )

