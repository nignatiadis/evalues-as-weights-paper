library(testthat)
library(dplyr)
ntest1 <- 500000
set.seed(1)
sim_6_df <- ttest_sim(ntest1, 0.9, 3.0, evalue_ncp=1.5, n_samples=6)
ev_6_fun <- attr(sim_6_df, "evalue_fun")

check_evalue <- dchisq(5*sim_6_df$sd^2, df=5, ncp=1.5)/dchisq(5*sim_6_df$sd^2, df=5, ncp=0)
expect_equal(check_evalue, sim_6_df$evalue)


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


# Check SIM Storey
sim_res_storey <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=TRUE)
sim_res_storey_detail <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=TRUE, detail=TRUE)
expect_gt  (sim_res_storey_detail$total_rjs, max(sim_res_storey_detail$rjs))
sim_res_storey_pi0 <- (1+sum(sim_res_storey_detail$p > 0.5))*2/20000
expect_equal(p.adjust(sim_res_storey_detail$p, method="BH") <= alpha/sim_res_storey_pi0,  sim_res_storey)

sim_res_storey_pvals2unity <- sim_mtp(sim_df$pvalue, rep(1, 20000), alpha, Storey=TRUE)
sim_res_storey_pvals2unity_detail <- sim_mtp(sim_df$pvalue, rep(1, 20000), alpha, Storey=TRUE, detail=TRUE)
expect_equal(sim_res_storey_pvals2unity_detail$theta_opt,  0)
expect_equal(sim_res_storey_pvals2unity_detail$p,  sim_df$pvalue)
expect_equal(sim_res_storey_pvals2unity, p.adjust(sim_df$pvalue, method="BH") <= alpha / pi0  )


sim_res_storey_pvals1unity <- sim_mtp(rep(1, 20000), sim_df$pvalue, alpha, Storey=TRUE)
sim_res_storey_pvals1unity_detail <- sim_mtp(rep(1, 20000), sim_df$pvalue, alpha, Storey=TRUE, detail=TRUE)
expect_equal(sim_res_storey_pvals1unity_detail$theta_opt,  pi/2)
expect_equal(sim_res_storey_pvals1unity_detail$p,  sim_df$pvalue)
expect_equal(sim_res_storey_pvals1unity, p.adjust(sim_df$pvalue, method="BH") <= alpha / pi0 )


# Check EP Storey
ep_res_storey <- ep_BH(sim_df$pvalue, sim_df$evalue, alpha, Storey=TRUE)
ep_vals <- sim_df$pvalue / sim_df$evalue * pi0
ep_adjp_naive_rjs <- p.adjust(ep_vals, method="BH") <= alpha

ep_vals[sim_df$pvalue > 0.5] <- Inf
ep_adjp_rjs <- p.adjust(ep_vals, method="BH") <= alpha
expect_equal(ep_adjp_rjs,  ep_res_storey)

# Check Storey
normalized_res_storey <- normalized_wBH(sim_df$pvalue, sim_df$evalue, alpha, Storey=TRUE)
ws_normalized_storey <- sim_df$evalue / mean(sim_df$evalue)
pi0_weighted <- 2*(max(ws_normalized_storey) + sum(ws_normalized_storey * (sim_df$pvalue > 0.5 )))/20000
wpbh_adj_naive_rjs <- p.adjust(sim_df$pvalue/ws_normalized_storey, method="BH") <= alpha/pi0_weighted
expect_lt( max(sim_df$pvalue[wpbh_adj_naive_rjs]), 0.5)
expect_equal(wpbh_adj_naive_rjs,  normalized_res_storey)

# Check Storey
IHW_res_storey <- ihw_nmeth_wrapper(sim_df$pvalue, sim_df$evalue, alpha, Storey=TRUE)
IHW_res_storey_detail <- ihw(sim_df$pvalue, sim_df$evalue, alpha, lambdas=Inf, null_proportion=TRUE)
expect_equal(IHW_res_storey, adj_pvalues(IHW_res_storey_detail) <= alpha)

all_methods_res <- apply_multiple_testing_methods(sim_df, alpha)

expect_equal(fdp_eval(sim_df$H, unweighted_res), select(filter(all_methods_res, method=="BH", weights=="Unweighted"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, unweighted_res_storey), select(filter(all_methods_res, method=="Storey", weights=="Unweighted"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, IHW_res), select(filter(all_methods_res, method=="BH", weights=="IHW"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, IHW_res_storey), select(filter(all_methods_res, method=="Storey", weights=="IHW"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, normalized_res), select(filter(all_methods_res, method=="BH", weights=="Norm. e-values"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, normalized_res_storey), select(filter(all_methods_res, method=="Storey", weights=="Norm. e-values"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, ep_res), select(filter(all_methods_res, method=="BH", weights=="e-values"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, ep_res_storey), select(filter(all_methods_res, method=="Storey", weights=="e-values"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, sim_res), select(filter(all_methods_res, method=="BH", weights=="SIM"), c(rjs, pow, FDP, FWER)))
expect_equal(fdp_eval(sim_df$H, sim_res_storey), select(filter(all_methods_res, method=="Storey", weights=="SIM"), c(rjs, pow, FDP, FWER)))
