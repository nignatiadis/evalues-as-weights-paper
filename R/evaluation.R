#' Run multiple testing procedures
#
#' @param sim_df `data.frame` returned from one of the simulation functions
#' @param alpha Nominal level at which to run the different multiple testing methods.
#' @importFrom dplyr bind_rows mutate
#' @export
apply_multiple_testing_methods <- function(sim_df, alpha){
  Hs <- sim_df$H

  Storey <- FALSE
  unweighted_res <- tau_weighted_bh(sim_df$pvalue, 1, Storey=Storey) <= alpha
  sim_res <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey)
  ep_res <- ep_BH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  normalized_res <- normalized_wBH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  IHW_res <- ihw_nmeth_wrapper(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  fisher_res <- fisher_BH(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey)

  # collect results
  res_BH <- bind_rows(
    mutate( fdp_eval(Hs, unweighted_res), weights="Unweighted"),
    mutate( fdp_eval(Hs, sim_res),   weights="SIM"),
    mutate( fdp_eval(Hs, ep_res), weights="e-values"),
    mutate( fdp_eval(Hs, normalized_res), weights="Norm. e-values"),
    mutate( fdp_eval(Hs, IHW_res),   weights="IHW"),
    mutate( fdp_eval(Hs, fisher_res),   weights="Fisher")
  )

  Storey <- TRUE
  unweighted_res <- tau_weighted_bh(sim_df$pvalue, 1, Storey=Storey) <= alpha
  sim_res <- sim_mtp(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey)
  ep_res <- ep_BH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  normalized_res <- normalized_wBH(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  IHW_res <- ihw_nmeth_wrapper(sim_df$pvalue, sim_df$evalue, alpha, Storey=Storey)
  fisher_res <- fisher_BH(sim_df$pvalue, sim_df$filter_pvalue, alpha, Storey=Storey)

  # collect results
  res_Storey <- bind_rows(mutate( fdp_eval(Hs, unweighted_res), weights="Unweighted"),
                          mutate( fdp_eval(Hs, sim_res),   weights="SIM"),
                          mutate( fdp_eval(Hs, ep_res), weights="e-values"),
                          mutate( fdp_eval(Hs, normalized_res), weights="Norm. e-values"),
                          mutate( fdp_eval(Hs, IHW_res),   weights="IHW"),
                          mutate( fdp_eval(Hs, fisher_res),   weights="Fisher")
  )

  res <- rbind( mutate(res_BH, method="BH"), mutate(res_Storey, method="Storey"))
}

#' Evaluate multiple testing procedure
#
#' @param Hs Vector with indicators of alternatives (1) and true nulls (0)
#' @param rjs Vector with indicator of rejected hypotheses
#' @return Data frame with columns `rjs` (total rejections), `pow` (Power), `FDP` (False discovery proportion)
#' @export
fdp_eval <- function(Hs, rjs){
  rjs_total <- sum(rjs)
  pow <- sum(rjs*Hs)/max(1,sum(Hs))
  FDP <- sum(rjs*(1-Hs))/max(1,rjs_total)
  FWER <- sum( (1-Hs)*rjs) > 0
  data.frame(rjs=rjs_total, pow=pow, FDP=FDP, FWER=FWER)
}
