#' Wrapper for the IHW procedure
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Es   Numeric vector of E-values
#' @param alpha   Significance level at which to apply method
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @export
ihw_nmeth_wrapper <- function(Ps, Es, alpha, Storey = FALSE) {
  if (Storey) {
    ihw_nmeth_fit <- IHW::ihw(
      Ps,
      Es,
      alpha,
      lambdas = Inf,
      null_proportion = TRUE,
      null_proportion_level = 0.5
    )
  } else {
    ihw_nmeth_fit <- IHW::ihw(Ps, Es, alpha, lambdas = Inf)
  }
  IHW::rejected_hypotheses(ihw_nmeth_fit)
}


divide_Ps_by_Ws <- function(Ps, Ws) {
  pmin(ifelse(Ws > 0, Ps / Ws, 1), 1)
}


#' The tau-weighted BH multiple testing procedure
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Ws   Numeric vector of multiple testing weights
#' @param tau  Numeric (default = 1), the level at which tau-censoring is applied. Forced to be <= tau_pi0 if Storey=TRUE.
#' @param tau_pi0  Numeric (default = 0.5), the threshold at which Storey is applied
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Vector of adjusted p-values
#' @export
tau_weighted_bh <-
  function(Ps,
           Ws,
           tau = 1,
           tau_pi0 = 0.5,
           Storey = FALSE) {
    if (length(Ws) == 1) {
      Ws <- rep(Ws, length(Ps))
    }
    if (Storey) {
      pi0_hat <-  weighted_storey_pi0(Ps, Ws, tau = tau_pi0)
      Ws <- Ws / pi0_hat
      tau <- min(tau, tau_pi0)
    }
    weighted_pvals <- divide_Ps_by_Ws(Ps, Ws)
    weighted_pvals[Ps > tau] <- Inf
    adj_p <- stats::p.adjust(weighted_pvals, method = "BH")
    adj_p
  }

#' Storey's pi0 estimator for weighted multiple testinng
#'
#' @param pvalues   Numeric vector of unadjusted p-values.
#' @param weights   Numeric vector of multiple testing weights
#' @param tau       Numeric (default = 0.5), the level at which tau-censoring is applied.
#' @param m         Total number of tests (default: `length(pvalues)`)
#'
#' @return Estimated null proportion
#' @export
weighted_storey_pi0 <-
  function(pvalues,
           weights = 1,
           tau = 0.5,
           m = length(pvalues)) {
    w_max <- max(weights)
    num <- w_max + sum(weights * (pvalues > tau))
    num / m / (1 - tau)
  }

sim_combiner <- function(z1, z2, theta) {
  if (theta == 0) {
    z_comb <- z1
  } else if (theta == pi / 2) {
    z_comb <- z2
  } else {
    z_comb <- cos(theta) * z1 + sin(theta) * z2
  }
  stats::pnorm(z_comb)
}

#' Single index modulated multiple testing
#'
#' @param pvalues   Numeric vector of unadjusted p-values.
#' @param pvalues2   Numeric vector of secondary unadjusted p-values.
#' @param alpha   Significance level at which to apply method
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#' @param grid_size Integer (default: 40): How many gridpoints to choose to discretize the interval [0, pi/2] (corresponding to the possible angles between the two statistics)
#' @param details Bool (default: FALSE): Whether to return a detailed list of the output of the multiple testing procedure (e.g., including the optimal angle theta_opt)
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#'
#' @importFrom stats p.adjust pnorm qnorm
#' @export
sim_mtp <-
  function(pvalues,
           pvalues2,
           alpha,
           grid_size = 40,
           details = FALSE,
           Storey = FALSE) {
    is_na_secondary_idx <- is.na(pvalues2)

    comb_grid <- seq(from = 0,
                     to = pi / 2,
                     length.out = grid_size)
    z1 <- stats::qnorm(pvalues)
    z2 <- stats::qnorm(pvalues2)
    rjs <- rep(NA, grid_size)
    for (i in 1:grid_size) {
      theta <- comb_grid[i]
      combined_pv <- sim_combiner(z1, z2, theta)
      combined_pv[is_na_secondary_idx] <- pvalues[is_na_secondary_idx]
      rjs[i] <-
        sum(stats::p.adjust(combined_pv, method = "BH") <= alpha)
    }
    idx_max <- which.max(rjs)

    theta_opt <- comb_grid[idx_max]
    p <- sim_combiner(z1, z2, theta_opt)
    p[is_na_secondary_idx] <- pvalues[is_na_secondary_idx]

    adj_p <- tau_weighted_bh(p, 1, Storey = Storey)

    if (!details) {
      return(adj_p <= alpha)
    } else {
      return(
        list(
          total_rjs = sum(adj_p <= alpha),
          theta_opt = theta_opt,
          rjs = rjs,
          Storey = Storey,
          p = p,
          adj_p = adj_p
        )
      )
    }
  }



#' ep-BH
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Es   Numeric vector of e-values
#' @param alpha   Significance level at which to apply method
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
ep_BH <- function(Ps, Es, alpha, Storey = FALSE) {
  if (Storey) {
    pi0_hat <- weighted_storey_pi0(Ps, 1, tau = 0.5)
    tau_cens <- 0.5
  } else {
    pi0_hat <- 1
    tau_cens <- 1
  }
  Es <- Es / pi0_hat
  adj_p <- tau_weighted_bh(Ps, Es, tau = tau_cens, Storey = FALSE)
  adj_p <= alpha
}


#' Normalized wBH
#'
#' @param Ps   Numeric vector of unadjusted p-values.
#' @param Es   Numeric vector of e-values
#' @param alpha   Significance level at which to apply method
#' @param Storey  Bool (default: FALSE): is the procedure pi0 adaptive or not?
#'
#' @return Binary vector of rejected/non-rejected hypotheses.
#' @export
normalized_wBH <- function(Ps, Es, alpha, Storey = FALSE) {
  normalized_weights <- Es / mean(Es)
  adj_p <- tau_weighted_bh(Ps, normalized_weights, Storey = Storey)
  adj_p <= alpha
}
