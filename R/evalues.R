#' LIMMA (Replicated Microarray) E-value
#'
#' @param tstat  Moderated t-statistic.
#' @param d   Likelihood degrees of freedom.
#' @param d0    Prior degrees of freedom.
#' @param gamma   Prior effect size variance scaling
#'
#' @return  E-value
#'
#' @export
limma_evalue <- function( tstat, d, d0, gamma){
  ratio <- d0*tstat^2/(gamma+d0)/(d + d0 + tstat^2)
  exponent <- - (d + d0 + 1)/2
  sqrt(gamma)/sqrt(gamma + d0)* (1-ratio)^exponent
}


#' Results table for LIMMA analysis that also includes e-values
#'
#' @param fit  Object of class `MArrayLM` on which the `eBayes` method has already been applied to.
#'             Typically this will be the output of the function `lmFit` followed by `eBayes`.
#' @param coef   column name specifying which coefficient or contrast of the linear model is of interest
#'
#' @return  A data.frame with the following columns:
#'    * t: moderated t-statistic
#'    * P.Value
#'    * B: the log odds ratio
#'    * evalue
#'
#' @importFrom limma eBayes
#' @export
limma_results_table_with_evalues <- function(fit, coef){
    coef_idx <- which(colnames(fit$t) == coef)

    fit$coefficients <- as.matrix(fit$coefficients)
    rn <- rownames(fit$coefficients)

    #	Check coef is length 1
    if(length(coef)!=1) {
      coef <- coef[1]
      warning("Only single coefficients handled")
    }


    #	Extract statistics for table
    M <- fit$coefficients[,coef]

    ebcols <- c("t","p.value","lods")


    tstat <- as.matrix(fit$t)[,coef]
    P.Value <- as.matrix(fit$p.value)[,coef]
    B <- as.matrix(fit$lods)[,coef]
    data_dof <- fit$df.residual
    prior_dof <- fit$df.prior
    gamma_hat_init <-  fit$var.prior[coef_idx]

    refit_eb <- limma::eBayes(fit, proportion=1)
    gamma_hat_prop1 <- refit_eb$var.prior[coef_idx]

    tab <- data.frame(t=tstat,P.Value=P.Value, B=B,
                      data_dof = data_dof)
    rownames(tab) <- rn

    evals <- mapply(limma_evalue, tab$t, tab$data_dof, prior_dof, gamma_hat_prop1)
    tab$evalue <- evals

    attr(tab, "prior_dof") <- prior_dof
    attr(tab, "gamma_hat_init") <- gamma_hat_init
    attr(tab, "gamma_hat_prop1") <- gamma_hat_prop1

    tab
}
