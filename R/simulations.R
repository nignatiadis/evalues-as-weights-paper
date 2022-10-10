#' t-test simulation
#'
#' Small t-test simulation
#'
#' @importFrom genefilter rowttests
#' @importFrom genefilter rowSds
#' @importFrom stats rnorm dchisq pchisq
#' @export
ttest_sim <- function(m, pi0, effect_size, evalue_ncp=2.5, n_samples=10){
  m0 <- ceiling(m*pi0)
  false_nulls <- sample(1:m, m-m0)
  z_table <- matrix(stats::rnorm(n_samples*m), ncol=n_samples)
  z_table[false_nulls,(n_samples/2+1):n_samples] <- matrix(stats::rnorm(n_samples/2*(m-m0), effect_size), ncol=n_samples/2)
  H <- rep(0,m)
  H[false_nulls] <- 1
  ttests <- rowttests(z_table, factor(rep(1:2,each=n_samples/2)))
  sds <- rowSds(z_table)
  chisq_stat <- (n_samples-1)*sds^2
  chisq_dof <- n_samples-1
  chisq_ncp <- n_samples*effect_size^2/4

  filter_pvalue <- 1-stats::pchisq(chisq_stat, chisq_dof)

  evalue <- stats::dchisq(chisq_stat, chisq_dof, ncp = evalue_ncp) / stats::dchisq(chisq_stat, chisq_dof)
  oracle_evalue <- stats::dchisq(chisq_stat, chisq_dof, ncp = chisq_ncp) / stats::dchisq(chisq_stat, chisq_dof)

  simDf <- data.frame(H=H, pvalue= ttests$p.value, sd= sds,
                      filter_pvalue=filter_pvalue,
                      evalue=evalue,
                      oracle_evalue=oracle_evalue)
  simDf
}
