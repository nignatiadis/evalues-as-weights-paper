#' t-test simulation
#'
#' Generates data of two sample t-test simulation
#'
#' @param m   Total number of hypotheses
#' @param pi0   Proportion of nulls
#' @param effect_size Mean of alternative hypotheses
#' @param evalue_ncp Noncentrality parameter used for calculation of e-values based on the likelihood ratio test
#' @param n_samples Total number samples used for each t-test. Each group has n_samples/2 samples.
#'
#' @importFrom genefilter rowttests
#' @importFrom genefilter rowSds
#' @importFrom stats rnorm dchisq pchisq
#' @export
ttest_sim <- function(m, pi0, effect_size, evalue_ncp=2.5, n_samples=10){
  stopifnot(n_samples >= 2)
  stopifnot(n_samples %% 2 == 0)
  stopifnot( (pi0 < 1) & (pi0 > 0))

  m0 <- ceiling(m*pi0)
  false_nulls <- sample(1:m, m-m0)
  z_table <- matrix(stats::rnorm(n_samples*m), ncol=n_samples)
  z_table[false_nulls,(n_samples/2+1):n_samples] <- z_table[false_nulls,(n_samples/2+1):n_samples] + effect_size
  H <- rep(0,m)
  H[false_nulls] <- 1
  ttests <- rowttests(z_table, factor(rep(1:2,each=n_samples/2)))
  sds <- rowSds(z_table)
  chisq_stat <- (n_samples-1)*sds^2
  chisq_dof <- n_samples-1

  filter_pvalue <- stats::pchisq(chisq_stat, chisq_dof, lower.tail=FALSE)

  evalue_fun <- function(x){
    stats::dchisq(x, chisq_dof, ncp = evalue_ncp) / stats::dchisq(x, chisq_dof, ncp=0)
  }

  evalue <- sapply(chisq_stat, evalue_fun)

  simDf <- data.frame(H=H, pvalue= ttests$p.value, sd= sds,
                      filter_pvalue=filter_pvalue,
                      evalue=evalue)
  attr(simDf, "evalue_fun") <- evalue_fun
  simDf
}

#' RNA-Seq & Microarray simulation
#'
#' A two-sample comparison based on both sources of information.
#'
#' @param K Integer, number of hypotheses/genes
#' @param m Total number samples (which is assumed to be the same for RNA-Seq and microarray data).
#' Each group of the two-group comparison has m/2 samples.
#' @param prop_signal_microarray Number in [0, 1] (default: 0.5):What proportion
#' of genes differentially expressed in the RNA-Seq synthetic data are also
#' differentially expressed in the synthetic microarray data?
#' @param es Number, binary logarithmic fold change for differentially expressed genes in RNASeq data
#' @param gamma Number, determines strength of effect sizes for alternative genes (for synthetic microarray data)
#' @param prior_df Number, Degrees of freedom for the prior of the inverse variances (for synthetic microarray data)
#' @param prior_s2 Number, Prior variance
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma lmFit eBayes
#' @importFrom stats rnorm rnbinom model.matrix rchisq rbinom
#' @export
rnaseq_microarray_sim <- function(K, m, mean_disp_pairs, prop_signal_microarray=0.5,
                                  es=0.5, gamma=1.0, prior_dof = 3.64, prior_s2 = 0.0144){

  stopifnot(m >= 4)
  stopifnot(m %% 2 == 0)
  stopifnot(K >= 20)

  K_alt <- floor(K * 8/10)

  mhalf <- m/2
  condition <- factor(rep(c("A","B"), each = mhalf))
  condition_df <- data.frame(condition=condition)
  x <- model.matrix(~condition)

  # Count data simulation and analysis via DESeq2
  beta <- c(rep(0, K_alt), sample(c(-es,es), K - K_alt, TRUE))
  H <- beta != 0
  sampled_rows <- sample(1:nrow(mean_disp_pairs), K, replace=TRUE)
  mu0 <- mean_disp_pairs[sampled_rows,1]
  disp <- mean_disp_pairs[sampled_rows,2]
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  mat <- matrix(rnbinom(K*m, mu=mu, size=1/disp), ncol=m)
  mode(mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(mat, condition_df, ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- as.data.frame(results(dds))

  mu0_o <- order(mu0, decreasing=TRUE) ### sort

  # Replicated microarray data simulation and analysis via DESeq2

  sigmai <- sqrt(prior_s2*prior_dof/rchisq(K,df=prior_dof))
  sigmai <- sort(sigmai, decreasing=TRUE) ###

  y_microarray <- matrix(rnorm(K*m,sd=sigmai), nrow=K)
  rownames(y_microarray) <- paste("Gene",1:nrow(y_microarray))

  n_signal_microarray <- floor(prop_signal_microarray*K*2/10)
  if (n_signal_microarray > 0){
    alternatives <- H[mu0_o] & (rbinom(K, 1, prob=prop_signal_microarray) == 1)
    y_microarray[alternatives,(mhalf+1):m] <-  y_microarray[alternatives,(mhalf+1):m] + rnorm(sum(alternatives), sd=sqrt(gamma)*sigmai[alternatives])
  }

  affy_lm_fit <- lmFit(y_microarray, x)
  affy_ebayes_lm_fit <- eBayes(affy_lm_fit)
  affy_res <- limma_results_table_with_evalues(affy_ebayes_lm_fit, coef="conditionB")

  res_df <- data.frame(pvalue = res$pvalue[mu0_o],
                       evalue = affy_res$evalue,
                       filter_pvalue = affy_res$P.Value,
                       H = H[mu0_o],
                       mu0 = mu0[mu0_o],
                       sigmai = sigmai)
  res_df[!is.na(res_df$pvalue), ]
}
