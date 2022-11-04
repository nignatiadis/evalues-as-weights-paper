library(testthat)
library(DESeq2)
library(limma)
library(dplyr)

K <- 50
m <- 10
es <- 0.5
gamma <- 1.0
prior_dof <- 3.64
prior_s2 <- 0.0144
prop_signal_microarray <- 0.5
data("meanDispPairs")
mean_disp_pairs <- filter(meanDispPairs, mean >= 1)

expect_error(rnaseq_microarray_sim(K=K, m=5, mean_disp_pairs=mean_disp_pairs))
expect_error(rnaseq_microarray_sim(K=9, m=m, mean_disp_pairs=mean_disp_pairs))

set.seed(1)
rnaseq_sim_df <- rnaseq_microarray_sim(K=K, m=m, mean_disp_pairs=mean_disp_pairs)
expect_false(is.unsorted(-rnaseq_sim_df$mu0))
expect_false(is.unsorted(-rnaseq_sim_df$sigmai))
expect_equal(sum(rnaseq_sim_df$H), 0.2*50 )


rnaseq_sim_df_as_sim <- rnaseq_microarray_sim(K=10000, m=20, mean_disp_pairs=mean_disp_pairs)

rnaseq_sim_df2 <- rnaseq_microarray_sim(K=1000, m=m, mean_disp_pairs=mean_disp_pairs, prop_signal_microarray=1, gamma = 2)

expect_lt(cor.test(rnaseq_sim_df2$pvalue, rnaseq_sim_df2$evalue)$p.value, 0.2)

rnaseq_sim_df3 <- rnaseq_microarray_sim(K=1000, m=m, mean_disp_pairs=mean_disp_pairs, prop_signal_microarray=0)
expect_gt(cor.test(rnaseq_sim_df3$pvalue, rnaseq_sim_df3$evalue)$p.value, 0.2)

# now run the simulation function and then compare to output when we implement it manually

set.seed(1)

mhalf <- m/2
condition <- factor(rep(c("A","B"), each = mhalf))
condition_df <- data.frame(condition=condition)
x <- model.matrix(~condition)

# Count data simulation and analysis via DESeq2
beta <- c(rep(0, K * 8/10), sample(c(-es,es), K * 2/10, TRUE))
H <- beta != 0

sampled_rows <- sample(1:nrow(mean_disp_pairs), K, replace=TRUE)
mu0 <- mean_disp_pairs[sampled_rows,1]
disp <- mean_disp_pairs[sampled_rows,2]
betafull <- cbind(log2(mu0), beta)
mu <- 2^(betafull %*% t(x))

mu_brute <- cbind(mu0,mu0,mu0,mu0,mu0,mu0*2^(beta),mu0*2^(beta),mu0*2^(beta),mu0*2^(beta),mu0*2^(beta))
colnames(mu_brute) <- colnames(mu)
expect_equal( mu, mu_brute)

mat <- matrix(rnbinom(K*m, mu=mu, size=1/disp), ncol=m)
mode(mat) <- "integer"


dds <- DESeqDataSetFromMatrix(mat, condition_df, ~ condition)
dds <- DESeq(dds,quiet=TRUE)
res <- as.data.frame(results(dds))

mu0_o <- order(mu0, decreasing=TRUE) ### sort
expect_false(is.unsorted(-mu[mu0_o]))
expect_equal(mu0[mu0_o], sort(mu0, decreasing=TRUE))

# Replicated microarray data simulation and analysis via DESeq2

sigmai <- sqrt(prior_s2*prior_dof/rchisq(K,df=prior_dof))
sigmai <- sort(sigmai, decreasing=TRUE) ###

y_microarray <- matrix(rnorm(K*m,sd=sigmai),nrow=K)
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
                     H = H[mu0_o])
res_df <- res_df[!is.na(res_df$pvalue), ]

expect_equal(res_df$pvalue, rnaseq_sim_df$pvalue)
expect_equal(res_df$evalue, rnaseq_sim_df$evalue)
expect_equal(res_df$filter_pvalue, rnaseq_sim_df$filter_pvalue)
expect_equal(res_df$H, rnaseq_sim_df$H)

affy_ebayes_lm_fit_prop05 <- eBayes(affy_lm_fit, proportion = 0.5)
expect_equal(setNames(affy_ebayes_lm_fit_prop05$lods[,2],NULL), log(res_df$evalue))

# now check the limma_results_table_with_evalues function
fit <- affy_ebayes_lm_fit
coef <- "conditionB"

fit$coefficients <- as.matrix(fit$coefficients)
rn <- rownames(fit$coefficients)




#	Extract statistics for table
M <- fit$coefficients[,coef]
expect_equal(M, fit$coefficients[,2])

coef_idx <- which(colnames(fit$t) == coef)
expect_equal(coef_idx, 2)
tstat <- as.matrix(fit$t)[,coef]
P.Value <- as.matrix(fit$p.value)[,coef]
B <- as.matrix(fit$lods)[,coef]
data_dof <- fit$df.residual
expect_equal(data_dof, rep(m-2, K))
prior_dof <- fit$df.prior
prior_s2 <- fit$s2.prior
gamma_hat_init <-  fit$var.prior[coef_idx]

refit_eb <- limma::eBayes(fit, proportion=0.5)
expect_equal(refit_eb$t, fit$t)
gamma_hat_prop1 <- refit_eb$var.prior[coef_idx]

tab <- data.frame(t=tstat,P.Value=P.Value, B=B,
                  data_dof = data_dof)
rownames(tab) <- rn

gamma_adjusted <- gamma_hat_prop1 / fit$stdev.unscaled[1,coef_idx]^2
evals <- mapply(limma_evalue, tab$t, tab$data_dof, prior_dof, gamma_adjusted)

expect_equal(integrate(function(t) {limma_evalue(t, tab$data_dof[1], prior_dof, gamma_adjusted)*dt(t, df= tab$data_dof[1]+prior_dof)}, -Inf, +Inf)$value, 1)



refit_eb2 <- limma::eBayes(fit, proportion=0.5)

refit_eb_from_lm <- limma::eBayes(affy_lm_fit, proportion=0.5)
lods2 <- as.matrix(refit_eb2$lods)[,coef]
lods3 <- as.matrix(refit_eb_from_lm$lods)[,coef]
expect_true( max(abs(evals - exp(lods3))) < 1e-13)


tmp_res <- apply_multiple_testing_methods(rnaseq_sim_df_as_sim, 0.05)

bh_rjs <- p.adjust(rnaseq_sim_df_as_sim$pvalue, method="BH") <= 0.05

expect_equal(fdp_eval(rnaseq_sim_df_as_sim$H, bh_rjs), dplyr::select(filter(tmp_res, method=="BH", weights=="Unweighted"), c(rjs, pow, FDP, FWER)))

eval_se<- sd(rnaseq_sim_df_as_sim$evalue[rnaseq_sim_df_as_sim$H == 0])/sqrt(10000)
eval_mean <- mean(rnaseq_sim_df_as_sim$evalue[rnaseq_sim_df_as_sim$H == 0])

expect_true( abs(eval_mean - 1 ) < 4*eval_se)
