library(tidyverse)
library(ggrepel)
library(cowplot)


stylized_sim <- function(m = 20,
                         n = 100 - m,
                         N = 10000,
                         sig.ratio = 1,
                         eps = 0.05,
                         b = 0.5) {
  t <- seq(0, b, b / 50)

  computep = function(del, X, Y) {
    #### Method 1: likelihood ratio

    delX = del * sig.ratio
    X = X + delX
    Y = Y + del
    Z = c(X, Y)
    TZ = (sum(X) * sig.ratio + sum(Y)) / sqrt(n * (sig.ratio) ^ 2 + m)
    PZ = 1 - pnorm(TZ)

    #### Method 2: P/E

    EX = exp(sum(X) * delX - n * (delX) ^ 2 / 2)
    PX = 1 - pnorm(sum(X) / sqrt(n))
    PY = 1 - pnorm(sum(Y) / sqrt(m))
    PE = PY / EX


    #### Method 3: Fisher

    PF = 1 - pchisq(-2 * log(PX * PY), 4)

    #### Return the p-values compared with threholds
    return(c(PZ <= eps, PE <= eps, PF <= eps))
  }

  ALLP <- matrix(0, 3, 51)

  for (i in 1:N) {
    X = rnorm(n)
    Y = rnorm(m)
    rej = sapply(t, function(t_single) computep(t_single, X, Y))
    ALLP  <- ALLP + rej
  }

  POW = ALLP / N
  #note code uses switched notation m<->n compared to manuscript
  label = paste0("m=",n,", n=",m)
  bind_rows(
    tibble(delta=t, power=POW[1,], method="Likelihood ratio", n=n, m=m, label=label),
    tibble(delta=t, power=POW[2,], method="P/E (Q-combiner)", n=n, m=m, label=label),
    tibble(delta=t, power=POW[3,], method="Fisher",n=n, m=m, label=label),
  )
}


set.seed(100)
sim_5 <- stylized_sim(m=95)

set.seed(100)
sim_20 <- stylized_sim(m=80)

set.seed(100)
sim_50 <- stylized_sim(m=50)

all_sims <- bind_rows(sim_5, sim_20, sim_50)

# http://tsitsul.in/blog/coloropt/
methods_list <- bind_rows(
  tibble(method = "Likelihood ratio", color="#8c9fb7"),
  tibble(method = "P/E (Q-combiner)",   color="black"), #color="#4053d3"),
  tibble(method = "Fisher", color="#b80058" )
) %>%
  mutate(method=factor(method, levels=method))

method_colors <- with(methods_list, setNames(color, method))

all_sims$method <- factor(all_sims$method, levels=methods_list$method)
all_sims$label <- factor(all_sims$label, levels=c("m=5, n=95", "m=20, n=80", "m=50, n=50"))

sim_plot <- ggplot(all_sims, aes(x=delta, y=power, linetype=method, color=method)) +
  geom_line(alpha=0.8) +
  facet_grid(~label) +
  ylab("Power") +
  xlab(expression(delta)) +
  scale_color_manual(values=method_colors) +
  theme_cowplot() +
  theme(legend.title=element_blank(), legend.position="top", legend.justification = "center")
sim_plot
save_plot("stylized_metaanalysis_power.pdf",sim_plot, base_height=3, base_width=8)





