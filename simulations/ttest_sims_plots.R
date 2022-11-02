library(tidyverse)
library(ggrepel)
library(cowplot)

r <- readRDS(file="ttest_simulation_results_with_storey.Rds")


res_r <-  bind_rows(r) %>%
  group_by(weights, method, effect_size, alpha, evalue_ncp, m, pi0) %>%
  dplyr::summarize(FDR = mean(FDP), Power=mean(pow), n=n(), sd_FDP = sd(FDP)) %>%
  arrange(effect_size, desc(Power)) %>%
  ungroup() %>%
  mutate(weights = case_when(
    weights == "e-values" ~ "ep-BH",
    weights == "Norm. e-values" ~ "wBH",
    TRUE ~ weights)
  )


merged_tbl <- res_r %>% pivot_wider(names_from=method, values_from=c(FDR, Power)) %>%
              mutate(power_ratio = Power_Storey/Power_BH, FDR_ratio = FDR_Storey/FDR_BH)


# http://tsitsul.in/blog/coloropt/
methods_list <- bind_rows(
  tibble(weights = "Unweighted", color="#00beff"),
  tibble(weights = "ep-BH",   color="#4053d3"), #4053d3
  tibble(weights = "wBH", color="#ddb310"),
  tibble(weights = "IHW", color="#8c9fb7"),
  tibble(weights = "SIM", color="#ff9287")
) %>%
  mutate(weights=factor(weights, levels=weights))

method_colors <- with(methods_list, setNames(color, weights))


single_panel <- function(res_sub, yaxis, ylabel=yaxis){
  sim_panel <- ggplot(res_sub,
                      aes(x=effect_size, y= !! sym(yaxis), color=weights, shape=weights)) +
    geom_point(alpha=0.75, size=0.8) +
    geom_text_repel(data=dplyr::filter(res_sub,
                                       effect_size == 3),
                    aes(x = effect_size,
                        y = !! sym(yaxis),
                        col = weights,
                        label = weights,
                        segment.square  = TRUE,
                        segment.inflect = TRUE),
                    segment.colour="darkgrey",
                    force = 1,
                    nudge_x           = 0.8,
                    direction         = "y",
                    segment.size      = 0.2,
                    segment.curvature = -0.002,
                    max.overlaps=Inf,
                    min.segment.length = 0) +
    geom_line(size=0.7, alpha=0.75) +
    scale_color_manual(values=method_colors) +
    scale_x_continuous(breaks = seq(1, to=3, by=0.5), limits=c(1.0, 3.6)) +
    ylab(ylabel) +
    xlab(expression(xi)) +
    theme_cowplot() +
    theme(legend.position="none",legend.title=element_blank())
  if (yaxis == "FDR"){
    sim_panel <- sim_panel +
      geom_segment(x=1.0, xend=3.0, y=0.1, yend = 0.1, linetype="dashed", color="black", alpha=0.3)
  }
  sim_panel
}



main_fig_plot <- plot_grid(single_panel(dplyr::filter(res_r, method=="BH"), "FDR", "FDR(BH)"),
                           single_panel(dplyr::filter(res_r, method=="Storey"), "FDR", "FDR(Storey)"),
                           single_panel(merged_tbl, "FDR_ratio", "FDR(Storey)/FDR(BH)"),
                           single_panel(dplyr::filter(res_r, method=="BH"), "Power", "Power(BH)"),
                           single_panel(dplyr::filter(res_r, method=="Storey"), "Power", "Power(Storey)"),
                           single_panel(merged_tbl, "power_ratio", "Power(Storey)/Power(BH)"),
  nrow = 2,
  labels=c("A","B", "C", "D", "E", "F"))


save_plot("ttest_sim_main.pdf", main_fig_plot, base_height=7, base_width=15)


k <- 9
evalue_plot <- function(chisq_stat){
  stats::dchisq(chisq_stat, k, ncp = 10) / stats::dchisq(chisq_stat, k, ncp=0)
}
chisq_stat <- seq(0.001, to=20, length=1000)
#sanity check
integrate(function(s) {evalue_plot(s)*dchisq(s, k, ncp=0)}, 0, 1000)

e_vs_s_plot <- ggplot(tibble(x=chisq_stat, E=sapply(chisq_stat, evalue_plot)), aes(x=x, y=E)) +
  geom_line() +
  xlab(expression(S[k])) +
  ylab(expression(E[k])) +
  geom_vline(xintercept =  qchisq(0.5, k), color="gray", linetype="dashed") +
  theme_cowplot()

e_vs_s_plot
save_plot("E_vs_S.pdf", e_vs_s_plot, base_height=3, base_width=6)

