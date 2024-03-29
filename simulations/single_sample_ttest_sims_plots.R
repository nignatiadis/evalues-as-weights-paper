library(tidyverse)
library(ggrepel)
library(cowplot)

r <- readRDS(file="single_sample_ttest_simulation_results.Rds")


res_r_all <-  bind_rows(r) %>%
  group_by(weights, method, effect_size, alpha, evalue_ncp, L, var_halfrange, m, pi0) %>%
  dplyr::summarize(FDR = mean(FDP), Power=mean(pow), n=n()) %>%
  arrange(effect_size, desc(Power)) %>%
  ungroup() %>%
  mutate(weights = case_when(
    weights == "e-values" ~ "ep-BH",
    weights == "Norm. e-values" ~ "wBH",
    TRUE ~ weights)
  )

res_r <- filter(res_r_all, var_halfrange == 0)

merged_tbl <- res_r %>% pivot_wider(names_from=method, values_from=c(FDR, Power)) %>%
  mutate(power_ratio = Power_Storey/Power_BH, FDR_ratio = FDR_Storey/FDR_BH)


# http://tsitsul.in/blog/coloropt/
methods_list <- bind_rows(
  tibble(weights = "Unweighted", color="#00beff"),
  tibble(weights = "ep-BH",   color="#4053d3"),
  tibble(weights = "wBH", color="#ddb310"),
  tibble(weights = "IHW", color="#8c9fb7"),
  tibble(weights = "SIM", color="#ff9287"),
  tibble(weights = "Fisher", color="#b80058" )
) %>%
  mutate(weights=factor(weights, levels=weights))

method_colors <- with(methods_list, setNames(color, weights))


single_panel <- function(res_sub, yaxis, ylabel=yaxis, ylim_upper_fdr=1){
  sim_panel <- ggplot(res_sub,
                      aes(x=effect_size, y= !! sym(yaxis), color=weights, shape=weights)) +
    geom_point(alpha=0.75, size=0.8) +
    geom_text_repel(data=dplyr::filter(res_sub,
                                       effect_size == 1.5),
                    aes(x = effect_size,
                        y = !! sym(yaxis),
                        col = weights,
                        label = weights,
                        segment.square  = TRUE,
                        segment.inflect = TRUE),
                    segment.colour="darkgrey",
                    force = 1,
                    nudge_x           = 0.3,
                    direction         = "y",
                    segment.size      = 0.1,
                    segment.curvature = -0.001,
                    max.overlaps=Inf,
                    min.segment.length = 0) +
    geom_line(linewidth=0.7, alpha=0.75) +
    scale_color_manual(values=method_colors) +
    scale_x_continuous(breaks = seq(0.5, to=1.5, by=0.25), limits=c(0.5, 1.9)) +
    ylab(ylabel) +
    xlab(expression(xi)) +
    theme_cowplot() +
    theme(legend.position="none",legend.title=element_blank())
  if (yaxis == "FDR"){
    sim_panel <- sim_panel +
      geom_segment(x=0.5, xend=1.5, y=0.1, yend = 0.1, linetype="dashed", color="black", alpha=0.3) +
      ylim(0, ylim_upper_fdr)
  }
  sim_panel
}



main_fig_plot <- plot_grid(single_panel(dplyr::filter(res_r, method=="BH"), "FDR", "FDR(BH)", ylim_upper=0.161),
                           single_panel(dplyr::filter(res_r, method=="Storey"), "FDR", "FDR(Storey)", ylim_upper=0.161),
                           single_panel(merged_tbl, "FDR_ratio", "FDR(Storey)/FDR(BH)"),
                           single_panel(dplyr::filter(res_r, method=="BH"), "Power", "Power(BH)"),
                           single_panel(dplyr::filter(res_r, method=="Storey"), "Power", "Power(Storey)"),
                           single_panel(merged_tbl, "power_ratio", "Power(Storey)/Power(BH)"),
                           nrow = 2,
                           labels=c("A","C", "E", "B", "D", "F"))

save_plot("ttest_sim_main.pdf", main_fig_plot, base_height=7, base_width=15)


res_var <- filter(res_r_all, var_halfrange > 0)
second_fig <- plot_grid(
single_panel(dplyr::filter(res_var, method=="BH"),  "FDR", "FDR(BH)") + facet_grid(.~var_halfrange, labeller = label_bquote(cols=tau==.(var_halfrange))),
single_panel(dplyr::filter(res_var, method=="BH"),  "Power", "Power(BH)") + facet_grid(.~var_halfrange, labeller = label_bquote(cols=tau==.(var_halfrange))),
nrow=2,
labels = c("A", "B"))

save_plot("ttest_sim_heterogeneous.pdf", second_fig, base_height=7, base_width=15)
