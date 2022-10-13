library(tidyverse)
library(ggrepel)
library(cowplot)

r <- readRDS(file="rnaseq_microarray_simulation_results.Rds")


res_r <-  bind_rows(r) %>%
  group_by(weights, method, es, prop_signal_microarray) %>%
  dplyr::summarize(FDR = mean(FDP), Power=mean(pow)) %>%
  ungroup() %>%
  mutate(prop_signal_microarray = factor(prop_signal_microarray, levels=c(1.0, 0.5, 0.0))) %>%
  mutate(weights = case_when(
    weights == "e-values" ~ "ep-BH",
    weights == "Norm. e-values" ~ "wBH",
    TRUE ~ weights)
  )

# http://tsitsul.in/blog/coloropt/
methods_list <- bind_rows(
  tibble(weights = "Unweighted", color="#00beff"),
  tibble(weights = "ep-BH",   color="#4053d3"), #4053d3
  tibble(weights = "wBH", color="#ddb310"),
  tibble(weights = "IHW", color="#8c9fb7"),
  tibble(weights = "SIM", color="#ff9287")  #color="#00bbad")
) %>%
  mutate(weights=factor(weights, levels=weights))

method_colors <- with(methods_list, setNames(color, weights))


single_panel <- function(res_sub, yaxis, ylabel=yaxis){
  sim_panel <- ggplot(res_sub,
                      aes(x=es, y= !! sym(yaxis), color=weights, shape=weights)) +
    geom_point(alpha=0.75, size=0.8) +
    geom_text_repel(data=dplyr::filter(res_sub,
                                       es == 0.9),
                    aes(x = es,
                        y=!! sym(yaxis),
                        col=weights,
                        label=weights,
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
    scale_x_continuous(breaks = seq(0.3, to=0.9, by=0.2), limits=c(0.3, 1.1)) +
    facet_grid(.~prop_signal_microarray, labeller=label_bquote(cols = pi[M]*' = '*.(as.character(prop_signal_microarray))))   +
    ylab(ylabel) +
    xlab(expression(xi)) +
    theme_cowplot() +
    theme(legend.position="none",legend.title=element_blank())
  if (yaxis == "FDR"){
    sim_panel <- sim_panel +
      geom_segment(x=0.3, xend=0.9, y=0.1, yend = 0.1, linetype= "dashed", color="black", alpha=0.3)
  }
  sim_panel
}


power_bh_methods <- single_panel(filter(res_r, method=="BH"), "Power")
fdr_bh_methods <- single_panel(filter(res_r, method=="BH"), "FDR")

save_plot("rnaseq_sim_main.pdf", power_bh_methods, base_height=4, base_width=15)
save_plot("rnaseq_sim_main.pdf", fdr_bh_methods, base_height=4, base_width=15)


rna_seq_main_fig_plot <- plot_grid(single_panel(filter(res_r, method=="BH"), "FDR"),
                                   single_panel(filter(res_r, method=="BH"), "Power"),
                                   nrow = 2,
                                   labels=c("A","B"))
save_plot("rnaseq_sim_main.pdf", rna_seq_main_fig_plot, base_height=8, base_width=15)


rna_seq_suppl_fig_plot <- plot_grid(single_panel(filter(res_r, method=="Storey"), "FDR"),
                                   single_panel(filter(res_r, method=="Storey"), "Power"),
                                   nrow = 2,
                                   labels=c("A","B"))
save_plot("rnaseq_sim_suppl.pdf", rna_seq_suppl_fig_plot, base_height=8, base_width=15)





