#####################################################
## First part is same as "Figure2_macsBehaviour.r" ##
## Here we only add one additional plot            ##
#####################################################

## For supplementary Figure 4
peaksAnno.all.df.space$group <- cut(
  peaksAnno.all.df.space$width,
  breaks = c(0,400,700,1200,2500,4000,Inf),
  labels = c(
    "very_narrow\n(0-400)",
    "narrow\n(400-700)",
    "medium\n(700-1200)",
    "broad\n(1200-2500)",
    "very_broad\n(2500-4000)",
    "extreme\n(>4000)"
  ),
  include.lowest = TRUE
)

plot_df <- peaksAnno.all.df.space %>% count(set_num, group)

#line-graph version of plot
p4 <-ggplot(plot_df, aes(x = set_num, y = n, group=group, color = group)) +
  geom_point() +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(
    breaks = sort(unique(plot_df$set_num)),
    labels = scales::label_number( suffix = "M")
  ) +
  scale_y_continuous(
    labels = scales::label_number( suffix = "M",scale = 1e-4)
  )+
  labs(
    x = "Number of reads",
    y = "Number of peaks", title = "B", color = "Peak width\ngroup"
  ) +
  pub_theme +
  theme(
    legend.position = "bottom")

design <- "AAADDD
           AAADDD
           BBBCCC
           BBBCCC"
combined_fig2 <- p1 + p2 + p3+ p4 +  plot_layout(design = design)
ggsave(
  filename = "SupplementaryFig_4_dataset1.pdf",
  plot = combined_fig2,
  width = 15,
  height = 12
)
