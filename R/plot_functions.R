volcano_plot_local <- function(d_FC, p_adjust, sig_FC_min, sig_p_value_min, symmetric_x, x_min, x_max, point_size, point_transparency,
                               scale_factor, scale_factor_species_label, squared = TRUE,  axistick =1, highlight_fdr = FALSE){

  if (p_adjust){
    d_FC$p_value_mod <- p.adjust(d_FC$p_value, method = "BH")
    y_lab = " FDR-adj. "}
  else {
    d_FC$p_value_mod <- d_FC$p_value
    y_lab = "Unadj. "}



  d_FC <- d_FC |> mutate(
    below_fdr = case_when(
      FDR <= 0.05 ~ "<= 5%",
      FDR <= 0.1 ~ "<= 10%",
      FDR > 0.1 ~ "> 10%",
      TRUE ~ NA
        ))

  d_FC$below_fdr <- factor(d_FC$below_fdr, c("<= 5%", "<= 10%", "> 10%"))


  d_FC$Significant <- ifelse((d_FC$p_value_mod < sig_p_value_min) & (abs(d_FC$log2FC) > log2(sig_FC_min)), "sign", "ns")

  FC_max_positive <- abs(max(d_FC$log2FC, na.rm = T))
  FC_max_negative <- abs(min(d_FC$log2FC, na.rm = T))
  FC_abs_max <- max(FC_max_negative, FC_max_positive)

  FC_seq <- seq(from = -(floor(FC_abs_max)), to=floor(FC_abs_max), by=1)

  y_max <- -log10(min(d_FC$p_value_mod))
  y_break_max = ceiling(y_max*1.2)
  y_breaks <- seq(0, y_break_max)
  y_labels <- 10^(-seq(0, y_break_max))


  p <- ggplot(d_FC, aes(x = log2FC, y = -log10(p_value_mod))) +
    #ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
    #geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed", size=0.1) +
    #geom_hline(yintercept = -log10(0.01), colour="#990000", linetype="dashed", size=0.1) +
    geom_hline(yintercept = -log10(0.05), colour="#d99898", linetype="dashed", linewidth=0.3) +
    geom_vline(xintercept = -log2(sig_FC_min), colour="#d99898", linetype="dashed", linewidth=0.3) +
    geom_vline(xintercept = log2(sig_FC_min), colour="#d99898", linetype="dashed", linewidth=0.3) +
    geom_vline(xintercept = 0, colour="#8fb0c4", size=.3) +
    geom_point(size = point_size, alpha = point_transparency, na.rm = T, stroke =.5, aes(color = Significant, fill = below_fdr), shape = 21) +
    scale_color_manual(guide = NULL, values = c("ns" = "grey", "sign" = "red3")) +
    scale_fill_manual(guide = NULL, name = "FDR", values = c("<= 5%" = "#f72e1e","<= 10%" = "#fab0aa", "> 10%" = "white"), labels = c(expression(""<=" 5%"), expression(""<="10%"), expression("">"10%"))) +
    #scale_shape_manual()
    theme_bw(base_size = 14) +
    labs(
      y = bquote(.(y_lab) ~ italic('P') ~ ' value (log scale)'),
      x = bquote(log[2] ~ 'FC')
    ) +
    scale_y_continuous(
      limits = c(0, y_max *1.2),
      breaks = y_breaks,
      labels = y_labels) +
    guides(fill = guide_legend(override.aes = list(color = "red"))) +


    #scale_y_continuous(
    #    breaks = c()
    #   breaks = c()
    #   #trans  = scales::compose_trans("log10", "reverse"),
    #   #labels = scales::label_log(base = 10),
    #   #breaks  = log_breaks_125
    # ) +
    #ggplot2::annotation_logticks(sides = 'tb', outside = TRUE, size = .2) +
    #scale_y_continuous(trans = "log1p", limits=c(0,5), breaks = c(1, 2, 3,4,5),  expand = c(0,0)) +
    ggrepel::geom_text_repel(data = d_FC %>% filter(Significant =="sign"),
                             aes(label = Compound),max.overlaps = 20,
                             size = 3 * scale_factor_species_label,
                             box.padding = unit(0.15, "lines"),
                             point.padding = unit(0.1, "lines")) +
    theme(axis.text = element_text(size=6 * scale_factor),
          axis.title = element_text(size=7 * scale_factor,face = "bold"),
          plot.title = element_text(size=8, face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.13,.88),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6, face = "bold"),
          legend.key.size = unit(2, 'mm'),
          legend.spacing.x = unit(0.01, "mm"),
          legend.spacing = unit(c(0,0,0,0), unit="mm"),
          legend.background = element_rect(fill=alpha('white', 0.0)))

  if(symmetric_x) {
    FC_seq <- seq(from = -(ceiling(FC_abs_max)), to=ceiling(FC_abs_max), by=axistick)
    #print(FC_seq)
    p <- p + scale_x_continuous( limits=c(-(FC_abs_max*1.1), FC_abs_max*1.1), breaks = FC_seq,  expand = c(0.02,0.02))
  } else
  {
    min_x <- ifelse(!is.na(x_min), x_min,-(FC_max_negative))
    max_x <- ifelse(!is.na(x_max), x_max,(FC_max_positive))
    FC_seq <- seq(from = min_x, to=max_x, by=axistick)
    #print(FC_seq)
    p <- p + scale_x_continuous(limits=c(min_x*1.1, max_x*1.1), breaks = FC_seq,  expand = c(0.02,0.02), oob = scales::oob_squish)
  }

  # p <- p + theme(
  #   panel.grid.major =  element_line(size=0.3),
  #   panel.grid.minor =  element_line(size=0.15))

  if (squared) p <- p + theme(aspect.ratio=1)

  p
}
