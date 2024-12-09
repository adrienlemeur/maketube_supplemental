{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(scales)
} #start

{
  genotube_perf <- read.table("maketube_genomes_genotube_performance_metrics_by_regions.tsv", sep = "\t", header = F)[,c(1, 2, 5, 7:13)]
  tbprofiler_perf <- read.table("maketube_genomes_TBprofiler_performance_metrics_by_regions.tsv", sep = "\t", header = F)[,c(1, 2, 5, 7:13)]
  all_pipe_perf <- na.omit(rbind(genotube_perf, tbprofiler_perf))
  
  colnames(all_pipe_perf) <- c("pipeline", "strain", "region", "variant_class", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")
  all_pipe_perf <- subset(all_pipe_perf, region != "duplicata_region")

  my_colors = c("firebrick", "chocolate1") ; names(my_colors) <- levels(c("snpmutator", "genotube"))

  structural_variants <- subset(all_pipe_perf,
                                  region %in% c("ANTIBIOTIC_GENE", "DR_scar_site", "duplicata_region", "duplicated_region", "insertion", "IS_flanking_region", "IS_scar_site"),
                                  TP+FN>1
                                )
  structural_variants$region <- factor(structural_variants$region, levels = c("ANTIBIOTIC_GENE", "insertion", "IS_flanking_region", "DR_scar_site", "IS_scar_site", "duplicated_region"), ordered = T)
  
  H37Rv_region <- subset(all_pipe_perf,
                                region %in% c("ANTIBIOTIC_GENE", "BLINDSPOTS", "PEPPE", "REPETITIVE", "RLC"),
                                TP+FN>1
  )
  
  H37Rv_region$region <- factor(H37Rv_region$region, levels = c("ANTIBIOTIC_GENE", "BLINDSPOTS", "PEPPE", "REPETITIVE", "RLC"), ordered = T)

} #importing data

{
  PRECISION_by_SV <- ggplot(subset(structural_variants, variant_class == "snp")) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = PRECISION, fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = PRECISION, fill = pipeline, alpha = 0.5),
               shape = 21,
               position = position_jitterdodge(jitter.width = 0.15)
    ) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(text = element_text(NA),
          strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none"
    ) +
    facet_wrap(~ region, nrow = 1, scales = "free_x") +
    scale_color_manual(name = "source", values = my_colors)

  RECALL_by_SV <- ggplot(subset(structural_variants, variant_class == "snp") ) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = RECALL, fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = RECALL, fill = pipeline, alpha = 0.5),
                shape = 21,
                position = position_jitterdodge(jitter.width = 0.15)
              ) +
    facet_wrap(~ region, nrow = 1, scales = "free_x") +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(text = element_text(NA),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()
    ) +
    scale_color_manual(name = "source", values = my_colors)
} #Structural variants plot

svg("precision_and_recall_SV.svg", width = 10, height = 7)
grid.arrange(PRECISION_by_SV, RECALL_by_SV, nrow = 2, ncol = 1)
dev.off()


{
  PRECISION_by_region <- ggplot(subset(H37Rv_region, variant_class == "snp")) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = PRECISION, fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = PRECISION, fill = pipeline, alpha = 0.5),
               shape = 21,
               position = position_jitterdodge(jitter.width = 0.15)
    ) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(text = element_text(NA),
          strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none"
    ) +
    facet_wrap(~ region, nrow = 1, scales = "free_x") +
    scale_color_manual(name = "source", values = my_colors)
  
  RECALL_by_region <- ggplot(subset(H37Rv_region, variant_class == "snp") ) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = RECALL, fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = RECALL, fill = pipeline, alpha = 0.5),
               shape = 21,
               position = position_jitterdodge(jitter.width = 0.15)
    ) +
    facet_wrap(~ region, nrow = 1, scales = "free_x") +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(text = element_text(NA),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()
    ) +
    scale_color_manual(name = "source", values = my_colors)
}

svg("precision_and_recall_region.svg", width = 10, height = 7)
grid.arrange(PRECISION_by_region, RECALL_by_region, nrow = 2, ncol = 1)
dev.off()
