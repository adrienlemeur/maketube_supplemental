{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(patchwork)
  library(gghalves)
  library(ggpp)
  library(ggpubr)

  all_strains_info <- read.table("a_database_to_rule_them_all.tsv", header = F, sep = "\t")[, 1:2]
  colnames(all_strains_info) <- c("strain", "source")
  all_strains_info$source <- factor(all_strains_info$source, ordered = T)
  my_colors = c(maketube = "firebrick", snpmutator = "chocolate1")
  
  custom_facet_labels <- c("nucmer" = "Nucmer", "PAF" = "Minimap2")
  custom_facet_labeller <- as_labeller(custom_facet_labels)

  my_PCH = c(1, 2, 4, 3)
  names(my_PCH) <- levels(all_strains_info$source)
}

maketube_results <- na.omit(read.table("nucmer_minimap2_results.tsv", sep = "\t", header = T))
subset(maketube_results)
colnames(maketube_results) <- c("sample", "pipeline", "structural_variant", "N", "TOTAL", "TOTAL2", "VARIANT_CLASS", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")

maketube_results_wo_SV <- na.omit(read.table("nucmer_paf_perfect2.tsv", sep = "\t", header = T))
colnames(maketube_results_wo_SV) <- c("sample", "pipeline", "FILTER", "TOTAL", "TOTAL2", "VARIANT_CLASS", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")


SNP <- subset(maketube_results,
                   (pipeline == "maketube" |
                      pipeline == "snpmutator") &
                     VARIANT_CLASS == "snp" &
                     PRECISION != 0 & RECALL != 0)


SNP_noSV <- subset(maketube_results_wo_SV,
                          (pipeline == "maketube" |
                          pipeline == "snpmutator") &
                          VARIANT_CLASS == "snp" &
                          PRECISION != 0 & RECALL != 0)

{
PRECISION_noSV <- ggplot(SNP_noSV, aes(pipeline, PRECISION, fill = pipeline)) +
  geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.01), linetype = "dotted") +
  geom_half_boxplot(center = TRUE, errorbar.draw = FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
  geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
  geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
                     ) +
  scale_fill_manual(name = "source", values = my_colors) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
        ) + ylab("PRECISION")

RECALL_noSV <- ggplot(SNP_noSV, aes(pipeline, RECALL, fill = pipeline)) +
  geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.01), linetype = "dotted") +
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
  geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
  geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
  theme_classic() +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
        ) + ylab("RECALL")
}
PRECISION_noSV / RECALL_noSV

{
PRECISION <- ggplot(SNP, aes(pipeline, PRECISION, fill = pipeline)) +
  geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.01), linetype = "dotted") +
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
  geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
  geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  scale_fill_manual(name = "source", values = my_colors) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
  ) + ylab("PRECISION")

RECALL <- ggplot(SNP, aes(pipeline, RECALL, fill = pipeline)) +
  geom_abline(slope = 0, intercept = seq(0.75, 1, by = 0.01), linetype = "dotted") +
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
  geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
  geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
  theme_classic() +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
  ) + ylab("RECALL")
}
PRECISION / RECALL

summary(subset(SNP, pipeline == "maketube" & FILTER == 'nucmer')$RECALL)
summary(subset(SNP, pipeline == "maketube" & FILTER == 'PAF')$RECALL)

summary(subset(SNP_noSV, pipeline == "maketube" & FILTER == 'nucmer')$RECALL)
summary(subset(SNP_noSV, pipeline == "maketube" & FILTER == 'PAF')$RECALL)
svg("fig4_prec_recall_nucmer_paf.svg")
grid.arrange(PRECISION, RECALL, nrow = 2, ncol = 1)
dev.off()

svg("S_fig4_prec_recall_nucmer_paf_no_duplication.svg")
grid.arrange(PRECISION_noSV, RECALL_noSV, nrow = 2, ncol = 1)
dev.off()

SNP



library(gghalves)
library(ggpp)

n = 0.02
ggplot(SNP, aes(pipeline, as.numeric(PRECISION), fill = pipeline)) +
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, 
                    width=0.5, nudge=0.02) +
  geom_half_violin(side="r", nudge=0.02) +
  geom_point(
             shape = 21,
             position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)
  ) +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x")



{
  PRECISION <- ggplot(SNP, aes(pipeline, as.numeric(PRECISION), fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.01), linetype = "dotted") +
    geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02) +
    geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
    geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
    facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x") +
    theme_classic() +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                       limits = c(0.90,1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none"
    ) + ylab("PRECISION")

  RECALL <- ggplot(SNP, aes(pipeline, as.numeric(PRECISION), fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.9, 1, by = 0.01), linetype = "dotted") +
    geom_boxplot(aes(pipeline, RECALL, fill = pipeline), show.legend = FALSE, outlier.shape = NA) + 
    geom_point(aes(x = pipeline, y = RECALL, fill = pipeline, alpha = 0.5),
               shape = 21,
               position = position_jitterdodge(jitter.width = 0.15)
    ) +
    facet_wrap(~ FILTER, nrow = 1) +
    theme_classic() +
    scale_fill_manual(name = "source", values = my_colors) +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                       limits = c(0.90, 1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "none"
    ) + ylab("RECALL") + ylim(c(0.90, 1))
  }



maketube_results <- na.omit(read.table("paf.shorter_alignement.tsv", sep = "\t", header = F))
colnames(maketube_results) <- c("pipeline", "FILTER", "TOTAL", "TOTAL2", "VARIANT_CLASS", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")
View(maketube_results)
shorter_alignement <- subset(maketube_results, VARIANT_CLASS == "snp")
ggplot(shorter_alignement, aes(pipeline, PRECISION, fill = pipeline)) +
  geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.01), linetype = "dotted") +
  geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
  geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
  geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0.90, 1, by = 0.05),
                     limits = c(0.90,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  scale_fill_manual(name = "source", values = my_colors) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
  ) + ylab("PRECISION")
