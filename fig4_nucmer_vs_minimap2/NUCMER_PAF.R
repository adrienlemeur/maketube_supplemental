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
  library(ggsignif)

  my_colors = c(maketube = "firebrick", snpmutator = "chocolate1")
  
  custom_facet_labels <- c("nucmer" = "Nucmer", "minimap2" = "Minimap2")
  custom_facet_labeller <- as_labeller(custom_facet_labels)
}

maketube_results <- na.omit(read.table("nucmer_minimap2_results.tsv", sep = "\t", header = T))
colnames(maketube_results) <- c("sample", "pipeline", "structural_variant", "FILTER", "TOTAL", "TOTAL2", "VARIANT_CLASS", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")

SNP <- subset(maketube_results,
                   (pipeline == "maketube" |
                      pipeline == "snpmutator") &
                     VARIANT_CLASS == "snp" & structural_variant == "w_duplicated")

SNP_noSV <- subset(maketube_results,
              (pipeline == "maketube" |
                 pipeline == "snpmutator") &
                VARIANT_CLASS == "snp" & structural_variant == "wo_duplicated")

{
  PRECISION_noSV <- ggplot(SNP_noSV, aes(pipeline, PRECISION, fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.025), linetype = "dotted") +
    geom_half_boxplot(center = TRUE, errorbar.draw = FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
    geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
    geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
    facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.025),
                       limits = c(0.90,1),
                       expand = expansion(mult = c(0, 0.3))
                       ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(strip.text = element_text(face="bold", size=20),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none"
          ) + ylab("PRECISION")

  RECALL_noSV <- ggplot(SNP_noSV, aes(pipeline, RECALL, fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.025), linetype = "dotted") +
    geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
    geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
    geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
    facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
    theme_classic() +
    scale_fill_manual(name = "source", values = my_colors) +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.025),
                       limits = c(0.90,1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "none"
          ) + ylab("RECALL")
}

PRECISION_noSV / RECALL_noSV

{
  PRECISION <-
ggplot(SNP, aes(pipeline, PRECISION, fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.025), linetype = "dotted") +
    geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
    geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
    geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
    facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.025),
                       limits = c(0.90,1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    scale_fill_manual(name = "source", values = my_colors) +
    theme(strip.text = element_text(face="bold", size=20),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none"
    ) + ylab("PRECISION")
  
  RECALL <- ggplot(SNP, aes(pipeline, RECALL, fill = pipeline)) +
    geom_abline(slope = 0, intercept = seq(0.90, 1, by = 0.025), linetype = "dotted") +
    geom_half_boxplot(center=TRUE, errorbar.draw=FALSE, width = 0.5, nudge = 0.02, outlier.shape = NA) +
    geom_half_violin(side = "r", nudge = 0.02, weight = 1/2) +
    geom_point(alpha = 0.5, shape = 21, position = position_jitternudge(width = 0.1, nudge.from = "jittered", x = -0.15)) +
    facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x", labeller = custom_facet_labeller) +
    theme_classic() +
    scale_fill_manual(name = "source", values = my_colors) +
    scale_y_continuous(breaks = seq(0.90, 1, by = 0.025),
                       limits = c(0.90,1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "none"
    ) + ylab("RECALL")
}
PRECISION / RECALL

png("fig4_prec_recall_nucmer_paf.png", res = 300, width = 1200, height = 1200)
grid.arrange(PRECISION, RECALL, nrow = 2, ncol = 1)
dev.off()

png("S_fig4_prec_recall_nucmer_paf_no_duplication.png", res = 300, width = 1200, height = 1200)
grid.arrange(PRECISION_noSV, RECALL_noSV, nrow = 2, ncol = 1)
dev.off()

library(flipr)

source = "maketube"
#source = "snpmutator"
caller="minimap2"

a <- subset(SNP, FILTER == "minimap2" & pipeline == "snpmutator")$PRECISION
b <- subset(SNP, FILTER == "nucmer" & pipeline == "snpmutator")$PRECISION


wilcox.test(x = rep(1, 25), y = rep(1, 25), alternative = "two.sided")
two.wilcox.test <- wilcox.test(x = a, y = b, alternative = "two.sided"); two.wilcox.test; two.wilcox.test$statistic
two.perm.test <- two_sample_test(x = a, y = b, B = 100000, stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue

two.wilcox.test$method

SNP_noSV
source = "maketube"
#source = "snpmutator"

a <- subset(SNP, FILTER == "minimap2" & pipeline == "snpmutator")$PRECISION
b <- subset(SNP_noSV, FILTER == "minimap2" & pipeline == "maketube")$PRECISION

two.wilcox.test <- wilcox.test(x = a, y = b, alternative = "two.sided"); two.wilcox.test; two.wilcox.test$statistic
two.perm.test <- two_sample_test(x = a, y = b, B = 100000, stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue
