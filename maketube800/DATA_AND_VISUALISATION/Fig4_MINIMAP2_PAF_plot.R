{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(scales)
  
  all_strains_info <- read.table("a_database_to_rule_them_all.tsv", header = F, sep = "\t")[, 1:2]
  colnames(all_strains_info) <- c("strain", "source")
  all_strains_info$source <- factor(all_strains_info$source, ordered = T)
  
  my_colors = c("firebrick", "chocolate1") ; names(my_colors) <- levels(c("snpmutator", "genotube"))
  
  my_PCH = c(1, 2, 4, 3)
  names(my_PCH) <- levels(all_strains_info$source)

  my_filter_color = c("yellow", "yellow3", "orange2", "firebrick3")
  names(my_filter_color) <- levels(c("F0", "F1", "F2", "F3"))
}

maketube_results <- na.omit(read.table("nucmer_paf.tsv", sep = "\t", header = T))
colnames(maketube_results) <- c("pipeline", "sample", "FILTER", "N", "TOTAL", "TOTAL2", "VARIANT_CLASS", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")

SNP <- subset(maketube_results,
                          (pipeline == "maketube" |
                          pipeline == "snpmutator") &
                          VARIANT_CLASS == "snp" & 
                          PRECISION != 0 & RECALL != 0)
PRECISION <- ggplot(SNP) +
  geom_abline(slope = 0, intercept = seq(0.75, 1, by = 0.05), linetype = "dotted") +
  geom_boxplot(aes(pipeline, PRECISION, fill = pipeline), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(x = pipeline, y = PRECISION, fill = pipeline, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ FILTER, nrow = 1, strip.position = "top", scales = "free_x") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0.75, 1, by = 0.05),
                     limits = c(0.74,1),
                     expand = expansion(mult = c(0, 0.3))
                     ) +
  scale_fill_manual(name = "source", values = my_colors) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
        ) + ylab("PRECISION")

RECALL <- ggplot(SNP) +
  geom_abline(slope = 0, intercept = seq(0.75, 1, by = 0.05), linetype = "dotted") +
  geom_boxplot(aes(pipeline, RECALL, fill = pipeline), show.legend = FALSE, outlier.shape = NA) + 
  geom_point(aes(x = pipeline, y = RECALL, fill = pipeline, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ FILTER, nrow = 1) +
  theme_classic() +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(breaks = seq(0.75, 1, by = 0.05),
                     limits = c(0.74,1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
        ) + ylab("RECALL")
grid.arrange(PRECISION, RECALL, nrow = 2, ncol = 1)

svg("nucmer_paf_precision_recall.svg")
grid.arrange(PRECISION, RECALL, nrow = 2, ncol = 1)
dev.off()

library(flipr)
#SNP$p
a = subset(SNP, pipeline == "snpmutator" & FILTER == "nucmer")$PRECISION
b = subset(SNP, pipeline == "maketube" & FILTER == "nucmer")$PRECISION

two.wilcox.test <- wilcox.test(x = a, y = b, alternative = "two.sided") ; two.wilcox.test; two.wilcox.test$statistic
two.perm.test <- two_sample_test(x = a, y = b, B = 100000, stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue

test







INDELS <- subset(maketube_results,
              (pipeline == "maketube" |
                 pipeline == "snpmutator") &
                VARIANT_CLASS == "indel" & 
                PRECISION != 0 & RECALL != 0)

INDEL_precision <- ggplot(INDELS) +
  geom_abline(slope = 0, intercept = 1, linetype = "dotted") +
  geom_boxplot(aes(pipeline, PRECISION, fill = pipeline), show.legend = FALSE, outlier.shape = NA) + 
  geom_point(aes(x = pipeline, y = PRECISION, fill = pipeline, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ FILTER, nrow = 1) +
  theme_classic() +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(breaks = breaks_pretty(n=10, min.n=0, max.n = 1)) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none"
  ) + ylab("PRECISION")

INDEL_recall <- ggplot(INDELS) +
  geom_abline(slope = 0, intercept = 1, linetype = "dotted") +
  geom_boxplot(aes(pipeline, RECALL, fill = pipeline), show.legend = FALSE, outlier.shape = NA) + 
  geom_point(aes(x = pipeline, y = RECALL, fill = pipeline, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ FILTER, nrow = 1) +
  theme_classic() +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(breaks = breaks_pretty(n=10, min.n=0, max.n = 1)) +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.position = "none"
  ) + ylab("RECALL")

grid.arrange(INDEL_precision, INDEL_recall, nrow = 2, ncol = 1)
