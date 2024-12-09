{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(flipr)

  all_strains_info <- read.table("all_strains_info.tsv", header = T, sep = "\t")
  
  source_order <- c('other_lineages', 'lineage4', 'snpmutator', 'maketube')
  all_strains_info$source <- factor(all_strains_info$source, levels = source_order, ordered = T)

  my_colors = c("dodgerblue", "chartreuse", "firebrick1", "chocolate1")
  names(my_colors) <- levels(all_strains_info$source)
  
  my_PCH = c(1, 2, 4, 3)
  names(my_PCH) <- levels(all_strains_info$source)

  my_filter_color = c("yellow", "yellow3", "orange2", "firebrick3")
  names(my_filter_color) <- levels(c("F0", "F1", "F2", "F3"))
} #start

{
  maketube_results <- read.table("maketube_metrics.tsv",
                                  sep = "\t", header = T
                                ) %>% na.omit()
  maketube_results <- maketube_results[as.numeric(maketube_results$RECALL) != 0 & as.numeric(maketube_results$PRECISION) != 0,]
  maketube_results <- maketube_results[maketube_results$filter != "nucmer" & maketube_results$filter != "PAF", ]
  maketube_results$filter_caller <- paste(maketube_results$filter, maketube_results$variant_caller, sep = "_")
  
  maketube_results$pipeline <- factor(maketube_results$pipeline, levels = c("snpmutator", "maketube"), ordered = T)
  
  filter_caller_order <- c("RAW_GATK", "RAW_freebayes", "MTBseq_GATK", "TBprofiler_mpileup", "genotube_freebayes")
  filter_caller_colors <- c("cadetblue1", "deepskyblue2", "chartreuse3", "tomato3", "firebrick2")
  names(filter_caller_colors) <- filter_caller_order
  
  maketube_results$filter_caller <- factor(maketube_results$filter_caller, levels = filter_caller_order, ordered = T)
  only_snp <- maketube_results[maketube_results$variant_class == "snp",]
  only_indel <- maketube_results[maketube_results$variant_class == "indel",]
} #import csv

snp_precision <- ggplot(only_snp) +
  geom_abline(slope = 0,
              intercept = c(0.5, 0.6, 0.7, 0.8, 0.9, 1),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(filter_caller, PRECISION, fill = filter_caller), show.legend = FALSE) +
  geom_point(aes(x = filter_caller, y = PRECISION, fill = filter_caller, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  scale_fill_manual(name = NULL, values = filter_caller_colors) +
  theme_classic() +
  ylab("PRECISION") +
  scale_y_continuous(limits = c(0.5, 1), n.breaks = 10, breaks = c(5:10)/10,
                     expand = expansion(mult = c(0, 0.3))) +
  scale_x_discrete(name="new axis name", labels=c("RAW_GATK" = "GATK \n(RAW)", "MTBseq_GATK" = "MTBseq", "TBprofiler_mpileup" = "TBprofiler", "RAW_freebayes" = "freebayes \n(RAW)", "genotube_freebayes" = "genotube")) +
  theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none"
  ) + ylab("PRECISION")

snp_recall <- ggplot(only_snp) +
  geom_abline(slope = 0, 
              intercept = c(0.875, 0.9, 0.925, 0.950, 0.975, 1),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(filter_caller, RECALL, fill = filter_caller), show.legend = FALSE) +
  geom_point(aes(x = filter_caller, y = RECALL, fill = filter_caller, alpha = 0.5),
             shape = 21,
             position = position_jitterdodge(jitter.width = 0.15)
  ) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  theme_classic() +
  ylab("RECALL") +
  scale_y_continuous(limits = c(0.875, 1), n.breaks = 6, breaks = seq(from = 0.875, to = 1, by = 0.025),
                     expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(name = "filter_caller", values = filter_caller_colors) + 
  scale_x_discrete(name="new axis name", labels=c("RAW_GATK" = "GATK \n(RAW)", "MTBseq_GATK" = "MTBseq", "TBprofiler_mpileup" = "TBprofiler", "RAW_freebayes" = "freebayes \n(RAW)", "genotube_freebayes" = "genotube")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0),
        legend.position = "none",
        strip.text = element_blank(),
        strip.background = element_blank()
        ) + ylab("RECALL")

svg("three_VC_SNP_precision_recall.svg")
grid.arrange(snp_precision, snp_recall, nrow = 2)
dev.off()


{
  a <- subset(only_snp, filter == "MTBseq" & pipeline == "maketube")$RECALL
  b <- subset(only_snp, filter == "TBprofiler" & pipeline == "maketube")$RECALL
  c <- subset(only_snp, filter == "genotube" & pipeline == "maketube")$RECALL
  
  two.wilcox.test <- wilcox.test(x = a, y = c, paired = T, alternative = "two.sided") ; two.wilcox.test; two.wilcox.test$statistic
  two.perm.test <- two_sample_test(x = a, y = c, B = 10000, stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue
} #rinse and repeat
