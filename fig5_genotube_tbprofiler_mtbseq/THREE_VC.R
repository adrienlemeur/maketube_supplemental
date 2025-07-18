{
  rm(list = ls())
  gc()
  #graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(flipr)
  library(patchwork)
} #start

{
  maketube_results <- read.table("three_variant_caller_results.tsv",
                                  sep = "\t", header = F
                                ) %>% na.omit()
  colnames(maketube_results) <- c("pipeline", "filter", "variant_caller", "experience", "region", "region2", "variant_class", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")

  # reformat and conquer for visualisation
  maketube_results$filter_caller <- paste(maketube_results$variant_caller, maketube_results$filter, sep = "_")

  maketube_results$pipeline <- factor(maketube_results$pipeline, levels = c("snpmutator", "maketube"), ordered = T)

  filter_callers_labels = c("freebayes_raw_freebayes" = "freebayes \n(raw)", "MTBseq_samtools" = "MTBseq", "TBprofiler_freebayes" = "TBprofiler", "genotube_freebayes" = "genotube", "samtools_raw_samtools" = "samtools \n(raw)")

  filter_caller_order <- c("freebayes_raw_freebayes", "samtools_raw_samtools", "MTBseq_samtools", "TBprofiler_freebayes", "genotube_freebayes")
  filter_caller_colors <- c("cadetblue1", "khaki", "deepskyblue2", "chartreuse1", "coral")
  names(filter_caller_colors) <- filter_caller_order
  
  maketube_results$filter_caller <- factor(maketube_results$filter_caller, levels = filter_caller_order, ordered = T)
} #import csv

# result with duplication
maketube_results_dupli <- subset(maketube_results, variant_class == "snp" & experience == "dupli")
my_breaks <- c(5:20)/20
PRECISION_dupli <- ggplot(maketube_results_dupli, aes(filter_caller, PRECISION, fill = filter_caller)) +
  geom_abline(slope = 0, intercept = my_breaks, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5), linetype = "dotted") +
  geom_vline(xintercept = c(5.5), linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  scale_fill_manual(name = NULL, values = filter_caller_colors) +
  theme_classic() + ylab("PRECISION") +
  scale_y_continuous(limits = c(0.5, 1), n.breaks = 10, breaks = my_breaks,
                     expand = expansion(mult = c(0, 0.3))) +
  scale_x_discrete(name="new axis name", labels=c("raw" = "GATK \n(RAW)", "MTBseq_GATK" = "MTBseq", "TBprofiler_mpileup" = "TBprofiler", "RAW_freebayes" = "freebayes \n(RAW)", "genotube_freebayes" = "genotube")) +
  theme(text = element_text(NA),
        strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), #element_text(size = 10, face = "bold", angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.ticks.x = element_blank(), axis.line.x = element_blank()
  ) +
  ylab("PRECISION")

RECALL_dupli <- ggplot(maketube_results_dupli, aes(filter_caller, RECALL, fill = filter_caller)) +
  geom_abline(slope = 0, intercept = my_breaks, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5), linetype = "dotted") +
  geom_vline(xintercept = c(5.5), linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  theme_classic() + ylab("RECALL") +
  scale_y_continuous(limits = c(0.5, 1), n.breaks = 10, breaks = my_breaks,
                     expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(name = "Variant caller", values = filter_caller_colors, labels = filter_callers_labels) + 
  scale_x_discrete(name="new axis name", labels = filter_callers_labels) +
  theme(text = element_text(NA),
        strip.text = element_blank(), #element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(), legend.position = "bottom", strip.background = element_blank(), panel.spacing = unit(0, "lines"), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, hjust = 1)
  ) +
  ylab("RECALL") + scale_alpha(guide = 'none')
PRECISION_dupli / RECALL_dupli

svg("fig5_precision_recall_3VC_only_violin.svg", width = 10, height = 7)
PRECISION_dupli / RECALL_dupli
dev.off()

maketube_results_WO_DUPLI <- subset(maketube_results, experience == "wo_dupli")  

PRECISION_WO_DUPLI <- ggplot(maketube_results_WO_DUPLI, aes(filter_caller, PRECISION, fill = filter_caller)) +
  geom_abline(slope = 0, intercept = my_breaks, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5), linetype = "dotted") +
  geom_vline(xintercept = c(5.5), linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  scale_fill_manual(name = NULL, values = filter_caller_colors) +
  theme_classic() + ylab("PRECISION") +
  scale_y_continuous(limits = c(0.5, 1), n.breaks = 10, breaks = my_breaks,
                     expand = expansion(mult = c(0, 0.3))) +
  scale_x_discrete(name="new axis name", labels=c("RAW_GATK" = "GATK \n(RAW)", "MTBseq_GATK" = "MTBseq", "TBprofiler_mpileup" = "TBprofiler", "RAW_freebayes" = "freebayes \n(RAW)", "genotube_freebayes" = "genotube")) +
  theme(text = element_text(NA),
        strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), #element_text(size = 10, face = "bold", angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.ticks.x = element_blank(), axis.line.x = element_blank()
  ) +
  ylab("PRECISION")

RECALL_WO_DUPLI <- ggplot(maketube_results_WO_DUPLI, aes(filter_caller, RECALL, fill = filter_caller)) +
  geom_abline(slope = 0, intercept = my_breaks, linetype = "dotted", color = "grey") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5), linetype = "dotted") +
  geom_vline(xintercept = c(5.5), linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  facet_wrap(~ pipeline, nrow = 1, scales = "free_x") +
  theme_classic() + ylab("RECALL") +
  scale_y_continuous(limits = c(0.5, 1), n.breaks = 10, breaks = my_breaks,
                     expand = expansion(mult = c(0, 0.3))) +
  scale_fill_manual(name = "Variant caller", values = filter_caller_colors, labels = filter_callers_labels) + 
  scale_x_discrete(name="new axis name", labels = filter_callers_labels) +
  theme(text = element_text(NA),
        strip.text = element_blank(), #element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(), legend.position = "bottom", strip.background = element_blank(), panel.spacing = unit(0, "lines"), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, hjust = 1)
  ) +
  ylab("RECALL") + scale_alpha(guide = 'none')

svg("fig5_precision_recall_3VC_only_violin_wo_dupli.svg", width = 10, height = 7)
PRECISION_WO_DUPLI / RECALL_WO_DUPLI
dev.off()

#my_filter_caller <- as.character(unique(maketube_results_dupli$filter_caller))


my_pipeline <- "maketube"

mtbseq_precision <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "MTBseq_samtools")$PRECISION
tbprofiler_precision <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "TBprofiler_freebayes")$PRECISION
genotube_precision <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "genotube_freebayes")$PRECISION

mtbseq_recall <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "MTBseq_samtools")$RECALL
tbprofiler_recall <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "TBprofiler_freebayes")$RECALL
genotube_recall <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "genotube_freebayes")$RECALL

tmp <- subset(maketube_results_dupli, pipeline == my_pipeline & filter_caller == "TBprofiler_freebayes")

two.wilcox.test <- wilcox.test(x = a, y = b, alternative = "two.sided"); two.wilcox.test
two.perm.test <- two_sample_test(y = a, x = b, B = 10000, alternative = "two_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue
