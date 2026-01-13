{
  rm(list = ls())
  gc()
  graphics.off()

  library(ggplot2)
  library(gridExtra)
  library(scales)
  library(grDevices)
  library(flipr)
  library(patchwork)
  library(gghalves)
  library(ggpp)

  library(ggdist)
} #start

results_by_region <- read.table("results_by_regions_with_without_duplicated.every_at_10k.bed", sep = "\t", header = F)
colnames(results_by_region) <- c("pipeline", "strain", "experience", "region", "region2", "variant_class", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")

sum(subset(results_by_region, region2 == "IS_scar_regions" & experience == "without_dupli")$TP)

unique(results_by_region$region2)

{
  results_by_region$pipeline <- factor(results_by_region$pipeline, levels = c("genotube", "TBprofiler"), ordered = T)
  results_by_region$pipeline <- factor(results_by_region$pipeline, levels = c("genotube", "TBprofiler"), ordered = T)
  results_by_region$region2 <- factor(results_by_region$region2, levels = c("NEUTRAL", "DR_scar_regions", "insertion_flanking_regions", "IS_flanking_regions", "IS_scar_regions", "duplicated_region"), ordered = T)

  results_by_region <- subset(results_by_region, region2 != "TOTAL" & variant_class == "snp" & PRECISION > 0)
  #results_by_region <- results_by_region[results_by_region$variant_class == "snp" & results_by_region$region2 != "TOTAL", ]
  my_colors = c(genotube = "firebrick", TBprofiler = "chocolate1")
  custom_facet_labels <- c("NEUTRAL" = "Neutral", "TOTAL" = "TOTAL",
                           "IS_flanking_regions" = "IS flanking\n regions",
                           "insertion_flanking_regions"="Ancestral-like\nflanking regions",
                           "IS_scar_regions"="IS scar\n regions",
                           "DR_scar_regions"="Deletion scar\nregions",
                           "duplicated_region"="Duplicated\nregion")
  custom_facet_labeller <- as_labeller(custom_facet_labels)
  custom_facet_colors <- c("NEUTRAL" = "black",
                           "DR_scar_regions" = "cornflowerblue",
                           "insertion_flanking_regions" = "darkorchid4",
                           "IS_flanking_regions" = "darkolivegreen",
                           "IS_scar_regions" = "chartreuse4",
                           "duplicated_region" = "deeppink4")
  
  duplication_labeller <- as_labeller(c(with_dupli = "With duplication regions", without_dupli = "Without duplication region"))
} #structural variants table

PRECISION <- ggplot(results_by_region, aes(x = region2, y = as.numeric(PRECISION), fill = region2)) +
  geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5, 5.5), linetype = "dotted") +
  geom_vline(xintercept = 6.5, linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  theme_classic()  + ylab("PRECISION") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_fill_manual(name = element_blank(), values = custom_facet_colors) +
  scale_alpha(guide = 'none') +
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
  facet_wrap(~ experience, nrow = 1, scales = "free_x",labeller = duplication_labeller, strip.position = "top")

RECALL <- ggplot(results_by_region, aes(x = region2, y = as.numeric(RECALL), fill = region2)) +
  geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
  geom_vline(xintercept = c(1.5 , 2.5, 3.5, 4.5, 5.5), linetype = "dotted") +
  geom_vline(xintercept = 6.5, linetype = "solid") +
  geom_violin(weight = 1/2, position = position_dodge(), scale = "width") +
  geom_point(alpha = 0.5, shape = 21, position = position_jitterdodge(jitter.width = 0.5)) +
  theme_classic()  + ylab("RECALL") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                     expand = expansion(mult = c(0, 0.3))
  ) +
  scale_x_discrete(labels = custom_facet_labels) +
  coord_cartesian(ylim=c(0, 1)) +
  scale_fill_manual(name = element_blank(), values = custom_facet_colors) +
  scale_alpha(guide = 'none') +
  theme(text = element_text(NA),
        strip.text = element_blank(), #element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_blank(), legend.position = "bottom", strip.background = element_blank(), panel.spacing = unit(0, "lines"), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(colour = rep(custom_facet_colors, 2), size = 10, face = "bold", angle = 45, hjust = 1)
  ) +
  facet_wrap(~ experience, nrow = 1, scales = "free_x", labeller = duplication_labeller)

PRECISION / RECALL

svg("fig6_precision_and_recall_region_10k.svg", width = 10, height = 7)
PRECISION / RECALL
dev.off()

unique(results_by_region$region2)
mean(subset(results_by_region, experience == "without_dupli" & region2 == "insertion_flanking_regions")$FP)

precision_neutral_with_dupli <- subset(results_by_region, experience == "with_dupli" & region2 == "NEUTRAL")$PRECISION
recall_neutral_with_dupli <- subset(results_by_region, experience == "with_dupli" & region2 == "NEUTRAL")$RECALL

precision_neutral_without_dupli <- subset(results_by_region, experience == "without_dupli" & region2 == "NEUTRAL")$PRECISION
recall_neutral_without_dupli <- subset(results_by_region, experience == "without_dupli" & region2 == "NEUTRAL")$RECALL

my_region = 'duplicated_region'

test_precision_with_dupli <- subset(results_by_region, experience == "with_dupli" & region2 == my_region)$PRECISION
test_recall_with_dupli <- subset(results_by_region, experience == "with_dupli" & region2 == my_region)$RECALL

test_precision_without_dupli <- subset(results_by_region, experience == "without_dupli" & region2 == my_region)$PRECISION
test_recall_without_dupli <- subset(results_by_region, experience == "without_dupli" & region2 == my_region)$RECALL

#PRECISION WITH
two.wilcox.test <- wilcox.test(y = precision_neutral_with_dupli, x = test_precision_with_dupli, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(y = precision_neutral_with_dupli, x = test_precision_with_dupli, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue
?wilcox.test
#RECALL WITH
two.wilcox.test <- wilcox.test(y = recall_neutral_with_dupli, x = test_recall_with_dupli, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(y = recall_neutral_with_dupli, x = test_recall_with_dupli, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue

#PRECISION WITHOUT
two.wilcox.test <- wilcox.test(y = precision_neutral_without_dupli, x = test_precision_without_dupli, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(y = precision_neutral_without_dupli, x = test_precision_without_dupli, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue

#RECALL WITHOUT
two.wilcox.test <- wilcox.test(y = recall_neutral_without_dupli, x = test_recall_without_dupli, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(y = recall_neutral_without_dupli, x = test_recall_without_dupli, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue




a = as.numeric(subset(genotube_subset, region == my_region)$PRECISION)
b = as.numeric(subset(genotube_theo, region == my_region)$PRECISION)
which(c == 0)
View(subset(tbprofiler_subset, region == my_region & PRECISION == 0))
c[61]
c = as.numeric(subset(tbprofiler_subset, region == my_region)$PRECISION)
d = as.numeric(subset(tbprofiler_theo, region == my_region)$PRECISION)

two.wilcox.test <- wilcox.test(x = a, y = b, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(x = a, y = b, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue

two.wilcox.test <- wilcox.test(x = c, y = d, alternative = "less"); two.wilcox.test
two.perm.test <- two_sample_test(x = c, y = d, B = 10000, alternative = "left_tail", stats = list(stat_welch), type = "exact") ; two.perm.test$observed ; two.perm.test$pvalue
