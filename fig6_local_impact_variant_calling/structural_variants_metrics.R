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
} #start

{
  results_by_region <- read.table("metrics_near_structural_variants.tsv", sep = "\t", header = F)
  colnames(results_by_region) <- c("pipeline", "strain", "experiment", "region", "region2", "variant_class", "TP", "FP", "FN", "RECALL", "PRECISION")

  neutral_genotube <- subset(results_by_region, pipeline == "genotube" & region == 'ANTIBIOTIC_GENE' & variant_class == "snp")

  #reference_region_size <- 3842481
  reference_region_size <- 50000

  PV <- (sum(as.numeric(neutral_genotube$TP)) + sum(as.numeric(neutral_genotube$FN)))/(reference_region_size*200)
  PIV <- 1 - PV

  PTP_genotube <- sum(neutral_genotube$TP) / (sum(neutral_genotube$TP) + sum(neutral_genotube$FN))
  PFN_genotube <- 1-PTP_genotube
  PFP_genotube <- sum(neutral_genotube$FP)/reference_region_size*200

  precision_calculator <- function(size, replicate, PTP, PFP, PFN, tag1, tag2){
    TP <- rbinom(replicate, size, PTP)
    FP <- rbinom(replicate, size, PFP)
    FN <- rbinom(replicate, size, PFN)
    PRECISION <- TP / (TP + FP)
    RECALL <- TP / (TP + FN)
    output <- cbind(tag1, tag2, size, TP, FP, FN, RECALL, PRECISION)

    return(output)
  }
}

{
  interval_size = c("DR_scar_regions" = 1800, "insertion_flanking_regions" = 4800, "IS_flanking_regions" = 9600, "IS_scar_regions" = 9600, "ANTIBIOTIC_GENE" = 50000, "duplicated_region" = 150000)

  tmp <- lapply(names(interval_size), function(region_name) {
    precision_calculator(interval_size[region_name], 200, PTP_genotube*PV, 0, PFN_genotube*PV, "genotube_theoretical", region_name)
  })
  genotube_theoretical_distribution <- as.data.frame(do.call(rbind, tmp))
  colnames(genotube_theoretical_distribution) <- c("pipeline", "region", "strain", "TP", "FP", "FN", "RECALL", "PRECISION")
  
  tmp <- lapply(names(interval_size), function(region_name) precision_calculator(interval_size[region_name], 2000, PTP_tbprofiler*PV, PFP_tbprofiler*PIV, PFN_tbprofiler*PV, "tbprofiler_theoretical", region_name))
  tb_profiler_theoretical_distribution <- as.data.frame(do.call(rbind, tmp))
  colnames(tb_profiler_theoretical_distribution) <- c("pipeline", "region", "strain", "TP", "FP", "FN", "RECALL", "PRECISION")

  both_distributions <- rbind(genotube_theoretical_distribution, tb_profiler_theoretical_distribution)[, -3]
  both_distributions$TP <- as.numeric(both_distributions$TP)
  both_distributions$FP <- as.numeric(both_distributions$FP)
  both_distributions$FN <- as.numeric(both_distributions$FN)
  both_distributions$PRECISION <- as.numeric(both_distributions$PRECISION)
  both_distributions$RECALL <- as.numeric(both_distributions$RECALL)
  both_distributions <- both_distributions[!(is.na(both_distributions$RECALL) | is.na(both_distributions$PRECISION)),]
  #write.table(both_distributions, file = "theoretical_distributions.tsv", sep = "\t", row.names = F)
}

{
  results_by_region <- read.table("metrics_near_structural_variants.tsv", sep = "\t", header = F)[,c(1,3,5:9)]
  results_by_region <- subset(results_by_region, V6 == "snp")[, -4]

  colnames(results_by_region) <- c("pipeline", "experiment", "region", "TP", "FP", "FN")
  results_by_region$PRECISION = results_by_region$TP / (results_by_region$FP + results_by_region$TP)
  results_by_region$RECALL = results_by_region$TP / (results_by_region$FN + results_by_region$TP)

  
  void_vector <- as.data.frame(t(sapply(unique(results_by_region$region), function(x) c("void", x, -1, -1, -1, -1, -1), USE.NAMES = F)))
  colnames(void_vector) <- c("pipeline", "region", "TP", "FP", "FN", "RECALL", "PRECISION")

  structural_variants <- rbind(results_by_region, void_vector)

  structural_variants$pipeline <- factor(structural_variants$pipeline, levels = c("genotube_theoretical", "genotube", "void", "tbprofiler_theoretical", "TBprofiler"), ordered = T)
  structural_variants <- na.omit(structural_variants[(structural_variants$pipeline == "genotube" | structural_variants$pipeline == "TBprofiler") & structural_variants$region != "TOTAL" & structural_variants$region != "ANTIBIOTIC_GENE", ])
  my_colors = c(genotube = "firebrick", tbprofiler_theoretical = "orange4", TBprofiler = "chocolate1", genotube_theoretical = "pink3")

} #structural variants table

{
  custom_facet_labels <- c("IS_flanking_regions" = "IS flanking\n regions",
                           "insertion_flanking_regions"="insertion flanking\n regions",
                           "IS_scar_regions"="IS scar\n regions",
                           "DR_scar_regions"="DR scar\n regions",
                           "duplicated_region"="duplicated region")
  custom_facet_labeller <- as_labeller(custom_facet_labels)
}

{
  PRECISION <- ggplot(structural_variants) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = as.numeric(PRECISION), fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = as.numeric(PRECISION), fill = pipeline, alpha = 0.5),
               shape = 21,
               position = position_jitterdodge(jitter.width = 0.15)
    ) +
    theme_classic()  + ylab("PRECISION") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_fill_manual(name = "source", values = my_colors, labels = c("genotube_theoretical" = "genotube\nunder H0", "genotube" = "genotube", "tbprofiler_theoretical" = "TBprofiler\nunder H0", "TBprofiler" = "TBprofiler")) +
    scale_alpha(guide = 'none') +
    theme(text = element_text(NA),
          strip.text = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none"
    ) +
    facet_wrap(~ pipeline, nrow = 1, scales = "free_x", labeller = custom_facet_labeller)

  RECALL <- ggplot(structural_variants) +
    geom_abline(slope = 0, intercept = seq(0, 1, by = 0.1), linetype = "dotted") +
    geom_boxplot(aes(x = region, y = as.numeric(RECALL), fill = pipeline), colour = "black", outliers = F) +
    geom_point(aes(x = region, y = as.numeric(RECALL), fill = pipeline, alpha = 0.5),
                shape = 21,
                position = position_jitterdodge(jitter.width = 0.15)
              ) +
    theme_classic() + ylab("RECALL") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1),
                       expand = expansion(mult = c(0, 0.3))
    ) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_fill_manual(name = "source", values = my_colors, labels = c("genotube_theoretical" = "genotube\nunder H0", "genotube" = "genotube", "tbprofiler_theoretical" = "TBprofiler\nunder H0", "TBprofiler" = "TBprofiler")) +
    scale_alpha(guide = 'none') +
    theme(text = element_text(NA),
          axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_blank()
    ) + 
    facet_wrap(~ pipeline, nrow = 1, scales = "free_x", labeller = custom_facet_labeller)
} #plot structural variants

PRECISION / RECALL

svg("fig6_precision_and_recall_region_with_H0_distribution_modified.svg", width = 10, height = 7)
PRECISION / RECALL
dev.off()

png("precision_and_recall_region_with_H0_distribution_modified.png", width = 3000, height = 2400, res = 400)
PRECISION / RECALL
dev.off()

View(results_by_region)


genotube_subset <- subset(structural_variants, pipeline == "genotube")
genotube_theo <- subset(structural_variants, pipeline == "genotube_theoretical")

tbprofiler_subset <- subset(structural_variants, pipeline == "TBprofiler")
tbprofiler_theo <- subset(structural_variants, pipeline == "tbprofiler_theoretical")
