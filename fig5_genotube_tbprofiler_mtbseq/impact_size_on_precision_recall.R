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
  results_by_region <- read.table("three_variant_caller_results.tsv", sep = "\t", header = F)
  colnames(results_by_region) <- c("pipeline", "filter", "variant_caller", "experience", "region", "region2", "variant_class", "TP", "FP", "FN", "RECALL", "PRECISION", "F1")
  
  total_genotube <- subset(results_by_region, pipeline == "maketube" & variant_caller == "genotube" & experience == 'dupli')

  reference_region_size <- 3842481

  PV <- (sum(total_genotube$TP) + sum(total_genotube$FN))/(reference_region_size*nrow(total_genotube))
  PIV <- 1 - PV

  PTP_genotube <- sum(total_genotube$TP) / (sum(total_genotube$TP) + sum(total_genotube$FN))
  PFN_genotube <- 1-PTP_genotube
  PFP_genotube <- sum(total_genotube$FP)/reference_region_size*nrow(total_genotube)

  theoretical_metrics_calculator <- function(size, replicate, PTP, PFP, PFN){
    TP <- rbinom(replicate, size, PTP)
    FP <- rbinom(replicate, size, PFP)
    FN <- rbinom(replicate, size, PFN)
    PRECISION <- TP / (TP + FP)
    RECALL <- TP / (TP + FN)
    output <- cbind(size, TP, FP, FN, RECALL, PRECISION)

    return(output)
  }
}

{
  interval_size = c(500, 1000, 2000, 5000, 10000, 20000, 50000)

  tmp <- lapply(interval_size, function(x) {
    theoretical_metrics_calculator(x, 1000, PTP_genotube*PV, PFP_genotube*PIV, PFN_genotube*PV)
  })

  genotube_theoretical_distribution <- na.omit(as.data.frame(do.call(rbind, tmp)))
  colnames(genotube_theoretical_distribution) <- c("region", "TP", "FP", "FN", "PRECISION", "RECALL")

  genotube_theoretical_distribution$TP <- as.numeric(genotube_theoretical_distribution$TP)
  genotube_theoretical_distribution$FP <- as.numeric(genotube_theoretical_distribution$FP)
  genotube_theoretical_distribution$FN <- as.numeric(genotube_theoretical_distribution$FN)
  genotube_theoretical_distribution$PRECISION <- as.numeric(genotube_theoretical_distribution$PRECISION)
  genotube_theoretical_distribution$RECALL <- as.numeric(genotube_theoretical_distribution$RECALL)
  genotube_theoretical_distribution$region <- as.factor(genotube_theoretical_distribution$region)
}

library(patchwork)

A <- ggplot(genotube_theoretical_distribution) + geom_boxplot(aes(region, as.numeric(PRECISION), group = region))
B <- ggplot(genotube_theoretical_distribution) + geom_boxplot(aes(region, as.numeric(RECALL), group = region))



