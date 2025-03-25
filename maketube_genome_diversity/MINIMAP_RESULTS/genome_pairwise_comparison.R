{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggforce)
}

dnadiff_data <- read.table("distance_to_reference.tsv", header = F, sep = "\t")
colnames(dnadiff_data) <- c("reference", "reference_lineage", "sample", "sample_lineage", "ref_length", "ref_aligned", "ref_ratio", "sample_length", "sample_aligned", "sample_ratio")
dnadiff_data$same_lineage <- dnadiff_data$reference_lineage == dnadiff_data$sample_lineage
dnadiff_data$same_lineage[dnadiff_data$same_lineage] <- "NO"
dnadiff_data$same_lineage[dnadiff_data$same_lineage == F] <- "YES"
dnadiff_data$same_lineage[dnadiff_data$sample_lineage == "snpmutator"] <- "snpmutator"
dnadiff_data$same_lineage[dnadiff_data$sample_lineage == "maketube"] <- "maketube"

ggplot(dnadiff_data) +
  geom_point(aes(ref_ratio, sample_ratio, col = same_lineage, pch = sample_lineage), size = 1) +
  scale_shape_manual(values=c(1:7)) + ylab("Aligned to reference") + xlab("Aligned to sample")
