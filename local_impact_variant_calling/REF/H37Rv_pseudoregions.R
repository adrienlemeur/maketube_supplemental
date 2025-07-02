{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(tidyr)
  library(dplyr)
}

all_genome_regions <- read.table("H37Rv_all_genomic_regions.pseudobed", header = T)
regions <- c("RLC", "PEPPE", "SUPPLY", "BLINDSPOTS", "REPETITIVE", "NEUTRAL") ; names(regions) <- c(1:5, "none")

replacing_regions <- function(x, regions) {
  region_name_vector <- sapply(strsplit(x, ","), function(z) regions[z])

  region_name_vector <-  region_name_vector[ ! region_name_vector %in% "SUPPLY" ]
  region_name_vector <- sort(unique(region_name_vector))
  paste0(region_name_vector, collapse = ";")
}

all_genome_regions$list <- sapply(all_genome_regions$list, function(x) replacing_regions(x, regions), USE.NAMES = F)
all_genome_regions$list[all_genome_regions$list == ""] <- "NEUTRAL"
independant_bed <- all_genome_regions[, c(1:3, 5)]

independant_unique_regions <- subset(independant_bed, list %in% regions)
independant_unique_regions$list <- paste0(independant_unique_regions$list, "_independant")
bed <- independant_bed %>%
  mutate(list = strsplit(as.character(list), ";")) %>% 
  unnest(list) %>% as.data.frame

write.table(independant_unique_regions, "H37Rv_independant_unique_region.bed", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(bed, "H37Rv_non_independant_region.bed", quote = F, sep = "\t", row.names = F, col.names = F)
independant_bed$list <- paste0(independant_bed$list, ";independant")
write.table(independant_bed, "H37Rv_independant_region.bed", quote = F, sep = "\t", row.names = F, col.names = F)
