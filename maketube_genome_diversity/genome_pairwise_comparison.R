{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggforce)
}

lineages_color <- c(L1 = "darkorchid1", L2 = "dodgerblue2", L3 = "blueviolet",L4 = "firebrick2",L5 =  "darkorange4" ,L6 = "darkgreen",L7 = "gold1", bovis = "black")

dnadiff_data <- read.table("distance_to_reference.tsv", header = T, sep = "\t")
colnames(dnadiff_data) <- c("reference", "reference_lineage", "reference_source", "sample", "sample_lineage", "sample_source", "ref_length", "ref_aligned", "sample_length", "sample_aligned")

dnadiff_data$relationship <- dnadiff_data$sample_lineage == dnadiff_data$reference_lineage
dnadiff_data$relationship[dnadiff_data$relationship == TRUE] <- "intralineage"
dnadiff_data$relationship[dnadiff_data$relationship == FALSE] <- "interlineage"

dnadiff_data$sample_source <- factor(dnadiff_data$sample_source, levels = c("natural", "snpmutator", "maketube"), ordered = T)
dnadiff_data$lineage_source <- paste(dnadiff_data$sample_lineage, dnadiff_data$sample_source, sep = ", ")

subset1 <- subset(dnadiff_data, !(sample_source == "snpmutator" & sample_lineage == "L2"))
subset_dnadiff <- subset(subset1,
                          (reference == "H37Rv" & sample_source == "natural") |
                          (reference_lineage == sample_lineage & sample_source == "natural") |
                          (sample_source != "natural" & reference_lineage == sample_lineage)
                        )

svg("fig3_A_pairwise_nucleotide_distance.svg", width = 10, height = 10)
ggplot(subset_dnadiff) +
  geom_hline(yintercept = seq(0, 0.013, by = 0.001),
              linetype = "dotted", color = "grey", ) +
  geom_vline(
              xintercept = seq(0, 0.013, by = 0.001),
              linetype = "dotted", color = "grey", ) +
  geom_point(aes(1-ref_aligned/ref_length, 1-sample_aligned/sample_length, fill = sample_lineage, pch = sample_source), size = 2) +
  geom_mark_hull(aes(x = 1-ref_aligned/ref_length, y = 1-sample_aligned/sample_length, col = sample_source, label = sample_source),
                 label.hjust = 0, fill = "grey",
                 con.type = "elbow", concavity = 10, expand = 0.03, radius = 0.03,
                 label.fontface = "bold", label.colour = "inherit", label.fontsize = 13,
                 con.colour = "inherit", label.fill = NA
  ) +
  scale_shape_manual(values = c(natural = 21, snpmutator = 24, maketube = 25), name = "Genome category :") +
  scale_fill_manual(values = lineages_color, name = "Strain lineage") +
  scale_color_manual(values = c("maketube" = "firebrick2", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  guides(fill = guide_legend(name = "test", override.aes = list(pch = 21)), 
         pch = guide_legend(name = "test2", order = 1)) +
  xlim(c(-0.001, 0.013)) + ylim(c(-0.001, 0.013)) +
  ylab("Distance to the reference") +
  xlab("Distance to the sample") +
  theme_light()
dev.off()

svg("fig3_B_pairwise_nucleotide_distribution.svg")
ggplot(subset(dnadiff_data, relationship == "intralineage")) +
  geom_vline(
              xintercept = seq(0, 0.013, by = 0.001),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(sample_source, ((1-ref_aligned/ref_length) + (1-sample_aligned/sample_length))/2, fill = sample_source)) +
  scale_fill_manual(values = c("maketube" = "firebrick2", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  xlab("") + ylab("Average pairwise \ndistance to H37Rv")
dev.off()

ggplot(test) +
  geom_hline(yintercept = seq(0, 0.013, by = 0.001),
             linetype = "dotted", color = "grey", ) +
  geom_vline(
    xintercept = seq(0, 0.013, by = 0.001),
    linetype = "dotted", color = "grey", ) +
  geom_point(aes(1-ref_aligned/ref_length, 1-sample_aligned/sample_length, fill = relationship, pch = sample_source), size = 2) +
  geom_mark_hull(aes(x = 1-ref_aligned/ref_length, y = 1-sample_aligned/sample_length, col = sample_source, label = sample_source),
                 label.hjust = 0, fill = "grey",
                 con.type = "elbow", concavity = 10, expand = 0.03, radius = 0.03,
                 label.fontface = "bold", label.colour = "inherit", label.fontsize = 13,
                 con.colour = "inherit", label.fill = NA
  ) +
  scale_shape_manual(values = c(natural = 21, snpmutator = 24, maketube = 25), name = "Genome category :") +
  guides(fill = guide_legend(name = "test", override.aes = list(pch = 21)), 
         pch = guide_legend(name = "test2", order = 1)) +
  xlim(c(-0.001, 0.013)) + ylim(c(-0.001, 0.013)) +
  ylab("Distance to the reference") +
  xlab("Distance to the sample") +
  theme_light()
