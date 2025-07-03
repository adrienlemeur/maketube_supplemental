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


# is the reference sequence from the same lineage as the sample considered ?
dnadiff_data$relationship <- dnadiff_data$sample_lineage == dnadiff_data$reference_lineage
dnadiff_data$relationship[dnadiff_data$relationship == TRUE] <- "intralineage"
dnadiff_data$relationship[dnadiff_data$relationship == FALSE] <- "interlineage"

# is the sample considered natural, created by snpmutator or created by maketube
dnadiff_data$sample_source <- factor(dnadiff_data$sample_source, levels = c("natural", "snpmutator", "maketube"), ordered = T)
dnadiff_data$lineage_source <- paste(dnadiff_data$sample_lineage, dnadiff_data$sample_source, sep = ", ")

# The sequence we used for snpmutator contained NA, resulting in unaligned bits.
# The % aligned on the reference was still 100% so we removed them for clarity.
sub1 <- subset(dnadiff_data, !(sample_source == "snpmutator" & sample_lineage == "L2"))

subset_dnadiff <- subset(sub1,
                          (reference == "H37Rv" & sample_source == "natural") | # only distance to H37Rv if artificial genome
#                          (reference_lineage == sample_lineage & sample_source == "natural") |
                          (sample_source != "natural" & reference_lineage == sample_lineage) # only intralineage distance for natural genome
                        )



#visualisation
fig3_A <- ggplot(data = subset_dnadiff) +
  geom_hline(yintercept = 0.001,
             linetype = "dotted", color = "red", ) +
  geom_hline(yintercept = seq(0, 0.013, by = 0.001),
             linetype = "dotted", color = "grey", ) +
  geom_vline(
    xintercept = seq(0, 0.013, by = 0.001),
    linetype = "dotted", color = "grey", ) +
  geom_point(aes(x = as.numeric(1-ref_aligned/ref_length), y = as.numeric(1-sample_aligned/sample_length), fill = sample_lineage, pch = sample_source), size = 2) +
  geom_mark_hull(aes(x = 1-ref_aligned/ref_length, y = 1-sample_aligned/sample_length, col = sample_source, label = sample_source),
                 label.hjust = 0, fill = "grey",
                 con.type = "elbow", concavity = 10, expand = 0.03, radius = 0.03,
                 label.fontface = "bold", label.fontsize = 13,
                 con.colour = "inherit", label.fill = NA
  ) +
  scale_shape_manual(values = c(natural = 21, snpmutator = 24, maketube = 25), name = "Genome category :") +
  scale_fill_manual(values = lineages_color, name = "Strain lineage") +
  scale_color_manual(values = c("maketube" = "firebrick", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(pch = 21)), 
         pch = guide_legend(order = 1)) +
  xlim(c(-0.001, 0.013)) + ylim(c(-0.001, 0.013)) +
  ylab("Distance to the reference") +
  xlab("Distance to the sample") +
  theme_light() +
  theme(axis.title.y = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
fig3_A

svg("fig3_A_pairwise_nucleotide_distance.svg", width = 10, height = 10)
fig3_A
dev.off()

svg("fig3_B_pairwise_nucleotide_distribution.svg")
ggplot(subset(subset_dnadiff, relationship == "intralineage")) +
  geom_vline(
              xintercept = seq(0, 0.013, by = 0.001),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(sample_source, ((1-ref_aligned/ref_length) + (1-sample_aligned/sample_length))/2, fill = sample_source)) +
  scale_fill_manual(values = c("maketube" = "firebrick", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  xlab("") + ylab("Average pairwise \ndistance to H37Rv") +
  theme(axis.title.y = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 14, hjust = 0.5),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = "bold"))
dev.off()
subset(dnadiff_data)
