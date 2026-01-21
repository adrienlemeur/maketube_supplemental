{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggforce)
  library(concaveman)
  ###removed twice invoked library: ggforce
}

lineages_color <- c(L1 = "darkorchid1", L2 = "dodgerblue2", L3 = "blueviolet",L4 = "firebrick2",L5 =  "darkorange4" ,L6 = "darkgreen",L7 = "gold1", bovis = "black")
###Table does not have headers, chaged to header = F
dnadiff_data <- read.table("distance_to_reference.tsv", header = F, sep = "\t") 
colnames(dnadiff_data) <- c("reference", "reference_lineage", "reference_source", "sample", "sample_lineage", "sample_source", "ref_length", "ref_aligned", "sample_length", "sample_aligned")

# is the reference sequence from the same lineage as the sample considered ?
dnadiff_data$relationship <- dnadiff_data$sample_lineage == dnadiff_data$reference_lineage
dnadiff_data$relationship[dnadiff_data$relationship == TRUE] <- "intralineage"
dnadiff_data$relationship[dnadiff_data$relationship == FALSE] <- "interlineage"

###What is the use of lineage_source?
# is the sample considered natural, created by snpmutator or created by maketube
dnadiff_data$sample_source <- factor(dnadiff_data$sample_source, levels = c("natural", "snpmutator", "maketube"), ordered = T)
dnadiff_data$lineage_source <- paste(dnadiff_data$sample_lineage, dnadiff_data$sample_source, sep = ", ")

# The sequence we used for snpmutator contained NA, resulting in unaligned bits.

###I don't think you need this, it just complicate your script and does not change anything
# The % aligned on the reference was still 100% so we removed them for clarity.
sub1 <- subset(dnadiff_data, !(sample_source == "snpmutator" & sample_lineage == "L2"))

###The following left 174 rows
subset_dnadiff <- subset(sub1,
  (reference_lineage == sample_lineage & sample_source == "natural") | # only distance to H37Rv if artificial genome ### What do you mean?
  (sample_source != "natural" & reference_lineage == sample_lineage) | # only intralineage distance for natural genome ###What fo you mean
  (sample_lineage == "L5" & reference_lineage == "bovis") |
  (sample_lineage == "L1" & reference_lineage == "L4")
)


#visualisation
#reformating, not interesting
subset_dnadiff$sample_lineage <- as.factor(subset_dnadiff$sample_lineage) ; subset_dnadiff$sample_source <- as.factor(subset_dnadiff$sample_source)

subset_dnadiff$x <- as.numeric(1-subset_dnadiff$ref_aligned/subset_dnadiff$ref_length)
subset_dnadiff$y <- as.numeric(1-subset_dnadiff$sample_aligned/subset_dnadiff$sample_length)

###Panels are switched, here below is for panel B
##fig3_A <- 

fig3_A <- ggplot(data = subset_dnadiff) +
  geom_hline(yintercept = seq(0, 0.013, by = 0.001),
             linetype = "dotted", color = "grey", ) +
  geom_vline(
    xintercept = seq(0, 0.013, by = 0.001),
    linetype = "dotted", color = "grey", ) +
  geom_point(aes(x = as.numeric(1-ref_aligned/ref_length), y = as.numeric(1-sample_aligned/sample_length), fill = sample_lineage, pch = sample_source), size = 2) +
  geom_mark_hull(data = subset_dnadiff, aes(x = x,
                                            y = y,
                                            col = sample_source, label = sample_source),
                 label.hjust = 0, fill = "grey",
                 con.type = "elbow", concavity = 10, expand = 0.03, radius = 0.03,
                 label.fontface = "bold", label.fontsize = 13
  ) +
  scale_shape_manual(values = c(natural = 21, snpmutator = 24, maketube = 25), name = "Genome category :") +
  scale_fill_manual(values = lineages_color, name = "Strain lineage") +
  scale_color_manual(values = c("maketube" = "firebrick", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(pch = 21)), 
         pch = guide_legend(order = 1)) +
  ylab("Sample not covered by reference") +  ###Changed label
  xlab("Reference not covered by sample") +  ###Changed label
  scale_y_continuous(breaks = seq(0, 1, by = 0.002), labels = scales::percent) + scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.002))+
  theme_light() +
  theme(axis.title.y = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = 18, hjust = 0.5),
        axis.text.y = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 20), legend.title = element_text(size = 18, face = "bold"))
fig3_A

png("fig3_A_pairwise_nucleotide_distance.png", width = 800, height = 800)
fig3_A
dev.off()

tmp <- subset(subset_dnadiff, sample_source == "maketube")
summary( 1 - (tmp$ref_aligned / tmp$ref_length)    )*100
summary( 1 - (tmp$sample_aligned / tmp$sample_length)    )*100

summary((as.numeric(1-tmp$ref_aligned/tmp$ref_length) + as.numeric(1-tmp$sample_aligned/tmp$sample_length))/2)

###Here below is for panel A
png("fig3_B_pairwise_nucleotide_distribution.png")
ggplot(subset_dnadiff) +
  geom_vline(
              xintercept = seq(0, 0.013, by = 0.001),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(sample_source, ((1-ref_aligned/ref_length)+(1-sample_aligned/sample_length))/2, fill = sample_source)) +
  scale_fill_manual(values = c("maketube" = "firebrick", "snpmutator" = "orange", "natural" = "darkgreen"), guide = "none") +
  xlab("") + ylab("Average pairwise ndistance") + ###Changed label
  scale_y_continuous(breaks = seq(0, 1, by = 0.002), labels = scales::percent) +
  theme(axis.title.y = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 18), legend.title = element_text(size = 20, face = "bold"))
dev.off()
