{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(ggforce)
}

dnadiff_data <- read.table("dnadiff_metrics.tsv", header = T, sep = "\t")
all_strains_info <- read.table("all_strains_info.tsv", header = T, sep = "\t")

source_order <- c('other_lineages', 'lineage4', 'maketube', 'snpmutator')
all_strains_info$source <- factor(all_strains_info$source, levels = source_order, ordered = T)

my_colors = c("dodgerblue", "green3", "firebrick1", "chocolate1")
names(my_colors) <- levels(all_strains_info$source)

my_lighter_colors <- sapply(my_colors, function(x) adjustcolor(x, alpha.f = 0.1))

all_info_merged <- merge(dnadiff_data, all_strains_info, by = 1)

all_info_merged$sample <- sapply(all_info_merged$sample, function(x) gsub(x, pattern = "H37Rv_mutated_", replacement = ""))
all_info_merged$sample <- sapply(all_info_merged$sample, function(x) gsub(x, pattern = "pop1_", replacement = ""))

my_PCH = c(1, 2, 4, 3)
names(my_PCH) <- levels(all_strains_info$source)

new_y_axis_breaks <- seq(from = 97, to = 100, by = 0.2)
new_x_axis_breaks <- seq(from = 97, to = 100, by = 0.2)


average_identity_mean <- ggplot(all_info_merged) +
  geom_abline(slope = 0,
              intercept = seq(98, 100, by = 0.125),
              linetype = "dotted", color = "grey") +
  geom_boxplot(aes(source, (sample_to_ref+ref_to_sample)/2, fill = source)) +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_y_continuous(
    breaks = seq(from = 98, to = 100, by = 0.125),
    limits = c(98, 100)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold", hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "none"
  ) + ylab("Average pairwise distance\nto H37Rv")

svg("average_pairwise_nucleotide.svg")
average_identity_mean
dev.off()



all_info_merged$light_fill_col <- sapply(all_info_merged$source, function(x) my_lighter_colors[x])
all_info_merged$sample[all_info_merged$source == "snpmutator" | all_info_merged$source == "maketube"] <- NA

svg("pairwise_nucleotidic_distance_to_reference.svg")
ggplot(all_info_merged) +
  geom_mark_hull(aes(x = sample_to_ref, y = ref_to_sample, fill = source, col = source, label = source),
                 label.hjust = 2,
                 con.type = "elbow", concavity = 10, expand = 0.03, radius = 0.03,
                 label.fontface = "bold", label.colour = "inherit", label.fontsize = 13,
                 con.colour = "inherit", label.fill = NA
                 ) +
  geom_point(aes(sample_to_ref, ref_to_sample, col = source, pch = source),
             size = 3, stroke = 1.5) +
  scale_fill_manual(name = "source", values = my_colors) +
  scale_color_manual(name = "source", values = my_colors) +
  scale_shape_manual(name = "source", values = my_PCH) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate(geom = "text", x = 97+0.1, y = 97, label ="x = y", cex = 5) +
  scale_x_continuous(breaks = new_x_axis_breaks, limits = c(97, 100)) +
  scale_y_continuous(breaks = new_y_axis_breaks, limits = c(97, 100)) +
  xlab("% of ref. aligned on sample") + ylab("% of sample aligned on ref.") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text=element_text(size=10),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  guides(pch = guide_legend(title = "Sample source :"),
         col = guide_legend(title = "Sample source :", fill = "Sample source :", label = "Sample source :")) +
  geom_label_repel(aes(sample_to_ref, ref_to_sample, label = sample, col = source),
                   nudge_y = -0.05,
                   size = 3, fill = NA, box.padding = 0.1,
                   point.padding = 0.5, segment.color = NA,
                   force_pull = 0.15, label.padding = 0,
                   label.size = 0, max.overlaps = 50
                  )
dev.off()

