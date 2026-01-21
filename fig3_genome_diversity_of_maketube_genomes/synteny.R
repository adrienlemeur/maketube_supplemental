#!/usr/bin/env

{
  rm(list = ls())
  gc()
  graphics.off()
  library(ggplot2)
  library(gridExtra)
  library(pafr)
  library(cowplot)
  library(patchwork)
  library(ggtext)
}

{
  dnadiff_data <- read.table("distance_to_reference.tsv", header = T, sep = "\t")
  colnames(dnadiff_data) <- c("reference", "reference_lineage", "reference_source", "sample", "sample_lineage", "sample_source", "ref_length", "ref_aligned", "sample_length", "sample_aligned")

  subset_dnadiff <- subset(dnadiff_data,
                           (reference_lineage == sample_lineage & sample_source == "natural") | # only distance to H37Rv if artificial genome
                             (sample_source != "natural" & reference_lineage == sample_lineage) | # only intralineage distance for natural genome
                             (sample_lineage == "L5" & reference_lineage == "bovis") |
                             (sample_lineage == "L1" & reference_lineage == "L4")
  )
  subset_dnadiff$relationship <- subset_dnadiff$sample_lineage == subset_dnadiff$reference_lineage
  subset_dnadiff$relationship[dnadiff_data$relationship == TRUE] <- "intralineage"
  subset_dnadiff$relationship[dnadiff_data$relationship == FALSE] <- "interlineage"
} #input_data


radical <- function(string) {unlist(strsplit(as.character(string), split = "-|_"))[1]}

paf2synteny <- function(pattern_input, target_append = NULL, col1 = "violet", col2 = "purple", my_subtitle = NA, reference = NA){
  PAF <- read_paf(list.files("PAF", pattern = pattern_input, full.names = T))
  if(!is.null(target_append)){
    PAF$qname <- paste(unique(PAF$qname)[1], target_append, sep = "\n")
  }
  if(is.na(my_subtitle)){
    title_format <- element_blank()
  } else {
    title_format <- element_markdown(
      face = "bold",
      hjust = 0.5, vjust = 0,
      box.color = "black",
      linewidth = 0.75, linetype = 1,
      padding = margin(10, 30, 10, 30),
      margin = margin(0)
    )
  }

  if(!is.na(reference) & (radical(unique(PAF$tname)[1]) == radical(reference))){
    tmp <- col1
    col1 <- col2
    col2 <- tmp
  }
  if(PAF$tlen[1] >= PAF$qlen[1]){
    paf_swapped <- PAF
    paf_swapped$qname  <- PAF$tname
    paf_swapped$tname  <- PAF$qname
    PAF <- paf_swapped

    plot <- plot_synteny(PAF, q_chrom = unique(PAF$qname)[1], t_chrom = unique(PAF$tname)[1], rc = F, centre = T) +
      ggtitle(label = my_subtitle) +
      theme(plot.title = title_format,
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
      scale_y_discrete(labels = c(PAF$qname[1], PAF$tname[1]))
    plot_object <- ggplot_build(plot)
    plot_object$data[[2]][1, "fill"] <- col1
    plot_object$data[[2]][2, "fill"] <- col2
    plot_recolored <- wrap_ggplot_grob(ggplot_gtable(plot_object))

  } else {
    plot <- plot_synteny(PAF, q_chrom = unique(PAF$qname)[1], t_chrom = unique(PAF$tname)[1], rc = F, centre = T) +
      ggtitle(label = my_subtitle) +
      theme(plot.title = title_format,
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
    plot_object <- ggplot_build(plot)
    plot_object$data[[2]][2, "fill"] <- col1
    plot_object$data[[2]][1, "fill"] <- col2
    plot_recolored <- wrap_ggplot_grob(ggplot_gtable(plot_object))
  }

  return(plot_recolored)
}
#PAF.plot <- plot_synteny(PAF, q_chrom = unique(PAF$qname)[1], t_chrom = unique(PAF$tname)[1], rc = F, centre = T)
#PAF.build <- ggplot_build(PAF.plot)

PAF <- read_paf(list.files("PAF", pattern = "H37Rv_mutated_9_L4_H37Rv", full.names = T))
plot_synteny(PAF, q_chrom = unique(PAF$qname)[1], t_chrom = unique(PAF$tname)[1])

{
  natA <- paf2synteny("CDC1551_L4_H37Rv", col1 = "darkgreen", col2 = "darkgreen", my_subtitle = "Natural genome /<br> reference")
  natB <- paf2synteny("W-148_L2_18b.paf", col1 = "darkgreen", col2 = "darkgreen")
  natC <- paf2synteny("BCG_Pasteur_bovis_AF2122-97", col1 = "darkgreen", col2 = "darkgreen")
  
  snpmutatorA <- paf2synteny("H37Rv_mutated_9_L4_H37Rv", "(snpmutator)", col1 = "darkgreen", col2 = "orange", my_subtitle = "snpmutator genome /<br> reference")
  snpmutatorB <- paf2synteny("18b_mutated_9_L2_18b", "(snpmutator)", col1 = "darkgreen", col2 = "orange")
  snpmutatorC <- paf2synteny("AF2122-97_mutated_1_bovis_AF2122-97", "(snpmutator)", col1 = "orange", col2 = "darkgreen")

  maketubeA <- paf2synteny("SV2_pop1_H9_L4_H37Rv", "(maketube)", col1 = "darkgreen", col2 = "firebrick", my_subtitle = "maketube genome /<br> reference")
  maketubeB <- paf2synteny("SV1_pop1_H1_L2_18b", "(maketube)", col1 = "darkgreen", col2 = "firebrick")
  maketubeC <- paf2synteny("SV3_pop1_H10_bovis_AF2122", "(maketube)", col1 = "darkgreen", col2 = "firebrick")
}

big_plot <- (natA / natB / natC) | (snpmutatorA / snpmutatorB / snpmutatorC) | (maketubeA / maketubeB / maketubeC)
big_plot

ggsave(big_plot, device = "png", filename = "fig3_C_syntenyplot_exemples.png", width = 4000, height = 2000, units = "px")

subset_dnadiff <- subset_dnadiff[order(subset_dnadiff$sample_source, subset_dnadiff$reference, decreasing = F),]

gigantic_list_of_plots <- apply(subset_dnadiff, 1, function(x) {
        switch (x["sample_source"],
          "natural" = {target_color = "darkgreen"; suffix = " "},
          "snpmutator" = {target_color = "orange"; suffix = "(snpmutator)"},
          "maketube" = {target_color = "firebrick"; suffix = "(maketube)"},
        )

      my_file <- list.files("PAF", pattern = paste(x["sample"], x["reference"], sep ="_"), full.names = F)
      paf2synteny(my_file, target_append = suffix, col1 = target_color, col2 = "darkgreen", reference = x["reference"])
})

ggsave(wrap_plots(gigantic_list_of_plots[c(1:25)], ncol = 5, nrow = 5), device = "pdf", filename = "A", width = 4000, height = 2000, units = "px")

ggsave(wrap_plots(gigantic_list_of_plots[c(26:50)], ncol = 5, nrow = 5), device = "pdf", filename = "B", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(51:75)], ncol = 5, nrow = 5), device = "pdf", filename = "C", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(76:100)], ncol = 5, nrow = 5), device = "pdf", filename = "D", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(101:125)], ncol = 5, nrow = 5), device = "pdf", filename = "E", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(126:150)], ncol = 5, nrow = 5), device = "pdf", filename = "F", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(151:170)], ncol = 5, nrow = 5), device = "pdf", filename = "G", width = 4000, height = 2000, units = "px")
ggsave(wrap_plots(gigantic_list_of_plots[c(171:173)], ncol = 5, nrow = 5), device = "pdf", filename = "H", width = 4000, height = 2000, units = "px")



