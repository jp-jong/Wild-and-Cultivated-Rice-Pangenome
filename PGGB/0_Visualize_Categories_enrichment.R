library(ggplot2)
library(dplyr)

# ============================= #
# Create Visualization function #
# ============================= #
plot_GO_over <- function(go_mutated,title, top_n = 10, plot_name) {
  go_mutated <- read.csv(go_mutated, sep = "\t", header = TRUE)
  colnames(go_mutated) <- c("PANTHER GO-Slim Terms", "REFLIST", "Query","Expected",
                            "Over_Under", "Fold_enrichment", "Raw_P", "FDR", "Sample")
  go_mutated <- go_mutated[!(go_mutated[[1]] %in% 
                               c("Unclassified (UNCLASSIFIED)", 
                                 "biological_process (GO:0008150)",
                                 "protein class (PC00000)",
                                 "molecular_function (GO:0003674)"
                                 )), ] # to exclude GO terms
  # Ensure FDR is numeric and log-transform it
  go_mutated$FDR <- as.numeric(go_mutated$FDR)
  go_mutated$log_FDR <- -log10(go_mutated$FDR)
  
  first_column_name <- names(go_mutated)[1]
  Sample_order <- unique(go_mutated$Sample)
  
  # filter tops
  go_filtered <- go_mutated %>%
    group_by(Sample) %>%
    slice_min(order_by = FDR, n = top_n) %>%
    ungroup()
  
  go_filtered$Sample <- factor(go_filtered$Sample, levels = Sample_order)
  
  # plot
  
  plt <- ggplot(go_filtered, aes(x = Sample, y = .data[[first_column_name]])) + # edit Y 
    geom_point(aes(size = Query, color = log_FDR)) +  
    # scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    scale_size_continuous(range = c(1, 9)) +
    scale_color_gradient(low = "#92B9C6", high = "#B44C32", name = "-log10(FDR)") +
    theme_classic2() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
          ) +
    labs(x = "Group/Sample", y = "Gene Ontology Terms", 
         size = "Number of genes", 
         title = title) +
    guides(color = guide_colorbar(barwidth = 1, barheight = 3))
  
  print(plt)
  # Save as PDF
  ggsave(
    filename = paste0(plot_name, ".jpeg"),
    plot = plt,
    width = 11,
    height = 9,
    device = "jpeg",
    dpi = 300
  )
}

#### PAN CATEGORIES #####

setwd("/Users/jongpaduhilao/Desktop/LAB_Files/pggb/sample_output_all_O_mer_p60/fasta_annot/tsv/panther/results")
pan_bp <- "pan_bp.txt"
pan_mf <- "pan_mf.txt"
pan_pc <- "pan_pc.txt"
read.csv(pan_bp)

plot_GO_over(pan_bp, "Top 10 Biological Process Enrichment", 10, "top10_pan_BP")
plot_GO_over(pan_mf, "Top 10 Molecular Function Enrichment", 10, "top10_pan_MF")
plot_GO_over(pan_pc, "Top 10 Protein Class Enrichment", 10, "top10_pan_PC")

#### DISPENSABLE CATEGORIES #####
disp_mf <- "variable_mf.txt"
disp_pc <- "variable_pc.txt"

plot_GO_over(disp_mf, "Top 10 Molecular Function Enrichment (Variable nodes)", 10, "top10_dsp_mf")
plot_GO_over(disp_pc, "Top 10 Protein Class Enrichment (Variable nodes)", 10, "top10_dsp_pc")

