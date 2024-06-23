# Volcano plots

# Load required libraries
library(ggplot2)
library(dplyr)

# Load the data
load_data <- function(file) {
  read.delim(file, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(log10pvalue = -log10(pvalue))
}

dea_data_1h <- load_data("DEA_Nutlin1h_vs_DMSO.txt")
dea_data_10h <- load_data("DEA_Nutlin10h_vs_DMSO.txt")

# Function to count upregulated and downregulated genes
count_genes <- function(data) {
  upregulated <- subset(data, padj < 0.05 & log2FoldChange >= 2)
  downregulated <- subset(data, padj < 0.05 & log2FoldChange <= -2)
  list(up = nrow(upregulated), down = nrow(downregulated))
}

counts_1h <- count_genes(dea_data_1h)
counts_10h <- count_genes(dea_data_10h)

print(paste("Number of up-regulated genes at 1h:", counts_1h$up))
print(paste("Number of down-regulated genes at 1h:", counts_1h$down))
print(paste("Number of up-regulated genes at 10h:", counts_10h$up))
print(paste("Number of down-regulated genes at 10h:", counts_10h$down))

# Define colors for up, down, and non-significant genes
colors <- c("upregulated" = "#37b578ff", "downregulated" = "#fde725ff", "not significant" = "grey")

# Function to create regulation column and volcano plot
create_volcano_plot <- function(data, title, up_color, down_color) {
  data <- data %>%
    mutate(regulation = ifelse(log2FoldChange >= 2 & padj < 0.05, "upregulated",
                               ifelse(log2FoldChange <= -2 & padj < 0.05, "downregulated", "not significant")))

  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("upregulated" = up_color, "downregulated" = down_color, "not significant" = "grey")) +
    geom_vline(xintercept = c(0, -2, 2), linetype = "dashed", color = c("black", "grey", "grey"), size = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted p-value", title = title) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))
}

# Create and display volcano plots
create_volcano_plot(dea_data_1h, "Volcano Plot of Differential Gene Expression at 1h", "#37b578ff", "#fde725ff")
create_volcano_plot(dea_data_10h, "Volcano Plot of Differential Gene Expression at 10h", "#43377fff", "#fde725ff")
