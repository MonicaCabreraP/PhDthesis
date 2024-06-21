# PCA and Volcano plot figures

```r
# Load required libraries
library(ggplot2)
library(pheatmap)

# Regularized log transformation for PCA
rld <- rlog(dds, blind = TRUE)

# Perform PCA and plot
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq data") +
  theme_minimal()

# Create a volcano plot
volcanoData <- data.frame(
  log2FoldChange = res$log2FoldChange,
  -log10padj = -log10(res$padj),
  gene = rownames(res)
)

# Label the up- and downregulated genes
volcanoData$regulation <- "Not significant"
volcanoData$regulation[volcanoData$log2FoldChange > 2 & volcanoData$`-log10padj` > -log10(0.05)] <- "Upregulated"
volcanoData$regulation[volcanoData$log2FoldChange < -2 & volcanoData$`-log10padj` > -log10(0.05)] <- "Downregulated"

ggplot(volcanoData, aes(x = log2FoldChange, y = -log10padj, color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes") +
  theme_minimal()
```
