# PCA and Volcano plot figures

```r
# Perform PCA and plot
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq data") +
  theme_minimal()

# Create a volcano plot for Nut10h vs Nut0h
volcanoData_Nut10h_vs_Nut0h <- data.frame(
  log2FoldChange = res_Nut10h_vs_Nut0h$log2FoldChange,
  -log10padj = -log10(res_Nut10h_vs_Nut0h$padj),
  gene = rownames(res_Nut10h_vs_Nut0h)
)

# Label the up- and downregulated genes
volcanoData_Nut10h_vs_Nut0h$regulation <- "Not significant"
volcanoData_Nut10h_vs_Nut0h$regulation[volcanoData_Nut10h_vs_Nut0h$log2FoldChange > 2 & volcanoData_Nut10h_vs_Nut0h$`-log10padj` > -log10(0.05)] <- "Upregulated"
volcanoData_Nut10h_vs_Nut0h$regulation[volcanoData_Nut10h_vs_Nut0h$log2FoldChange < -2 & volcanoData_Nut10h_vs_Nut0h$`-log10padj` > -log10(0.05)] <- "Downregulated"

ggplot(volcanoData_Nut10h_vs_Nut0h, aes(x = log2FoldChange, y = -log10padj, color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes (Nut10h vs Nut0h)") +
  theme_minimal()

# Create a volcano plot for Nut1h vs Nut0h
volcanoData_Nut1h_vs_Nut0h <- data.frame(
  log2FoldChange = res_Nut1h_vs_Nut0h$log2FoldChange,
  -log10padj = -log10(res_Nut1h_vs_Nut0h$padj),
  gene = rownames(res_Nut1h_vs_Nut0h)
)

# Label the up- and downregulated genes
volcanoData_Nut1h_vs_Nut0h$regulation <- "Not significant"
volcanoData_Nut1h_vs_Nut0h$regulation[volcanoData_Nut1h_vs_Nut0h$log2FoldChange > 2 & volcanoData_Nut1h_vs_Nut0h$`-log10padj` > -log10(0.05)] <- "Upregulated"
volcanoData_Nut1h_vs_Nut0h$regulation[volcanoData_Nut1h_vs_Nut0h$log2FoldChange < -2 & volcanoData_Nut1h_vs_Nut0h$`-log10padj` > -log10(0.05)] <- "Downregulated"

ggplot(volcanoData_Nut1h_vs_Nut0h, aes(x = log2FoldChange, y = -log10padj, color = regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-2, 2), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle("Volcano Plot of Differentially Expressed Genes (Nut1h vs Nut0h)") +
  theme_minimal()
```
