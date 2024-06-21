vsd <- vst(dds, blind = FALSE)

# Extract the normalized expression matrix
expr_matrix <- assay(vsd)

# Perform PCA
pca_res <- prcomp(t(expr_matrix), scale. = TRUE)

# Plot PCA
pca_data <- as.data.frame(pca_res$x)
pca_data$sample <- rownames(pca_data)

ggplot(pca_data, aes(x = PC1, y = PC2, label = sample)) +
  geom_point() +
  geom_text(vjust = -0.5) +
  labs(title = "PCA of RNA-seq Data", x = "Principal Component 1", y = "Principal Component 2")
