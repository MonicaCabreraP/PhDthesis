Explanation script

# PCA Plot

- **rlog Transformation**: Regularized log transformation is performed on the count data to stabilize the variance across the mean, which is useful for PCA.

- **PCA Data Preparation**: The plotPCA function from DESeq2 is used to calculate the PCA, and the results are extracted for plotting.

- **Plotting**: ggplot2 is used to create the PCA plot, where the percentage variance explained by each principal component is displayed on the axes.

# Volcano Plot:

- **Data Preparation**: A data frame is created with log2 fold change and -log10 adjusted p-values.

- **Labeling**: Genes are labeled as "Upregulated", "Downregulated", or "Not significant" based on their log2 fold change and adjusted p-value.

- **Plotting**: ggplot2 is used to create the volcano plot, with different colors for upregulated and downregulated genes. Dashed lines indicate the thresholds for significance.
