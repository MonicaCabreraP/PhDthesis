options(stringsAsFactors = F)
####
library(cluster)
library(dendextend)
library(factoextra)
library(ggfortify)
library(dplyr)       # basic data manipulation and plotting
library(ggplot2)     # data visualization
################################################################################################
# Importing the PeakMatrix created by the chicagotool makePeakmatrix that contains 
## the chicago score (CS) per interaction and sample if in at least one sample the CS is >=5
#################################################################################################
peak_matrix_downsampled <- read.table("/home/monica/Desktop/projects/p53/Step2_CHiCAGO/NUTWT_recalibration/peakMatrix/peakMatrix_downsampled.txt", header = T)
peak_matrix_downsampled <- peak_matrix_downsampled[12:19] #taking only the columns of the samples

all <- na.omit(peak_matrix_downsampled) # Esto simplemente elimina cualquier fila que contenga NA
all <- scale(all)
all_t <- data.frame(t(all)) #To transform the table as required

##################################################################################################
# Principal Component Analysis (PCA)
##################################################################################################
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# https://bradleyboehmke.github.io/HOML/pca.html

pdf("/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/plots/PCA_downsampled.pdf", width = 15, height = 10)

pca <- prcomp(all_t) 
summary(all_t)

fviz_eig(pca) #Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.

PCA <- fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#39558CFF", "#FDE725FF"),
             repel = TRUE     # Avoid text overlapping
)
 
print(PCA)
dev.off()

##################################################################################################
# Hierarchical clustering
##################################################################################################

#Theory Hierarchical clustering
#################################
# It does not require to pre-specify the number of clusters to be generated and it uses pairwise distance matrix between observations as clustering criteria.
# Can be divided into two main types: 
#   - 1) **Agglomerative clustering**: 
#   It’s also known as AGNES (Agglomerative Nesting). It works in a bottom-up manner.That is, each object is initially considered as a single-element cluster (leaf). At each step of the algorithm, the two clusters that are the most similar are combined into a new bigger cluster (nodes). This procedure is iterated until all points are member of just one single big cluster (root) (see figure below). The result is a tree which can be plotted as a dendrogram.
#   *Agglomerative clustering is good at identifying small clusters.*
#   
#   The most common agglomeration methods are:
#   - **Maximum or complete linkage clustering** : It computes all pairwise dissimilarities
#     between the elements in cluster 1 and the elements in cluster 2, and considers 
#     the largest value (i.e., maximum value) of these dissimilarities as the distance 
#     between the two clusters. It tends to produce more compact clusters.
#   - **Minimum or single linkage clustering**: It computes all pairwise dissimilarities between the 
#     elements in cluster 1 and the elements in cluster 2, and considers the smallest of 
#     these dissimilarities as a linkage criterion. It tends to produce long, “loose” clusters. 
#   - **Mean or average linkage clustering**: It computes all pairwise dissimilarities between 
#     the elements in cluster 1 and the elements in cluster 2, and considers the average of 
#     these dissimilarities as the distance between the two clusters.
#   - **Ward’s minimum variance method**: It minimizes the total within-cluster variance. 
#     At each step the pair of clusters with minimum between-cluster distance are merged.

#   - 2) **Divisive hierarchical clustering**: 
#   It’s also known as DIANA (Divise Analysis) and it works in a top-down manner. The algorithm is an inverse order of AGNES. It begins with the root, in which all objects are included in a single cluster. At each step of iteration, the most heterogeneous cluster is divided into two. The process is iterated until all objects are in their own cluster (see figure below).
#   *Divisive hierarchical clustering is good at identifying large clusters.*
#   
#   ➙ **R FUNCTIONS**
#   The commonly used functions are:
#   - For agglomerative hierarchical clustering (HC): 
#   - hclust() [in stats package]
#   - agnes() [in cluster package]
#   - For divisive HC: 
#   - diana() [in cluster package]
#####
pdf("/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/plots/Hierarchical_cluster_downsampled.pdf", width = 15, height = 10)

# AGLOMERATIVE
################
# Euclidean distance
all_dist_Euclidean <- dist(all_t, method = "euclidean") #Here we use the Euclidean distance

## hclust()
all_hclust_complete_E <- hclust(all_dist_Euclidean, method = "complete")
#plot(all_hclust_complete_E, main = "Complete linkage") 

#all_hclust_single_E <- hclust(all_dist_Euclidean, method = "single")
#plot(all_hclust_single_E, main = "Simple linkage")

#all_hclust_average_E <- hclust(all_dist_Euclidean, method = "average") 
#plot(all_hclust_average_E, main = "Average linkage")

#all_hclust_ward_E <- hclust(all_dist_Euclidean, method = "ward.D") 
#plot(all_hclust_ward_E, main = "Ward linkage")

## agnes()
# all_agnes_complete <- agnes(all_dist_Euclidean, method = "complete")

# Manhattan distance
# all_dist_Manhattan <- dist(all_t, method = "manhattan")#This metric is less affected by outliers than the Euclidean, is more robust
# all_hclust_complete_M <- hclust(all_dist_Manhattan, method = "complete")
# plot(all_hclust_complete_M, main = "Complete linkage")
####################
# Divisive Hierarchical Clustering
####################
# all_divisive_E <- diana(all_dist_Euclidean) # To compute divisive hierarchical clustering
# all_divisive_M <- diana(all_dist_Manhattan) # To compute divisive hierarchical clustering
# all_divisive_E %>% pltree(main="Euclidean")
####################
# For viridis colors: 
## to know the rgb: library(scales)
                  # show_col(viridis_pal()(20))
mycolors <- c("#FDE725FF", "#32648EFF")
col_dend <-  all_hclust_complete_E %>%
  as.dendrogram() %>%
  set("branches_k_color", mycolors, k = 2)

plot(col_dend, main = "Hierarchical clustering \n hclust (complete)")

dev.off()
