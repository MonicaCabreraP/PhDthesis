options(stringsAsFactors = F)
####
library(reshape2)
library(ggplot2)
library(corrplot)
#####

###############################################################
# PRE-CHICAGO: Creating a fake PeakMatrix with the downsampled chinputs 
## that will contain the number of reads x sample x interaction
###############################################################
# setwd("/home/monica/Desktop/projects/data/")
# chinputs <- "3.CHiCAGO_chinputs.txt"
# chinputs <- file(chinputs,open="r")
# chinputs <-readLines(chinputs)
# 
# df <- as.data.frame(matrix(ncol = 2))
# colnames(df) <- c("baitID", "otherEndID")
# df <- df[-1,]
# 
# for (i in 1:length(chinputs)) 
# {
#   a <- read.table(chinputs[i],header = T)
#   a <- a[,1:3]
#   colnames(a)[3] <- gsub(pattern = "_.*",replacement = "",basename(chinputs[i]))
#   df <- merge(df, a,by = c("baitID", "otherEndID"), all = T)
# }
# 
# df[is.na(df)] <- 0

##################################################################################################
# Correlation matrix showing correlation coefficients between biological replicates
## The default method used is pearson correlation coefficient - measures the linear dependence between two variables.  
## A value of 1 implies that a linear equation describes the relationship between X and Y perfectly,
## A negative value implies oposite proprocional.
## Near 0 implies low correlation between these two variables
### kendall and spearman correlation methods are non-parametric rank-based correlation test.##################################################################################################

# correlation <- cor(df[,-c(1:2)],method = "pearson")
# write.table(correlation, file ="/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/correlation_PCA_dendogram/prechicago_correlation.txt", row.names = T, col.names = T, quote = F, sep = "\t")
##################################################################################################
# Pre-CHiCAGO correlation matrix
##################################################################################################
pdf("/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/plots/Correlation_pre_chicago_downsampled.pdf", width = 15, height = 10)

prechicago <- read.table("/home/monica/Desktop/projects/p53/Step2_CHiCAGO/downsampled/NUTWT_recalibration/correlation_PCA_dendogram/prechicago_correlation.txt", header = T)
cormat <- as.matrix(prechicago)

# To identify the hidden pattern in the matrix it is usefull to reorder the matrix.

#To know the default orther: dimnames(...) )
cormat <- cormat[,c(8,7,6,5,4,3,2,1)] # en ese orden sale el mismo plot que utilizando la funcion reorder pero con la diagonal invertida 

# # # hclust for hierarchical clustering order is used in the example below.
#    reorder_cormat <- function(cormat){
# # #  dd <- as.dist((1-cormat)/2)
#     dd <- dist(cormat, method = "euclidean") # same as using above distance 
#     hc <- hclust(dd, method = "complete") #default cluster method: complete
# # #   # method: the agglomeration method to be used. 
# # #   ## This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#     cormat <- cormat[hc$order, hc$order]
#   }
# # # 
# # # # Reorder the correlation matrix
#  cormat <- reorder_cormat(cormat)

# As a correlation matrix has redundant information. We can use the functions below to set half matrix.
# Get lower triangle of the correlation matrix
# get_lower_triangle <- function(cormat){
#   cormat[upper.tri(cormat)] <- NA
#   return(cormat)
# }
# 
# Get upper triangle of the correlation matrix
# get_upper_triangle <- function(cormat){
#   cormat[lower.tri(cormat)] <- NA
#   return(cormat)
# }

# cormat <- get_upper_tri(cormat)
# cormat <- get_lower_tri(cormat)

# To change diagonal orientation 
# cormat <- cormat[,c(8:1)]

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)

# Create a ggheatmap
prechicago_correlation <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "black", mid = "grey", 
                       midpoint = median(melted_cormat$value), limit = c(min(melted_cormat$value), max(melted_cormat$value)), space = "Lab", 
                       name="Pearson Correlation") +
  theme_minimal()+ labs(title=paste("Pre-CHiCAGO correlation matrix")
                        #, subtitle = ""
                        )+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4, label=sprintf("%0.3f", round(melted_cormat$value, digits = 3))) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0, 0),
    legend.position = c(1, 0.9),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+
  coord_fixed()

# Print the heatmap
print(prechicago_correlation)
dev.off()

###############################################################################################
# POST-CHICAGO: Importing the PeakMatrix created by the chicagotool makePeakmatrix that contains 
## the chicago score (CS) per interaction and sample if in at least one sample the CS is >=5
# And calculating Pearson correlation matrix
###############################################################################################
#BR
pdf("/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/plots/Correlation_post_chicago_BR_downsampled.pdf", width = 15, height = 10)
peak_matrix_downsampled <- read.table("/home/monica/Desktop/projects/p53/Step2_CHiCAGO/downsampled/NUTWT_recalibration/correlation_PCA_dendogram/peakMatrix_downsampled.txt", header = T)
peak_matrix_downsampled <- peak_matrix_downsampled[12:19] #taking only the columns of the samples

cormat <- cor(as.matrix(peak_matrix_downsampled), method = c("pearson"))

#Merged
pdf("/home/monica/Desktop/projects/p53/Step3_analysis/NUTWT_recalibration/plots/Correlation_post_chicago_merged_downsampled.pdf", width = 15, height = 10)
peak_matrix_downsampled_merged <- read.table("/home/monica/Desktop/projects/p53/Step2_CHiCAGO/downsampled/NUTWT_recalibration/correlation_PCA_dendogram/peakMatrix_downsampled_merged.txt", header = T)
peak_matrix_downsampled_merged <- peak_matrix_downsampled_merged[12:15] #taking only the columns of the samples

cormat <- cor(as.matrix(peak_matrix_downsampled_merged), method = c("pearson"))

##################################################################################################
# Post-CHiCAGO correlation matrix
##################################################################################################
# To identify the hidden pattern in the matrix it is usefull to reorder the matrix.

#To know the default orther: dimnames(...) )
dimnames(cormat)
#BR
cormat <- cormat[,c(8,7,6,5,4,3,2,1)] # en ese orden sale el mismo plot que utilizando la funcion reorder pero con la diagonal invertida 
#Merged
cormat <- cormat[,c(4,3,2,1)] # en ese orden sale el mismo plot que utilizando la funcion reorder pero con la diagonal invertida 

## hclust for hierarchical clustering order is used in the example below.
reorder_cormat <- function(cormat){
  #dd <- as.dist((1-cormat)/2)
  dd <- dist(cormat, method = "euclidean") # same as using above distance 
  hc <- hclust(dd, method = "complete") #default cluster agglomaration method: complete
  # method: the agglomeration method to be used. 
  ## This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
  cormat <- cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)

# change orientation diagonal 
cormat <- cormat[,c(8:1)] # for all biological replicates
cormat <- cormat[,c(4:1)] #for merged

# As a correlation matrix has redundant information. We’ll use the functions below to set half matrix.
# Get lower triangle of the correlation matrix
# get_lower_triangle <- function(cormat){
#   cormat[upper.tri(cormat)] <- NA
#   return(cormat)
# }
# 
# Get upper triangle of the correlation matrix
# get_upper_triangle <- function(cormat){
#   cormat[lower.tri(cormat)] <- NA
#   return(cormat)
# }

# cormat <- get_upper_tri(cormat)
# cormat <- get_lower_tri(cormat)

#
# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)

# Create a ggheatmap
postchicago_correlation <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2( low = "white",
                        mid = "grey97",
                        high = "black",
                        #midpoint = 0.47, #for BR
                        midpoint = 0.6, #for merged
                        space = "Lab",
                        name="Pearson Correlation") +
  theme_minimal()+ labs(title=("Post-CHiCAGO correlation matrix"), sub= "caca"
                        #, subtitle = ""
  )+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4, label=sprintf("%0.3f", round(melted_cormat$value, digits = 3))) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0, 0.3),
    legend.position = c(1, 0.9),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+
  coord_fixed()

# Print the heatmap
print(postchicago_correlation)

dev.off()
##################################################################################################

