# Fig 4. Heatmap ----
pdf(file="Fig4.Heatmap.pdf",width = 12, height = 8)
## All compartments ----
Number_dynamic_compartments_Freq_ordered$order <- seq(1:nrow(Number_dynamic_compartments_Freq_ordered))

dynamic_compartments_categoryinfo <- Compartments[(! grepl("NULL", Compartments$Category_compartments)) & Compartments$state_compartments == "dynamic",]
dynamic_compartments_categoryinfo$order = Number_dynamic_compartments_Freq_ordered$order[match(dynamic_compartments_categoryinfo$Category_compartments,Number_dynamic_compartments_Freq_ordered$Category_compartments)]
dynamic_compartments_categoryinfo$Freq = Number_dynamic_compartments_Freq_ordered$Percentage_switch[match(dynamic_compartments_categoryinfo$Category_compartments,Number_dynamic_compartments_Freq_ordered$Category_compartments)]
dynamic_compartments_categoryinfo_ordered <- dynamic_compartments_categoryinfo %>% arrange(-Freq)

dynamic_compartments_categoryinfo_top1x100 <- dynamic_compartments_categoryinfo_ordered[dynamic_compartments_categoryinfo_ordered$Freq >=1,]
dynamic_compartments_categoryinfo_top1x100 <- dynamic_compartments_categoryinfo_top1x100[order(dynamic_compartments_categoryinfo_top1x100$order,decreasing=FALSE),]

mycols <- colorRamp2(breaks = c(-0.5,0, 0.5), colors = c("blue3", "snow1","red3"))

### Basic heatmap ----
Heatmap(as.matrix(dynamic_compartments_categoryinfo_top1x100[,c(4:9)]),column_title = "DON'T FORGET TO CHANGE",  cluster_rows = F, cluster_columns = F, show_row_names = F, col = mycols, clustering_distance_rows = "pearson")

### Clustered by kmeans ----
#Elbow Method for finding the optimal number of clusters
set.seed(103)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- dynamic_compartments_categoryinfo_top1x100[,c(4:9)]
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

Heatmap(as.matrix(data),column_title = "DON'T FORGET TO CHANGE",  cluster_rows = F, cluster_columns = T, show_row_names = F, col = mycols, clustering_distance_rows = "pearson",column_km = 5, row_km=5)
#heatmap.2(as.matrix(dynamic_compartments_categoryinfo_top1x100[,c(4,7:11)]), scale = "none", hclustfun = hclust, col = bluered(100),tracecol="grey50",density.info="none" )

### Splitted by chr ----
Heatmap(as.matrix(dynamic_compartments_categoryinfo_top1x100[,c(4:9)]), name = "DON'T FORGET TO CHANGE", cluster_row_slices = FALSE, show_row_names = F, cluster_columns = F,  row_split = factor(as.numeric(dynamic_compartments_categoryinfo_top1x100$chrom)),row_title = c("chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18","chr19", "chr20","chr21", "chr22","chrX"), row_title_rot = 0,row_title_gp = gpar(col = c("black"),fontsize=6),row_gap = unit(1, "mm"))


## Comparmtents with p53 binding ----
### Basic heatmap ----
Heatmap(as.matrix(dynamic_compartments_categoryinfo_top1x100[dynamic_compartments_categoryinfo_top1x100$p53 == "p53_binding",][,c(4:9)]),column_title = "DON'T FORGET TO CHANGE - P53", cluster_rows = F, cluster_columns = F, show_row_names = F, col = mycols, clustering_distance_rows = "pearson")

### Clustered by kmeans ----
#Elbow Method for finding the optimal number of clusters
set.seed(103)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
p53_data <- dynamic_compartments_categoryinfo_top1x100[dynamic_compartments_categoryinfo_top1x100$p53 == "p53_binding",][,c(4:9)]
wss <- sapply(1:k.max, 
              function(k){kmeans(p53_data, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

Heatmap(as.matrix(p53_data),column_title = "DON'T FORGET TO CHANGE - P53", cluster_rows = F, cluster_columns = T, show_row_names = F, col = mycols, clustering_distance_rows = "pearson",column_km = 3, row_km=4)

### Splitted by chr ----
Heatmap(as.matrix(dynamic_compartments_categoryinfo_top1x100[dynamic_compartments_categoryinfo_top1x100$p53 == "p53_binding",][,c(4:9)]), name = "DON'T FORGET TO CHANGE", cluster_row_slices = F, show_row_names = F, cluster_columns = F,  row_split = factor(as.numeric(dynamic_compartments_categoryinfo_top1x100[dynamic_compartments_categoryinfo_top1x100$p53 == "p53_binding",][,1])), row_title = c("chr1", "chr2", "chr3", "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18","chr19", "chr20","chr21", "chr22"), row_title_rot = 0,row_title_gp = gpar(col = c("black"),fontsize=6),row_gap = unit(1, "mm"))

dev.off()
#####
