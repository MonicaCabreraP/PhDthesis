## Script name: Load_HiC_data.R               
## Purpose: Load all data for HiC analysis (Compartments, TADs and epigenomic data to overlap)      
## Last modification: 06 Sep 2022  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com

# Notes ### ###
## When merging compartments with TADs some compartments might contain up to 3 Insulation scores, therefore, I repeat the compartment lines adding the different IS values in the macrotable (from X compartments to X compartments)
### ### ### ###

samples=c("WT.DMSO","WT.Nutlin.1h", "WT.Nutlin.4h", "WT.Nutlin.7h", "WT.Nutlin.10h", "WT.Nutlin.24h")

# Packages
packages <- c("reshape", "GenomicRanges", "dplyr","ggplot2", "reshape", "hrbrthemes", "stringr", "ComplexHeatmap", "gplots", "tidyverse", "circlize","RColorBrewer")
invisible(lapply(packages, library, character.only = TRUE))

# set the working directory where the data is
#wd <- "/home/mcabrera/Desktop/MN/projects/DATA/"
wd <- "/home/monica/Desktop/MN/"
setwd(paste0(wd,"projects/DATA/HCT116/HiC/"))

### ### ### ### ### 
# Compartments ----

Compartments <- read.table("p53_compartments_v37.tsv", sep="\t", header=T)
Compartments <- Compartments[,-c(1,6,7)] # Removing id column and KOs data (DMSO and Nutlin1h)
colnames(Compartments) <- c("chrom","start", "end","WT.DMSO","WT.Nutlin.1h", "WT.Nutlin.4h", "WT.Nutlin.7h", "WT.Nutlin.10h", "WT.Nutlin.24h")
Compartments[is.na(Compartments)] <- 0 #To not delete all rows with one NA at some time point, I convert the NA to 0 
Compartments <- Compartments[rowSums(Compartments[,c(4:9)])!=0,] #removing all rows that have 0 (NA) at all time points  

## Categorizing compartments between A and B with the compartment values
for (i in samples) { Compartments[[paste0("Category_",i)]] <- as.character(ifelse(Compartments[i] < 0, 'B', ifelse(Compartments[i] > 0, 'A', 'NULL'))) }
Compartments$Category_compartments <- do.call(paste, c(Compartments[,c(10:15)], sep="-"))
Compartments$state_compartments <- as.character(ifelse(Compartments$Category_compartments == "A-A-A-A-A-A" , 'static in A', ifelse(Compartments$Category_compartments == "B-B-B-B-B-B", 'static in B', 'dynamic')))

# 1) How many compartments do we have for the analysis? ----
print(paste0("There are ", nrow(Compartments), " compartments with information at some time-point")) 
# "There are 28.084 compartments with value at some points of the time-course" ----
print(paste0("There are ", nrow(Compartments[- grep("NULL", Compartments$Category_compartments),]), " compartments with information at all time-points")) 
# "There are 26.976 compartments with value at all points of the time-course" ----


# 2) Is the A/B distribution equal in all samples? ----
table_AB <- data.frame(apply(Compartments[- grep("NULL", Compartments$Category_compartments),][,c(4:9)], 2, function(x) { table(sign(x)) }))
rownames(table_AB) <- c("B", "A")
# WT.DMSO WT.Nutlin.1h WT.Nutlin.4h WT.Nutlin.7h WT.Nutlin.10h WT.Nutlin.24h compartment
# B   11182        10315        10595        10809         12320         11802           B
# A   15794        16661        16381        16167         14656         15174           A
table_AB$compartment <- rownames(table_AB)
table_AB_melt <- melt(table_AB, id="compartment")

# 3) How many compartments are static and how many compartments change residency? ----
static_dynamic_compartments <- data.frame(sort(table(Compartments[- grep("NULL", Compartments$Category_compartments),]$state_compartments), decreasing = T))
### static in A: 13.838 ----
### static in B: 9.574  ----
### dynamic: 3.564  ----
colnames(static_dynamic_compartments) <- c("Category", "Total_number")

# 4) How is the compartment residency happening? ----
dynamic_compartments_Freq <- data.frame(sort(table(Compartments$Category_compartments),decreasing = T))
# A-A-A-A-A-A                B-B-B-B-B-B             A-A-A-A-A-NULL                A-A-A-A-B-A 
# 13838                       9574                        799                        575 
# A-A-A-A-A-B                A-A-A-A-B-B                B-A-A-A-B-B                B-A-B-B-B-B 
# 417                        389                        317                        313 
# B-A-A-B-B-B                B-A-A-A-B-A                B-B-B-B-B-A             B-B-B-B-B-NULL 
# 155                        153                        153                        104 
# B-B-A-B-B-B       NULL-A-A-A-NULL-NULL                B-A-B-A-B-B                A-B-B-B-B-B 
# 86                         69                         65                         62 
# B-A-A-A-A-A                A-A-A-B-B-B                A-B-A-A-B-B                B-A-A-B-B-A 
# 62                         49                         49                         46 
# B-A-B-B-B-A                A-B-A-B-B-B                B-B-B-B-A-B                A-A-B-A-B-B 
# 43                         39                         35                         30 
# A-B-A-A-A-A                A-B-B-A-B-B                B-B-B-B-A-A                A-A-B-A-A-B 
# 30                         27                         24                         23 
# A-A-B-B-B-B                B-B-B-A-B-B                B-B-A-A-B-B                A-A-B-A-A-A 
# 23                         23                         21                         20 
# A-B-A-A-A-B                A-B-B-B-A-B          A-A-A-A-NULL-NULL                A-B-B-B-A-A 
# 20                         20                         19                         19 
# B-B-A-B-B-A                A-A-A-B-B-A                A-A-A-B-A-A                A-B-B-B-B-A 
# 19                         18                         17                         17 
# A-NULL-A-A-A-NULL                A-B-B-A-A-A                B-A-B-A-B-A                A-A-B-B-A-A 
# 15                         14                         14                         12 
# A-A-B-B-A-B                B-B-B-A-B-A                A-A-B-A-B-A                A-B-A-B-B-A 
# 12                         12                         11                         11 
# A-A-A-B-A-B                B-A-A-A-A-B                A-B-A-B-A-A                A-B-B-A-A-B 
# 10                         10                          9                          9 
# NULL-NULL-A-A-NULL-NULL                B-A-B-B-A-A                B-A-B-B-A-B             A-A-A-A-B-NULL 
# 9                          8                          8                          7 
# A-B-A-A-B-A                B-A-A-B-A-A                B-B-A-A-B-A                A-A-B-B-B-A 
# 7                          7                          7                          6 
# A-B-B-A-B-A    A-NULL-NULL-NULL-A-NULL             B-B-B-B-A-NULL          NULL-A-A-A-A-NULL 
# 6                          6                          6                          6 
# A-B-A-B-A-B                B-A-A-B-A-B          B-NULL-B-B-B-NULL NULL-NULL-NULL-A-NULL-NULL 
# 5                          5                          5                          5 
# B-A-B-A-A-A                B-B-A-A-A-A                B-B-A-B-A-A                B-B-B-A-A-A 
# 4                          4                          4                          4 
# NULL-NULL-A-NULL-NULL-NULL             A-B-A-A-A-NULL             B-A-A-A-B-NULL             B-A-B-B-B-NULL 
# 4                          3                          3                          3 
# B-B-B-A-A-B             B-B-B-A-B-NULL       NULL-NULL-B-B-B-NULL NULL-NULL-NULL-B-NULL-NULL 
# 3                          3                          3                          3 
# NULL-NULL-NULL-NULL-A-NULL             A-B-A-B-B-NULL       A-NULL-NULL-A-A-NULL             B-A-A-A-A-NULL 
# 3                          2                          2                          2 
# B-B-A-A-A-B    NULL-NULL-NULL-B-B-NULL NULL-NULL-NULL-NULL-B-NULL          A-B-A-NULL-A-NULL 
# 2                          2                          2                          1 
# A-B-B-B-A-NULL             A-B-B-B-B-NULL       A-NULL-A-A-NULL-NULL       A-NULL-NULL-A-B-NULL 
# 1                          1                          1                          1 
# B-A-B-A-A-B             B-B-A-A-B-NULL             B-B-A-B-B-NULL             B-B-B-A-A-NULL 
# 1                          1                          1                          1 
# B-B-B-B-NULL-NULL          B-NULL-A-B-B-NULL             B-NULL-B-B-B-B       B-NULL-NULL-B-B-NULL 
# 1                          1                          1                          1 
# B-NULL-NULL-NULL-A-NULL    B-NULL-NULL-NULL-B-NULL    NULL-A-A-NULL-NULL-NULL             NULL-B-B-B-B-B 
# 1                          1                          1                          1 
# NULL-B-B-B-B-NULL       NULL-B-B-B-NULL-NULL       NULL-NULL-A-A-A-NULL       NULL-NULL-B-A-B-NULL 
# 1                          1                          1                          1 
# NULL-NULL-B-B-NULL-NULL    NULL-NULL-NULL-A-A-NULL       NULL-NULL-NULL-B-B-B 
# 1                          1                          1 

## Adding the category type id into the compartment matrix:
dynamic_compartments_Freq$id <- seq(1:nrow(dynamic_compartments_Freq))
colnames(dynamic_compartments_Freq) <- c("Category_compartments", "Freq", "id")
Compartments$compartment_category_id = dynamic_compartments_Freq$id[match(Compartments$Category_compartments,dynamic_compartments_Freq$Category_compartments)]

Compartments$ID <- seq(1:nrow(Compartments)) #this will be used when adding TAD value information as rows will be duplicated or triplicated according to the number of TADs found in each compartment row 

## Creating GRanges
Compartments.gr <- makeGRangesFromDataFrame(Compartments,seqnames.field = "chrom", start.field = "start",end.field = "end", keep.extra.columns = T)

###############

### ### ### ### ###
# TADs ----
TADs <- read.table("p53_tads_v37.tsv", sep="\t", header=T)
TADs <- TADs[,-c(1)] # Removing rowID and KO's information
TADs <- TADs[,c(1:4,8,12,16,20,24)] ## Selecting Insulation Score
TADs[is.na(TADs)] <- 0

TADs <- TADs[rowSums(TADs[,c(4:9)])>0,] #removing all rows with 0 

TADs$individual_TADstate <- data.frame(apply(TADs[,c(4:9)], 2, function(x) { ifelse(x == 0, 'NoTAD', 'TAD') }))
TADs$Category_tads <- do.call(paste, c(TADs[,c(10)], sep="-"))

TADs$state_tads <- as.character(ifelse(TADs$Category_tads == "TAD-TAD-TAD-TAD-TAD-TAD" , 'static', 'stage specific'))

# 5) How many regions with TAD insulation score do we have for the analysis? ----
print(paste0("There are ", nrow(TADs), " TAD boundaries")) 
# "There are 5.832 TAD boundaries with information at least in one sample" ----

# 6) How many regions are static and stage specific? ----
static_stage_specific <- data.frame(sort(table(TADs$state_tads), decreasing = T))
# "There are 4.397 TAD boundaries stage specific " ----
# "There are 1.435 TAD boundaries static" ----

# 7) Which is the average IS of each sample? ----
IS_average_value <- data.frame(apply(TADs[,c(4:9)], 2, function(x) { sum(x)}))
# score_WT.DMSO_IS                                     2230.841
# score_WT.Nutlin.1h_IS                                2009.344
# score_WT.Nutlin.4h_IS                                2030.871
# score_WT.Nutlin.7h_IS                                2114.337
# score_WT.Nutlin.10h_IS                               2334.787
# score_WT.Nutlin.24h_IS                               2130.720

# 8) How many regions with IS values are in each sample? (X out of 5.832 meaning the rest are NA values) ----
data.frame(apply(TADs[,c(4:9)], 2, function(x) { nrow(TADs[x>0,])}))
# score_WT.DMSO_IS                                         3776
# score_WT.Nutlin.1h_IS                                    3542
# score_WT.Nutlin.4h_IS                                    3580
# score_WT.Nutlin.7h_IS                                    3655
# score_WT.Nutlin.10h_IS                                   3785
# score_WT.Nutlin.24h_IS                                   3450

# 9) How is the TAD distribution? Are TADs appearing, disappearing or static through the time-course? ----
sort(table(TADs$Category_tads), decreasing = T)
# TAD-TAD-TAD-TAD-TAD-TAD         TAD-TAD-TAD-TAD-TAD-NoTAD NoTAD-NoTAD-NoTAD-NoTAD-NoTAD-TAD 
# 1435                               280                               240 
# NoTAD-TAD-NoTAD-NoTAD-NoTAD-NoTAD         TAD-NoTAD-TAD-TAD-TAD-TAD   NoTAD-NoTAD-NoTAD-NoTAD-TAD-TAD 
# 165                               157                               145 
# NoTAD-NoTAD-TAD-NoTAD-NoTAD-NoTAD         TAD-TAD-NoTAD-TAD-TAD-TAD     TAD-NoTAD-NoTAD-NoTAD-TAD-TAD 
# 140                               137                               121 
# TAD-TAD-TAD-TAD-NoTAD-NoTAD NoTAD-NoTAD-NoTAD-TAD-NoTAD-NoTAD       TAD-NoTAD-NoTAD-TAD-TAD-TAD 
# 112                               110                                99 
# NoTAD-TAD-TAD-TAD-TAD-TAD   TAD-NoTAD-TAD-NoTAD-NoTAD-NoTAD         TAD-TAD-TAD-NoTAD-TAD-TAD 
# 93                                93                                91 
# TAD-TAD-NoTAD-NoTAD-NoTAD-NoTAD       TAD-NoTAD-TAD-TAD-TAD-NoTAD     NoTAD-TAD-TAD-TAD-NoTAD-NoTAD 
# 88                                85                                83 
# TAD-TAD-NoTAD-TAD-TAD-NoTAD TAD-NoTAD-NoTAD-NoTAD-NoTAD-NoTAD   TAD-NoTAD-NoTAD-NoTAD-NoTAD-TAD 
# 82                                78                                78 
# TAD-NoTAD-NoTAD-TAD-TAD-NoTAD   NoTAD-TAD-TAD-NoTAD-NoTAD-NoTAD       TAD-TAD-TAD-NoTAD-TAD-NoTAD 
# 74                                68                                68 
# TAD-NoTAD-NoTAD-TAD-NoTAD-NoTAD     TAD-NoTAD-TAD-NoTAD-TAD-NoTAD         TAD-TAD-TAD-TAD-NoTAD-TAD 
# 66                                64                                64 
# NoTAD-TAD-TAD-TAD-TAD-NoTAD   NoTAD-NoTAD-NoTAD-TAD-TAD-NoTAD   NoTAD-TAD-NoTAD-TAD-NoTAD-NoTAD 
# 63                                61                                58 
# TAD-NoTAD-TAD-NoTAD-TAD-TAD   NoTAD-NoTAD-TAD-TAD-NoTAD-NoTAD       NoTAD-NoTAD-TAD-TAD-TAD-TAD 
# 58                                56                                56 
# NoTAD-TAD-NoTAD-TAD-TAD-TAD       TAD-TAD-NoTAD-NoTAD-TAD-TAD   NoTAD-TAD-NoTAD-NoTAD-TAD-NoTAD 
# 56                                54                                53 
# TAD-TAD-NoTAD-NoTAD-TAD-NoTAD     TAD-TAD-NoTAD-TAD-NoTAD-NoTAD   NoTAD-NoTAD-TAD-NoTAD-TAD-NoTAD 
# 53                                53                                51 
# NoTAD-NoTAD-NoTAD-NoTAD-TAD-NoTAD   NoTAD-NoTAD-NoTAD-TAD-NoTAD-TAD     NoTAD-NoTAD-TAD-NoTAD-TAD-TAD 
# 50                                45                                45 
# NoTAD-NoTAD-TAD-NoTAD-NoTAD-TAD     TAD-NoTAD-TAD-TAD-NoTAD-NoTAD     TAD-TAD-TAD-NoTAD-NoTAD-NoTAD 
# 44                                44                                44 
# NoTAD-NoTAD-NoTAD-TAD-TAD-TAD       NoTAD-TAD-TAD-TAD-NoTAD-TAD   NoTAD-TAD-NoTAD-NoTAD-NoTAD-TAD 
# 43                                43                                40 
# NoTAD-NoTAD-TAD-TAD-TAD-NoTAD   TAD-NoTAD-NoTAD-NoTAD-TAD-NoTAD     NoTAD-TAD-NoTAD-NoTAD-TAD-TAD 
# 39                                38                                36 
# NoTAD-TAD-TAD-NoTAD-TAD-NoTAD       NoTAD-TAD-TAD-NoTAD-TAD-TAD     TAD-NoTAD-TAD-NoTAD-NoTAD-TAD 
# 35                                35                                30 
# NoTAD-NoTAD-TAD-TAD-NoTAD-TAD       TAD-NoTAD-TAD-TAD-NoTAD-TAD     NoTAD-TAD-NoTAD-TAD-TAD-NoTAD 
# 29                                29                                28 
# TAD-TAD-NoTAD-TAD-NoTAD-TAD     NoTAD-TAD-NoTAD-TAD-NoTAD-TAD     TAD-TAD-NoTAD-NoTAD-NoTAD-TAD 
# 27                                26                                26 
# TAD-TAD-TAD-NoTAD-NoTAD-TAD     TAD-NoTAD-NoTAD-TAD-NoTAD-TAD     NoTAD-TAD-TAD-NoTAD-NoTAD-TAD 
# 26                                22                                20 

## Creating GRanges
TADs.gr <- makeGRangesFromDataFrame(TADs,seqnames.field = "chrom", start.field = "start",end.field = "end", keep.extra.columns = T)
###################

### ### ### ### ###
# p53 binding ----
p53 <- read.table(paste0(wd,"/HCT116/ChIPseq_TF/p53/Nutlin3a/WT/HCT116_ChIP_TF_Nutlin3a_WT.peaks.txt"))
p53 <- p53[,c("V2","V3","V4")]
colnames(p53) <- c("chr", "start", "end")
p53 <- data.frame(lapply(p53, gsub, pattern='chr', replacement=''))

# 10) How many p53 bindings do we have? ----
print(paste0("There are ", nrow(p53), " p53 bindings")) 
# "There are 5.667 p53 bindings" ----

## Creating GRanges
p53.gr <- makeGRangesFromDataFrame(p53,seqnames.field = "chr", start.field = "start",end.field = "end", keep.extra.columns = T)
###################

# Creating macrotables ----

## Adding to the compartments table the p53 binding information
overlap_Compartments_p53 <- countOverlaps(Compartments.gr,p53.gr)
Compartments$p53 <- overlap_Compartments_p53
Compartments$p53 <- as.character(ifelse(Compartments$p53 == "0" , 'no_p53_binding', ifelse(Compartments$p53 >= "1", 'p53_binding', 'nothing')))

for (i in samples) { Compartments[[paste0("p53_Category_",i)]] <- as.character(ifelse(Compartments[paste0("Category_",i)] == "A" & Compartments$p53 == "p53_binding", 'A_p53', 
                                                                                                ifelse(Compartments[paste0("Category_",i)] == "A" & Compartments$p53 == "no_p53_binding", 'A',
                                                                                                       ifelse(Compartments[paste0("Category_",i)] == "B" & Compartments$p53 == "p53_binding", 'B_p53', 
                                                                                                              ifelse(Compartments[paste0("Category_",i)] == "B" & Compartments$p53 == "no_p53_binding", 'B', 'NULL'))))) }

# 11) How many compartments have p53 binding? ----
table_compartments_p53 <- sort(table(as.data.frame(Compartments[- grep("NULL", Compartments$Category_compartments),]$p53)), decreasing = T)
### no_p53_binding: 22.325 ----
### p53_binding : 4.651 ----

# 12) To which compartment is p53 binding? ----
table_AB_p53_info <- data.frame(apply(Compartments[- grep("NULL", Compartments$Category_compartments),][,c(21:26)], 2, function(x) { table(x)}))
colnames(table_AB_p53_info) <-  sub("p53_Category_", "", colnames(table_AB_p53_info))
table_AB_p53_info$compartment <- rownames(table_AB_p53_info)

# WT.DMSO WT.Nutlin.1h WT.Nutlin.4h WT.Nutlin.7h WT.Nutlin.10h WT.Nutlin.24h compartment
# A       12286        13006        12776        12573         11267         11712           A
# A_p53    3508         3655         3605         3594          3389          3462       A_p53
# B       10039         9319         9549         9752         11058         10613           B
# B_p53    1143          996         1046         1057          1262          1189       B_p53

table_AB_p53_info_melt <- melt(table_AB_p53_info, id="compartment")

# 13) How are the compartments with p53 binding (static or dynamic)? ----
compartments_p53_binding <- Compartments[-grep("NULL", Compartments$Category_compartments),]
compartments_p53_binding <- compartments_p53_binding[-grep("no_p53_binding", compartments_p53_binding$p53),]
table(compartments_p53_binding$state_compartments)
# dynamic: 543 ----
# static in A: 3.240 ----
# static in B: 868 ----

# 14) How is the compartment residency happening in regions with p53 binded? ----
dynamic_compartments_p53_Freq <- data.frame(sort(table(compartments_p53_binding$Category_compartments),decreasing = T))
# A-A-A-A-A-A A-A-A-A-A-B A-A-A-A-B-A A-A-A-A-B-B A-A-A-B-A-A A-A-A-B-A-B A-A-A-B-B-A A-A-A-B-B-B A-A-B-A-A-A A-A-B-A-A-B 
# 3240          64          82          49           4           1           3           1           4           7 
# A-A-B-A-B-A A-A-B-A-B-B A-A-B-B-A-A A-A-B-B-A-B A-A-B-B-B-A A-A-B-B-B-B A-B-A-A-A-A A-B-A-A-A-B A-B-A-A-B-A A-B-A-A-B-B 
# 2           4           4           1           1           1           6           3           2           2 
# A-B-A-B-A-A A-B-A-B-B-B A-B-B-A-A-B A-B-B-A-B-A A-B-B-A-B-B A-B-B-B-A-A A-B-B-B-A-B A-B-B-B-B-A A-B-B-B-B-B B-A-A-A-A-A 
# 2           4           1           1           2           2           2           1          12          15 
# B-A-A-A-A-B B-A-A-A-B-A B-A-A-A-B-B B-A-A-B-A-A B-A-A-B-A-B B-A-A-B-B-A B-A-A-B-B-B B-A-B-A-A-A B-A-B-A-B-B B-A-B-B-A-A 
# 2          22          51           2           1           7          13           2          16           1 
# B-A-B-B-A-B B-A-B-B-B-A B-A-B-B-B-B B-B-A-A-A-A B-B-A-A-A-B B-B-A-A-B-A B-B-A-A-B-B B-B-A-B-A-A B-B-A-B-B-A B-B-A-B-B-B 
# 1           7          47           2           1           1           3           1           2          19 
# B-B-B-A-A-A B-B-B-A-A-B B-B-B-A-B-A B-B-B-A-B-B B-B-B-B-A-A B-B-B-B-A-B B-B-B-B-B-A B-B-B-B-B-B 
# 1           1           4           4          10           8          31         868 

colnames(dynamic_compartments_p53_Freq) <- c("Category_compartments", "Freq")

# 15) How many compartments have insulation score values? ----
## Adding to the compartments table the TAD information
overlap_Compartment_TADs <- countOverlaps(query = Compartments.gr , subject = TADs.gr)
Compartments$tads <- overlap_Compartment_TADs

## Creating GRanges with TADs information
Compartments.gr <- makeGRangesFromDataFrame(Compartments,seqnames.field = "chrom", start.field = "start",end.field = "end", keep.extra.columns = T)

## Adding the Insulation score of the TADs in the compartments with TAD borders information
overlap_Compartment_TADs <- findOverlaps(query = Compartments.gr , subject = TADs.gr)

# Features from gr1 with overlaps in gr2
# Note: The same feature from gr1 can overlap with multiple features from gr2
Compartments_TADs.gr <- Compartments.gr[queryHits(overlap_Compartment_TADs)];

# Add the metadata from gr2
mcols(Compartments_TADs.gr) <- cbind.data.frame(
  mcols(Compartments_TADs.gr),
  mcols(TADs.gr[subjectHits(overlap_Compartment_TADs)]));

Compartments_TADs <- data.frame(Compartments_TADs.gr)

macrotable_Compartments <- dplyr::full_join(Compartments, Compartments_TADs, by="ID")

table(Compartments$tads)
## 20.768 compartments with no IS value ----
## 6.030 with 1 IS value ----
## 1.261 with 2 IS values ---- 
## 25 with 3 IS values. 
#Therefore, there are 1.311 (25*2 + 1261) compartments lines repeated

# 16) How are the compartments with IS values? ----
sort(table(Compartments_TADs[-grep("NULL", Compartments_TADs$Category_compartments),]$state_compartments))
# dynamic: 1.122 ---- 
# static in B: 2.404 ----
# static in A: 4.787 ----

table_categories_compartments_with_TADs <- sort(table(Compartments_TADs$Category_compartments), decreasing = T)
# A-A-A-A-A-A                B-B-B-B-B-B             A-A-A-A-A-NULL                A-A-A-A-B-A                A-A-A-A-A-B 
# 4787                       2404                        240                        184                        144 
# A-A-A-A-B-B                B-A-A-A-B-B                B-A-B-B-B-B                B-A-A-A-B-A                B-A-A-B-B-B 
# 111                        110                        101                         49                         47 
# B-B-B-B-B-A             B-B-B-B-B-NULL                B-A-A-A-A-A                A-B-B-B-B-B                B-A-B-A-B-B 
# 42                         34                         27                         26                         23 
# B-B-A-B-B-B                A-A-A-B-B-B                B-A-A-B-B-A                B-A-B-B-B-A                A-A-B-B-B-B 
# 23                         19                         15                         14                         11 
# B-B-B-B-A-B                A-A-A-B-B-A                A-B-A-A-A-A                A-B-A-A-B-B                A-B-B-B-A-B 
# 9                          8                          8                          8                          8 
# A-NULL-A-A-A-NULL       NULL-A-A-A-NULL-NULL                A-A-A-B-A-A                A-A-B-A-B-B                A-B-B-A-B-B 
# 8                          8                          7                          7                          7 
# B-B-A-B-B-A                B-B-B-B-A-A                A-A-A-B-A-B                A-A-B-A-A-A                A-B-B-A-B-A 
# 7                          7                          6                          6                          6 
# A-B-B-B-A-A                A-B-B-B-B-A                B-B-A-A-B-B                A-A-B-B-A-A                A-B-A-A-B-A 
# 6                          6                          6                          5                          5 
# A-B-A-B-B-B                A-B-B-A-A-A    NULL-NULL-A-A-NULL-NULL             A-A-A-A-B-NULL                A-B-A-A-A-B 
# 5                          5                          5                          4                          4 
# A-B-A-B-B-A                B-A-A-A-A-B                B-B-A-B-A-A                B-B-B-A-B-B                B-A-A-B-A-A 
# 4                          4                          4                          4                          3 
# B-A-B-B-A-B                A-A-B-A-A-B                A-A-B-A-B-A                B-A-A-B-A-B                B-A-B-A-A-A 
# 3                          2                          2                          2                          2 
# B-A-B-B-A-A          B-NULL-A-B-B-NULL NULL-NULL-A-NULL-NULL-NULL          A-A-A-A-NULL-NULL                A-A-B-B-A-B 
# 2                          2                          2                          1                          1 
# A-A-B-B-B-A             A-B-A-A-A-NULL                A-B-A-B-A-A                A-B-B-A-A-B    A-NULL-NULL-NULL-A-NULL 
# 1                          1                          1                          1                          1 
# B-A-A-A-A-NULL                B-A-B-A-B-A                B-B-B-A-A-A                B-B-B-A-A-B                B-B-B-A-B-A 
# 1                          1                          1                          1                          1 
# B-B-B-B-A-NULL          B-B-B-B-NULL-NULL          B-NULL-B-B-B-NULL          NULL-A-A-A-A-NULL          NULL-B-B-B-B-NULL 
# 1                          1                          1                          1                          1 
# NULL-NULL-A-A-A-NULL NULL-NULL-NULL-NULL-A-NULL 
# 1                          1 

# 17) How are the compartments with p53 binding and TAD IS value? ----
table_compartments_categories_with_TADs_p53 <- sort(table(Compartments_TADs[-grep("no_p53_binding", Compartments_TADs$p53),]$Category_compartments), decreasing = T)
# A-A-A-A-A-A          B-B-B-B-B-B       A-A-A-A-A-NULL          A-A-A-A-B-A          A-A-A-A-A-B          A-A-A-A-B-B          B-A-A-A-B-B 
# 1125                  246                   67                   28                   27                   17                   14 
# B-A-B-B-B-B          B-A-A-A-B-A          B-A-B-A-B-B          B-B-A-B-B-B NULL-A-A-A-NULL-NULL          B-B-B-B-B-A       B-B-B-B-B-NULL 
# 11                    9                    8                    6                    6                    5                    5 
# A-B-B-B-B-B          A-B-A-A-B-A          B-A-A-A-A-A          A-A-B-A-A-A          A-A-B-A-B-B          A-A-B-B-A-A          A-B-B-A-B-A 
# 4                    3                    3                    2                    2                    2                    2 
# A-B-B-B-A-B    A-NULL-A-A-A-NULL          B-A-A-B-A-A          B-A-A-B-B-B          B-A-B-A-A-A          B-A-B-B-B-A          B-B-A-B-A-A 
# 2                    2                    2                    2                    2                    2                    2 
# B-B-B-A-B-B          B-B-B-B-A-B       A-A-A-A-B-NULL          A-A-A-B-A-B          A-A-A-B-B-B          A-A-B-A-A-B          A-A-B-B-B-B 
# 2                    2                    1                    1                    1                    1                    1 
# A-B-A-A-A-A          A-B-A-A-A-B          A-B-A-B-B-B          A-B-B-A-A-B          A-B-B-A-B-B          B-A-A-A-A-B          B-B-A-A-B-B 
# 1                    1                    1                    1                    1                    1                    1 
# B-B-A-B-B-A          B-B-B-A-B-A          B-B-B-B-A-A    NULL-B-B-B-B-NULL 
# 1                    1                    1                    1 

table_static_dynamic_compartments_with_TADs_p53 <- sort(table(Compartments_TADs[-grep("no_p53_binding", Compartments_TADs$p53),]$state), decreasing = T)
#dynamic: 173
# static in A: 1125
# static in B: 246


# 
# 

# ### 2) Pairwise comparison against DMSO ----
# Nutlin_timecourse=c("WT.Nutlin.1h", "WT.Nutlin.4h", "WT.Nutlin.7h", "WT.Nutlin.10h","WT.Nutlin.24h")
# 
# pairwise <- list()
# statistics_Ncompartments <- list()
# statistics_Category_compartments <- list()
# statistics_Dynamism_compartments <- list()
# pairwise_tads <- list()
# count_AB_compartments <- list()
# count_AB_compartments_melt <- list()
# compartments_with_tads <- list()
# count_AB_compartments_with_tads <- list()
# 
# for hh
# {
#   pairwise[[i]] <- Compartments[,c("chrom","start", "end","WT.DMSO",i,"Category_WT.DMSO", paste0("Category_",i),"state", "p53", "p53_Category_WT.DMSO", paste0("p53_Category_",i),"tads")]
#   pairwise[[i]]$Category <- do.call(paste, c(pairwise[[i]][,c(6,7)], sep="-"))
#   pairwise[[i]] <- pairwise[[i]][- grep ("NULL", pairwise[[i]]$Category),]
#   
#   statistics_Ncompartments[[i]] <- nrow(pairwise[[i]])
#   statistics_Category_compartments[[i]] <- table(as.data.frame(pairwise[[i]]$Category))
#   statistics_Dynamism_compartments[[i]] <- table(as.data.frame(pairwise[[i]]$state))
# 
#   count_AB_compartments[[i]] <- data.frame(apply(pairwise[[i]][,grep(pattern="^Category_", colnames(pairwise[[i]]))], 2, function(x) { table(x) }))
#   count_AB_compartments[[i]]$compartments <- rownames(count_AB_compartments[[i]])
#   count_AB_compartments_melt[[i]] <- melt(count_AB_compartments[[i]], id="compartments")
#   
#   pairwise[[i]] <- makeGRangesFromDataFrame(pairwise[[i]],seqnames.field = "chrom", start.field = "start",end.field = "end", keep.extra.columns = T)
# 
#   ## Adding the Insulation score of the TADs in the compartments with TAD borders information
#   overlap_Compartment_TADs <- findOverlaps(query = pairwise[[i]] , subject = TADs.gr)
# 
#   # Features from gr1 with overlaps in gr2
#   # Note: The same feature from gr1 can overlap with multiple features from gr2
#   pairwise_tads[[i]] <- pairwise[[i]][queryHits(overlap_Compartment_TADs)];
# 
#   # Add the metadata from gr2
#   mcols(pairwise_tads[[i]]) <- cbind.data.frame(
#     mcols(pairwise_tads[[i]]),
#     mcols(pairwise_tads[[i]][subjectHits(overlap_Compartment_TADs)]));
# 
#   count_AB_compartments_with_tads[[i]] <- table(pairwise_tads[[i]]$Category)
# 
# }
# 
# # Saving environment as HiC.RData in projects

