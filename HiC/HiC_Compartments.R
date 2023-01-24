# Script information ----
## Script name: HiC_Compartments.R                 
## Purpose: Clean the compartments file obtained from TADbit and prepare it for visualization 
## Last modification: 24 Jan 2023  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Notes ----
## Instead of removing all rows with NA values (na.omit), I convert them into 0 and I only delete the rows with all samples containing a 0,
## The aim of this strategy is to not lose information when I do the comparision of sometime points as the majority of missing values are at 24h 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Settings ---- 
# wd <- "/home/mcabrera/"
wd <- "/home/monica/"
setwd(paste0(wd,"Desktop/MN/projects/p53/data/HCT116/HiC/")) #path where the file is located after processing the HiC data with TADbit 
samples <- c("DMSO","p53_1h","p53_4h","p53_7h","p53_10h", "p53_24h")

packages <- c("colormap","imputeTS")
invisible(lapply(packages, library, character.only = TRUE))

## Selecting 6 colors from the 'viridis' palette for representing the samples time-course
colors_samples=c("#fde725ff", "#79d051ff","","#26a784ff","#2a768eff","#404284ff","440154ff") # Obtained with: scales::show_col(colormap(colormap = colormaps$viridis,nshades=6), labels = T)
## Selecting 2 colors from the 'RdBu' palette for representing the A-B compartments
color_A_B_Compartments=c("#050aacff","#b20a1cff") # Obtained with: scales::show_col(colormap(colormap = colormaps$RdBu,nshades=2))
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Loading and cleaning the dataset ----
Compartments <- read_tsv("aligned_Compartments_TADbit.tsv")
Compartments <- Compartments[,-c(1,6:7)] # ***Specific step from my data*** Deleting the column with rownames and the KO samples

colnames(Compartments) <- c("Chromosome","Start","End",samples) #Re-naming the columns

summary(Compartments[,c(samples)])
# DMSO             p53_1h            p53_4h            p53_7h           p53_10h           p53_24h       
# Min.   :-1.4124   Min.   :-1.4253   Min.   :-1.3498   Min.   :-1.2753   Min.   :-1.3727   Min.   :-1.2951  
# 1st Qu.:-0.4242   1st Qu.:-0.3710   1st Qu.:-0.3979   1st Qu.:-0.4359   1st Qu.:-0.4959   1st Qu.:-0.4355  
# Median : 0.1308   Median : 0.1689   Median : 0.1530   Median : 0.1615   Median : 0.0971   Median : 0.0760  
# Mean   :-0.0281   Mean   :-0.0114   Mean   :-0.0184   Mean   :-0.0174   Mean   :-0.0426   Mean   :-0.0526  
# 3rd Qu.: 0.3433   3rd Qu.: 0.3465   3rd Qu.: 0.3424   3rd Qu.: 0.3674   3rd Qu.: 0.3896   3rd Qu.: 0.3342  
# Max.   : 1.1798   Max.   : 1.2438   Max.   : 1.2995   Max.   : 1.2622   Max.   : 1.0353   Max.   : 1.1042  
# NA's   :2407      NA's   :2332      NA's   :2296      NA's   :2261      NA's   :2406      NA's   :3128     

# As it is observed above, there are missing values (NA's), so the way I chose to handle them is replacing the NA by 0.
Compartments_NA20 <- na_replace(Compartments,0) # To not delete all rows with one NA at some time point, I first convert the NA to 0
Compartments_NA20 <- Compartments_NA20[rowSums(Compartments_NA20[,c(samples)])!=0,] # Then remove all rows that have 0 (NA) at all time points (if they sum 0 they will be removed), otherwise, it means it has a compartment value in at least one point and it can be useful if we want to subset this table

summary(Compartments_NA20[,c(samples)]) #Visualizing that NA are removed
# DMSO              p53_1h             p53_4h             p53_7h            p53_10h            p53_24h       
# Min.   :-1.41238   Min.   :-1.42527   Min.   :-1.34980   Min.   :-1.27530   Min.   :-1.37269   Min.   :-1.2951  
# 1st Qu.:-0.41926   1st Qu.:-0.36843   1st Qu.:-0.39575   1st Qu.:-0.43577   1st Qu.:-0.49214   1st Qu.:-0.4095  
# Median : 0.12788   Median : 0.16730   Median : 0.15254   Median : 0.16134   Median : 0.09144   Median : 0.0561  
# Mean   :-0.02794   Mean   :-0.01135   Mean   :-0.01836   Mean   :-0.01739   Mean   :-0.04238   Mean   :-0.0510  
# 3rd Qu.: 0.34255   3rd Qu.: 0.34621   3rd Qu.: 0.34218   3rd Qu.: 0.36743   3rd Qu.: 0.38835   3rd Qu.: 0.3273  
# Max.   : 1.17977   Max.   : 1.24379   Max.   : 1.29946   Max.   : 1.26223   Max.   : 1.03526   Max.   : 1.1042  

## Let's categorize the Compartments values in A and B according to their sign (A compartment has positive eigenvector values and B negative)
for (i in samples) { Compartments_NA20[[paste0("Category_",i)]] <- as.character(ifelse(Compartments_NA20[i] < 0, 'B', ifelse(Compartments_NA20[i] > 0, 'A', 'NULL'))) }
Compartments_NA20$combinations <- do.call(paste, c(Compartments_NA20[,grep("Category", names(Compartments_NA20), value = TRUE)], sep="-"))
Compartments_NA20$state_Compartments_NA20 <- as.character(ifelse(Compartments_NA20$combinations == paste(replicate(length(samples), "A"), collapse = "-") , 'static in A', ifelse(Compartments_NA20$combinations == paste(replicate(length(samples), "B"), collapse = "-") , 'static in B', 'dynamic')))

head(Compartments_NA20,2)

## 1) How many Compartments_NA20 do we have for the analysis? ----
nrow(Compartments_NA20) # Compartments with information at some time-point
# 28120
nrow(Compartments_NA20[- grep("NULL", Compartments_NA20$combinations),]) #Compartments_NA20 with information at all time-points
# 27245

Compartments.gr <- makeGRangesFromDataFrame(Compartments_NA20, keep.extra.columns = T)
write.table(Compartments_NA20, file="clean_aligned_Compartments_TADbit",sep = "\t", quote = F)
save.image(file='Compartments.RData')
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #