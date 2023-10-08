# Script information ----
## Script name: TADs.R                 
## Purpose: Clean the TADs file obtained from TADbit and manual alignment and prepare it for visualization 
## Last modification: 06 Feb 2023  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Notes ----
# Handleing NA values: If we remove the missing values (NA), we end up with: 2.264 TADborders (2.091 Invariant and 171 Variables)
# 44 TAD-TAD-TAD-TAD-TAD-NoTAD, 27 TAD-NoTAD-TAD-TAD-TAD-TAD, 24 TAD-TAD-NoTAD-TAD-TAD-TAD, 19 NoTAD-TAD-TAD-TAD-TAD-TAD, 11 TAD-TAD-TAD-NoTAD-TAD-TAD, 9 TAD-NoTAD-NoTAD-TAD-TAD-TAD, 9 TAD-TAD-TAD-TAD-NoTAD-TAD, 5 TAD-TAD-NoTAD-NoTAD-TAD-TAD  --> The rest of combinations have less frequency than 5 times
# Of the 27 TAD-NoTAD-TAD-TAD-TAD-TAD: 4 chr7 , 3 chr3, 3 chr5, 3 chrX, 2 chr12, 2 chr2, 2 chr20, 1 chr1, 1 chr10, 1 chr13, 1 chr15, 1 chr16, 1 chr17, 1 chr19, 1 chr8 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Settings ----
wd <- "/home/mcabrera/Desktop/MN/p53/results/HCT116/HiC/"
samples <- c("Nut.0h","Nut.1h","Nut.4h","Nut.7h","Nut.10h", "Nut.24h") #Add here the name of the samples to analyze as they appear in the columns of the TADbit file

packages <- c("readr","imputeTS","GenomicRanges")
invisible(lapply(packages, library, character.only = TRUE))
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Loading and cleaning the dataset ----
TADs_TADbitScore <- read_tsv(paste0(wd,"aligned_TADborders_TADbit_score.tsv"),col_names = T) # Loading the TAD borders scores obtained from TADbit and manually aliged by François Serra
head(TADs_TADbitScore) #Visualizing that the data is properly loaded

TADs_TADbitScore <- TADs_TADbitScore %>% select(Chromosome,position,starts_with("score_")) # Taking only the columns that contain the TADbit score from the samples of interest
colnames(TADs_TADbitScore) <- c("Chromosome","TADborder",samples) #Re-naming the columns

chromosome_order<-c(paste("chr",1:22,sep=""),"chrX")
TADs_TADbitScore$Chromosome<-factor(TADs_TADbitScore$Chromosome, levels=chromosome_order)
# TADs_TADbitScore$Chromosome <- sapply(TADs_TADbitScore$Chromosome, gsub, pattern='chr', replacement='') # Remove the chr string (if needed) for following steps 

head(TADs_TADbitScore) #Visualizing that what I did above to the dataset happened in the way I wanted

summary(TADs_TADbitScore)
# DMSO             p53_1h           p53_4h           p53_7h          p53_10h         p53_24h      
# Min.   : 0.000  Min.   : 0.000   Min.   : 0.000   Min.   : 0.000   Min.   : 0.000 Min.   : 0.00 
# 1st Qu.: 7.000  1st Qu.: 7.000   1st Qu.: 7.000   1st Qu.: 8.000   1st Qu.: 8.000 1st Qu.: 7.00 
# Median :10.000  Median :10.000   Median :10.000   Median :10.000   Median :10.000 Median :10.00 
# Mean   : 8.574  Mean   : 8.355   Mean   : 8.466   Mean   : 8.765   Mean   : 8.783 Mean   : 8.41  
# 3rd Qu.:10.000  3rd Qu.:10.000   3rd Qu.:10.000   3rd Qu.:10.000   3rd Qu.:10.000 3rd Qu.:10.00  
# Max.   :10.000  Max.   :10.000   Max.   :10.000   Max.   :10.000   Max.   :10.000 Max.   :10.00  
# NA's   :622     NA's   :1422     NA's   :1535     NA's   :1588     NA's   :1201   NA's   :1812    

# As it is observed, there are missing values (NA's), the way I chose to handle them is replacing them by 0.
TADs_TADbitScore_NA20 <- na_replace(TADs_TADbitScore,0) 
TADs_TADbitScore_NA20 <- TADs_TADbitScore_NA20[rowSums(TADs_TADbitScore_NA20[,c(samples)])!=0,] # Remove all rows that have 0 (NA) at all time points (if they sum 0 they will be removed), otherwise, it means it has a compartment value in at least one point and it can be useful if we want to subset this table

summary(TADs_TADbitScore_NA20) #Visualizing that NA are removed
# DMSO             p53_1h           p53_4h        p53_7h          p53_10h         p53_24h      
# Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  
# 1st Qu.: 6.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  
# Median :10.000  Median : 7.000  Median : 7.000  Median : 8.000  Median : 9.000  Median : 6.000  
# Mean   : 7.458  Mean   : 5.869  Mean   : 5.747  Mean   : 5.852  Mean   : 6.576  Mean   : 5.221  
# 3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000   
# Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000 

TADs_TADbitScore_NA20$individual_TADstate <- data.frame(apply(TADs_TADbitScore_NA20[,-c(1,2)], 2, function(x) { ifelse(x <= 4 , 'NoTAD', 'TAD') }))
TADs_TADbitScore_NA20$Category_tads <- do.call(paste, c(TADs_TADbitScore_NA20$individual_TADstate, sep="-"))

TADs_TADbitScore_NA20$state_tads <-as.character(ifelse(TADs_TADbitScore_NA20$Category_tads == paste(replicate(length(samples), "TAD"), collapse = "-") , 'Invariant', 
                                                       ifelse(a$Category_tads == paste(replicate(length(samples), "NoTAD"), collapse = "-") , 'NoTAD', 'Variable')))

## 1) How many TADs in total do we have for the analysis? ----
nrow(TADs_TADbitScore_NA20) 
# 4.777 TADs

TADs.gr <- makeGRangesFromDataFrame(TADs_TADbitScore_NA20, seqnames.field = "Chromosome", start.field = "TADborder", end.field = "TADborder", keep.extra.columns = T)
write.table(TADs_TADbitScore_NA20, file="aligned_TADborders_TADbit_clean",sep = "\t", quote = F)
save.image(file='TADs.RData')
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #
