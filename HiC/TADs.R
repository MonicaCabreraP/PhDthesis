## Script information ----
## Script name: TADs.R                 
## Purpose: Clean the TADs file and prepare ir for visualization 
## Last modification: 16 Jan 2023  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com

# Notes ----
# If we remove the missing values (NA), we end up with: 2.264 TADborders (2.091 Invariant and 171 Variables)
# 44 TAD-TAD-TAD-TAD-TAD-NoTAD, 27 TAD-NoTAD-TAD-TAD-TAD-TAD, 24 TAD-TAD-NoTAD-TAD-TAD-TAD, 19 NoTAD-TAD-TAD-TAD-TAD-TAD, 11 TAD-TAD-TAD-NoTAD-TAD-TAD, 9 TAD-NoTAD-NoTAD-TAD-TAD-TAD, 9 TAD-TAD-TAD-TAD-NoTAD-TAD, 5 TAD-TAD-NoTAD-NoTAD-TAD-TAD  --> The rest of combinations have less frequency than 5 times
# Of the 27 TAD-NoTAD-TAD-TAD-TAD-TAD: 4 chr7 , 3 chr3, 3 chr5, 3 chrX, 2 chr12, 2 chr2, 2 chr20, 1 chr1, 1 chr10, 1 chr13, 1 chr15, 1 chr16, 1 chr17, 1 chr19, 1 chr8 

## Settings ----
wd <- "/home/mcabrera/"
#wd <- "/home/monica/"
setwd(paste0(wd,"Desktop/MN/projects/p53/data/HCT116/HiC/")) #path where the file is located after processing the HiC data with TADbit 

packages <- c("imputeTS")
invisible(lapply(packages, library, character.only = TRUE))

colors_samples=c("#FDE725FF", "#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF")

## Working with the dataset ----
TADs_TADbitScore <- read_tsv("aligned_TADborders_TADbit",col_names = T) # Loading the TAD borders scores obtained from TADbit and manually aliged by François Serra

head(TADs_TADbitScore) #Visualizing the data loaded

TADs_TADbitScore <- TADs_TADbitScore[,grep(pattern="score_|Chromosome|position", x=colnames(TADs_TADbitScore))] # Taking only the columns that contain the TADbit score from my samples
TADs_TADbitScore <- TADs_TADbitScore[order(TADs_TADbitScore$Chromosome),] # Sorting by chromosome as is disorganized
# TADs_TADbitScore$Chromosome <- sapply(TADs_TADbitScore$Chromosome, gsub, pattern='chr', replacement='') #To remove the chr string if needed

samples <- c("DMSO","p53_1h","p53_4h","p53_7h","p53_10h", "p53_24h")
colnames(TADs_TADbitScore) <- c("Chromosome","TADborder",samples) #Re-naming the columns

head(TADs_TADbitScore) #Visualizing all what I did to the dataset happened

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
head(TADs_TADbitScore_NA20,2) #Visualizing all what I did to the dataset happened

summary(TADs_TADbitScore_NA20)
# DMSO             p53_1h           p53_4h        p53_7h          p53_10h         p53_24h      
# Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  Min.   : 0.000  
# 1st Qu.: 6.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  1st Qu.: 0.000  
# Median :10.000  Median : 7.000  Median : 7.000  Median : 8.000  Median : 9.000  Median : 6.000  
# Mean   : 7.458  Mean   : 5.869  Mean   : 5.747  Mean   : 5.852  Mean   : 6.576  Mean   : 5.221  
# 3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000  3rd Qu.:10.000   
# Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000  Max.   :10.000 
  
TADs_TADbitScore_NA20$individual_TADstate <- data.frame(apply(TADs_TADbitScore_NA20[,-c(1,2)], 2, function(x) { ifelse(x <= 4 , 'NoTAD', 'TAD') }))
TADs_TADbitScore_NA20$Category_tads <- do.call(paste, c(TADs_TADbitScore_NA20$individual_TADstate, sep="-"))

TADs_TADbitScore_NA20$state_tads <-as.character(ifelse(TADs_TADbitScore_NA20$Category_tads == paste(replicate(length(samples), "TAD"), collapse = "-") , 'Invariant', ifelse(a$Category_tads == paste(replicate(length(samples), "NoTAD"), collapse = "-") , 'NoTAD', 'Variable')))

write.table(TADs_TADbitScore_NA20, file="clean_aligned_TADborders_TADbit",sep = "\t", quote = F)
save.image(file='TADs.RData')
