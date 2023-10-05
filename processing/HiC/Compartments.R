# Script information ----
## Script name: HiC_Compartments.R                 
## Purpose: Clean the compartments file obtained from TADbit and prepare it for visualization 
## Last modification: 06 Feb 2023  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Notes ----
## Saved files as R Data to use in next scripts:
### Compartments.gr --> File with the compartment values + categories + state (contains NA/(NULL values)) + genes and ENG located in each compartment collapsed in commas
### samples --> File with the samples name to be analyzed to maintain same name in all the analysis
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Settings ---- 
wd <- "/home/mcabrera/Desktop/MN/p53/results/HCT116/HiC/" # Set working directory
samples <- c("Nut.0h","Nut.1h","Nut.4h","Nut.7h","Nut.10h", "Nut.24h") #Add here the name of the samples to analyze as they appear in the columns of the TADbit file

packages <- c("readr","dplyr","imputeTS","GenomicRanges") # Needed libraries
invisible(lapply(packages, library, character.only = TRUE)) # Load needed libraries
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #

# Load the compartments data ----
compartments_file <- paste0(wd, "aligned_Compartments_TADbit.tsv")

if (!file.exists(compartments_file)) {
  stop("Compartments file does not exist.")
}

Compartments <- read_tsv(compartments_file)
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

# Replace missing values with zeros
Compartments_NA20 <- Compartments %>% 
  mutate(across(.cols = all_of(samples), ~if_else(is.na(.), 0, .)))

# Remove rows with all zeros
Compartments_NA20 <- Compartments_NA20 %>% 
  filter_at(vars(all_of(samples)), any_vars(. != 0))

summary(Compartments_NA20[,c(samples)]) #Visualizing that NA are removed
# Nut.0h             Nut.1h             Nut.4h            Nut.7h           Nut.10h            Nut.24h       
# Min.   :-1.41238   Min.   :-1.42527   Min.   :-1.34980   Min.   :-1.27530   Min.   :-1.37269   Min.   :-1.2951  
# 1st Qu.:-0.41926   1st Qu.:-0.36843   1st Qu.:-0.39575   1st Qu.:-0.43577   1st Qu.:-0.49214   1st Qu.:-0.4095  
# Median : 0.12788   Median : 0.16730   Median : 0.15254   Median : 0.16134   Median : 0.09144   Median : 0.0561  
# Mean   :-0.02794   Mean   :-0.01135   Mean   :-0.01836   Mean   :-0.01739   Mean   :-0.04238   Mean   :-0.0510  
# 3rd Qu.: 0.34255   3rd Qu.: 0.34621   3rd Qu.: 0.34218   3rd Qu.: 0.36743   3rd Qu.: 0.38835   3rd Qu.: 0.3273  
# Max.   : 1.17977   Max.   : 1.24379   Max.   : 1.29946   Max.   : 1.26223   Max.   : 1.03526   Max.   : 1.1042  

## Let's categorize the Compartments values in A and B according to their sign (A compartment has positive eigenvector values and B negative)
for (i in samples) { Compartments_NA20[[paste0(i,"_compartment_category")]] <- as.character(ifelse(Compartments_NA20[i] < 0, 'B', 
                                                                                       ifelse(Compartments_NA20[i] > 0, 'A', 'NULL'))) }

Compartments_NA20$all_compartments_combinations <- do.call(paste, c(Compartments_NA20[,grep("category", names(Compartments_NA20), value = TRUE)], sep="-"))
Compartments_NA20$all_compartments_state <- as.character(ifelse(Compartments_NA20$all_compartments_combinations == paste(replicate(length(samples), "A"), collapse = "-") , 'static in A', ifelse(Compartments_NA20$all_compartments_combinations == paste(replicate(length(samples), "B"), collapse = "-") , 'static in B', 'dynamic')))

Compartments.gr <- makeGRangesFromDataFrame(Compartments_NA20, keep.extra.columns = T)
Compartments.gr$location_compartments <- paste0(Compartments.gr@seqnames,":",Compartments.gr@ranges)

## 1) How many Compartments do we have for the analysis? ----
nrow(data.frame(Compartments.gr)) 
### 28120 Compartments with information at some time-point ----

nrow(data.frame(Compartments.gr)[- grep("NULL", Compartments.gr$all_compartments_combinations),])
### 27245 Compartments_NA20 with information at all time-points ----

Compartments.df <- data.frame(Compartments.gr)

# Calculate the number of compartments that correspond to 5% of the total
# n <- nrow(Compartments.df)
# n_top <- ceiling(n * 0.05)
# n_bottom <- floor(n * 0.05)

# Compartments <- list()
# df_top <- list()
# df_bottom <- list()
# Compartments_strong <- list()

for (sample in samples[-1]) 
  {
  Compartments.df[[paste0(sample,"_Nut.0h_diff")]] <- Compartments.df[[sample]] - Compartments.df$Nut.0h

    for (threshold in c("0.1","0.15","0.2"))
    {
  # # Select the top 5% of rows by distribution
  # df_top[[sample]] <- head(Compartments.df[order(-Compartments.df[[paste0(sample,"_Nut.0h_diff")]]), ], n_top)
  #   
  # # Select the bottom 5% of rows by distribution
  # df_bottom[[sample]] <- tail(Compartments.df[order(-Compartments.df[[paste0(sample,"_Nut.0h_diff")]]), ], n_bottom)
  #   
  # Compartments[[sample]] <- rbind(df_top[[sample]],df_bottom[[sample]])
  #   
  # Compartments[[sample]][[paste0(sample,"_Nut.0h_strenght_compartment")]] <- as.character(ifelse(Compartments[[sample]][[paste0(sample,"_Nut.0h_diff")]] > 0, paste0('Top5%_',sample), ifelse(Compartments[[sample]][[paste0(sample,"_Nut.0h_diff")]] < 0, paste0('Bottom5%_',sample), paste0(sample,"_No_strong"))))
  #   
  # Compartments_strong[[sample]] <- as.data.frame(do.call(cbind, Compartments[[sample]]))
  #   
  # Compartments.df[[paste0(sample,"_Nut.0h_strenght_compartment")]] = Compartments_strong[[sample]][[paste0(sample,"_Nut.0h_strenght_compartment")]][match(Compartments.df$location_compartments,Compartments_strong[[sample]]$location_compartments)] #we add to the compartments table the frequency number
  #   
  # Compartments.df[[paste0(sample,"_Nut.0h_strenght_compartment")]] <- Compartments.df[[paste0(sample,"_Nut.0h_strenght_compartment")]] %>% replace_na('none')
  #   
  # Compartments.df[[paste0(sample,"_Nut.0h_strenght_compartment_",threshold)]] <- as.character(ifelse(Compartments.df[[paste0(sample,"_Nut.0h_diff")]] > threshold, paste0('Strong activation ',sample),
                                                                                   # ifelse(Compartments.df[[paste0(sample,"_Nut.0h_diff")]] < -(threshold), paste0('Strong inactivation ',sample),"None")))
  Compartments.df[[paste0(sample,"_Nut.0h_strenght_compartment_",threshold)]] <- as.character(ifelse(Compartments.df[[paste0(sample,"_Nut.0h_diff")]] > threshold, paste0('Strong activation ',sample),
                                                                                                     ifelse(Compartments.df[[paste0(sample,"_Nut.0h_diff")]] < paste0("-",threshold), paste0('Strong inactivation ',sample),"None")))
  
}
}
Compartments.gr <- makeGRangesFromDataFrame(Compartments.df, keep.extra.columns = T)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #
write.table(data.frame(Compartments_NA20), file=paste0(wd,"aligned_Compartments_TADbit_clean.tsv"),sep = "\t", quote = F)
save(Compartments.gr, file = paste0(wd,"Compartments.RData"),compress = T)
saveRDS(samples, paste0(wd,"samples.rds"))
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- #
