options(stringsAsFactors = F)
##############################
library(JavierreLabR)
library(GenomicInteractions)
library(ggplot2)
library(gridExtra)
library(data.table) 
library(UpSetR)
##############################

ibeds <- list.files("/home/monica/Desktop/scratch/p53/Step2_CHiCAGO/","*_5.ibed",recursive = T, full.names = T)
#ibeds <- ibeds[-length(ibeds)]

#samples <- file(ibeds,open="r")
#samples <-readLines(samples)

##ATENTION! If warning: In readLines(merged) : incomplete final line found on 'ibeds_merged.txt'
# means you need to put a empty line (click enter in the last row) 

interactions_list <- list()
setwd("/home/monica/Desktop/scratch/p53/Step2_CHiCAGO/")

#for (i in 1:length(samples)){
for (i in ibeds)
{
  print(i)
  name <- basename(i)
  #names processing for BR
  name <- data.frame(do.call('rbind', strsplit(as.character(name),'.',fixed=TRUE)))
  #names processing for merged
  #name <- data.frame(do.call('rbind', strsplit(as.character(name),'_',fixed=TRUE)))
  
  #name <- basename(dirname(i))
  
  #########################################################
  #Step 1: Removing duplicated interactions and changing anotation to transcript id anotation instead of genes
  #########################################################
  #a <- load_interactions(samples[i])
  a <- load_interactions(i)
  a <- annotate_interactions(a, annotation = "/home/monica/Desktop/scratch/Softwares/CHiCAGO/designDir/digest_and_probes_homo_sapiens_hg19_updated_16_02_2020.txt", header =T)
  interactions2ibed(a, paste0(name$X1,"_clean")) #, over.write = T
  
  #########################################################
  #Step 2: Dividing the interactions in cis and trans and saving the cis
  #########################################################
  a_cis <- a[is.cis(a)]
  interactions2ibed(a_cis, paste0(name$X1,"_clean_cis")) #, over.write = T
  
  interactions_list[[paste0(name$X1)]] <- a_cis
}

