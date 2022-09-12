## Script name: Compartment_figures.R = Fig1.AB_distribution.pdf               
## Last modification: 09 September 2022  		 
## Author: Monica Cabrera-Pasadas
## Email: monica.cabrera.pasadas@gmail.com

# Packages ----
packages <- c("ggplot2", "reshape", "hrbrthemes", "GenomicRanges", "dplyr","ComplexHeatmap","remotes", "gplots", "tidyverse","circlize", "viridis","dplyr", "ggforce","gplots", "ggpubr","RColorBrewer","data.table", "ggrepel","cluster","factoextra")
invisible(lapply(packages, library, character.only = TRUE))

# wd <- "/home/monica/Desktop/MN/"
wd <- "/home/monica/Desktop/MN/"
setwd(paste0(wd,"projects/p53/HiC/Figures"))
load(paste0(wd,"/projects/DATA/HCT116/HiC/HiC.RData"))

# Fig 1.AB_distribution.pdf ----

pdf(file="Fig1.Stacked_AB_general_distribuion.pdf",width = 12, height = 8)

compartment_order <- factor(table_AB_melt$compartment, levels=c("A","B"))

## Stacked barplot all compartments A/B genome distribution ----

### Horizontal ----
sample_order_horizontal <- factor(table_AB_melt$variable, levels=c("WT.Nutlin.24h", "WT.Nutlin.10h", "WT.Nutlin.7h", "WT.Nutlin.4h", "WT.Nutlin.1h", "WT.DMSO"))

print(ggplot(data = table_AB_melt, aes(x = sample_order_horizontal, y = value, fill = compartment_order)) +
        geom_col() +
        geom_text(aes(label = paste0(value," compartments: ", round(value/nrow(Compartments)*100,digits = 2),"%")),color="white", position = position_stack(vjust = .5), size=4)+
        scale_fill_manual(values=c("red","blue"))+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,26973), expand = c(0,0))+     
        theme_classic()+ 
        xlab("")+
        ylab(paste0("Distribution of the total number of compartments: ", nrow(Compartments[- grep("NULL", Compartments$Category),]), " compartments"))+
        ggtitle("A/B genome distribution of all compartments")+
        theme(axis.text.y=element_text(size=rel(1), angle=0),
              axis.text.x=element_text(size=rel(1)),
              axis.title=element_text(size=rel(1),face="bold"))+
        geom_hline(yintercept=table_AB$WT.DMSO[1], color="white")+
        coord_flip())

### Vertical ----
sample_order_vertical <- factor(table_AB_melt$variable, levels=c("WT.DMSO","WT.Nutlin.1h","WT.Nutlin.4h","WT.Nutlin.7h","WT.Nutlin.10h","WT.Nutlin.24h"))

print(ggplot(data = table_AB_melt, aes(x = sample_order_vertical, y = value, fill = compartment_order)) +
        geom_col() +
        geom_text(aes(label = paste0(value,"\n","(", round(value/nrow(Compartments)*100,digits = 2),"% )")),position = position_stack(vjust = .5),color="white", size=8)+
        scale_fill_manual(values=c("red","blue"))+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,25000,26973), expand = c(0,0))+     
        theme_classic()+ 
        ylab("")+
        xlab("")+
        ggtitle(paste0("A/B genome distribution of the total number of compartments: ", nrow(Compartments[- grep("NULL", Compartments$Category_compartments),]), " compartments"))+
        theme(axis.text.y=element_text(size=rel(2), angle=0),
              axis.text.x=element_text(size=rel(2)),
              axis.title=element_text(size=rel(2),face="bold"))+
        geom_hline(yintercept=table_AB$WT.DMSO[1], color="white"))


## Stacked barplot compartments with p53 A/B genome distribution ----
table_AB_p53 <- table_AB_p53_info[grep("A_p53|B_p53", table_AB_p53_info$compartment),]
table_AB_p53_melt <- melt(table_AB_p53, id="compartment")

compartment_order_p53 <- factor(table_AB_p53_melt$compartment, levels=c("A_p53","B_p53"))
### Horizontal ----
print(ggplot(data = table_AB_p53_melt, aes(x = sample_order_horizontal, y = value, fill = compartment_order_p53)) +
        geom_col() +
        geom_text(aes(label = paste0(value,"\n", round(value/nrow(compartments_p53_binding[-grep("no_p53_binding", compartments_p53_binding$p53),])*100,digits = 2),"%")),color="white", position = position_stack(vjust = .5), size=6.5)+
        scale_fill_manual(values=c("red","blue"))+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,5000), expand = c(0,0))+     
        theme_classic()+ 
        xlab("")+
        ylab("")+
        ggtitle(paste("A/B genome distribution of compartments with p53: ","\n", nrow(compartments_p53_binding[-grep("no_p53_binding", compartments_p53_binding$p53),])))+
        theme(axis.text.y=element_text(size=rel(1), angle=0),
              axis.text.x=element_text(size=rel(1)),
              axis.title=element_text(size=rel(1),face="bold"))+
        geom_hline(yintercept=table_AB_p53_melt$WT.DMSO[1], color="white")+
        coord_flip())

### Vertical ----
print(ggplot(data = table_AB_p53_melt, aes(x = sample_order_vertical, y = value, fill = compartment_order_p53)) +
        geom_col()+
        geom_text(aes(label = paste0(value,"\n", round(value/nrow(compartments_p53_binding[-grep("no_p53_binding", compartments_p53_binding$p53),])*100,digits = 2),"%")),color="white", position = position_stack(vjust = .5), size=7)+
        scale_fill_manual(values=c("red","blue"))+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,5000), expand = c(0,0))+ 
        theme_classic()+
        xlab("")+
        ylab("")+
        ggtitle(paste0("A/B genome distribution of the total number of compartments: ", nrow(Compartments[- grep("NULL", Compartments$Category_compartments),]), " compartments"))+
        theme(axis.text.y=element_text(size=rel(2), angle=0),
              axis.text.x=element_text(size=rel(2)),
              axis.title=element_text(size=rel(2),face="bold"))+
        geom_hline(yintercept=table_AB_p53_melt$WT.DMSO[1], color="white"))

# Stacked barplot compartments with and without p53 A/B compartment genome distribution ----
compartment_order_p53_info <- factor(table_AB_p53_info_melt$compartment, levels=c("A","A_p53","B","B_p53"))
### Horizontal ----
sample_order_horizontal_p53_info <- factor(table_AB_p53_info_melt$variable, levels=c("WT.Nutlin.24h", "WT.Nutlin.10h", "WT.Nutlin.7h", "WT.Nutlin.4h", "WT.Nutlin.1h", "WT.DMSO"))

print(ggplot(data = table_AB_p53_info_melt, aes(x = sample_order_horizontal_p53_info, y = value, fill = compartment_order_p53_info)) +
        geom_col() +
        geom_text(aes(label = paste0(value,"\n ", round(value/nrow(Compartments)*100,digits = 2),"%")),color="white", position = position_stack(vjust = .5), size=3)+
        scale_fill_manual(values=c("red","indianred1", "blue", "lightskyblue"))+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,26973), expand = c(0,0))+
        theme_classic()+
        xlab("")+
        ylab(paste0("Distribution of the total number of compartments with p53 information associated: \n compartments without p53 binding ", nrow(Compartments[Compartments$p53 == "no_p53_binding",]), " or with p53 binding ", nrow(Compartments[Compartments$p53 == "p53_binding",]), ")"))+
        ggtitle("A/B genome distribution of all compartments with and without p53 binding")+
        theme(axis.text.y=element_text(size=rel(1), angle=0),
              axis.text.x=element_text(size=rel(1)),
              axis.title=element_text(size=rel(1),face="bold"))+
        geom_hline(yintercept=sum(table_AB_p53_info$WT.DMSO[3]+table_AB_p53_info$WT.DMSO[4]), color="white")+
        coord_flip())

### Vertical ----
sample_order_vertical <- factor(table_AB_p53_info_melt$variable, levels=c("WT.DMSO","WT.Nutlin.1h","WT.Nutlin.4h","WT.Nutlin.7h","WT.Nutlin.10h","WT.Nutlin.24h"))

print(ggplot(data = table_AB_p53_info_melt, aes(x = sample_order_vertical, y = value, fill = compartment_order_p53_info)) +
        geom_col() +
        geom_text(aes(label = paste0(value,"\n ", round(value/nrow(Compartments)*100,digits = 2),"%")),color="white", position = position_stack(vjust = .5), size=3)+
        scale_fill_manual(values=c("red","indianred1", "blue", "lightskyblue"))+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,26973), expand = c(0,0))+
        theme_classic()+
        xlab("")+
        ylab(paste0("Distribution of the total number of compartments with p53 information associated: \n compartments without p53 binding ", nrow(Compartments[Compartments$p53 == "no_p53_binding",]), " or with p53 binding ", nrow(Compartments[Compartments$p53 == "p53_binding",])))+
        ggtitle("A/B genome distribution of all compartments with and without p53 binding")+
        theme(axis.text.y=element_text(size=rel(1), angle=0),
              axis.text.x=element_text(size=rel(1)),
              axis.title=element_text(size=rel(1),face="bold"))+
        geom_hline(yintercept=sum(table_AB_p53_info$WT.DMSO[3]+table_AB_p53_info$WT.DMSO[4]), color="white"))
        
dev.off()
#####
