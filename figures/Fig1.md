---
title: "Figure1_HiC_Compartments "
output:
  pdf_document: default
  date: "`r format(Sys.time(), '%d %B, %Y')`"
  author: "Mónica Cabrera-Pasadas"
---

# Fig 1. Compartment dynamics across TP53 activation

```{r setup-chunk, message=FALSE}
wd <- "/home/mcabrera/Desktop/MN/p53/"
samples <- readRDS(paste0(wd,"samples.rds"))
results <- paste0(wd,"results/HCT116/")  #path where the clean files are

knitr::opts_chunk$set(dev = "pdf",
                      dpi = 300,
                      # echo = FALSE, #to print code or not
                      cache = TRUE,
                      root.dir = paste0(results))
getwd()
packages <- c("ggforce","ggfortify","dendextend","colormap","paletteer", "ggthemes","scico","reshape","ggplot2","GenomicRanges","hrbrthemes","dplyr","tidyverse","ComplexHeatmap","circlize","factoextra", "corrplot")
invisible(lapply(packages, library, character.only = TRUE))

colors_samples=c("#fde725ff", "#37b578ff","#21908dff","#31668dff","#43377fff","#440154ff") 
color_A_B_Compartments=c("#881010","#23447f") 
color_heatmap_compartments <- colorRamp2(breaks = c(-0.5,0,0.5), colors = c("#23447f", "snow1","#881010"))
```

```{r loading_HiC_Compartments_dataset}
load(paste0(results,"HiC/Compartments.RData"))
```

```{r Fig1D_HiC_Compartments_HiC_Compartments_stable_dynamic,fig.cap = " "}
static_dynamic_compartments <- data.frame(sort(table(data.frame(Compartments.gr)[- grep("NULL", Compartments.gr$all_compartments_combinations),]$all_compartments_state), decreasing = T))

colnames(static_dynamic_compartments) <- c("Category", "Total_number")

static_dynamic_compartments$Category <- as.character(static_dynamic_compartments$Category)
static_dynamic_compartments[nrow(static_dynamic_compartments) + 1,] = c('static', sum(static_dynamic_compartments[1:2,2]))
static_dynamic_compartments$Category <- as.factor(static_dynamic_compartments$Category)
static_dynamic_compartments$Category <- factor(static_dynamic_compartments$Category, levels = c("static in A", "dynamic", "static in B"))
static_dynamic_compartments$Total_number <- as.integer(static_dynamic_compartments$Total_number)

ggplot(static_dynamic_compartments[1:3,], aes(x="", y = Total_number, fill=Category)) +
        geom_bar(position="stack",stat="identity",width = 0.2 )+
        scale_fill_manual(values=c("#881010","black", "#23447f"))+
        ggtitle(paste0("Out of the total ",nrow(data.frame(Compartments.gr)[- grep("NULL", Compartments.gr$combinations),])," compartments: \n"))+
  labs(title = "Genomic distribution of the compartments", subtitle = "Divided by stable and switching")+
        xlab("")+
        ylab("Number of compartments")+
        theme_ipsum() +
        theme_classic()+
        theme_minimal()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        geom_text(aes(label = paste0(Category,"\n",Total_number," \n (",round((Total_number/nrow(data.frame(Compartments.gr)[- grep("NULL", Compartments.gr$all_compartments_combinations),])*100),digits=2)," %)")),position = position_stack(vjust = .5), color="white", size=2.5)
```

```{r Fig1D_HiC_Compartments_HiC_Compartments_Freq_dynamism}
compartments_Freq <- data.frame(sort(table(data.frame(Compartments.gr)$all_compartments_combinations),decreasing = T)) 
compartments_Freq$order_combinations <- seq(1:nrow(compartments_Freq)) # we add the frequency rank as ID to match it with the Compartment table with the pourpose of ordering the heatmap by Frequency.

Compartments.gr$Frequency_combinations = compartments_Freq$Freq[match(as.character(data.frame(Compartments.gr)$all_compartments_combinations),compartments_Freq$Var1)] #we add to the compartments table the frequency number

# Compartments.gr$order_combinations = compartments_Freq$order_combinations[match(data.frame(Compartments.gr)$Frequency_combinations,compartments_Freq$Var1)] #we add to the compartments table the frequency order ID

dynamic_compartments_Freq <- data.frame(sort(table(data.frame(Compartments.gr)[! grepl("NULL", Compartments.gr$all_compartments_combinations) & Compartments.gr$all_compartments_state=="dynamic",]$all_compartments_combinations),decreasing = T)) #We check now only the frequency of the dynamic regions

dynamic_compartments_Freq$Percentage <- round(dynamic_compartments_Freq$Freq/nrow(data.frame(Compartments.gr)[(! grepl("NULL", Compartments.gr$all_compartments_combinations)) & Compartments.gr$all_compartments_state == "dynamic",])*100, digits = 2)

dynamic_compartments_Freq$colors <- as.character(ifelse(grepl('^A', dynamic_compartments_Freq$Var1) & dynamic_compartments_Freq$Percentage >=1, '#23447f',ifelse(grepl('^B', dynamic_compartments_Freq$Var1) & dynamic_compartments_Freq$Percentage >=1,'#881010','grey50')))

ggplot(dynamic_compartments_Freq, aes(x="",y=Percentage,col=labels,fill=Var1)) +
        geom_bar(position="stack",stat="identity",width = 1 ,color="black") +
        scale_fill_manual(values = dynamic_compartments_Freq$colors) +
        labs(title = "Switching compartments distribution of all compartments",subtitle = "Divided by all possibilities")+
        xlab("")+
        ylab("% of switching compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(axis.ticks.x=element_blank())+
        theme(aspect.ratio = 2.5)+
        scale_y_continuous(expand = c(0,0))+
        geom_text(aes(label = paste0(Percentage," %")),position = position_stack(vjust = .5), color="white", size=5)+
        # geom_text(data = dynamic_compartments_Freq[dynamic_compartments_Freq$Percentage >= 1,], color = "white") +
        ggtitle("")
```
```{r Chi-Square}
# Chi-square test examines whether rows and columns of a contingency table are statistically significantly associated.
## rows are the different variables, values are the frequencies of the variables done by each category
A_B_distribution <- sapply(data.frame(Compartments.gr)[- grep("NULL", Compartments.gr$all_compartments_combinations),grepl("compartment_category", colnames(data.frame(Compartments.gr)))], function(x) {table(x)})

colnames(A_B_distribution) <- samples

chisq <- chisq.test(t(A_B_distribution))
print(chisq)
corrplot(chisq$residuals, is.cor = FALSE,tl.col="black", method="circle",col=colorRampPalette(c("white","black"))(100), sig.level = 0.005, insig = "blank") # method=pie
```

```{r Fig1E_HiC_Compartments_Heatmap_dynamic_compartments}
dynamics_Freq_ordered <- data.frame(Compartments.gr)[Compartments.gr$all_compartments_state=="dynamic" & !grepl("NULL",Compartments.gr$all_compartments_combinations),] %>% arrange(-Frequency_combinations)

dynamics_Freq_ordered$Percentage_Freq = dynamic_compartments_Freq$Percentage[match(dynamics_Freq_ordered$all_compartments_combinations,dynamic_compartments_Freq$Var1)] #we add to the compartments table the frequency number

dynamics_Freq_ordered$types <-as.character(

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-B-B-B-B-B", 'DMSO specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-A-A-A-A-A", 'p53 specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-A-B-B-B-B", 'p53 1h specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-B-A-A-A-A", 'p53 1h specific inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-B-A-B-B-B", 'p53 4h specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-A-B-A-A-A", 'p53 4h specific inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-B-B-A-B-B", 'p53 7h specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-A-A-B-A-A", 'p53 4h specific inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-B-B-B-A-B", 'p53 10h specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-A-A-A-B-A", 'p53 10h specific inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-B-B-B-B-A", 'p53 24h specific activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-A-A-A-A-B", 'p53 24h specific inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-A-A-A-B-B" | dynamics_Freq_ordered$all_compartments_combinations == "B-A-A-B-B-B", 'Early activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-B-B-B-A-A" | dynamics_Freq_ordered$all_compartments_combinations == "A-B-B-A-A-A", 'Early inactivation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "B-B-B-B-A-A", 'Late activation',

    ifelse(dynamics_Freq_ordered$all_compartments_combinations == "A-A-A-A-B-B", 'Late inactivation','Highly dynamic')))))))))))))))))

# dynamics_Freq_ordered$types[dynamics_Freq_ordered$Percentage_Freq <9] <- "<10%"

Compartments_types_freq <- data.frame(sort(table(dynamics_Freq_ordered$types),decreasing = T))
Compartments_types_freq$type_freq <- Compartments_types_freq$Freq*100/nrow(dynamics_Freq_ordered)
Compartments_types_freq$name_freq <- as.character(Compartments_types_freq$Var1)
Compartments_types_freq$name_freq[Compartments_types_freq$type_freq <5] <- "<5%"

dynamics_Freq_ordered$heatmap = Compartments_types_freq$type_freq[match(dynamics_Freq_ordered$types,Compartments_types_freq$Var1)] #we add to the compartments table the frequency number
dynamics_Freq_ordered$heatmap_types = Compartments_types_freq$name_freq[match(dynamics_Freq_ordered$types,Compartments_types_freq$Var1)]

dynamics_Freq_ordered <- dynamics_Freq_ordered[order(dynamics_Freq_ordered$heatmap, decreasing = T),]

dynamics_Freq_ordered$caca <- paste0(dynamics_Freq_ordered$heatmap_types,"\n",round(dynamics_Freq_ordered$heatmap, digits = 2),"%")

Heatmap(as.matrix(dynamics_Freq_ordered[,c(samples)]),column_title = "Heatmap dynamic regions", column_title_gp = gpar(fontsize = 8, fontface = "bold"),
cluster_rows = F, show_row_names = F,row_split=factor(dynamics_Freq_ordered$types, levels = c("Early activation","p53 1h specific activation","p53 10h specific inactivation","p53 24h specific inactivation","Late inactivation","Highly dynamic","<5%")),
row_title_rot = 0,row_title_gp = gpar(col = c("black"),fontsize=5),row_gap = unit(1.5, "mm")
,border = FALSE,border_gp = gpar(col = "black", lwd = 2),
cluster_columns=F,col = color_heatmap_compartments)
```

```{r Fig1F_HiC_Compartments_Region_chr12_compartment_dynamism, warning=FALSE,fig.cap = "caca"}
region <- melt(data.frame(Compartments.gr)[! grepl("NULL", Compartments.gr$all_compartments_combinations) & data.frame(Compartments.gr)$seqnames == "12" & data.frame(Compartments.gr)$start >= 55000000 & data.frame(Compartments.gr)$end < 95000000,][,c(samples,"start")], id.vars = c("start"))
region$colour <- ifelse(region$value < 0,"negative","positive")

ggplot(region,aes(start,value))+
  geom_bar(stat="identity",position="identity",aes(fill = colour),size=10)+
  scale_fill_manual(values=c(positive="#881010",negative="#23447f"))+
  theme_classic()+
  ylim(-1,1)+
  xlab("Genomic position (in bp) of the chromosome 12")+
  facet_wrap_paginate(. ~ variable, nrow = 6, ncol = 1)

```

```{r Fig1F_HiC_Compartments_Diffregion_chr12_compartment_dynamism}
diff_region <- melt(data.frame(Compartments.gr)[! grepl("NULL", Compartments.gr$all_compartments_combinations) & data.frame(Compartments.gr)$seqnames == "12" & data.frame(Compartments.gr)$start >= 55000000 & data.frame(Compartments.gr)$end < 95000000,grep(pattern = "_diff|start", x=colnames(data.frame(Compartments.gr)))], id.vars = c("start"))

diff_region$colour <- ifelse(diff_region$value < 0,"negative","positive")

ggplot(diff_region,aes(start,value))+
  geom_bar(stat="identity",position="identity",aes(fill = colour),size=10)+
  scale_fill_manual(values=c(positive="#881010",negative="#23447f"))+
  theme_classic()+
  ylim(-1,1)+
  xlab("Genomic position (in bp) of the chromosome 12")+
  facet_wrap_paginate(. ~ variable, nrow = 6, ncol = 1)
```
