
# Fig 3. Compartment categories distribution----

pdf(file="Fig3.Stacked_dynamic_compartment_distribuion.pdf",width = 12, height = 8)

myColors <- c("red","blue","darkblue","blue","lightskyblue","red","tomato","indianred1","brown3","tomato4","red3","darkred","blue","red","lightskyblue","dodgerblue","salmon","red","blue","red","blue","blue","blue","red","blue","blue","red","red","blue","blue","blue","blue","red","blue","blue","blue","blue","red","blue","blue","red","blue","blue","blue","red","blue","blue","red","red","blue","red", "red","blue","blue","blue","red","red","red","red","red","red","red","red")
Number_dynamic_compartments_Freq_ordered <- as_tibble(dynamic_compartments_Freq[- grep("NULL", dynamic_compartments_Freq$Category_compartments),])
Number_dynamic_compartments_Freq_ordered$Category_compartments <- unfactor(Number_dynamic_compartments_Freq_ordered$Category_compartments)

category_colors <- data.frame(Category_compartments=Number_dynamic_compartments_Freq_ordered$Category_compartments,colors=myColors)

## All compartments ----
Number_dynamic_compartments_Freq_ordered <- merge(Number_dynamic_compartments_Freq_ordered[-c(1:2),],category_colors[-c(1:2),], by="Category_compartments")
Number_dynamic_compartments_Freq_ordered <- Number_dynamic_compartments_Freq_ordered %>% arrange(-Freq)

Number_dynamic_compartments_Freq_ordered$Percentage_switch <- round(Number_dynamic_compartments_Freq_ordered$Freq/nrow(Compartments[(! grepl("NULL", Compartments$Category_compartments)) & Compartments$state_compartments == "dynamic",])*100, digits = 2)
Number_dynamic_compartments_Freq_ordered <- as_tibble(Number_dynamic_compartments_Freq_ordered) %>% arrange(-Freq)
Number_dynamic_compartments_Freq_ordered$Category_compartments <- factor(Number_dynamic_compartments_Freq_ordered$Category_compartments, levels = unique(Number_dynamic_compartments_Freq_ordered$Category_compartments))

print(ggplot(Number_dynamic_compartments_Freq_ordered, aes(x="",y=Percentage_switch,fill=Category_compartments)) + 
  geom_bar(position="stack",stat="identity",width = 0.9 ,color="black") + 
  scale_fill_manual(values = Number_dynamic_compartments_Freq_ordered$colors, 
                    breaks=Number_dynamic_compartments_Freq_ordered$Category_compartments[Number_dynamic_compartments_Freq_ordered$Percentage_switch >= 1])+
  labs(title = "Switching compartments distribution of all compartments",subtitle = "Divided by all possibilities (64)")+
  xlab("")+
  ylab("% of switching compartments")+
  theme_ipsum() +
  theme_classic()+
  theme(legend.position="right",legend.justification="right",
        legend.margin=margin(100,100,100,100),
        legend.box.margin=margin(10,10,10,10))+
  theme(axis.ticks.x=element_blank())+
  theme(aspect.ratio = 2.5)+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(label = paste0(Percentage_switch," %")),position = position_stack(vjust = .5), color="white", size=8)+
  ggtitle(""))

## Compartments with p53 ----
Number_dynamic_compartments_p53_Freq_ordered <- as_tibble(dynamic_compartments_p53_Freq)
Number_dynamic_compartments_p53_Freq_ordered <- merge(Number_dynamic_compartments_p53_Freq_ordered[-c(1:2),],category_colors[-c(1:2),], by="Category_compartments")
Number_dynamic_compartments_p53_Freq_ordered <- Number_dynamic_compartments_p53_Freq_ordered %>% arrange(-Freq)

Number_dynamic_compartments_p53_Freq_ordered$Percentage_switch <- round((Number_dynamic_compartments_p53_Freq_ordered$Freq/nrow(compartments_p53_binding[compartments_p53_binding$state_compartments=="dynamic",])*100),digits = 2)
Number_dynamic_compartments_p53_Freq_ordered <- as_tibble(Number_dynamic_compartments_p53_Freq_ordered) %>% arrange(-Freq)
Number_dynamic_compartments_p53_Freq_ordered$Category_compartments <- factor(Number_dynamic_compartments_p53_Freq_ordered$Category_compartments, levels = unique(Number_dynamic_compartments_p53_Freq_ordered$Category_compartments))

print(ggplot(Number_dynamic_compartments_p53_Freq_ordered, aes(x="",y=Percentage_switch,fill=Category_compartments)) + 
        geom_bar(position="stack",stat="identity",width = 0.9 ,color="black") + 
        scale_fill_manual(values = Number_dynamic_compartments_p53_Freq_ordered$colors, 
                          breaks=Number_dynamic_compartments_p53_Freq_ordered$Category_compartments[Number_dynamic_compartments_p53_Freq_ordered$Percentage_switch >= 1])+
        labs(title = "Switching compartments distribution of all compartments",subtitle = "Divided by all possibilities (64)")+
        xlab("")+
        ylab("% of switching compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="right",legend.justification="right",
              legend.margin=margin(100,100,100,100),
              legend.box.margin=margin(10,10,10,10))+
        theme(axis.ticks.x=element_blank())+
        theme(aspect.ratio = 2.5)+
        scale_y_continuous(expand = c(0,0))+
        geom_text(aes(label = paste0(Percentage_switch," %")),position = position_stack(vjust = .5), color="white", size=8)+
        ggtitle(""))


dev.off()
#####
