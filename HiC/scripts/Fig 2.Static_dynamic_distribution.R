# Fig 2.Static_dynamic_distribution.pdf ----

pdf(file="Fig2.Stacked_static_dynamic_distribuion.pdf",width = 12, height = 8)

## AB distributions separated by static or dynamic compartments ----
static_dynamic_compartments$Category <- as.character(static_dynamic_compartments$Category)
static_dynamic_compartments[nrow(static_dynamic_compartments) + 1,] = c('static', sum(static_dynamic_compartments[2:3,2]))
static_dynamic_compartments$Category <- as.factor(static_dynamic_compartments$Category)
static_dynamic_compartments$Total_number <- as.integer(static_dynamic_compartments$Total_number)

print(ggplot(static_dynamic_compartments[c(3,4),], aes(x="", y = Total_number, fill=Category)) +     
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        # labs(title = "Genomic distribution of the compartments", subtitle = "Divided by stable and switching")+
        xlab("")+
        ylab("Number of compartments")+
        ggtitle(paste0("Out of the total ",nrow(Compartments[- grep("NULL", Compartments$Category_compartments),])," compartments with values in all time points: \n"))+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,23500,26973), expand = c(0,0))+
        theme(aspect.ratio = 3.5)+
        geom_text(aes(label = paste0(static_dynamic_compartments[c(3,4),]$Total_number," \n",round((static_dynamic_compartments[c(3,4),]$Total_number/nrow(Compartments[- grep("NULL", Compartments$Category_compartments),])*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=8))

## AB distributions separated by static in A, static in B or dynamic compartments ----
order_dynamism <- factor(static_dynamic_compartments[1:3,]$Category, levels=c("dynamic", "static in A", "static in B"))

print(ggplot(static_dynamic_compartments[1:3,], aes(x="", y = Total_number, fill=order_dynamism)) +     
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        # labs(title = "Genomic distribution of the compartments", subtitle = "Divided by stable  and switching")+
        xlab("")+
        ylab("Number of compartments")+
        ggtitle(paste0("Out of the total ",nrow(Compartments[- grep("NULL", Compartments$Category_compartments),])," compartments with values in all time points: \n"))+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,5000,10000, 15000,20000,25000,26973), expand = c(0,0))+
        theme(aspect.ratio = 2.5)+
        geom_text(aes(label = paste0(Total_number," \n",round((Total_number/nrow(Compartments[- grep("NULL", Compartments$Category_compartments),])*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=5))


### compartments with p53 binding ----
#### AB distributions separated by static or dynamic compartments ----
p53_state <- data.frame(table(Compartments[- grep("NULL", Compartments$Category_compartments),]$state, Compartments[- grep("NULL", Compartments$Category_compartments),]$p53))
colnames(p53_state) <- c("Category", "binding","Total_number")

p53_state$Category <- as.character(p53_state$Category)
p53_state$binding <- as.character(p53_state$binding)
p53_state[nrow(p53_state) + 1,] = c('static', 'no_p53_binding' , sum(p53_state[2:3,3]))

p53_state$Total_number <- as.integer(p53_state$Total_number)
p53_state[nrow(p53_state) + 1,] = c('static', 'p53_binding' , sum(p53_state[5:6,3]))
p53_state$Total_number <- as.integer(p53_state$Total_number)
p53_state$binding <- as.factor(p53_state$binding)
p53_state$Category <- as.factor(p53_state$Category)

order_dynamism_p53 <- factor(p53_state[p53_state$binding == "p53_binding",]$Category, levels=c("dynamic", "static in A", "static in B"))

print(ggplot(p53_state[p53_state$binding == "p53_binding",][c(1,4),], aes(x="", y = Total_number, fill=Category)) +
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        labs(title = "Genomic distribution of the dynamic compartments", subtitle = "Divided by p53 binding")+
        xlab("")+
        ylab("Number of compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,sum(p53_state[p53_state$binding == "p53_binding",][c(1,4),]$Total_number)), expand = c(0,0))+
        theme(aspect.ratio = 2.5)+
        geom_text(aes(label = paste0(p53_state[p53_state$binding == "p53_binding",][c(1,4),]$Total_number," \n",round((p53_state[p53_state$binding == "p53_binding",][c(1,4),]$Total_number/sum(p53_state[p53_state$binding == "p53_binding",][c(1,4),]$Total_number)*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=7)+
        ggtitle(paste0("Out of the total ",sum(p53_state[p53_state$binding == "p53_binding",][c(1,4),]$Total_number)," compartments with p53 binding: \n")))

#### AB distributions separated by static static in A, static in B or dynamic compartments ----
print(ggplot(p53_state[p53_state$binding == "p53_binding",][c(1:3),], aes(x="", y = Total_number, fill=Category)) +
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        labs(title = "Genomic distribution of the dynamic compartments", subtitle = "Divided by p53 binding")+
        xlab("")+
        ylab("Number of compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,sum(p53_state[p53_state$binding == "p53_binding",][c(1:3),]$Total_number)), expand = c(0,0))+
        theme(aspect.ratio = 2.5)+
        geom_text(aes(label = paste0(p53_state[p53_state$binding == "p53_binding",][c(1:3),]$Total_number," \n",round((p53_state[p53_state$binding == "p53_binding",][c(1:3),]$Total_number/sum(p53_state[p53_state$binding == "p53_binding",][c(1:3),]$Total_number)*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=7)+
        ggtitle(paste0("Out of the total ",sum(p53_state[p53_state$binding == "p53_binding",][c(1:3),]$Total_number)," compartments with p53 binding: \n")))

### compartments without p53 binding  ----
#### AB distributions separated by static or dynamic compartments ----

order_dynamism_p53 <- factor(p53_state[p53_state$binding == "p53_binding",]$Category, levels=c("dynamic", "static in A", "static in B"))

print(ggplot(p53_state[p53_state$binding == "no_p53_binding",][c(1,4),], aes(x="", y = Total_number, fill=Category)) +
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        labs(title = "Genomic distribution of the dynamic compartments", subtitle = "Divided by p53 binding")+
        xlab("")+
        ylab("Number of compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,sum(p53_state[p53_state$binding == "no_p53_binding",][c(1,4),]$Total_number)), expand = c(0,0))+
        theme(aspect.ratio = 2.5)+
        geom_text(aes(label = paste0(p53_state[p53_state$binding == "no_p53_binding",][c(1,4),]$Total_number," \n",round((p53_state[p53_state$binding == "no_p53_binding",][c(1,4),]$Total_number/sum(p53_state[p53_state$binding == "no_p53_binding",][c(1,4),]$Total_number)*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=7)+
        ggtitle(paste0("Out of the total ",sum(p53_state[p53_state$binding == "no_p53_binding",][c(1,4),]$Total_number)," compartments with p53 binding: \n")))

#### AB distributions separated by static static in A, static in B or dynamic compartments ----
print(ggplot(p53_state[p53_state$binding == "no_p53_binding",][c(1:3),], aes(x="", y = Total_number, fill=Category)) +
        geom_bar(position="stack",stat="identity",width = 0.9 )+
        scale_fill_manual(values=c("black","grey60", "grey80"))+
        labs(title = "Genomic distribution of the dynamic compartments", subtitle = "Divided by p53 binding")+
        xlab("")+
        ylab("Number of compartments")+
        theme_ipsum() +
        theme_classic()+
        theme(legend.position="bottom")+
        theme(axis.ticks.x=element_blank())+
        scale_y_continuous(breaks = c(0,1000,2000,3000,4000,sum(p53_state[p53_state$binding == "no_p53_binding",][c(1:3),]$Total_number)), expand = c(0,0))+
        theme(aspect.ratio = 2.5)+
        geom_text(aes(label = paste0(p53_state[p53_state$binding == "no_p53_binding",][c(1:3),]$Total_number," \n",round((p53_state[p53_state$binding == "no_p53_binding",][c(1:3),]$Total_number/sum(p53_state[p53_state$binding == "no_p53_binding",][c(1:3),]$Total_number)*100),digits=2)," %")),position = position_stack(vjust = .5), color="white", size=7)+
        ggtitle(paste0("Out of the total ",sum(p53_state[p53_state$binding == "no_p53_binding",][c(1:3),]$Total_number)," compartments with p53 binding: \n")))

dev.off()
#####
