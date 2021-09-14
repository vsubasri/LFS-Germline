library(ggplot2)
library(ggpubr)
library(extrafont)

###################################### FIGURE 3B ####################################### 

## Cohort comparison for LFS, KiCS, 1000G ##
pfreq <- read.csv('~/research/lfs_wt/PClassificationComparison_1000G_KiCS.csv')
pfreq_stattest <- read.csv('~/research/lfs_wt/PClassificationComparison_1000G_KiCS_stattest.csv') %>% 
  mutate(y.position = c(0.5, 0.6, 0.7, 0.8))

pfreq_plot <- ggplot(pfreq) +
  geom_bar(aes(Cohort,Class.3,fill=Type),stat="identity",color="black",width = 0.9,position="dodge",show.legend = FALSE) +
  scale_fill_manual(values = c("#DE3533","#0047AB","#006644","#006644")) +
  labs(x = "Cohort", y = "Proportion of Patients with Class 3 Mutations")+
  theme_classic() +
  theme(text = element_text(size=16,family="Helvetica Neue")) +
  stat_pvalue_manual(pfreq_stattest, label = "p",size=3) +
  scale_y_continuous(limits = c(0, 1))
  
pfreq_plot

  ###################################### FIGURE 3D ####################################### 

## proprotion of vairants explained ##
ce <- read.csv('~/research/lfs_wt/cancer_explained_rate.csv')
ce$Category <- factor(ce$Category,levels=c("Class 1","Class 2","Class 3","Class 4","Class 5"))
ce_plot <- ggplot(ce,aes(x=TP53_Status,y=Proportion,fill=Category),width = 0.9) +
  geom_bar(stat="identity",position="dodge",color="black") + 
  theme_classic() +
  labs(x = "TP53 Status", y = "Proportion of Patients")+
  theme(text = element_text(size=16,family="Helvetica Neue")) +
  scale_y_continuous(limits = c(0, 1))
ce_plot

###################################### FIGURE 3E ####################################### 

## Actionable Variants ##
av <- read.csv('~/research/lfs_wt/actionable_variants.csv')
av_plot <- ggplot(av,aes(Therapy,Number.of.Patients)) +
  geom_bar(stat="identity",position="dodge",width = 0.9,fill="grey",colour="black") +
  theme_classic()+
  labs(x = "Therapy", y = "Number of Patients")+
  theme(text = element_text(size=16,family="Helvetica Neue")) 
av_plot
