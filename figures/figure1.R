library(ggplot2)
library(ggpubr)
library(dplyr)  

###########################################################################################
######################################## FUNCTIONS ########################################
###########################################################################################

remove_outliers <- function(pc,n) {
  pc$PC1 <- scale(pc$PC1)
  pc$PC2 <- scale(pc$PC2)
  u_thres_pc1 <- mean(pc$PC1) + n*sd(pc$PC1)
  l_thres_pc1 <- mean(pc$PC1) - n*sd(pc$PC1)
  u_thres_pc2 <- mean(pc$PC2) + n*sd(pc$PC2)
  l_thres_pc2 <- mean(pc$PC2) - n*sd(pc$PC2)
  pc2 <- pc[pc$PC1 < u_thres_pc1 & pc$PC1 > l_thres_pc1, ]
  pc2 <- pc2[pc2$PC2 < u_thres_pc2 & pc2$PC2 > l_thres_pc2, ]
  outlier_id <- as.character(pc$SentrixID[!pc$SentrixID %in% pc2$SentrixID])
  cat(paste0(length(outlier_id)," outliers removed","\n"))
  keep_id <- as.character(pc2$SentrixID)
  return(keep_id)
}

div_age <- function(clin.dat) {
  clin.plot <- clin.dat[!is.na(clin.dat$ageofonset),]
  mean_age = mean(clin.plot$ageofonset)
  clin.plot$type <- ifelse(clin.plot$ageofonset < mean_age, "below", "above")
  clin.plot <- clin.plot[order(clin.plot$ageofonset), ]
  clin.plot$sample <- factor(clin.plot$sample, levels = clin.plot$sample) 
  ggplot(clin.plot, aes(x=sample, y=scale(ageofonset,scale=F), label=scale(ageofonset,scale=F))) + 
    geom_bar(position = 'identity',stat='identity', aes(fill=type), width=.5)  +
    scale_fill_manual(name="Age of onset", 
                      labels = c("Above Average", "Below Average"), 
                      values = c("above"="#DE3533", "below"="#0047AB")) + 
    labs(subtitle=paste("Mean age =", mean(clin.plot$ageofonset), "(n = ",length(clin.plot$ageofonset) ,")"), 
         title= "Distribution of Age of Onset") + 
    coord_flip() + 
    geom_hline(aes(yintercept = 0 ), linetype="dashed")+
    theme_bw() 
}

plot_age <- function(clin.plot) {
  plt <- ggplot(clin.plot, aes(x=ageofonset, label=ageofonset)) + 
    geom_density(aes(fill=cancer_diagnosis), alpha = 0.5) +
    theme_classic() +
    labs(x="Age of Diagnosis (months)", y = "Proportion of individuals") + 
    theme(text = element_text(size=16,family="Helvetica Neue")) + 
    guides(fill=guide_legend("Cancer Type"))
  ## + labs(subtitle=paste0("n = ",length(clin.plot$ageofonset)),title= "Distribution of Age of Diagnosis",)
  return(plt)
}

plot_tt <- function(clin.plot) {
  plt <- ggplot(clin.plot[clin.plot$gender %in% c("M","F"),]) + 
    geom_bar(aes(x=cancer_diagnosis, ..count..,fill = gender),color="black",width = 0.9, position = "dodge",alpha=0.8) + 
    scale_fill_manual(values = c("#DE3533","#0047AB")) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size=16,family="Helvetica Neue")) + 
    labs(x="Tumor Type", y = "Number of individuals") 
  ## + labs(subtitle=paste("n =",length(clin.plot$cancer_diagnosis)), title= "Distribution of Cancer Diagnosis",)
  return(plt)
}

### summary of families, grouped by cancer status and p53 mut 
plot_family <- function(clin.all) {
  clin.all$cancerstatus <- ifelse(clin.all$cancer_diagnosis == "Unaffected","Unaffected","Cancer")
  plt <- ggplot(clin.all, aes(x=family, ..count..)) + 
    scale_fill_manual(values = c("#DE3533","#0047AB","#006644","#10C25B")) +
    geom_bar(aes(fill =interaction(cancerstatus,p53)), width = 0.9,position = "dodge",alpha=0.8) + 
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          text = element_text(size=16,family="Helvetica Neue")) + 
    labs(x="Family", y = "Number of individuals") +
    guides(fill=guide_legend(title="Cancer and p53 Status")) 
  ## + labs(subtitle=paste("n =",length(unique(clin.all$family))), title= "Distribution of Families",)
  return(plt)
}

###########################################################################################
######################## PLOT DISTRIBUTION OF AGE, TT AND FAMILY ##########################
###########################################################################################

## read in and clean methylation clinical data
lfs.clin <- read.csv('~/research/resources/lfs_mut_clinical_comprehensive.csv',stringsAsFactors = F)

## identify samples to remove
cog_ids <- lfs.clin$sample[lfs.clin$Hospital %in% c("COG")]

## load and clean pca data for 450k methylation 
pc_450 <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_450k_PEER_PCA_apcsgt.csv',stringsAsFactors = F)
pc_450 <- pc_450[!pc_450$ids %in% cog_ids,]
pc_450 <- pc_450[!duplicated(pc_450$tm_donor),]
ids_keep <- remove_outliers(pc_450,3)
pc_450 <- pc_450[pc_450$SentrixID %in% ids_keep,]

## load and clean pca data for 850k methylation 
pc_850 <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_850k_PEER_PCA_apcsgt.csv',stringsAsFactors = F)
pc_850 <- pc_850[!pc_850$ids %in% cog_ids,]
pc_850 <- pc_850[!duplicated(pc_850$tm_donor),]
ids_keep <- remove_outliers(pc_850,3)
pc_850 <- pc_850[pc_850$SentrixID %in% ids_keep,]

## load and clean pca data for nci methylation 
pc_nci <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_nci_PEER_PCA_apcsgt.csv',stringsAsFactors = F)
pc_nci <- pc_nci[!pc_nci$ids %in% cog_ids,]
pc_nci <- pc_nci[!duplicated(pc_nci$tm_donor),]
ids_keep <- remove_outliers(pc_nci,3)
pc_nci <- pc_nci[pc_nci$SentrixID %in% ids_keep,]
nci_fam <- read.csv('~/research/resources/NIH_LFS_clinical.csv')
pc_nci$family <- paste0("NCI_",nci_fam$family[match(pc_nci$ids,nci_fam$Malkin.ID)])

## combine
clin.all <- do.call(rbind,list(pc_450[1:44],pc_850[1:44],pc_nci[1:44]))
clin.all <- clin.all[!duplicated(clin.all$tm_donor),]

### summary of families, grouped by cancer status and p53 mut 
clin.all$family <- ifelse(is.na(clin.all$family),clin.all$tm_donor,clin.all$family)
clin.all <- clin.all[!clin.all$ids %in% cog_ids,]
clin.tt <- clin.all[clin.all$cancer_diagnosis %in% c("ACC","Breast","CPC","OS","RMS","Unaffected","Glioma","Leukemia","Skin"),]
colnames(clin.all)[1] <- "sample"
clin.all$cancer_diagnosis <- ifelse(clin.all$cancer_diagnosis == "High Grade Sarcoma","Sarcoma",clin.all$cancer_diagnosis) 

## plots methylation summary data
meth.age <- plot_age(clin.tt)
meth.tt <- plot_tt(clin.all)
meth.fam <- plot_family(clin.all)

## read in and clean wgs clinical data
lfs.wgs.mut <- lfs.clin[lfs.clin$WGS %in% c("wgs"),]
lfs.wgs.mut <- lfs.wgs.mut[!lfs.wgs.mut$Hospital %in% c("COG"),]
lfs.wgs.wt <- read.csv('~/research/lfs_wt/lfswt_clin.csv',stringsAsFactors = F)
clincols <- c("sample","p53","family","gender","tm_donor","cancer_diagnosis","ageofonset")
lfs.wgs <- rbind(lfs.wgs.mut[clincols],lfs.wgs.wt[clincols])
lfs.wgs$cancer_diagnosis <- ifelse(grepl("Breast Cancer|BRCA|DCIS|IDC",lfs.wgs$cancer_diagnosis),"Breast",lfs.wgs$cancer_diagnosis)
lfs.wgs$cancer_diagnosis <- ifelse(grepl("RMS",lfs.wgs$cancer_diagnosis),"RMS",lfs.wgs$cancer_diagnosis)
lfs.wgs$cancer_diagnosis <- ifelse(grepl("PNET|LGG|DA",lfs.wgs$cancer_diagnosis),"Glioma",lfs.wgs$cancer_diagnosis)
lfs.wgs$cancer_diagnosis <- ifelse(grepl("ALL|AML",lfs.wgs$cancer_diagnosis),"Leukemia",lfs.wgs$cancer_diagnosis)
lfs.wgs$cancer_diagnosis <- ifelse(grepl("Basal|Melanoma",lfs.wgs$cancer_diagnosis),"Skin",lfs.wgs$cancer_diagnosis)
lfs.wgs$cancerstatus <- ifelse(lfs.wgs$cancer_diagnosis == "Unaffected","Unaffected","Cancer")
lfs.wgs.tt <- lfs.wgs[lfs.wgs$cancer_diagnosis %in% c("ACC","Breast","CPC","OS","RMS","Unaffected","Glioma","Leukemia","Skin"),]

## plots wgs summary data
wgs.age <- plot_age(lfs.wgs.tt)
wgs.tt <- plot_tt(lfs.wgs)
wgs.fam <- plot_family(lfs.wgs)

###########################################################################################
###################################### COMBINE PLOTS ######################################
###########################################################################################

ggarrange(wgs.fam,meth.fam,wgs.age,meth.age,wgs.tt,meth.tt,
          nrow=3,ncol=2,common.legend = TRUE, align = "v",
          labels=c("B","C","D","E","F","G"),
          legend = "bottom")


ggarrange(wgs.fam, meth.fam,
          nrow=1,ncol=2,common.legend = TRUE,
          labels=c("B","C"),
          font.label=list(color="black",size=18),
          legend = "bottom")

ggarrange(wgs.age,meth.age,
          nrow=1,ncol=2,common.legend = TRUE,
          labels=c("D","E"),
          font.label=list(color="black",size=18),
          legend = "bottom")

ggarrange(wgs.tt,meth.tt,
          nrow=1,ncol=2,common.legend = TRUE, align = "h",
          labels=c("F","G"),
          font.label=list(color="black",size=18),
          legend = "bottom")
