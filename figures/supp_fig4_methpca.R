library(ggplot2)
library(ggpubr)
library(RColorBrewer)

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

pc <- read.csv('~/research/lfs_wt/meth_pcas/Noob_beta_PCA.csv')
pc <- pc[pc$Meth != "Problem",]
pc <- pc[!is.na(pc$ids),]
cols <- colorRampPalette(brewer.pal(n = 9, 'Set1'))(length(unique(pc$Project)))

pc_450_before <- pc[pc$array == "450",]
ids_keep <- remove_outliers(pc_450_before,3)
pc_450_before <- pc_450_before[pc_450_before$SentrixID %in% ids_keep,]
b450 <- ggplot(pc_450_before,aes(PC1,PC2,color=Project)) + 
    geom_point(size = 3, alpha = 0.7) +
    xlab('PC1') + ylab('PC2') +
    scale_color_manual(name = '', values = cols) + 
    geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
    theme_minimal(base_size = 18)


pc_nci_before <- pc[pc$ids %in% seq(6208,6369,1),]
ids_keep <- remove_outliers(pc_nci_before,3)
pc_nci_before <- pc_nci_before[pc_nci_before$SentrixID %in% ids_keep,]
bnci <- ggplot(pc_nci_before,aes(PC1,PC2,color=Project)) + 
  geom_point(size = 3, alpha = 0.7) +
  xlab('PC1') + ylab('PC2') +
  scale_color_manual(name = '', values = cols) + 
  geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme_minimal(base_size = 18)

pc_850_before <- pc[!pc$ids %in% c(as.character(pc_450_before$ids),as.character(pc_nci_before$ids)),]
ids_keep <- remove_outliers(pc_850_before,2)
pc_850_before <- pc_850_before[pc_850_before$SentrixID %in% ids_keep,]
b850 <- ggplot(pc_850_before,aes(PC1,PC2,color=Project)) + 
  geom_point(size = 3, alpha = 0.7) +
  xlab('PC1') + ylab('PC2') +
  scale_color_manual(name = '', values = cols) + 
  geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme_minimal(base_size = 18)

lfspath <- read.csv('~/research/resources/lfs_mut_clinical_comprehensive.csv',stringsAsFactors = F)
pc_450 <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_450k_PEER_PCA_apcsgt.csv')
ids_keep <- remove_outliers(pc_450,3)
pc_450 <- pc_450[pc_450$SentrixID %in% ids_keep,]
pc_450 <- pc_450[!pc_450$ids %in% unique(lfspath$sample[lfspath$Hospital == "COG"]),]
pc_450 <- pc_450[!duplicated(pc_450$tm_donor),]

a450 <- ggplot(pc_450,aes(PC1,PC2,color=Project)) + 
  geom_point(size = 3, alpha = 0.7) +
  xlab('PC1') + ylab('PC2') +
  scale_color_manual(name = '', values = cols) + 
  geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme_minimal(base_size = 18)

pc_850 <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_850k_PEER_PCA_apcsgt.csv')
ids_keep <- remove_outliers(pc_850,3)
pc_850 <- pc_850[pc_850$SentrixID %in% ids_keep,]
pc_850 <- pc_850[!pc_850$ids %in% unique(lfspath$sample[lfspath$Hospital == "COG"]),]
pc_850 <- pc_850[!duplicated(pc_850$tm_donor),]

a850 <- ggplot(pc_850,aes(PC1,PC2,color=Project)) + 
  geom_point(size = 3, alpha = 0.7) +
  xlab('PC1') + ylab('PC2') +
  scale_color_manual(name = '', values = cols) + 
  geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme_minimal(base_size = 18)


pc_nci <- read.csv('~/research/lfs_wt/meth_pcas/Noob_ComBatProject_beta_nci_PEER_PCA_apcsgt.csv')
ids_keep <- remove_outliers(pc_nci,3)
pc_nci <- pc_nci[pc_nci$SentrixID %in% ids_keep,]
pc_nci <- pc_nci[!pc_nci$ids %in% unique(lfspath$sample[lfspath$Hospital == "COG"]),]
pc_nci <- pc_nci[!duplicated(pc_nci$tm_donor),]

anci <- ggplot(pc_nci,aes(PC1,PC2,color=Project)) + 
  geom_point(size = 3, alpha = 0.7) +
  xlab('PC1') + ylab('PC2') +
  scale_color_manual(name = '', values = cols) + 
  geom_hline(yintercept= 0, linetype="dashed",color = "grey", size=1) +
  geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
  theme_minimal(base_size = 18)


ggarrange(b450, b850, bnci, a450, a850, anci,
          labels = c("450 (before)", "850 (before)", "nci (before)","450 (after)","850 (after)","nci (after)"),
          ncol = 3, nrow = 2)

