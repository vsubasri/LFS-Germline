library(stringr)
library(qvalue)
library(bedr)
library(ggpubr)
library(reshape2)
library(qqman)
library(dplyr)
library(tidyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ggplot2)

############################# Get significant diagnostic probes ############################

discovery <- read.csv('~/research/EWAS_meQTL/EWAS_Noob_ComBatProject_beta_450k_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
discovery$AT_ADJP <- qvalue(discovery$AT_P)$qvalue
discovery <- discovery[discovery$AT_ADJP <0.1,]
validation <- read.csv('~/research/EWAS_meQTL/EWAS_Noob_ComBatProject_beta_850k_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
validation$AT_ADJP <- qvalue(validation$AT_P)$qvalue
validation <- validation[validation$AT_ADJP <0.1,]
extvalidation <- read.csv('~/research/EWAS_meQTL/EWAS_Noob_ComBatProject_beta_nci_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
extvalidation$AT_ADJP <- qvalue(extvalidation$AT_P)$qvalue
extvalidation <- extvalidation[extvalidation$AT_ADJP <0.1,]

######################### Get probes sig in discovery and validation ########################

sig_final <- discovery[discovery$PROBE %in% validation$PROBE,]
sig_final <- sig_final[sig_final$PROBE %in% extvalidation$PROBE,]
discovery <- discovery[discovery$PROBE %in% sig_final$PROBE,]
validation <- validation[validation$PROBE %in% sig_final$PROBE,]
extvalidation <- extvalidation[extvalidation$PROBE %in% sig_final$PROBE,]

########################## Remove probes with contradicting results #########################

disv_probes <- discovery$PROBE[discovery$AT_direction != validation$AT_direction]
disev_probes <- discovery$PROBE[discovery$AT_direction != extvalidation$AT_direction]
vev_probes <- validation$PROBE[validation$AT_direction != extvalidation$AT_direction]
discr_probes <- c(disv_probes,disev_probes,vev_probes)
discovery <- discovery[!discovery$PROBE %in% discr_probes,]
validation <- validation[!validation$PROBE %in% discr_probes,]
extvalidation <- extvalidation[!extvalidation$PROBE %in% discr_probes,]
sig_final <- sig_final[!sig_final$PROBE %in% discr_probes,]

########################## Get enrichment in discovery vs validation #########################

allp <- data.frame(probe=discovery$PROBE,d_ES = discovery$AT_ES, v_ES = validation$AT_ES, ev_ES = extvalidation$AT_ES)
allp$d_ES <- ifelse(discovery$AT_direction == 1,discovery$AT_ES*-1,discovery$AT_ES)
allp$v_ES <- ifelse(validation$AT_direction == 1,validation$AT_ES*-1,validation$AT_ES)
allp$ev_ES <- ifelse(extvalidation$AT_direction == 1,extvalidation$AT_ES*-1,extvalidation$AT_ES)

overlap_mm_top <- read.csv('~/research/EWAS_meQTL/cancerstatus_ComBatPEERbeta_nocog_931_cisCSCE.txt',sep='\t',stringsAsFactors = F)
overlap_mm_top$AT_P <- sig_final$AT_P[match(sig_final$PROBE,overlap_mm_top$probe)]
overlap_mm_top <- overlap_mm_top[overlap_mm_top$qval < 0.05,]

allp$cisCSCE <- ifelse(allp$probe %in% overlap_mm_top$probe,"sig cis-CSCE","non sig cis-CSCE")

################## Plot enrichment score in the discovery vs validation set ##################

v_ES <- ggplot(allp,aes(x=d_ES,y=v_ES)) +
  geom_point(aes(color=cisCSCE)) +
  geom_point(data=subset(allp, cisCSCE == "sig cis-CSCE" ),aes(color=cisCSCE)) + 
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16,family="Helvetica Neue")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=0),linetype="dashed",color="grey") +
  geom_abline(aes(slope=1, intercept=0),linetype="dashed",color="grey") +
  xlab("Discovery Enrichment Score") +
  ylab("Validation Enrichment Score")

ev_ES <- ggplot(allp,aes(x=d_ES,y=ev_ES)) +
  geom_point(aes(color=cisCSCE)) +
  geom_point(data=subset(allp, cisCSCE == "sig cis-CSCE" ),aes(color=cisCSCE)) + 
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16,family="Helvetica Neue")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept=0),linetype="dashed",color="grey") +
  geom_abline(aes(slope=1, intercept=0),linetype="dashed",color="grey") +
  xlab("Discovery Enrichment Score") +
  ylab("External Validation Enrichment Score")

ggarrange(v_ES,ev_ES,
          ncol = 2,nrow=1) 

################################# Plot significant cis-CSCE #################################

data <- readRDS("~/research/EWAS_meQTL/Noob_ComBatProject_beta_PEER_apcsgt.rds")

####################### Format SNV data, remove ./. and reformat #######################

snvs <- read.csv('~/research/EWAS_meQTL/cancerstatus_meQTL_disc_val_allgenes_931_probesX10kbsnvs.input',sep='\t',check.names = F,stringsAsFactors = F)
colnames(snvs)[8:length(snvs)] <- do.call(rbind,str_split(colnames(snvs)[8:length(snvs)],"\\."))[,1]
snvs[snvs=="0"] = NA

#########################################################################################

clin <- read.csv('~/research/EWAS_meQTL/lfs_mut_clinical_comprehensive.csv',stringsAsFactors = F)
cog_ids <- unique(na.omit(clin$tm_donor[clin$Hospital == "COG"]))
patientids <- clin$tm_donor[match(colnames(snvs)[8:length(snvs)],clin$sample)]
patientids <- patientids[!patientids %in% cog_ids]
data <- data[data$tm_donor %in% patientids,]
data <- data %>%distinct(tm_donor, .keep_all = TRUE) %>%distinct(ids, .keep_all = TRUE)
clin <- clin[clin$sample %in% colnames(snvs)[8:length(snvs)],]
data$ids <- clin$sample[match(data$tm_donor,clin$tm_donor)]
data <- data[data$ids %in% colnames(snvs),]
snvs <- snvs[c(colnames(snvs)[1:7],colnames(snvs)[8:length(snvs)][colnames(snvs)[8:length(snvs)] %in% data$ids])]
snvs$snv_id <- paste0("chr",snvs$CHROM, "_",snvs$POS)

############################ Intersect SNVs and probes  ################################
overlap_snvs <- data.frame(t(snvs[8:(length(snvs))]))
colnames(overlap_snvs) <- snvs$snv_id
overlap_snvs <-overlap_snvs[1:(dim(overlap_snvs)[1]-1),]
## Remove SNPs missing calls ## 
overlap_snvs <- overlap_snvs[, which(colMeans(!is.na(overlap_snvs)) > 0.5)]
## Remove SNPs with no variation ## 
overlap_snvs2 <- overlap_snvs <- overlap_snvs[apply(overlap_snvs, 2,function(x) length(unique(x)) > 1)]
overlap_snvs_plot <-overlap_snvs <- as.matrix(overlap_snvs[order(match(rownames(overlap_snvs),data$ids)),])
overlap_mm_top <- read.csv('~/research/EWAS_meQTL/cancerstatus_ComBatPEERbeta_nocog_931_cisCSCE.txt',sep='\t',stringsAsFactors = F)
overlap_mm_top <- overlap_mm_top[overlap_mm_top$qval < 0.05,]

#######################################################################################################
############################################# Plot meQTLs #############################################
#######################################################################################################

overlap_snvs_plot[is.na(overlap_snvs_plot)] = "No call" ; overlap_snvs_plot[overlap_snvs_plot == " 4"] = "A/A" ; overlap_snvs_plot[overlap_snvs_plot == " 3"] = "R/A" ;overlap_snvs_plot[overlap_snvs_plot == " 1"] = "R/R"
#overlap_snvs_plot[overlap_snvs_plot != 1] = "snv" ; overlap_snvs_plot[overlap_snvs_plot == 1] = "no snv"
overlap_snvs_plot <- data.frame(overlap_snvs_plot)
overlap_snvs_plot$ids <- rownames(overlap_snvs_plot)
final_overlap_plot <- merge(data[c("ids","cancerstatus","cancer_diagnosis",as.character(unique(overlap_mm_top$probe)))],overlap_snvs_plot,by="ids")
final_overlap_plot <- final_overlap_plot[final_overlap_plot$ids %in% colnames(snvs)[8:length(snvs)]]

########################### Top cisCSCE in cancer genes ######################################

ETV6 <- ggplot(final_overlap_plot[final_overlap_plot$chr12_11915932 != "No call",],aes(x=chr12_11915932,y=cg23738760)) +
  geom_boxplot() + 
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status"))

ASXL1 <- ggplot(final_overlap_plot,aes(x=chr20_30941503,y=cg12592359)) +
  geom_boxplot() + 
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status"))

TET3 <- ggplot(final_overlap_plot,aes(x=chr2_74280115,y=cg02956499)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue"))+
  guides(color=guide_legend(title="Cancer Status"))

MIR143 <- ggplot(final_overlap_plot,aes(x=chr5_148808390,y=cg04317047)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue"))+
  guides(color=guide_legend(title="Cancer Status"))

LRG1 <- ggplot(final_overlap_plot,aes(x=chr19_4540429,y=cg03882382)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status"))

PARVB <- ggplot(final_overlap_plot,aes(x=chr22_44419010,y=cg11553066)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus),width = 0.1) +
  theme_classic() +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue"))+
  guides(color=guide_legend(title="Cancer Status"))


ggarrange(ETV6,TET3,MIR143,ASXL1,PARVB,LRG1,
          labels = c("ETV6","TET3","MIR143", "ASXL1","PARVB","LRG1"),
          ncol = 2,nrow=3,
          common.legend=TRUE,
          legend = "bottom")

########################## FIGURE 5C ##########################

data <- readRDS("~/research/EWAS_meQTL/Noob_ComBatProject_beta_PEER_apcsgt.rds")
data <- data[!is.na(data$cancerstatus),]
data <- data[!duplicated(data$tm_donor),]
cpdb_genes <- fread("~/research/pathway_analysis/CPDB_pathways_genes.tab", sep='\t') 
wnt_genes <- unlist(str_split(cpdb_genes$hgnc_symbol_ids[cpdb_genes$pathway == "Wnt signaling pathway - Homo sapiens (human)"],','))
wnt_csce <- overlap_mm_top[overlap_mm_top$Gene %in% wnt_genes,]

SMAD3 <- ggplot(data,aes(x=cancerstatus,y=cg17184477)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status"))+
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  labs(x="Cancer Status")

SERPINF <- ggplot(data,aes(x=cancerstatus,y=cg22242539)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  guides(color=guide_legend(title="Cancer Status"))+
  labs(x="Cancer Status")

SFRP4 <- ggplot(data,aes(x=cancerstatus,y=cg09594069)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  guides(color=guide_legend(title="Cancer Status"))+
  labs(x="Cancer Status")

WNT4 <- ggplot(data,aes(x=cancerstatus,y=cg10606730)) +
  geom_boxplot() + 
  geom_jitter(aes(color=cancerstatus)) +
  theme_classic()+
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status"))+
  labs(x="Cancer Status")

ggarrange(WNT4,SMAD3,SFRP4, SERPINF,
          labels = c("WNT4","SMAD3", "SFRP4","SERPINF"),
          ncol = 4,nrow=1,
          common.legend = TRUE,
          legend="bottom")

########################## FIGURE 5D ##########################

LEF1 <- ggplot(data,aes(x=cancerstatus,y=cg03041109)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=cancerstatus)) +
  theme_classic()+
  scale_color_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue")) +
  guides(color=guide_legend(title="Cancer Status")) +
  labs(x="Cancer Status")

fam <- read.csv('~/research/resources/NIH_LFS_clinical.csv')
fam$family_code <- paste0("NIH",fam$family_code)
data$family <- ifelse(data$ids %in% fam$Malkin.ID,fam$family_code[match(data$ids,fam$Malkin.ID)],data$family) 
data_fam <- data[!duplicated(data$tm_donor),]
data_fam <- data_fam[c(colnames(data_fam)[1:44],"cg03041109")]
data_fam_unaff <- data_fam[data_fam$cancerstatus != "Unaffected",]
data_fam_cancer <- data_fam[data_fam$cancerstatus != "Cancer",]
fam_unaff_cancer_ids <- data_fam_unaff$family[data_fam_unaff$family %in% data_fam_cancer$family]
data_fam <- data_fam[data_fam$family %in% fam_unaff_cancer_ids,]
data_fam <- data_fam[data_fam$family %in% names(table(data_fam$family)[table(data_fam$family)>1]),]
data_fam$family <- factor(data_fam$family)

lef1_fam <- ggplot(data_fam,aes(x=family,y=cg03041109,fill=cancerstatus)) +
  geom_bar(stat="identity",position=position_dodge2(width = 1, preserve = "single")) +
  theme_classic()+
  scale_fill_manual(values = c("#DE3533","#0047AB")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue"))+
  guides(color=guide_legend(title="Cancer Status")) +
  labs(x="Family") 

ggarrange(LEF1,lef1_fam,
          ncol=2,nrow=1,
          align="hv",widths=c(1,4),
          common.legend = TRUE,
          legend="bottom")


#######################################################################################################
################################### Calculate GWAS using top meQTLs ###################################
#######################################################################################################

rsSNPs <- read.csv('~/research/EWAS_meQTL/rsSNPs_literature/gwas_catalog_v1.0-associations_e96_r2019-07-30.tsv',sep='\t')
cancer_rsSNPs <- rsSNPs[grepl('osteosarcoma|Osteosarcoma|cancer|Cancer|glioma|Glioma|leukemia|carcinoma|malignancies|lymphoma|Lymphoma|Glioblastoma|Neuroblastoma|tumor',rsSNPs$DISEASE.TRAIT),]
rsSNP_snvids <- overlap_mm_top$ID[overlap_mm_top$ID %in% cancer_rsSNPs$SNPS]

overlap_mm_denovo <- overlap_mm_top[!(overlap_mm_top$ID %in% rsSNP_snvids),] 
overlap_mm_lit <- overlap_mm_top[(overlap_mm_top$ID %in% rsSNP_snvids),] 

#######################################################################################################
################################### Feature enrichment plot ###################################
#######################################################################################################

library(effsize)

scores <- read.csv('~/research/EWAS_meQTL/259_final_scores.txt',sep='\t')
scores$FEATURE <- do.call(rbind,str_split(scores$SAMPLE,"-"))[,2]
scores$SAMPLE <- do.call(rbind,str_split(scores$SAMPLE,"-"))[,1]
predictive <- scores[scores$GROUP =="predictive",]
sampled <- scores[scores$GROUP =="sampled",]

d_all <- list() ; p_all <- list()
for (feat in unique(scores$FEATURE))  {
  print(feat)
  d <- cohen.d(predictive$SCORE[predictive$FEATURE == feat],sampled$SCORE[sampled$FEATURE == feat])$estimate
  p <- wilcox.test(predictive$SCORE[predictive$FEATURE == feat],sampled$SCORE[sampled$FEATURE == feat])$p.value
  p_all[feat] <- p
  d_all[feat] <- d
}

score_stats <- melt(data.frame(feat=unique(scores$FEATURE),
                               d=format(unlist(d_all),digits=2,format="f"),
                               label=formatC(unlist(p_all),digits=2,format="e"),
                               p=-log(unlist(p_all))),id=c("feat","label"))
score_stats$label <- ifelse(score_stats$variable == "p",as.character(score_stats$label),score_stats$value)
score_stats %>% 
  group_by(variable) %>%
  mutate(rescale = scales::rescale(as.numeric(value))) %>%
  ggplot(., aes(x = variable, y =factor(feat))) +
  geom_tile(aes(alpha = rescale, fill = variable,  width=0.95, height=0.95), color = "white") +
  geom_text(aes(label=label,size=16,family="Helvetica Neue")) +
  scale_alpha(range = c(0.1, 1)) +
  scale_fill_manual(values = c("#DE3533","#0047AB")) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        text = element_text(size=16,family="Helvetica Neue"),
        axis.text.x = element_text(angle = 180, vjust = 0.5, hjust=1))



