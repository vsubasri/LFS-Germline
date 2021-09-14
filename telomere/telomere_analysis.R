library(ggplot2)
library(stringr)
library(ggpubr)
library(ggpmisc)
library(plyr)

mutclin <- read.csv("~/research/resources/lfs_mut_clinical_comprehensive.csv",header=TRUE,stringsAsFactors=FALSE)
mutclin <- mutclin[mutclin$Hospital != "COG",]
mutclin <- mutclin[!is.na(mutclin$WGS),]
wtclin <- read.csv("~/research/lfs_wt/lfswt_clin.csv", header=TRUE,stringsAsFactors=FALSE)
common_col_names <- intersect(names(mutclin), names(wtclin))
clin <- merge(mutclin, wtclin, by=common_col_names, all.x=TRUE, all.y = TRUE)

files <- list.files("~/research/telomere/final_results", "telseq", full.names = T)

all_samples <- list()
for (i in 1:length(files)) {
  print(files[i])
  telseq_sample <- read.csv(files[i], sep='\t')
  telseq_sample$ID <- gsub("\\.telseq","",basename(files[i]))
  all_samples[[i]] <- telseq_sample
}

telseq_all <- do.call(rbind, all_samples)


telseq_all$LENGTH_ESTIMATE <- as.numeric(telseq_all$LENGTH_ESTIMATE)
telseq_all$tumortype <- clin$cancer_diagnosis[match(telseq_all$ID, clin$sample)]
telseq_all$agesamplecollection <- as.numeric(clin$agesamplecollection[match(telseq_all$ID, clin$sample)])
telseq_all$ageofonset <- as.numeric(clin$ageofonset[match(telseq_all$ID, clin$sample)])

telseq_all$agestatus <- ifelse((telseq_all$agesamplecollection < 18*12), "Children", "Adult")
telseq_all$cancerstatus <- ifelse((telseq_all$tumortype == "Unaffected"), "Unaffected", "Cancer")
telseq_all$p53status <- clin$p53[match(telseq_all$ID, clin$sample)]
telseq_all$status <- ifelse((telseq_all$p53status == "Mut"), 
                            ifelse((telseq_all$cancerstatus == "Unaffected"), "MUT_Unaffected", "MUT_Cancer"),ifelse((telseq_all$cancerstatus == "Unaffected"), "WT_Unaffected", "WT_Cancer"))
telseq_all$overallstatus <- ifelse((telseq_all$p53status == "Mut"), 
                                   ifelse((telseq_all$agestatus == "Children"), ifelse((telseq_all$cancerstatus == "Unaffected"), "C_MUT_Unaffected", "C_MUT_Cancer"), 
                                          ifelse((telseq_all$cancerstatus == "Unaffected"), "A_MUT_Unaffected", "A_MUT_Cancer")), 
                                   ifelse((telseq_all$agestatus == "Children"), ifelse((telseq_all$cancerstatus == "Unaffected"), "C_WT_Unaffected", "C_WT_Cancer"), 
                                          ifelse((telseq_all$cancerstatus == "Unaffected"), "A_WT_Unaffected", "A_WT_Cancer")))

telseq_all <- telseq_all[!is.na(telseq_all$overallstatus),]

telseq_mut <- telseq_all[telseq_all$p53status == "Mut", ]
telseq_wt <- telseq_all[telseq_all$p53status == "Wt", ]
telseq_adult <- telseq_all[telseq_all$agestatus == "Adult", ]
telseq_child <- telseq_all[telseq_all$agestatus == "Children", ]

status_mean <- ddply(telseq_all, .(overallstatus), summarise, med = median(agesamplecollection))
  
fig6a <- ggplot(telseq_all[!is.na(telseq_all$overallstatus),], aes(x=overallstatus, y=LENGTH_ESTIMATE,fill=cancerstatus)) + 
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(comparisons = list(c("A_WT_Cancer","A_WT_Unaffected"), c("C_WT_Cancer","C_WT_Unaffected"),
                                        c("C_MUT_Cancer","C_MUT_Unaffected"), c("A_MUT_Cancer","A_MUT_Unaffected"),
                                        c("C_MUT_Cancer","C_WT_Cancer"), c("A_MUT_Cancer","A_WT_Cancer"))) 

fig6b <-ggplot(telseq_all[!is.na(telseq_all$status),], aes(x=status, y=LENGTH_ESTIMATE,fill=cancerstatus)) + 
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(method="wilcox.test",comparisons = list(c("WT_Cancer","WT_Unaffected"),
                                                             c("MUT_Cancer","MUT_Unaffected"))) 
ggarrange(fig6a,fig6b,nrow=2,labels=c("A","B"))

ggplot(telseq_all, aes(x=agesamplecollection, y=LENGTH_ESTIMATE)) + 
  geom_point(aes(color=interaction(cancerstatus,p53status))) + 
  geom_smooth(method = "lm",se=FALSE,aes(color=interaction(cancerstatus,p53status))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title= "", x="Age of Sample Collection", y = "Mean Telomere Lengthn") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16)) +
  stat_fit_glance(method = 'lm', geom = 'text',
                  aes(label = paste("p = ", signif(..p.value.., digits = 4), sep = "")),label.x.npc = 'right', label.y.npc = 1, size = 3)

telseq_ao <- telseq_all[!is.na(telseq_all$ageofonset),]
fig5b <- ggplot(telseq_ao, aes(x=ageofonset, y=LENGTH_ESTIMATE,color=p53status)) + 
  geom_point() + 
  geom_smooth(method = "lm",se=FALSE,aes(color=p53status)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title= "", x="Age of Onset", y = "Mean Telomere Lengthn") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16)) +
  stat_fit_glance(method = 'lm', geom = 'text_npc',
                  mapping = aes(label = sprintf('r^2=%.3f, p=%.2g',
                                                stat(r.squared), stat(p.value))),
            #      aes(label = paste("p = ", signif(..p.value.., digits = 3), sep = "")),
                  label.x = "right", label.y = "top", size = 3,hstep=0,vstep=0.075)


ggplot(telseq_all, aes(x=agesamplecollection, y=LENGTH_ESTIMATE,fill=status)) + 
  geom_point(aes(color=status)) + 
  geom_smooth(method = "lm",aes(color=status),se=F) +
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(color = status, label = paste("P-value = ", signif(..p.value.., digits = 3), sep = "")),label.x.npc = 1, label.y.npc = 1, size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title= "", x="Age of Sample Collection", y = "Mean Telomere Lengthn") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16))

telseq_all$agestatus<- factor(telseq_all$agestatus, levels = c("Children", "Adult"))
telseq_all$p53status<- factor(telseq_all$p53status, levels = c("Mut", "Wt"))
telseq_all$cancerstatus <- factor(telseq_all$cancerstatus, levels = c("Cancer", "Unaffected"))
ggplot(telseq_all, aes(x=agestatus, y=LENGTH_ESTIMATE,color=interaction(p53status,cancerstatus))) + 
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_jitter(position=position_jitter(width=.1, height=0))

telseq_mut_tt <- telseq_mut[telseq_mut$tumortype %in% c("ACC","Breast ca","CPC","RMS","OS","Unaffected"),]
ggplot(telseq_mut_tt, aes(x=tumortype, y=LENGTH_ESTIMATE,color=agestatus)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_jitter(position=position_jitter(width=.1, height=0))

telseq_mut_tt$tumortype <- ordered(telseq_mut_tt$tumortype, levels = c("CPC", "RMS", "Breast ca","OS","ACC","Unaffected"))
ggscatter(telseq_mut_tt, x="agesamplecollection", y="LENGTH_ESTIMATE",
       color = "tumortype", palette = "jco",
       add = "reg.line", sort.by.groups = FALSE, repel = TRUE) + 
  stat_cor(aes(color=tumortype)) + 
  facet_wrap('tumortype') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(telseq_mut_tt, aes(x=tumortype, y=LENGTH_ESTIMATE,color=tumortype)) + 
  geom_boxplot()

go_telomere <- read.table("~/research/telomere/GO_telomere_maintenance.txt",sep='\t', NULL,header=T, row.names = NULL)

