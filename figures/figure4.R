library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)
library(msigdbr)
library(ggpubr)
library(ComplexHeatmap)
library(plyr)
library(circlize)
library(RColorBrewer)
library(ggpubr)

lfsmut <- read.csv('~/research/LFS Mut/lfs_mut_clinical_comprehensive.csv',stringsAsFactors = F)
lfsmut <- lfsmut[lfsmut$WGS %in% c("wgs"),]
lfsmut <- lfsmut[!lfsmut$Hospital %in% c("COG"),]
lfsmut$sample <- str_replace_all(lfsmut$sample,"-","_")
lfsmut <- lfsmut[!is.na(lfsmut$WGS),]
lfswt <- read.csv('~/research/lfs_wt/lfswt_clin.csv',stringsAsFactors = F)
clincols <- c("sample","p53","tm_donor","family","gender","cancer_diagnosis","ageofonset")

lfspath <- rbind.fill(lfsmut,lfswt)
lfspath <- lfspath[order(lfspath$p53,lfspath$cancer_diagnosis),]
lfspath$Class1 <- str_replace_all(paste0(lfspath$P_SNV_AD),"NA;|NA","")
lfspath$Class2 <- str_replace_all(paste0(lfspath$P_SNV_AD, lfspath$P_SNV_AR),"NA;|NA","")
lfspath$Class3 <- str_replace_all(paste0(lfspath$Class2,lfspath$P_SNV_OtherCancerGene,lfspath$P_SNV_Kinase,lfspath$P_SNV_TS,lfspath$SV_CPG),"NA;|NA","")
lfspath$Class4 <- str_replace_all(paste0(lfspath$Class3,lfspath$CancerOnly_Candidate_P_nonCPG),"NA;|NA","")
lfspath$Class5 <- str_replace_all(paste0(lfspath$Class4,lfspath$CancerOnly_VUS_AD,lfspath$CancerOnly_VUS_AR,
                                         lfspath$CancerOnly_VUS_Kinase,lfspath$CancerOnly_VUS_OtherCancerGene,
                                         lfspath$CancerOnly_VUS_TS),"NA;|NA","")
lfspath$Class3 <- str_replace_all(lfspath$Class4,"LEF1;|WNT;|AXIN1;","")
lfspath$Class4 <- str_replace_all(lfspath$Class4,"LEF1;|WNT;|AXIN1;","")
lfspath$Class5 <- str_replace_all(lfspath$Class5,"LEF1;|WNT;|AXIN1;","")
lfspath$All <- str_replace_all(paste0(lfspath$Class3,lfspath$P_SNV_nonCPG,lfspath$SV_nonCPG),"NA;|NA","")
lfspath$n_Class1 <- lengths(regmatches(lfspath$Class1, gregexpr(";", lfspath$Class1)))
lfspath$n_Class2 <- lengths(regmatches(lfspath$Class2, gregexpr(";", lfspath$Class2)))
lfspath$n_Class3 <- lengths(regmatches(lfspath$Class3, gregexpr(";", lfspath$Class3)))
lfspath$n_Class4 <- lengths(regmatches(lfspath$Class4, gregexpr(";", lfspath$Class4)))
lfspath$n_Class5 <- lengths(regmatches(lfspath$Class5, gregexpr(";", lfspath$Class5)))
lfspath$n_All <- lengths(regmatches(lfspath$All, gregexpr(";", lfspath$All)))

lfspath$cancerstatus <- ifelse(lfspath$cancer_diagnosis == "Unaffected","Unaffected","Cancer")
lfspath$cancer_diagnosis <- ifelse(grepl("Breast Cancer|BRCA|DCIS|IDC",lfspath$cancer_diagnosis),"Breast",lfspath$cancer_diagnosis)
lfspath$cancer_diagnosis <- ifelse(grepl("RMS",lfspath$cancer_diagnosis),"RMS",lfspath$cancer_diagnosis)
lfspath$cancer_diagnosis <- ifelse(grepl("PNET|LGG|DA",lfspath$cancer_diagnosis),"Glioma",lfspath$cancer_diagnosis)
lfspath$cancer_diagnosis <- ifelse(grepl("ALL|AML",lfspath$cancer_diagnosis),"Leukemia",lfspath$cancer_diagnosis)
lfspath$cancer_diagnosis <- ifelse(grepl("Basal|Melanoma",lfspath$cancer_diagnosis),"Skin",lfspath$cancer_diagnosis)
lfspath$is_Class3 <- ifelse(lfspath$n_Class3 < 1,"No Class3","Class3")
plotmuts <- melt(lfspath[c("sample","p53","gender","cancerstatus","cancer_diagnosis","is_Class3","ageofonset",colnames(lfspath)[grepl("^n_Class",colnames(lfspath))])],
                 id=c("sample","p53","gender","cancerstatus","cancer_diagnosis","is_Class3","ageofonset"))
lfspath_cancer <- lfspath[lfspath$cancer_diagnosis != "Unaffected",]

############### cancer status vs number of P variants ###############

cancerstatus_plot <- ggplot(plotmuts,aes(x=cancerstatus,y=value)) +
  geom_boxplot(position=position_dodge()) + 
  geom_point(aes(color=cancerstatus),position = position_jitterdodge()) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=16,family="Helvetica Neue")) + 
  labs(x = "Cancer Status", y = "Number of Variants per Patient")+
  facet_grid(~variable+p53) + 
  stat_compare_means(method="wilcox.test",method.args = list(alternative="less"))

############### mutation status vs presence of P variants in cancer genes ###############
mutstatus_plot <- ggplot(plotmuts_tt,aes(x=p53,y=value)) +
  geom_boxplot(position=position_dodge()) + 
  geom_point(aes(color=p53),position = position_jitterdodge()) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=16,family="Helvetica Neue")) + 
  labs(x = "TP53 Status", y = "Number of Variants per Patient")+
  facet_grid(~variable+cancerstatus) + 
  stat_compare_means(method="wilcox.test",method.args = list(alternative="greater"))


############### pathogenic variants in oncogenic gene sets ###############
m_df = data.frame(msigdbr(species = "Homo sapiens",category = "H"))  %>% 
  group_by(gs_name) %>% 
  dplyr::summarise(genes = paste0(human_gene_symbol, collapse = ";|;"),n=n()) %>% 
  dplyr::select(gs_name,genes,n) %>% 
  unique()

m_df$genes <- paste0(";",m_df$genes,';')
m_df[44,]$genes <- paste0(m_df[44,]$genes,"|;SMAD4;"); m_df[40,]$genes <- paste0(m_df[40,]$genes,"|;MTOR;");

############### plot landscape of pathogenic variants ###############

lfspath$P_SNV_Class3 <- str_replace_all(paste0(lfspath$Class2,lfspath$P_SNV_TS,lfspath$P_SNV_Kinase,lfspath$P_SNV_OtherCancerGene),"NA;|NA","")
oncoprint <- unique(melt(cbind(lfspath['sample'],data.frame(do.call(rbind,str_split(lfspath$P_SNV_Class3,"\\;")))),id="sample")[c("sample","value")])
oncoprint <- oncoprint[oncoprint$value != "",]
oncoprint$mut <- "SNV_INDEL"
oncoprint_del <- unique(melt(cbind(lfspath['sample'],data.frame(do.call(rbind,str_split(str_replace_all(lfspath$SV_CPG,"\\(del\\)",""),"\\;")))),id="sample")[c("sample","value")])
oncoprint_del <- oncoprint_del[!is.na(oncoprint_del$value),]
oncoprint_del <- oncoprint_del[oncoprint_del$value != "",]
oncoprint_del$mut <- "DEL"
oncoprint_del <- oncoprint_del[!grepl("dup",oncoprint_del$value),]
oncoprint_dup <- unique(melt(cbind(lfspath['sample'],data.frame(do.call(rbind,str_split(str_replace_all(lfspath$SV_CPG,"\\(CNV dup\\)|\\(cnvdup\\)",""),"\\;")))),id="sample")[c("sample","value")])
oncoprint_dup <- oncoprint_dup[!is.na(oncoprint_dup$value),]
oncoprint_dup <- oncoprint_dup[oncoprint_dup$value != "",]
oncoprint_dup$mut <- "DUP"
oncoprint_dup <- oncoprint_dup[!grepl("del",oncoprint_dup$value),]

oncoprint <- reshape(rbind(oncoprint,oncoprint_del,oncoprint_dup), idvar="sample", timevar = "value",direction = "wide")
oncoprint_clin <- lfspath[c("sample","p53","family","cancer_diagnosis")]
oncoprint_clin <- oncoprint_clin[oncoprint_clin$sample %in% oncoprint$sample,]
rownames(oncoprint) <- oncoprint$sample 
oncoprint_clin <- oncoprint_clin[match(rownames(oncoprint),oncoprint_clin$sample),]
oncoprint$sample <- oncoprint_clin$sample <- NULL
colnames(oncoprint) <- str_replace(colnames(oncoprint),"mut.","")

ha <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                        Family = str_replace_all(oncoprint_clin$family,"family_",""),
                        CancerType = oncoprint_clin$cancer_diagnosis,
                        p53Status=oncoprint_clin$p53,
                        col = list(
                          p53Status = c("Mut" = "#000000", "Wt" = "#E0E0E0"),
                          CancerStatus = c("Cancer"="#000000","Unaffected"="#E0E0E0"),
                          CancerType=c("ACC"="#fa3c5a",
                                       "AA"="#f97171",
                                       "Breast"="#ea6042",
                                       "CPC"="#c89137",
                                       "Colon"="#ffde00",
                                       "CSA"="#003830",
                                       "Glioma"="#81be41",
                                       "High_Grade_Sarcoma"="#00b5b5",
                                       "Leukemia"="#bae3ff",
                                       "LMS" = "#1ad1e5",
                                       "MFH"="#bfbee7",
                                       "MFS"="#84508d",
                                       "MPH"="#fcf1e1",
                                       "OS"="#432156",
                                       "RMS"="#010134",
                                       "Sarcoma"="#ba7979",
                                       "Skin"="#794044",
                                       "Unaffected"="#8b7d7b")))

col = c(SNV_INDEL = "red", DEL = "blue", DUP = "purple")
oncoPrint(t(oncoprint),
          alter_fun = list(
            background = function(x, y, w, h) 
              grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
            SNV_INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["SNV_INDEL"], col = NA)),
            DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                   gp = gpar(fill = col["DUP"], col = NA)),
            DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                gp = gpar(fill = col["DEL"], col = NA))
          ), col = col,top_annotation = ha,
          row_names_gp = gpar(fontsize = 12)
)

