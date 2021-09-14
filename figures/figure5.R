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

############### pathogenic variants in oncogenic gene sets ###############
m_df = data.frame(msigdbr(species = "Homo sapiens",category = "H"))  %>% 
  group_by(gs_name) %>% 
  dplyr::summarise(genes = paste0(human_gene_symbol, collapse = ";|;"),n=n()) %>% 
  dplyr::select(gs_name,genes,n) %>% 
  unique()

m_df$genes <- paste0(";",m_df$genes,';')
m_df[44,]$genes <- paste0(m_df[44,]$genes,"|;SMAD4;"); m_df[40,]$genes <- paste0(m_df[40,]$genes,"|;MTOR;");

lfspath$search <- paste0(";",str_replace_all(lfspath$All,"\\(del\\)|\\(dup\\)|\\(cnvdup\\)|\\(CNV dup\\)",""))
m_df_hits <- data.frame(t(do.call(rbind,lapply(m_df$genes, function(i) lengths(regmatches(lfspath$search,gregexpr(i,lfspath$search)))))))
colnames(m_df_hits) <- m_df$gs_name
m_df_hits <- m_df_hits[colSums(m_df_hits) > 0]
msigdb_muts <- cbind(lfspath,m_df_hits)
plot_pathmuts <- melt(msigdb_muts[c("sample","p53","cancerstatus","cancer_diagnosis",colnames(m_df_hits))])
plot_pathmuts$variable <- str_replace_all(plot_pathmuts$variable,"HALLMARK_","")

############### wilcoxon test between burden of P variants in unaff vs cancer ###############
unaffpath <- lfspath[lfspath$cancer_diagnosis == "Unaffected",]
cancerpath <- lfspath[!lfspath$cancer_diagnosis == "Unaffected",]

m_df_unaffhits <- data.frame(t(do.call(rbind,lapply(m_df$genes, function(i) lengths(regmatches(unaffpath$search,gregexpr(i,unaffpath$search)))))))
colnames(m_df_unaffhits) <- m_df$gs_name
msigdb_unaffmuts <- cbind(lfspath[lfspath$cancer_diagnosis == "Unaffected",],m_df_unaffhits)

m_df_cancerhits <- data.frame(t(do.call(rbind,lapply(m_df$genes, function(i) lengths(regmatches(cancerpath$search,gregexpr(i,cancerpath$search)))))))
colnames(m_df_cancerhits) <- m_df$gs_name
msigdb_cancermuts <- cbind(lfspath[!lfspath$cancer_diagnosis == "Unaffected",],m_df_cancerhits)
msigdb_cancermuts_wt <- msigdb_cancermuts[msigdb_cancermuts$p53 == "Wt",]
msigdb_unaffmuts_wt <- msigdb_unaffmuts[msigdb_unaffmuts$p53 == "Wt",]
msigdb_cancermuts_mut <- msigdb_cancermuts[msigdb_cancermuts$p53 == "Mut",]
msigdb_unaffmuts_mut <- msigdb_unaffmuts[msigdb_unaffmuts$p53 == "Mut",]

msigdb_mut <- rbind(msigdb_unaffmuts_mut,msigdb_cancermuts_mut)
msigdb_wt <- rbind(msigdb_unaffmuts_wt,msigdb_cancermuts_wt)

p <- list() ; p_mut <- list() ; p_wt <- list() ; p_mutwt <- list()
for (i in 1:nrow(m_df)) {
  col <- as.character(factor(m_df$gs_name)[i])
  print(col)
  p[[i]] <- wilcox.test(msigdb_cancermuts[,col],msigdb_unaffmuts[,col])$p.value
  p_mut[[i]] <- wilcox.test(msigdb_cancermuts_mut[,col],msigdb_unaffmuts_mut[,col])$p.value
  p_wt[[i]] <- wilcox.test(msigdb_cancermuts_wt[,col],msigdb_unaffmuts_wt[,col])$p.value
  p_mutwt[[i]] <- wilcox.test(msigdb_mut[,col],msigdb_wt[,col])$p.value
}

pathway_wt <- data.frame(pathway=m_df$gs_name,
                         p_cu=qvalue(do.call(rbind,p))$qvalue,
                         p_cu_mut=qvalue(do.call(rbind,p_mut))$qvalue,
                         p_cu_wt=qvalue(do.call(rbind,p_wt))$qvalue,
                         p_mutwt=qvalue(do.call(rbind,p_mutwt))$qvalue)

############### pathway analysis ###############

source('~/research/pathway_analysis/msigdb_pathway_analysis.R')

unaffgenes <- unlist(str_split(unaffpath$search,";"))
unaffgenes <- unaffgenes[unaffgenes != ""]
unaff_pathways <- fisher_test(unaffgenes)
unaff_pathways$CancerStatus <- "Unaffected"
colnames(unaff_pathways)[-1] <- paste0("unaff_",colnames(unaff_pathways)[-1])

cancergenes <- unlist(str_split(cancerpath$search,";"))
cancergenes <- cancergenes[cancergenes != ""]
cancer_pathways <- fisher_test(cancergenes)
cancer_pathways$CancerStatus <- "Cancer"
colnames(cancer_pathways)[-1] <- paste0("cancer_",colnames(cancer_pathways)[-1])

lfs_pathways <- merge(unaff_pathways,cancer_pathways, by="gs_name")

mutpath <- lfspath[lfspath$p53 == "Mut",]
wtpath <- lfspath[lfspath$p53 == "Wt",]

mutgenes <- unlist(str_split(mutpath$search,";"))
mutgenes <- mutgenes[mutgenes != ""]
mut_pathways <- fisher_test(mutgenes)
mut_pathways$CancerStatus <- "Mut"
colnames(mut_pathways)[-1] <- paste0("mut_",colnames(mut_pathways)[-1])

wtgenes <- unlist(str_split(wtpath$search,";"))
wtgenes <- wtgenes[wtgenes != ""]
wt_pathways <- fisher_test(wtgenes)
wt_pathways$CancerStatus <- "Wt"
colnames(wt_pathways)[-1] <- paste0("wt_",colnames(wt_pathways)[-1])

lfs_pathways2 <- merge(mut_pathways,wt_pathways, by="gs_name")
lfs_pathways2$WT_q <- p.adjust(lfs_pathways2$WT_p,method="fdr")


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

############### plot landscape of pathogenic variants by pathway + mutant vs wt ###############
mut_path <- data.frame(t(do.call(rbind,lapply(m_df$genes, function(i) regmatches(lfspath$search,gregexpr(i,lfspath$search))))))
colnames(mut_path) <- str_replace_all(m_df$gs_name,"HALLMARK_","")
rownames(mut_path) <- lfspath$sample
delsearch = unique(str_replace_all(unlist(str_split(c(lfspath$SV_CPG,lfspath$SV_nonCPG),';')),"\\(del\\)",""))
delsearch = paste0(delsearch[!delsearch %in% c(NA,"")],collapse = "|")
mut_path[apply(mut_path,2,function(x) grepl(delsearch,x))] <- "DEL"
dupsearch = unique(str_replace_all(unlist(str_split(c(lfspath$SV_CPG,lfspath$SV_nonCPG),';')),"\\(dup\\)|\\(cnvdup\\)| \\(CNV dup\\)",""))
dupsearch = paste0(dupsearch[!dupsearch %in% c(NA,"")],collapse = "|")
mut_path[apply(mut_path,2,function(x) grepl(dupsearch,x))] <- "DUP"
mut_path[apply(mut_path,2,function(x) grepl(";",x))] <- "SNV_INDEL"
mut_path[] = apply(mut_path, 2, function(x) as.character(x))
mut_path[mut_path == "character(0)"] = NA

mutids = lfspath$sample[lfspath$p53 == "Mut"]
wtids = lfspath$sample[lfspath$p53 == "Wt"]
mut_plot <- mut_path[rownames(mut_path) %in% mutids,]
wt_plot <- mut_path[rownames(mut_path) %in% wtids,]
relevant_paths <- unique(c(colnames(mut_plot[sapply(mut_plot, function(x) sum(!is.na(x))) > 0]),colnames(wt_plot[sapply(wt_plot, function(x) sum(!is.na(x))) > 0])))
mut_plot <- mut_plot[relevant_paths]
wt_plot <- wt_plot[relevant_paths]

lfs_pathways2 <- lfs_pathways2[lfs_pathways2$gs_name %in% paste0("HALLMARK_",colnames(mut_plot)),]

ha_mutbypath <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
      #                  Family = str_replace(lfspath$family[lfspath$sample %in% mutids],"family_",""),
                        CancerType = lfspath$cancer_diagnosis[lfspath$sample %in% mutids],
                        CancerStatus = lfspath$cancerstatus[lfspath$sample %in% mutids],
                        p53Status=lfspath$p53[lfspath$sample %in% mutids],
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
col = c(SNV_INDEL = "red", DEL = "blue",DUP="purple",qval="green")
omut <- oncoPrint(t(mut_plot),
          alter_fun = list(
            background = function(x, y, w, h) 
              grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
            SNV_INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["SNV_INDEL"], col = NA)),
            DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                gp = gpar(fill = col["DEL"], col = NA)),
            DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["DUP"], col = NA))
          ),
          remove_empty_columns = TRUE, 
          left_annotation = rowAnnotation(Pathway_qval = anno_simple(lfs_pathways2$mut_pval, col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
#                                          Wilcox_qval = anno_simple(lfs_pathways2$WT_q,col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
          col = col, 
          top_annotation = ha_mutbypath,
          row_names_gp = gpar(fontsize = 10),
          column_order = mutids)


ha_wtbypath <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
          #                        Family = str_replace(lfspath$family[lfspath$sample %in% wtids],"family_",""),
                                  CancerType = lfspath$cancer_diagnosis[lfspath$sample %in% wtids],
                                  CancerStatus = lfspath$cancerstatus[lfspath$sample %in% wtids],
                                  p53Status=lfspath$p53[lfspath$sample %in% wtids],
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
owt <- oncoPrint(t(wt_plot),
          alter_fun = list(
            background = function(x, y, w, h) 
              grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
            SNV_INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["SNV_INDEL"], col = NA)),
            DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                gp = gpar(fill = col["DEL"], col = NA)),
            DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                 gp = gpar(fill = col["DUP"], col = NA))
          ), 
          remove_empty_columns = TRUE, 
          left_annotation = rowAnnotation(Pathway_qval = anno_simple(lfs_pathways2$wt_pval, col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
          col = col, 
          top_annotation = ha_wtbypath,
          row_names_gp = gpar(fontsize = 10),
          column_order = wtids)


draw(omut+owt,annotation_legend_list=list(Legend(col_fun=colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu")),title="qvalue")))

############### plot landscape of pathogenic variants by pathway + mutant vs wt ###############

unaffids = lfspath$sample[lfspath$cancerstatus == "Unaffected"]
cancerids = lfspath$sample[lfspath$cancerstatus == "Cancer"]
cancer_plot <- mut_path[rownames(mut_path) %in% cancerids,]
unaff_plot <- mut_path[rownames(mut_path) %in% unaffids,]
relevant_paths <- unique(c(colnames(cancer_plot[sapply(cancer_plot, function(x) sum(!is.na(x))) > 0]),colnames(unaff_plot[sapply(unaff_plot, function(x) sum(!is.na(x))) > 0])))
cancer_plot <- cancer_plot[relevant_paths]
unaff_plot <- unaff_plot[relevant_paths]

lfs_pathways <- lfs_pathways[lfs_pathways$gs_name %in% paste0("HALLMARK_",colnames(cancer_plot)),]

ha_cancerbypath <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
      #                            Family = str_replace(lfspath$family[lfspath$sample %in% cancerids],"family_",""),
                                  CancerType = lfspath$cancer_diagnosis[lfspath$sample %in% cancerids],
                                  CancerStatus = lfspath$cancerstatus[lfspath$sample %in% cancerids],
                                  p53Status=lfspath$p53[lfspath$sample %in% cancerids],
                                  col = list(p53Status = c("Mut" = "#000000", "Wt" = "#E0E0E0"),
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
col = c(SNV_INDEL = "red", DEL = "blue",DUP="purple")
ocancer <- oncoPrint(t(cancer_plot),
                  alter_fun = list(
                    background = function(x, y, w, h) 
                      grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
                    SNV_INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                               gp = gpar(fill = col["SNV_INDEL"], col = NA)),
                    DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["DEL"], col = NA)),
                    DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                        gp = gpar(fill = col["DUP"], col = NA))
                  ), col = col, remove_empty_columns = TRUE, 
                  left_annotation = rowAnnotation(Pathway_qval = anno_simple(lfs_pathways$cancer_pval, col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
                                         #         Wilcox_qval = anno_simple(lfs_pathways$WT_q,col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
                  top_annotation = ha_cancerbypath,
                  row_names_gp = gpar(fontsize = 10),
                  column_order = cancerids)

ha_unaffbypath <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                #                 Family = str_replace(lfspath$family[lfspath$sample %in% unaffids],"family_",""),
                                 CancerType = lfspath$cancer_diagnosis[lfspath$sample %in% unaffids],
                                 CancerStatus = lfspath$cancerstatus[lfspath$sample %in% unaffids],
                                 p53Status=lfspath$p53[lfspath$sample %in% unaffids],
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
ounaff <- oncoPrint(t(unaff_plot),
                 alter_fun = list(
                   background = function(x, y, w, h) 
                     grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
                   SNV_INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                              gp = gpar(fill = col["SNV_INDEL"], col = NA)),
                   DEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["DEL"], col = NA)),
                   DUP = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                       gp = gpar(fill = col["DUP"], col = NA))
                 ), col = col, remove_empty_columns = TRUE, 
                 left_annotation = rowAnnotation(Pathway_qval = anno_simple(lfs_pathways$unaff_pval, col = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu"),transparency = 0.2))),
                 top_annotation = ha_unaffbypath,
                 row_names_gp = gpar(fontsize = 10),
                 column_order = unaffids)

olist = ocancer + ounaff
draw(olist,auto_adjust=FALSE,annotation_legend_list=list(Legend(col_fun=colorRamp2(c(1,0.8,0.6,0.4,0.2,0), brewer.pal(6, "BuPu")),title="qvalue")))

