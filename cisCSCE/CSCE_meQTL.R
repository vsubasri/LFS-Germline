#library(qqman)
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(aziztest))
suppressMessages(require(argparse))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))

parser <- ArgumentParser()
parser$add_argument("--value", action="store")

args <- parser$parse_args()
value <- args$value

setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')
source('Scripts/util_functions.R')
set.seed(123)

cat("[ Formatting data ]","\n")
data <- readRDS(paste0("rds/",value,".rds"))
data <- data %>%distinct(tm_donor, .keep_all = TRUE) %>%distinct(ids, .keep_all = TRUE)
data$cancer_diagnosis <- as.character(data$cancer_diagnosis)
data$cancer_diagnosis <- ifelse(grepl("Breast|DCIS|IDC|hylloides|mammary", data$cancer_diagnosis), "Breast ca", data$cancer_diagnosis)
data$cancer_diagnosis <- ifelse(grepl("ARMS|ERMS|RMS", data$cancer_diagnosis), "RMS", data$cancer_diagnosis)
data$cancer_diagnosis <- ifelse(grepl("strocytoma|pendymoma|glioma|GBM", data$cancer_diagnosis), "Glioma", data$cancer_diagnosis)
data$cancer_diagnosis <- ifelse(grepl("ACT", data$cancer_diagnosis), "ACC", data$cancer_diagnosis)

cat("[ Reading in annotation ]","\n")
## be careful about which samples are removed and what cancer they have 
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% data.frame()
anno <- anno[anno$Name %in% colnames(data),]
anno$gene <- do.call(rbind,str_split(anno$UCSC_RefGene_Name,";"))[,1]
beta <- data.frame(data[45:length(data)], row.names = data$ids)

cat("[ Performing association testing ]","\n")
spear_pval <- list() ; spear_rho <- list() ; wilc_pval <- list()
aziz_pval <- list() ; aziz_es <- list() ; aziz_oddcas <- list() ; aziz_oddratio <- list() ; aziz_direction <- list()
cancerstatus <- ifelse(data$cancerstatus == "Unaffected",0,1)
null <- calibrate_test(cancerstatus,rep=100000000)
cols <- colnames(data)[45:length(data)]
for (i in 1:length(cols)) {
  probe <- cols[i]
  es <- aziz.test(cancerstatus, data[,probe],rep=0)
  pval <- get_calibrated_pvalues(null,es$es)
  aziz_pval[[probe]] <- pval
  aziz_es[[probe]] <- es$es
  aziz_oddcas[[probe]] <- es$oddcas
  aziz_oddratio[[probe]] <- es$oddratio
  aziz_direction[[probe]] <- es$direction
  wilc_pval[[probe]] <- wilcox.test(data[,probe],cancerstatus)$p.value
  corr <- cor.test(cancerstatus, data[,probe],  method = "spearman")
  spear_pval[[probe]] <- corr$p.value
  spear_rho[[probe]] <- corr$estimate
  if (i %% 1000 == 0) {
    print(i)
  }
}

rsSNP_format <- data.frame(AT_P = do.call(rbind,aziz_pval),
                           AT_ES = do.call(rbind,aziz_es),
                           AT_oddcas = do.call(rbind,aziz_oddcas),
                           AT_oddratio = do.call(rbind,aziz_oddratio),
			   AT_direction =do.call(rbind,aziz_direction),
                           WT_P = do.call(rbind,wilc_pval),
                           S_P = do.call(rbind,spear_pval),
                           S_rho = do.call(rbind,spear_rho))
rsSNP_format$ADJP <- p.adjust(rsSNP_format$AT_P)
rsSNP_format$SNP <- rownames(rsSNP_format)
rsSNP_format$CHR <- anno$chr[match(rsSNP_format$SNP,anno$Name)]
rsSNP_format$CHR <- str_replace(rsSNP_format$CHR,"chrX","23")
rsSNP_format$CHR <- str_replace(rsSNP_format$CHR,"chrY","24")
rsSNP_format$CHR <- as.numeric(str_replace(rsSNP_format$CHR,"chr",""))
rsSNP_format$BP <- anno$pos[match(rsSNP_format$SNP,anno$Name)]
rsSNP_format$BP <- as.numeric(rsSNP_format$BP) ; rsSNP_format$ADJP <- as.numeric(rsSNP_format$ADJP)
manhattan(rsSNP_format,cex = 0.5, cex.axis = 0.8, suggestiveline=-log10(1e-2),genomewideline=-log10(1e-5))

colnames(rsSNP_format) <- c("AT_P","AT_ES","AT_oddcas","AT_oddratio","AT_direction","WT_P","S_P","S_rho","AT_ADJP","PROBE","CHR","BP")

rsSNP_format$GENE <- anno$UCSC_RefGene_Name[match(rsSNP_format$PROBE,anno$Name)]
rsSNP_format$GROUP <- anno$UCSC_RefGene_Group[match(rsSNP_format$PROBE,anno$Name)]
sig_rsSNP_format <- rsSNP_format[rsSNP_format$AT_ADJP < 0.1,]
write.csv(rsSNP_format,paste0('Output/EWAS_',value,'_cancerstatus_aziztest.csv'),quote=F,row.names = F)

