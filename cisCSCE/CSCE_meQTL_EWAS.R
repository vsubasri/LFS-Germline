require(stringr)
require(qvalue)
require(bedr)
require(reshape2)
require(dplyr)
require(tidyr)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
require(ggplot2)

setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/Scripts/LFSgermline_manuscript/')
set.seed(123)


methfile <- "/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/Noob_ComBatProject_beta_combined_PEER_apcsgt.rds"
bed <- '/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/CSCE_analysis/meQTL_EWAS/output_combatPEERbeta_azizFDRp0.1/cancerstatus_meQTL_disc_val_allgenes_931.bed'
snvfile <- '/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/CSCE_analysis/meQTL_EWAS/output_combatPEERbeta_azizFDRp0.1/probesX10kbsnvs.input'

######################### Read in methylation and clinical data ########################
cat("[ Read in data ]","\n")

data <- readRDS(methfile)
data <- data %>%distinct(tm_donor, .keep_all = TRUE) %>%distinct(ids, .keep_all = TRUE)

######################### Get list of cancer genes ########################

cpgenes <- read.csv('/hpf/largeprojects/adam/projects/lfs/resources/CPG_NEJM1508054.csv')

######################### Read in diagnostic methylation probes ########################

meQTL_probes <- read.csv(bed,stringsAsFactors = F,sep='\t',header= F)
meQTL_probes.sort1 <- bedr(engine = "bedtools", input = list(i = meQTL_probes), method = "sort", params = "");

####################### Format SNV data, remove ./. and reformat #######################

cat("[ Reformat and clean data ]","\n")

snvs <- read.csv(snvfile,sep='\t',check.names = F,stringsAsFactors = F)
colnames(snvs)[8:length(snvs)] <- do.call(rbind,str_split(colnames(snvs)[8:length(snvs)],"\\."))[,1]
snvs[snvs=="0"] = NA

############################ Impute with most common allele #############################
#maxSNP <- apply(snvs,1,function(x) names(which.max(table(x))))
#snvs <- data.frame(apply(snvs,2,function(x) ifelse(is.na(x),maxSNP,x)),check.names = F)
#########################################################################################

colnames(snvs)[8:length(snvs)] <- str_replace_all(colnames(snvs)[8:length(snvs)],"_2|A|B|_S1|","")
data <- data[data$ids %in% colnames(snvs),]
snvs <- snvs[c(colnames(snvs)[1:7],colnames(snvs)[8:length(snvs)][colnames(snvs)[8:length(snvs)] %in% data$ids])]
snvs$snv_id <- paste0("chr",snvs$CHROM, "_",snvs$POS)

snvs_coordinates <- data.frame(paste0("chr",str_replace(as.character(snvs$CHROM)," ","")),snvs$POS,snvs$POS,stringsAsFactors = F)
colnames(snvs_coordinates) <- c("chrom","start","end")
snvs_coordinates$start <- as.numeric(as.character(snvs_coordinates$start));snvs_coordinates$end <- as.numeric(as.character(snvs_coordinates$start))

snvs.sort1 <- bedr(engine = "bedtools", input = list(i = snvs_coordinates), method = "sort", params = "");

################################ Intersect SNVs and probes  ################################

cat("[ Get SNPs associated with each probe ]","\n")

snvs_meQTL <- bedr(engine="bedtools", input=list(a=meQTL_probes.sort1,b=snvs.sort1),method="window",params = "-w 10000")
snvs_meQTL <- unique(snvs_meQTL) 
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% data.frame()
snvs_meQTL$snv_id <- paste0(snvs_meQTL$V5,"_",snvs_meQTL$V6)
rm(snvs_coordinates,snvs.sort1,meQTL_probes.sort1)

overlap_snvs <- data.frame(t(snvs[8:(length(snvs))]))
colnames(overlap_snvs) <- snvs$snv_id
overlap_snvs <-overlap_snvs[1:(dim(overlap_snvs)[1]-1),]
overlap_snvs <- overlap_snvs[, which(colMeans(!is.na(overlap_snvs)) > 0.5)]

############################## Remove SNPs with no variation ##############################

overlap_snvs2 <- overlap_snvs <- overlap_snvs[apply(overlap_snvs, 2,function(x) length(unique(x)) > 1)]

########################### Match order of methylation and SNVs ###########################

overlap_snvs_plot <-overlap_snvs <- as.matrix(overlap_snvs[order(match(rownames(overlap_snvs),data$ids)),])

#################################### Calculate meQTLs #####################################

cat("[ Perform meQTL analysis ]","\n")

overlap_mm <- list() 
for (i in 1:dim(overlap_snvs)[2]) {
  probes <- snvs_meQTL[snvs_meQTL$snv_id==colnames(overlap_snvs)[i],]$V4
  if (length(na.omit(overlap_snvs[,i])) < 25) {
    next
  }
  if (i %% 1000 == 0) {
    print(i)
  }
  for (probe in probes) {
    rr <- cor.test(as.numeric(na.omit(overlap_snvs[,i])), data[which(!is.na(overlap_snvs[,i])),probe],  method = "spearman")
    id <- paste0(probe,"_",colnames(overlap_snvs)[i])
    overlap_mm[[id]] <- c(rr$p.value,rr$estimate,colnames(overlap_snvs)[i])
  }
} 

overlap_mm <- data.frame(do.call(rbind,overlap_mm),stringsAsFactors = F)
overlap_mm$probe <- do.call(rbind,str_split(rownames(overlap_mm),"_"))[,1]
overlap_mm <- overlap_mm[!is.na(overlap_mm$rho),]
overlap_mm$V1 <- as.numeric(overlap_mm$V1)
overlap_mm <- merge(overlap_mm,snvs[c(colnames(snvs)[1:7],"snv_id")],by.x="V3",by.y="snv_id")
overlap_mm_top <- overlap_mm %>% group_by(probe) %>% 
  dplyr::slice(which.min(V1))
overlap_mm_top$qval <- p.adjust(overlap_mm_top$V1, method="fdr")
overlap_mm_top <- merge(overlap_mm_top,anno,by.x="probe",by.y="Name")
overlap_mm_top$Gene <- do.call(rbind,str_split(overlap_mm_top$UCSC_RefGene_Name,';'))[,1]
overlap_mm_top$CPG <- ifelse(overlap_mm_top$Gene %in% cpgenes$Gene,"CPG",NA)
overlap_mm_top$TS <- ifelse(overlap_mm_top$Gene %in% as.character(cpgenes$Gene[grepl("Suppressor",cpgenes$Category)]),"TS",NA)
overlap_mm_top$ADR <- ifelse(overlap_mm_top$Gene %in% as.character(cpgenes$Gene[grepl("Autosomal",cpgenes$Category)]),"AD/R",NA)
write.table(overlap_mm_top,'Output/ComBatPEERbeta_931_cisCSCE.txt',sep='\t',quote=F,row.names=F)

