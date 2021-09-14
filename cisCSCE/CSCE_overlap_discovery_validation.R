library(stringr)
library(dplyr)
library(reshape2)
library(qvalue)

outdir = "output_combatPEERbeta_azizFDRp0.1/"
discovery <- read.csv('Output/EWAS_Noob_ComBatProject_beta_450k_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
discovery$AT_ADJP <- qvalue(discovery$AT_P)$qvalue
discovery <- discovery[discovery$AT_ADJP <0.1,]
validation <- read.csv('Output/EWAS_Noob_ComBatProject_beta_850k_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
validation$AT_ADJP <- qvalue(validation$AT_P)$qvalue
validation <- validation[validation$AT_ADJP <0.1,]
extvalidation <- read.csv('Output/EWAS_Noob_ComBatProject_beta_nci_PEER_apcsgt_cancerstatus_aziztest.csv',stringsAsFactors = F)
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

#################################### Write input bed file ###################################

sig_bed <- sig_final[c("CHR","BP","BP","PROBE")] ; sig_bed$CHR <- str_replace_all(sig_bed$CHR,"chr23","chrX")
write.table(sig_bed,paste0('/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/CSCE_analysis/meQTL_EWAS/',outdir,"tmp_probes.bed"),sep='\t',row.names=F,col.names=F, quote=F)

