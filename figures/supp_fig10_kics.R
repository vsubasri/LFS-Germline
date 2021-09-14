library(ggplot2)

## Kics supplementary plots ##
kics_hered <- c('333', '345', '167', '359', '220', '295', '220', '51', '167', '98', '222', '63', '83', '120', '141', '156', '171', '232', '219')
kics_p <- read.csv('~/research/KiCS/kics.hg19.pathogenic.cpg.txt',sep='\t')
kics_p <- kics_p[!kics_p$kics_id %in% kics_hered,]
kics_clin <- read.csv('~/research/KiCS/KiCS_clinical_20200501.txt',sep='\t')
kics_ids <- read.csv('~/research/KiCS/kics_ids.txt',sep='\t')
kics_ids <- kics_ids[kics_ids$vs_pipeline == "Y",]
kics_clin <- kics_clin[kics_clin$KiCS.ID %in% kics_ids$kics_id,]
kics_clin <- kics_clin[!kics_clin$KiCS.ID %in% kics_hered,]
kics_clin <- kics_clin[!duplicated(kics_clin$KiCS.ID),]
kics_clin$P_Germline <- ifelse(kics_clin$KiCS.ID %in% kics_p$kics_id,"P",NA)

ggplot(kics_clin, aes(x=abbreviation, ..count..)) + 
    geom_bar(aes(fill = P_Germline), alpha=0.8) + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(subtitle=paste("n =",length(kics_clin$abbreviation)), 
         title= "Distribution of Cancer Diagnosis", x="Tumor Type", y = "Number of individuals") 


