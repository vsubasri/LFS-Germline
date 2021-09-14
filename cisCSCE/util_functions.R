suppressMessages(require(dplyr))
suppressMessages(require(reshape2))

##########
# function that Loops through list, preprocesses, and convert to beta, m, and cn values
##########

preprocessMethod <- function(data, preprocess) {
  if (preprocess == 'raw') {
    Mset <- preprocessRaw(data)
  }
  if (preprocess == 'quan') {
    Mset <- preprocessQuantile(data, fixOutliers = TRUE,
                               removeBadSamples = TRUE, badSampleCutoff = 10.5,
                               quantileNormalize = TRUE, stratified = TRUE,
                               mergeManifest = FALSE, sex = NULL)
  }
  if (preprocess == 'illumina') {
    Mset  <- preprocessIllumina(data)
  }
  if (preprocess == 'swan') {
    Mset  <-preprocessSWAN(data)
  }
  if (preprocess == 'funnorm') {
    Mset <-preprocessFunnorm(data)
  }
  if (preprocess == 'noob') {
    Mset <-preprocessNoob(data, dyeMethod="single")
  }
  return(Mset)
}

processData <- function(directory,id_map_clin,method,valtype,arraytype) {
	rgCases <- read.metharray.exp(directory, recursive = T, force=TRUE)
	rgCases  <- convertArray(rgCases,outType = c("IlluminaHumanMethylation450k"))
	preprocessedCases <- preprocessMethod(rgCases,method)
	pheno <- id_map_clin[match(rownames(preprocessedCases@colData),id_map_clin$SentrixID),]
	pheno$cancerstatus <- ifelse(pheno$cancer_diagnosis != "Unaffected","Cancer","Unaffected")
	pheno$array <- arraytype
	GRset <- mapToGenome(preprocessedCases)
	estSex <- getSex(GRset)
	if (valtype == "beta") {
		meth <- data.frame(getBeta(GRset))
	} else if (valtype == "M") {
		meth <- data.frame(getM(GRset))
	} else {
		cat("Invalid methylation type value")
	}
	pheno$gender <- estSex[, "predictedSex"]
	data <- cbind(pheno,t(meth))
	return(data)
}

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

##########
# Function that gets malkin id from sample_name
##########
getIdName <- function(data) {
  column_split <- strsplit(as.character(data$ids), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  sub_ids <- gsub('-', '', sub_ids)
  data$ids <- sub_ids
  data$identifier <- NULL
  return(data)
}


getTrainTestVal <- function(data,nseed,nsamps) {
	set.seed(nseed)
	ind <- 44
	## create test set
	data_val <- data[data$ids %in% seq(6056,6147,1),]
	## create validation set
	data <- data[!data$ids %in% seq(6056,6147,1),]
	testids <- data[1:ind] %>%
        	filter(!ids %in% seq(6056,6147,1)) %>%
        	group_by(ageoutcome = ifelse(data$ageofonset > 12*6 | is.na(data$ageofonset),"GT6","LT6")) %>%
        	sample_n(nsamps) %>%
        	select(tm_donor)
	data_test <- data[data$tm_donor %in% testids$tm_donor,]
	## create training set
	data_train <- data[!data$tm_donor %in% c(data_test$tm_donor,data_val$tm_donor),]

	## remove outliers in training data
	p_train <- prcomp(data_train[(ind+1):length(data_train)],scale=T)
	p_train <- cbind(data_train[1:ind],p_train$x)
	keep_train <- remove_outliers(p_train,3)
	data_train <- data_train[data_train$SentrixID %in% keep_train,]
	return(list(data_train,data_test,data_val))
}

generate_pcsummary <- function(pc_clin, filename, outdir) {
  pc_clin_tt <- pc_clin[pc_clin$cancer_diagnosis %in% c("ACC","RMS","CPC","Breast","OS","Glioma","Unaffected"),]
  pcs <- colnames(pc_clin)[grepl( "PC" , names(pc_clin) )] ;
  anova_array <- list() ; anova_batch <- list() ; anova_cancerstatus <- list() ; anova_cancertype <- list()
  for (i in pcs) {
    anova_array[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$array))$p.value
    anova_batch[[i]] <- summary(aov(pc_clin[,i] ~ as.factor(pc_clin$Project)))[[1]][1,"Pr(>F)"]
    anova_cancerstatus[[i]] <- wilcox.test(pc_clin[,i] ~ as.factor(pc_clin$cancerstatus))$p.value
    anova_cancertype[[i]] <- summary(aov(pc_clin_tt[,i] ~ as.factor(pc_clin_tt$cancer_diagnosis)))[[1]][1,"Pr(>F)"]
  }

  anova <- data.frame(PC=pcs,
                      batch_p=as.numeric(do.call("cbind",anova_batch)),
                      array_p=as.numeric(do.call("cbind",anova_array)),
                      cancerstatus_p=as.numeric(do.call("cbind",anova_cancerstatus)),
		      cancertype_p=as.numeric(do.call("cbind",anova_cancertype)))
  anova <- cbind(anova,data.frame(t(summary(pc)$importance)))
  anova$batch_padj <- anova$batch_p/anova$Proportion.of.Variance
  anova$array_padj <- anova$array_p/anova$Proportion.of.Variance
  anova$cancerstatus_padj <- anova$cancerstatus_p/anova$Proportion.of.Variance
  anova$cancertype_padj <- anova$cancertype_p/anova$Proportion.of.Variance
  write.csv(anova,paste0(outdir,'Output/',filename),quote=F,row.names=F)
}

generate_pcplots <- function(pc_clin,output,outdir) {
  pdf(paste0(outdir,"Plots/",output,"_PCA_cancerstatus.pdf"),width=9,height=7)
  cancerstatusplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=cancerstatus)) +
    geom_point() +
    theme_bw()
  print(cancerstatusplot)
  suppressMessages(dev.off())

  pdf(paste0(outdir,"Plots/",output,"_PCA_gender.pdf"),width=9,height=7)
  genderplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=gender)) +
    geom_point() +
    theme_bw()
  print(genderplot)
  suppressMessages(dev.off())

  if (length(table(pc_clin$array)) > 1) {
  pdf(paste0(outdir,"Plots/",output,"_PCA_array.pdf"),width=9,height=7)
    arrayplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=array)) +
      geom_point() +
      theme_bw()
    print(arrayplot)
    suppressMessages(dev.off())
  }

  pdf(paste0(outdir,"Plots/",output,"_PCA_age.pdf"),width=9,height=7)
  ageplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=agesamplecollection)) +
    geom_point() +
    theme_bw()
  print(ageplot)
  suppressMessages(dev.off())

  pc_clin$array <- factor(pc_clin$array)
  pdf(paste0(outdir,"Plots/",output,"_PCA_confounders.pdf"),width=9,height=7)
  sourceplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=Project, shape=array)) +
      geom_point() +
      theme_bw()
  print(sourceplot)
  suppressMessages(dev.off())

  pdf(paste0(outdir,"Plots/",output,"_PCA_project.pdf"),width=9,height=7)
  projectplot <- ggplot(pc_clin,aes(x=PC1,y=PC2,color=Project)) +
    geom_point() +
    theme_bw()
  print(projectplot)
  suppressMessages(dev.off())

}

remove_SNPs_xRP <- function(GRset,pheno,output) {
  GRset_noSNPs <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)
  B_noSNPs <- data.frame(getBeta(GRset_noSNPs))
  B_noSNPs <- B_noSNPs[rowSums(is.na(B_noSNPs)) == 0,]
  data_B_noSNPs <- cbind(pheno,t(B_noSNPs))
  saveRDS(data_B_noSNPs,paste0("rds/",output,"_noSNPs.rds"))

  xReactiveProbes <- read.csv(file="Resources/cross_reactive_probes.csv", stringsAsFactors=FALSE)
  keep <- !(featureNames(GRset_noSNPs) %in% xReactiveProbes$TargetID)
  GRset_noSNPs_noxRP <- GRset_noSNPs[keep,]
  B_noSNPs_noxRP <- data.frame(getBeta(GRset_noSNPs_noxRP))
  B_noSNPs_noxRP <- B_noSNPs_noxRP[rowSums(is.na(B_noSNPs_noxRP)) == 0,]
  data_B_noSNPs_noxRP <- cbind(pheno,t(B_noSNPs_noxRP))
  saveRDS(data_B_noSNPs_noxRP,paste0("rds/",output,"_noSNPs_noxRP.rds"))
}

remove_sex <- function(cleaned,pheno,output) {
  sexprobes <- scan('Scripts/sexprobes.txt',what="character")
  nosex <- cleaned[, !colnames(cleaned) %in% sexprobes]
  data_nosex <- cbind(pheno,nosex)
  saveRDS(data_nosex,paste0("rds/",output,"_NoSex.rds"))
  return(nosex)
}

get_technicalreplicates <- function(data) {
  duplicated_ids_diffarray = data[c("tm_donor","array")] %>%
    filter(!is.na(tm_donor)) %>%
    group_by(tm_donor) %>%
    mutate(types = n_distinct(array)) %>%
    filter(types > 1) %>%
    distinct(tm_donor, .keep_all = TRUE) %>%
    pull(tm_donor)
  duplicated_data <- data[data$tm_donor %in% duplicated_ids_diffarray,] %>%
    distinct(tm_donor, array, .keep_all = TRUE)
  return(duplicated_data)
}

plot_concordance <- function(duplicated_data,value,outdir) {
  allprobes = colnames(duplicated_data)[grepl("cg",colnames(duplicated_data))]
  sample_probes = sample(allprobes,100000,replace=FALSE)
  sample_beta = duplicated_data[c("tm_donor","array",sample_probes)]
  plot_dups <- melt(sample_beta,measure.vars=sample_probes, id.vars= c("tm_donor","array")) %>%
  dcast(tm_donor + variable ~ array)
  tech_corr <- round(cor(plot_dups$`450`,plot_dups$`850`,method="spearman"),digits=3)
  plot_dups$diff <- abs(plot_dups$`450`-plot_dups$`850`)
  distsum <- sum(plot_dups[, 'diff'])

  pdf(paste0(outdir,'Plots/',value,'_450_850_concordance.pdf'),width=9,height=7)
  cplot <- ggplot(plot_dups,aes(x=`450`,y=`850`)) +
    geom_point(aes(colour = diff)) +
    theme_bw() +
    ggtitle(paste0("Correlation = ",tech_corr,", Distance = ",distsum)) +
    scale_colour_gradient2(low="yellow",mid="black",high="blue") +
    labs(x="450k Array",y="850k Array",color=value)
  print(cplot)
  suppressMessages(dev.off())
}
