library(stringr)
library(ggplot2)
library(lme4)
library(reshape2)
library(ggpubr)
library(parallel)
library(qvalue)
library(stringr)
library(dplyr)

mutclin <- read.csv("~/research/lfs_wt/PRS/lfsmut.fam", header=TRUE,stringsAsFactors=FALSE,sep='\t')
wtclin <- read.csv("~/research/lfs_wt/PRS/lfswt.fam", header=TRUE,stringsAsFactors=FALSE,sep='\t')
clin <- rbind(mutclin,wtclin)

## Filter rsSNPs from gwas catalogue
rsSNPs <- read.csv('~/research/EWAS_meQTL/rsSNPs_literature/gwas_catalog_v1.0-associations_e96_r2019-07-30.tsv',sep='\t',stringsAsFactors = F)
cancer_rsSNPs <- rsSNPs[grepl('osteosarcoma|Osteosarcoma|cancer|Cancer|glioma|Glioma|leukemia|carcinoma|malignancies|lymphoma|Lymphoma|Glioblastoma|Neuroblastoma|tumor|sarcoma|Sarcoma',rsSNPs$DISEASE.TRAIT),]
notvalid <-c("30738427","22318345","30738427","24025145","22293537","29432556","19176441","28090653","23478653",
             "23233662","26237429","22180457","30952644","27515689","30161160","22637743","30277654","28120872",
             "27670397","22142827","20350937","20592726","22664479","27723779","28779025","28429243","29266176",
             "23242368","26460308","26546620","25964295","25987655","25963547","19901119","25987655","21079520",
             "30648747","22872573","24974847","25710658","30299488","29278617","30410027","30305637","28817678",
             "25867717","22076464","29743610","25192705","27145994","24743840","25826619","29942513","23544012",
             "30152087","23817570","23354978","22232737","20729852","20686608","23103227","24836286","29551738",
             "23263487","25038754","27354352","30529582","30898391","30323354","26443449","24448986","29356057",
             "22383897","21242260","26732429","21866343","22158540","20676098","20676096","20871597","28274756",
             "20876614","27501781","22797724","22923026","24143190","19698717","29385134","29471430","18326623",
             "28163062","25105248","30557369","24700089","22951594","21499248","25281661","30281874","21908515",
             "19219042","19478819","19664746","20101243","21118971","21642993","21725308","22807686","28662289",
             "22960999","23023329","23341777","23644492","23749188","24658283","25129146","25134534","30412241",
             "25145502","26701879","27206850","27436580","29246937","30714141","21540461","22452962","23143601",
             "23144319","24737549","28703219","22923054","23468962","27393504","28171663","28295283","27207650",
             "26965516","26084801","28924153","29884837","26129866","26222057","26586795","25939597","29210060",
             "29059430","29698419","25824743","28081215")
cancer_rsSNPs <- cancer_rsSNPs[!cancer_rsSNPs$PUBMEDID %in% notvalid,]
cancer_rsSNPs <- cancer_rsSNPs[!is.na(cancer_rsSNPs$OR.or.BETA),]

lfs_wt_cancerSNPs <- read.csv('~/research/lfs_wt/PRS/lfswt_cancer_rsSNP.txt',sep='\t',check.names = F,stringsAsFactors = F)
lfs_mut_cancerSNPs <- read.csv('~/research/lfs_wt/PRS/lfsmut_cancer_rsSNP.txt',sep='\t',check.names = F,stringsAsFactors = F)

process_SNPs <- function(lfs_cancerSNPs,cancer_rsSNPs) {
    lfs_cancerSNPs$OR.or.BETA <- cancer_rsSNPs$OR.or.BETA[match(lfs_cancerSNPs$ID,cancer_rsSNPs$SNPS)]
    lfs_cancerSNPs$P.VALUE <- cancer_rsSNPs$P.VALUE[match(lfs_cancerSNPs$ID,cancer_rsSNPs$SNPS)]
    lfs_cancerSNPs$MAF <- cancer_rsSNPs$RISK.ALLELE.FREQUENCY[match(lfs_cancerSNPs$ID,cancer_rsSNPs$SNPS)]
    lfs_cancerSNPs$A1<- do.call(rbind,str_split( cancer_rsSNPs$STRONGEST.SNP.RISK.ALLELE[match(lfs_cancerSNPs$ID,cancer_rsSNPs$SNPS)],'-'))[,2]
    lfs_cancerSNPs <- lfs_cancerSNPs[lfs_cancerSNPs$A1 != "?",]
    lfs_cancerSNPs$A2 <- ifelse(lfs_cancerSNPs$A1 == lfs_cancerSNPs$REF,lfs_cancerSNPs$ALT,lfs_cancerSNPs$REF)
    lfs_cancerSNPs$Effect_Allele <- ifelse(lfs_cancerSNPs$A1 == lfs_cancerSNPs$REF,"REF","ALT")
    lfs_cancerSNPs <- lfs_cancerSNPs[lfs_cancerSNPs$ID %in% cancer_rsSNPs$SNPS,]
    return(lfs_cancerSNPs)
}

get_genotype <- function(lfs_cancerSNPs) {
  lastcol <- length(lfs_cancerSNPs)-7
  lfs_cancerSNPs_GT <- data.frame(apply(t(lfs_cancerSNPs[8:lastcol]),2,function(x) str_replace(x,"/"," ")))
  colnames(lfs_cancerSNPs_GT) <- lfs_cancerSNPs$ID; rownames(lfs_cancerSNPs_GT) <- colnames(lfs_cancerSNPs[8:lastcol])
  return(lfs_cancerSNPs_GT)
}

lfs_wt_cancerSNPs <- process_SNPs(lfs_wt_cancerSNPs,cancer_rsSNPs)
lfs_wt_cancerSNPs_GT <- get_genotype(lfs_wt_cancerSNPs)
lfs_mut_cancerSNPs <- process_SNPs(lfs_mut_cancerSNPs,cancer_rsSNPs)
lfs_mut_cancerSNPs_GT <- get_genotype(lfs_mut_cancerSNPs)

lfs_wt_cancerSNPs_GT <- lfs_wt_cancerSNPs_GT[colnames(lfs_wt_cancerSNPs_GT) %in% colnames(lfs_mut_cancerSNPs_GT)] 
lfs_mut_cancerSNPs_GT <- lfs_mut_cancerSNPs_GT[colnames(lfs_wt_cancerSNPs_GT)] 

lfs_cancerSNPs_GT <- rbind(lfs_wt_cancerSNPs_GT,lfs_mut_cancerSNPs_GT)
lfs_cancerSNPs_GT <- lfs_cancerSNPs_GT[, sapply(lfs_cancerSNPs_GT, function(col) length(unique(col))) < 5]


clin <- clin[clin$Individual_ID %in% rownames(lfs_cancerSNPs_GT),]
lfs_cancerSNPs_GT <- lfs_cancerSNPs_GT[match(clin$Individual_ID,rownames(lfs_cancerSNPs_GT)),]
## Filter rsSNPs from gwas catalogue
SNP.pvalue <- cancer_rsSNPs[c("SNPS","P.VALUE")]
SNP.pvalue <- SNP.pvalue[!duplicated(SNP.pvalue$SNPS),]
SNP.pvalue <- SNP.pvalue[match(colnames(lfs_cancerSNPs_GT),SNP.pvalue$SNPS),]
map_plink <- unique(rbind(lfs_wt_cancerSNPs[c("CHROM","ID","POS")],lfs_mut_cancerSNPs[c("CHROM","ID","POS")]))
map_plink$d <- 0 ; map_plink <- map_plink[c("CHROM","ID","d","POS")]
ped_plink <- cbind(clin,lfs_cancerSNPs_GT)
ped_plink <- Filter(function(x) mean(x == '. .') <= 0, ped_plink)
map_plink <- map_plink[map_plink$ID %in% colnames(ped_plink),]
map_plink <- map_plink[match(colnames(ped_plink),map_plink$ID),]

cancer_rsSNPs <- cancer_rsSNPs[cancer_rsSNPs$SNPS %in% colnames(ped_plink),]
os_rsSNPs <- cancer_rsSNPs[grepl('osteosarcoma|Osteosarcoma',cancer_rsSNPs$DISEASE.TRAIT),]
glioma_rsSNPs <- cancer_rsSNPs[grepl('glioma|Glioma',cancer_rsSNPs$DISEASE.TRAIT),]
bc_ids= c("25751625","29059683")
bc_rsSNPs <- cancer_rsSNPs[cancer_rsSNPs$PUBMEDID %in% bc_ids,]
nonmelanoma_skin_rsSNPs <- cancer_rsSNPs[grepl('Non-melanoma',cancer_rsSNPs$DISEASE.TRAIT),]
ALL_rsSNPs <- cancer_rsSNPs[grepl('Acute lymphoblastic leukemia',cancer_rsSNPs$DISEASE.TRAIT),]
leukemia_rsSNPs <- cancer_rsSNPs[grepl('leukemia|Leukemia',cancer_rsSNPs$DISEASE.TRAIT),]
relevant_cancer_rsSNPs <- unique(as.character(c(os_rsSNPs$SNPS,glioma_rsSNPs$SNPS,bc_rsSNPs$SNPS,
                                                nonmelanoma_skin_rsSNPs$SNPS,ALL_rsSNPs$SNPS,leukemia_rsSNPs$SNPS)))

effect_cancerSNPs <- cancer_rsSNPs[c("SNPS","OR.or.BETA","P.VALUE","RISK.ALLELE.FREQUENCY","X95..CI..TEXT.","PUBMEDID")]
N <-apply(do.call(rbind,str_split(str_replace_all(gsub("[^0-9., ]", "", gsub( "[-+*/]", " , ", cancer_rsSNPs$INITIAL.SAMPLE.SIZE) ),",","")," ")),2,as.numeric)
N[is.na(N)] <- 0
effect_cancerSNPs$N <- rowSums(N)
effect_cancerSNPs <- effect_cancerSNPs[!duplicated(effect_cancerSNPs$SNPS),]
effect_cancerSNPs$A1 <- lfs_wt_cancerSNPs$A1[match(effect_cancerSNPs$SNPS,lfs_wt_cancerSNPs$ID)]
colnames(effect_cancerSNPs) <- c("ID","OR","PVALUE","MAF","95CI","PUBMEDID","N","A1")
effect_cancerSNPs <- effect_cancerSNPs[c("ID","A1","OR","PVALUE","MAF","95CI","N","PUBMEDID")]
  
#write.table(ped_plink,'~/research/PRS/lfs.ped',sep=' ',quote=F,row.names = F,col.names = F)
#write.table(map_plink,'~/research/PRS/lfs.map',sep='\t',quote=F,row.names = F,col.names = F)
#write.table(map_plink$ID,'~/research/PRS/lfs.valid.snp',sep='\n',quote=F,row.names = F,col.names = F)
#write.table(effect_cancerSNPs,'~/research/PRS/gwas_catalogue_cancerSNPs',sep='\t',quote=F,row.names = F)

lfs_wt_cancerSNPs_input <- read.csv('~/research/lfs_wt/PRS/LFSwt_cancer_rsSNP.input.txt',sep='\t',check.names = F,stringsAsFactors = F)
lfs_mut_cancerSNPs_input <- read.csv('~/research/lfs_wt/PRS/LFSmut_cancer_rsSNP.input.txt',sep='\t',check.names = F,stringsAsFactors = F)
lfs_wt_cancerSNPs_input <- lfs_wt_cancerSNPs_input[lfs_wt_cancerSNPs_input$ID %in% colnames(lfs_cancerSNPs_GT),]
lfs_mut_cancerSNPs_input <- lfs_mut_cancerSNPs_input[lfs_mut_cancerSNPs_input$ID %in% colnames(lfs_cancerSNPs_GT),]
lfs_cancerSNPs_input <- cbind(lfs_mut_cancerSNPs_input,lfs_wt_cancerSNPs_input[8:length(lfs_wt_cancerSNPs_input)])
cancerstatus <- clin$Phenotype[match(colnames(lfs_cancerSNPs_input[8:length(lfs_cancerSNPs_input)]),clin$Individual_ID)]

p <- list() ; rho <- list()
for (i in 1:dim(lfs_cancerSNPs_input)[1]) {
  corr <- cor.test(cancerstatus, as.numeric(lfs_cancerSNPs_input[i,8:length(lfs_cancerSNPs_input)]),  method = "spearman")
  p[[i]] <- corr$p.value ; rho[[i]] <- corr$estimate
}

results <- data.frame(rsSNP=lfs_cancerSNPs_input$ID
                      ,p=do.call(rbind,p),rho=do.call(rbind,rho))
results$q <- qvalue(results$p)$qvalue
results$trait <- cancer_rsSNPs$DISEASE.TRAIT[match(results$rsSNP,cancer_rsSNPs$SNPS)]
results$study_p <- cancer_rsSNPs$P.VALUE[match(results$rsSNP,cancer_rsSNPs$SNPS)]
results$OR <- cancer_rsSNPs$OR.or.BETA[match(results$rsSNP,cancer_rsSNPs$SNPS)]


rsSNP_plot <- data.frame(t(lfs_cancerSNPs_input[8:length(lfs_cancerSNPs_input)]))
colnames(rsSNP_plot) <- lfs_cancerSNPs_input$ID
rsSNP_plot <- merge(clin,rsSNP_plot,by.x="Individual_ID",by.y=0)

ggplot(rsSNP_plot,aes(x=as.factor(Phenotype),y=rs1004030)) +
  geom_boxplot() + geom_point(aes(color=rs1004030),position = position_jitterdodge())


######### results ##############

mutclin <- read.csv("~/research/resources/lfs_mut_normals_wgs_clinical.txt",sep='\t',row.names=1)[c('cancer_diagnosis')]
wtclin <- read.csv("~/research/lfs_wt/lfswt_clin.csv",row.names = 2)[c('cancer_diagnosis')]
colnames(wtclin) <- c("cancer_diagnosis")
ttclin <- rbind(mutclin,wtclin) 

read_scores <- function(dir,tt,clin) {
  files <- list.files(dir,paste0(tt,".out?.*profile?"),full.names = T)
  all_scores <- list()
  for (file in files) {
    pvalthres <- str_replace(str_split(file,"out.")[[1]][2],".profile","")
    score <- as.numeric(paste0("0.",read.csv(file,sep='.')$FID......IID..PHENO....CNT...CNT2....SCORE))
    all_scores[[pvalthres]] <- score
  }
  all_scores <- data.frame(t(do.call(rbind,all_scores)))
#  all_scores <- all_scores[ , colSums(is.na(all_scores)) == 0]
  clin_scores <- cbind(clin,tt=as.character(ttclin[match(clin$Individual_ID,rownames(ttclin)),]),all_scores)
  ## remove cog samples
  clin_scores <- clin_scores[!clin_scores$Individual_ID %in% cog_ids,]
  cols = c("Individual_ID","tt",colnames(clin_scores)[8:length(clin_scores)])
  clin_plot <- melt(clin_scores[cols],idvar="Individual_ID")
  if (tt == "bc") {
    clin_plot$tt <- str_replace_all(clin_plot$tt,"Breast Cancer|IDC|BC|Breast","bc")
    clin_plot$cancerstatus <- ifelse(clin_plot$tt == tt,"Breast Cancer","Not Breast Cancer")
    clin_scores$tt_ofinterest <- factor(ifelse(grepl("Breast Cancer|IDC|BC|Breast",clin_scores$tt),1,0))
  } else if (tt == "os") {
    clin_plot$tt <- str_replace_all(clin_plot$tt,"OST|OS","os")
    clin_plot$cancerstatus <- ifelse(clin_plot$tt == tt,"OS","Not OS")
    clin_scores$tt_ofinterest <- factor(ifelse(grepl("OST|OS",clin_scores$tt),1,0))
  } else if (tt == "leukemia") {
    clin_plot$tt <- str_replace_all(clin_plot$tt,"AML|ALL","leukemia")
    clin_plot$cancerstatus <- ifelse(clin_plot$tt == tt,"Leukemia","Not Leukemia")
    clin_scores$tt_ofinterest <- factor(ifelse(grepl("AML|ALL",clin_scores$tt),1,0))
  } else if (tt == "glioma") {
      clin_plot$tt <- str_replace_all(clin_plot$tt,"LGG|DA|PNET|Glioma","glioma")
      clin_plot$cancerstatus <- ifelse(clin_plot$tt == tt,"Glioma","Not Glioma")
      clin_scores$tt_ofinterest <- factor(ifelse(grepl("LGG|DA|PNET|Glioma",clin_scores$tt),1,0))
  } else if (tt == "nonmelanoma_skin") {
    clin_plot$tt <- str_replace_all(clin_plot$tt,"Basal Cell Carcinoma","nonmelanoma_skin")
    clin_plot$cancerstatus <- ifelse(clin_plot$tt == tt,"Skin","Not Skin")
    clin_scores$tt_ofinterest <- factor(ifelse(grepl("Basal Cell Carcinoma",clin_scores$tt),1,0))
  } else if (tt == "lfs") {
    clin_plot$cancerstatus <- ifelse(clin_plot$tt == "Unaffected", "Unaffected","Cancer")
    clin_scores$tt_ofinterest <- factor(ifelse(grepl("Unaffected",clin_scores$tt),1,0))
  } else {
    print("Invalid cancer type")
    exit()
  }
  anno_df = compare_means(value ~ variable, group.by = "cancerstatus", data = clin_plot) %>%
    mutate(y_pos = 40)
  p <- ggplot(clin_plot,aes(x=cancerstatus,y=value)) +
    geom_boxplot(position=position_dodge()) + 
    geom_point(aes(color=cancerstatus),position = position_jitterdodge(),na.rm = TRUE) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          text = element_text(size=14,family="Helvetica Neue")) +   
    facet_wrap(~variable,nrow = 1) + 
    stat_compare_means(method="wilcox.test",method.args = list(alternative="greater")) + 
    scale_color_manual(values = c("#DE3533","#0047AB")) +
    guides(fill=guide_legend(title="Cancer Type"))
  return(list(clin_scores,p))
}

#https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/

########### account for family differences #############

lfs.clin <- read.csv('~/research/EWAS_meQTL/lfs_mut_clinical_comprehensive.csv',stringsAsFactors = F)
cog_ids <- unique(na.omit(lfs.clin$tm_donor[lfs.clin$Hospital == "COG"]))

mixed_model_covars <- function(mm_scores) {
  mm_scores <- mm_scores[complete.cases(mm_scores),]
# estimate the model and store results in m
  m <- glmer(tt_ofinterest ~ PRS + MutStatus + Sex + (1 | Family_ID), data = mm_scores, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
  fstat <- anova(m)[1,4]
  df2 <- dim(mm_scores)[1]-1
  pval <- pf(q=fstat, df1=1, df2=df2, lower.tail=FALSE)
  print(isSingular(m))
  
  
  # print the mod results without correlations among fixed effects
  print(m, corr = FALSE)
    
    
  ## 95% CI
  se <- sqrt(diag(vcov(m)))
  # table of estimates with 95% CI
  (tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 *
                  se))
  ## odds ratio
  exp(tab)
  
  sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
    cid <- unique(dat[, clustervar[1]])
    ncid <- length(cid)
    recid <- sample(cid, size = ncid * reps, replace = TRUE)
    if (replace) {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
                                        size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
      })
    } else {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
      })
    }
    dat <- as.data.frame(do.call(rbind, rid))
    dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
                                labels = FALSE))
    dat$NewID <- factor(dat$NewID)
    return(dat)
  }
  
  set.seed(20)
  tmp <- sampler(mm_scores, "Family_ID", reps = 100)
  
  bigdata <- cbind(tmp, mm_scores[tmp$RowID, ])
  
  f <- fixef(m)
  r <- getME(m, "theta")
  
  
  cl <- makeCluster(4)
  clusterExport(cl, c("bigdata", "f", "r"),envir=environment())
  clusterEvalQ(cl, require(lme4))
  
  myboot <- function(i) {
    object <- try(glmer( tt_ofinterest ~ PRS + MutStatus+ Sex + (1 | NewID), data = bigdata, subset = Replicate == i, family = binomial,
                        nAGQ = 1, start = list(fixef = f, theta = r)), silent = TRUE)
    if (class(object) == "try-error")
      return(object)
    c(fixef(object), getME(object, "theta"))
  }
  
  start <- proc.time()
  res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
  end <- proc.time()
  
  # shut down the cluster
  stopCluster(cl)
  
  # calculate proportion of models that successfully converged
  success <- sapply(res, is.numeric)
  mean(success)
  
  # combine successful results
  bigres <- do.call(cbind, res[success])
  
  # calculate 2.5th and 97.5th percentiles for 95% CI
  (ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))
  
  # All results
  finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres),
                      ci)
  # round and print
  print(round(finaltable, 3))
  
  # temporary data
  tmpdat <- mm_scores[, c("PRS","MutStatus","Sex","Family_ID")]
  
  summary(mm_scores$PRS)
  
  jvalues <- with(mm_scores, seq(from = min(PRS), to = max(PRS), length.out = 100))
  
  # calculate predicted probabilities and store in a list
  pp <- lapply(jvalues, function(j) {
    tmpdat$PRS <- j
    predict(m, newdata = tmpdat, type = "response")
  })
  
  # average marginal predicted probability across a few different PRS
  sapply(pp[c(1, 20, 40, 60, 80, 100)], mean)
  
  
  # get the means with lower and upper quartiles
  plotdat <- t(sapply(pp, function(x) {
    c(M = mean(x), quantile(x, c(0.25, 0.75)))
  }))
  
  # add in PRS values and convert to data frame
  plotdat <- as.data.frame(cbind(plotdat, jvalues))
  
  # better names and show the first few rows
  colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "PRS")
  head(plotdat)
  
  # plot average marginal predicted probabilities
  ggplot(plotdat, aes(x = PRS, y = PredictedProbability)) + geom_line() + geom_linerange(aes(ymin = Lower, ymax = Upper)) + geom_line(size = 2) + ylim(c(0, 1))
  
  # calculate predicted probabilities and store in a list
  sexprobs <- lapply(levels(mm_scores$Sex), function(sex) {
    tmpdat$Sex[] <- sex
    lapply(jvalues, function(j) {
      tmpdat$PRS <- j
      predict(m, newdata = tmpdat, type = "response")
    })
  })
  
  # get means and quartiles for all jvalues for each Sex
  plotsex <- lapply(sexprobs, function(X) {
    temp <- t(sapply(X, function(x) {
      c(M=mean(x), quantile(x, c(.25, .75)))
    }))
    temp <- as.data.frame(cbind(temp, jvalues))
    colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "PRS")
    return(temp)
  })
  
  # collapse to one data frame
  plotsex <- do.call(rbind, plotsex)
  
  # add cancer stage
  plotsex$Sex <- factor(rep(levels(mm_scores$Sex), each = length(jvalues)))
  
  p <- ggplot(plotsex, aes(x = PRS, y = PredictedProbability)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Sex), alpha = .15) +
    geom_line(aes(colour = Sex), size = 2) +
    ylim(c(0, 1)) + facet_wrap(~  Sex) + 
    scale_fill_manual(values = c("#DE3533","#0047AB")) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          text = element_text(size=14,family="Helvetica Neue")) 
  
  biprobs <- lapply(levels(mm_scores$MutStatus), function(stage) {
    tmpdat$MutStatus[] <- stage
    lapply(jvalues, function(j) {
      tmpdat$PRS <- j
      predict(m, newdata = tmpdat, type = "response")
    })
  })
  
  # get means and quartiles for all jvalues for each tumour type
  plotdat2 <- lapply(biprobs, function(X) {
    temp <- t(sapply(X, function(x) {
      c(M=mean(x), quantile(x, c(.25, .75)))
    }))
    temp <- as.data.frame(cbind(temp, jvalues))
    colnames(temp) <- c("PredictedProbability", "Lower", "Upper", "PRS")
    return(temp)
  })
  
  
  # collapse to one data frame
  plotdat2 <- do.call(rbind, plotdat2)
  
  # add mutation status
  plotdat2$MutStatus <- factor(rep(levels(mm_scores$MutStatus), each = length(jvalues)))

  p2 <- ggplot(plotdat2, aes(x = PRS, y = PredictedProbability)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = MutStatus), alpha = .15) +
    geom_line(aes(colour = MutStatus), size = 2) +
    ylim(c(0, 1)) + facet_wrap(~  MutStatus) + 
    theme_classic() + 
    scale_fill_manual(values = c("#DE3533","#0047AB")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          text = element_text(size=14,family="Helvetica Neue")) 
  return(list(p,p2))
}

mixed_model_nocovars <- function(mm_scores) {
  mm_scores <- mm_scores[complete.cases(mm_scores),]
  # estimate the model and store results in m
  m <- glmer(tt_ofinterest ~ PRS + (1 | Family_ID), data = mm_scores, family = binomial, control = glmerControl(optimizer = "bobyqa"),
             nAGQ = 10)
  fstat <- anova(m)[1,4]
  df2 <- dim(mm_scores)[1]-1
  pval <- pf(q=fstat, df1=1, df2=df2, lower.tail=FALSE)
  print(isSingular(m))
  
  
  # print the mod results without correlations among fixed effects
  print(m, corr = FALSE)
  
  
  ## 95% CI
  se <- sqrt(diag(vcov(m)))
  # table of estimates with 95% CI
  (tab <- cbind(Est = fixef(m), LL = fixef(m) - 1.96 * se, UL = fixef(m) + 1.96 *
                  se))
  ## odds ratio
  exp(tab)
  
  sampler <- function(dat, clustervar, replace = TRUE, reps = 1) {
    cid <- unique(dat[, clustervar[1]])
    ncid <- length(cid)
    recid <- sample(cid, size = ncid * reps, replace = TRUE)
    if (replace) {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = sample(which(dat[, clustervar] == recid[i]),
                                        size = length(which(dat[, clustervar] == recid[i])), replace = TRUE))
      })
    } else {
      rid <- lapply(seq_along(recid), function(i) {
        cbind(NewID = i, RowID = which(dat[, clustervar] == recid[i]))
      })
    }
    dat <- as.data.frame(do.call(rbind, rid))
    dat$Replicate <- factor(cut(dat$NewID, breaks = c(1, ncid * 1:reps), include.lowest = TRUE,
                                labels = FALSE))
    dat$NewID <- factor(dat$NewID)
    return(dat)
  }
  
  set.seed(20)
  tmp <- sampler(mm_scores, "Family_ID", reps = 100)
  
  bigdata <- cbind(tmp, mm_scores[tmp$RowID, ])
  
  f <- fixef(m)
  r <- getME(m, "theta")
  
  
  cl <- makeCluster(4)
  clusterExport(cl, c("bigdata", "f", "r"))
  clusterEvalQ(cl, require(lme4))
  
  myboot <- function(i) {
    object <- try(glmer( tt_ofinterest ~ PRS +  (1 | NewID), data = bigdata, subset = Replicate == i, family = binomial,
                         nAGQ = 1, start = list(fixef = f, theta = r)), silent = TRUE)
    if (class(object) == "try-error")
      return(object)
    c(fixef(object), getME(object, "theta"))
  }
  
  start <- proc.time()
  res <- parLapplyLB(cl, X = levels(bigdata$Replicate), fun = myboot)
  end <- proc.time()
  
  # shut down the cluster
  stopCluster(cl)
  
  # calculate proportion of models that successfully converged
  success <- sapply(res, is.numeric)
  mean(success)
  
  # combine successful results
  bigres <- do.call(cbind, res[success])
  
  # calculate 2.5th and 97.5th percentiles for 95% CI
  (ci <- t(apply(bigres, 1, quantile, probs = c(0.025, 0.975))))
  
  # All results
  finaltable <- cbind(Est = c(f, r), SE = c(se, NA), BootMean = rowMeans(bigres),
                      ci)
  # round and print
  print(round(finaltable, 3))
  
  # temporary data
  tmpdat <- mm_scores[, c("PRS","MutStatus","Sex","Family_ID")]
  
  summary(mm_scores$PRS)
  
  jvalues <- with(mm_scores, seq(from = min(PRS), to = max(PRS), length.out = 100))
  
  # calculate predicted probabilities and store in a list
  pp <- lapply(jvalues, function(j) {
    tmpdat$PRS <- j
    predict(m, newdata = tmpdat, type = "response")
  })
  
  # average marginal predicted probability across a few different Lengths of
  # Stay
  sapply(pp[c(1, 20, 40, 60, 80, 100)], mean)
  
  
  # get the means with lower and upper quartiles
  plotdat <- t(sapply(pp, function(x) {
    c(M = mean(x), quantile(x, c(0.25, 0.75)))
  }))
  
  # add in LengthofStay values and convert to data frame
  plotdat <- as.data.frame(cbind(plotdat, jvalues))
  
  # better names and show the first few rows
  colnames(plotdat) <- c("PredictedProbability", "Lower", "Upper", "PRS")
  head(plotdat)
  
  # plot average marginal predicted probabilities
  p <- ggplot(plotdat, aes(x = PRS, y = PredictedProbability)) + 
    geom_line() + 
    geom_linerange(aes(ymin = Lower, ymax = Upper)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = .15) +
    geom_line(size = 2) + ylim(c(0, 1)) +
    scale_fill_manual(values = c("#DE3533","#0047AB")) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size=14,family="Helvetica Neue")) 

  print(p)
    
}

dir="~/research/lfs_wt/PRS/PRS_results"

lfs_scores <- read_scores(dir,"lfs",clin)
leukemia_scores <- read_scores(dir,"leukemia",clin)
skin_scores <- read_scores(dir,"nonmelanoma_skin",clin)
glioma_scores <- read_scores(dir,"glioma",clin)
os_scores <- read_scores(dir,"os",clin)
breast_scores <- read_scores(dir,"bc",clin)

ggarrange(lfs_scores[[2]],leukemia_scores[[2]],breast_scores[[2]],
          nrow=3,ncol=1)
ggarrange(skin_scores[[2]],glioma_scores[[2]],os_scores[[2]],
          nrow=3,ncol=1)

ggarrange(lfs_scores[[2]],leukemia_scores[[2]],skin_scores[[2]],
          glioma_scores[[2]],os_scores[[2]],breast_scores[[2]],
          nrow=6,ncol=1)

get_fam_scores <- function(clin_scores,cutoff) {
    mm_scores <- clin_scores[c("Individual_ID","Family_ID","Sex","tt",cutoff,"tt_ofinterest")]
    mm_scores$Family_ID <- factor(mm_scores$Family_ID)
    mm_scores$MutStatus <- factor(str_replace_all(gsub('[[:digit:]]+', '', mm_scores$Family_ID),"MisFam",""))
    mm_scores$Sex <- factor(str_replace(str_replace(mm_scores$Sex,"2","Female"),"1","Male"))
    colnames(mm_scores) <-c("Individual_ID","Family_ID","Sex","tt","PRS","tt_ofinterest","MutStatus")
    model_plts <- mixed_model_covars(mm_scores)
    
    relevantfamilies <- mm_scores$Family_ID[mm_scores$tt_ofinterest == 1]
    tt_scores <- clin_scores[clin_scores$Family_ID %in% relevantfamilies,]
    Tumor_Label <- "Unaffected"
    tt_scores$tt_ofinterest <- ifelse(tt_scores$tt_ofinterest == 1,Tumor_Label,paste("Not",Tumor_Label))
    tt_scores <- tt_scores %>% 
      mutate_if(is.numeric, function(x){round(x,2)}) 
    plt <- ggplot(tt_scores,aes_string(x="Individual_ID", y=cutoff)) +
      facet_grid(~Family_ID, scales = "free_x", space = "free_x")+
      geom_bar( stat="identity", position="dodge",aes(fill=tt_ofinterest)) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            legend.position = "none",
            text = element_text(size=14,family="Helvetica Neue")) +   
      geom_text(aes_string(vjust=0.01, label = cutoff),size=2)+
      scale_fill_manual(values=c("grey94","red"))
    return(list(plt,model_plts[[1]],model_plts[[2]]))
}

lfs_prsmodel <- get_fam_scores(lfs_scores[[1]],"X0.0000000000000001")
leukemia_prsmodel <- get_fam_scores(leukemia_scores[[1]],"X0.0000000001")
skin_prsmodel <- get_fam_scores(skin_scores[[1]],"X0.0000000000000001")
glioma_prsmodel <- get_fam_scores(glioma_scores[[1]],"X0.0000000001")
os_prsmodel <- get_fam_scores(os_scores[[1]],"X0.0000001")
breast_prsmodel <- get_fam_scores(breast_scores[[1]],"X0.0000000000000001")


ggarrange(lfs_prsmodel[[1]],
          leukemia_prsmodel[[1]],
          breast_prsmodel[[1]],
          nrow=3,ncol=1,
          widths = c(1, 0.25,0.25))

ggarrange(skin_prsmodel[[1]],
          glioma_prsmodel[[1]],
          os_prsmodel[[1]],
          nrow=3,align="v")


ggarrange(leukemia_prsmodel[[2]],leukemia_prsmodel[[3]],
          breast_prsmodel[[2]],breast_prsmodel[[3]],
          skin_prsmodel[[2]],skin_prsmodel[[3]],
          nrow=3,ncol=2,align="h",labels = c("A","","B","","C",""))

