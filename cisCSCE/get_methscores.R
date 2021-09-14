library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)

data <- readRDS("/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/rds/NoobComBat_beta_PEER_residuals_covs_RF.rds")
meth <- data.frame(t(data[13:length(data)]))

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% data.frame()
meth$Relation_to_Island <- anno$Relation_to_Island[match(rownames(meth),anno$Name)]

get_sampscores <- function(meth,relation) {
  meth_relation <- meth[meth$Relation_to_Island == relation,]
  meth_relation$Relation_to_Island <- NULL
  samp_scores <- colSums(meth_relation)
  return(samp_scores)
}

all_scores <- list()
relations <- unique(anno$Relation_to_Island)
for (i in 1:length(relations)) {
  relation_scores <- get_sampscores(meth,relations[i])
  all_scores[[i]] <- relation_scores
}

all_scores <- do.call("rbind",all_scores)
all_scores <- data.frame(t(all_scores))
colnames(all_scores) <- relations
clin_scores <- cbind(data[1:12],all_scores)
write.table(clin_scores, 'Output/methScore_relationToIsland.txt',sep='\t',quote=F,row.names=F)

