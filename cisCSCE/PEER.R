library(peer)
library(ggplot2)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--value", action="store")

args <- parser$parse_args()

setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')
source('Scripts/util_functions.R')

cat("[ Read in data ]","\n")
data <- readRDS(paste0('rds/',args$value,'.rds'))
meth <- data[grepl("cg",colnames(data))]
covs <- data[c("array","Project","cancerstatus","agesamplecollection","gender","tm_donor")]
model = PEER()

cat("[ Run PEER model ]","\n")
PEER_setCovariates(model, as.matrix(covs))
PEER_setPhenoMean(model,as.matrix(meth))
PEER_setNk(model,100)

cat("[ Number of Hidden Factors] : ")
PEER_getNk(model)
PEER_update(model)

cat("[ Get output ]","\n")
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

covs <- "apcsgt"
cat("[ Plot precision ]","\n")
pdf(paste0('Plots/',args$value,'_PEER_precision_',covs,'.pdf'),height=7, width=9)
plot(precision)
dev.off()

cat("[ Save object ]","\n")
saveRDS(factors,paste0('rds/',args$value,'_PEER_factors_',covs,'.rds'))
saveRDS(weights,paste0('rds/',args$value,'_PEER_weights_',covs,'.rds'))
saveRDS(precision,paste0('rds/',args$value,'_PEER_precision_',covs,'.rds'))
saveRDS(residuals,paste0('rds/',args$value,'_PEER_residuals_',covs,'.rds'))

pc <- prcomp(residuals, scale = TRUE)
pc_clin <- cbind(data[!grepl("cg",colnames(data))],pc$x)
write.csv(pc_clin,paste0('Output/',args$value,'_PEER_PCA_',covs,'.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin, paste0(args$value,'_PEER_PCA_',covs,'_summary.csv'))
generate_pcplots(pc_clin, paste0(args$value,'_PEER_',covs))

