suppressMessages(require(sva))
suppressMessages(require(limma))
suppressMessages(require(dplyr))
suppressMessages(require(corrplot))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))
suppressMessages(require(gplots))
suppressMessages(require(ggfortify))
suppressMessages(require(argparse))
suppressMessages(library(umap))

parser <- ArgumentParser()
parser$add_argument("--value", action="store")
parser$add_argument("--infile", action="store")
parser$add_argument("--outdir", action="store")

args <- parser$parse_args()
value <- args$value
outdir <- args$outdir
setwd('/hpf/largeprojects/adam/projects/lfs/lfs_germline/methyl_data/')
source('Scripts/util_functions.R')

###################################################################################################################
# Combat batch correction
###################################################################################################################

cat("[ Reading in data ]","\n")

data <- readRDS(paste0('rds/',args$infile))
data <- data[data$Meth != "Problem",]
data <- data[!is.na(data$Project),]

pc <- prcomp(as.matrix(data[45:length(data)]), scale = TRUE)
pc_clin <- cbind(data[1:44], pc$x)
keep_ids <- remove_outliers(pc_clin,3)
data <- data[data$SentrixID %in% keep_ids,]
clin <- data[1:44] ; beta <- data.frame(data[45:length(data)])
beta <- as.matrix(t(sapply( beta, as.numeric )))

cat("[ Calling ComBat ]","\n")
modcombat <- model.matrix(~1, data = data)
batches <- as.numeric(data$Project)
beta_ComBat <- ComBat(dat=beta,batch=batches,mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
beta_ComBat <- data.frame(t(beta_ComBat))
data_ComBat <- cbind(clin,beta_ComBat)
saveRDS(data_ComBat,paste0('rds/',value,'.rds'))

cat("[ PCA after batch correction with ComBat ]","\n")
pc <- prcomp(as.matrix(beta_ComBat), scale = TRUE)
pc_clin <- cbind(clin, pc$x)
write.csv(pc_clin, paste0('Output/',value,'_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0(value,'_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,value,outdir)

###################################################################################################################
# Visualize differences before and after correction 
###################################################################################################################

cat("[ Plotting technical replicates ]","\n")

duplicated_data <- get_technicalreplicates(data_ComBat)
plot_concordance(duplicated_data,value)

pc_beta <- data.frame(duplicated_data[45:length(duplicated_data)], row.names = paste0(duplicated_data$ids," (",duplicated_data$array,")"))
pc <- prcomp(as.matrix(pc_beta), scale = TRUE)
pc_clin <- cbind(duplicated_data[1:44],pc$x)
write.csv(pc_clin,paste0('Output/',value,'_TechnicalReplicates_PCA.csv'),quote=F,row.names=F)

generate_pcsummary(pc_clin,paste0(value,'_TechnicalReplicates_PCA_summary.csv'),outdir)
generate_pcplots(pc_clin,paste0(value,'_TechnicalReplicates'),outdir)

u <- umap(beta_ComBat) ; ud <- cbind(data_ComBat[1:45],u$layout)
write.csv(ud,paste0('Output/Umap_',value,'.csv'))
