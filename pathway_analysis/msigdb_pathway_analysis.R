library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)
library(msigdbr)
library(ggpubr)

# calculate how many genes are involved in each pathway
format_msigdb <- function(db) {
  db$genes <- str_replace_all(db$genes,';',"")
  msigdb.sep <- db %>% separate(genes, paste0("Gene", 0:467), sep="\\|")
  msigdb.long <- melt(msigdb.sep, id.vars=c("gs_name","n"), measure.vars=colnames(msigdb.sep)[2:(length(msigdb.sep)-1)]) %>% na.omit() %>% as.data.table()
  msigdb.long$sample_frequency <- 0
  msigdb.long$gene_frequency <- 0
  return(msigdb.long)
}

inc <- function(x)
{
  eval.parent(substitute(x <- x + 1))
}

sample_frequency <- function(msigdb.long, genes) {
  explored_genes <- list()
  for (gene in genes) {
    if (!gene %in% explored_genes) {
      inc(msigdb.long[msigdb.long$value == gene]$gene_frequency)
    }
    inc(msigdb.long[msigdb.long$value == gene]$sample_frequency)
    explored_genes[[gene]] <- gene
  }
  msigdb.pathways <- msigdb.long %>% 
    dplyr::group_by(gs_name,n) %>% 
    dplyr::summarise(sample_frequency = sum(sample_frequency),gene_frequency = sum(gene_frequency),genes = paste(unique(value), collapse = ','))
  return(msigdb.pathways)
}

fisher_test <- function(genes) {
  universe_genes = unlist(unique(unlist(strsplit(paste0(str_replace_all(m_df$genes,";",""),collapse="|"), "\\|"))))
  p_values = list()
  msigdb.pathways <- sample_frequency(msigdb.long, genes)
  for (i in 1:nrow(msigdb.pathways)) {
    if (msigdb.pathways[i,]$gene_frequency > 1) {
      gene_nopath <- length(unique(genes)) - msigdb.pathways[i,]$gene_frequency
      nogene_inpath <- msigdb.pathways[i,]$n - msigdb.pathways[i,]$gene_frequency
      nogene_nopath <- 21306-length(unique(c(genes,unlist(str_split(msigdb.pathways[i,]$genes,',')))))-3
      x <- matrix(c(msigdb.pathways[i,]$gene_frequency, gene_nopath, nogene_inpath, nogene_nopath), byrow = TRUE, 2, 2)
      p <- fisher.test(x,alternative="greater",simulate.p.value=T)
      p <- p$p.value
      p_values[[i]] <- p
    } else {
      p_values[[i]] <- NA
    }
  }
  msigdb.pathways$pval <- unlist(p_values)
  qval <- p.adjust(msigdb.pathways$pval,method="fdr")
  msigdb.pathways$qval <- qval
  #normalized enrichment score = -log10(qval)/mean(es)
  msigdb.pathways$nes <- -log(msigdb.pathways$pval)*log(msigdb.pathways$sample_frequency)
  msigdb.pathways <- msigdb.pathways[order(msigdb.pathways$pval),]
  return(msigdb.pathways)
}

plot_pathways <- function(pathways_all) {
  if (dim(pathways_all)[1] < 20) {
    pathways <- pathways_all
  } else {
    pathways <- pathways_all[order(pathways_all$nes, decreasing = TRUE),][1:20,]
  }
  p <- ggplot(pathways, aes(nes, gs_name))
  p <- p + geom_point(aes(colour=qval, size=sample_frequency)) +
    scale_color_gradientn(colours=rainbow(10), limits=c(0, 1)) +
    geom_vline(xintercept=0, size=0.5, colour="gray50") +
    theme(panel.background=element_rect(fill="white", colour="gray95"),
          panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
          axis.title.y=element_blank()) +
    expand_limits(x=c(-1,10)) +
    scale_x_continuous(breaks=c(seq(-1,8))) +
    scale_y_discrete(limits=rev(pathways$gs_name))
  p
}

msigdb.long <- format_msigdb(m_df)
