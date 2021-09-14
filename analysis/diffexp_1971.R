library(fgsea) 
library(dplyr)
library(tibble)
library(ggplot2)


diffexp <- read.csv('~/research/lfs_wt/Affy_U133PLUS2_S-SM-genes_up&down_Valli.csv')[c("Genes","Median_FCh")]
vus <- read.csv('~/research/lfs_wt/family_2/1971.hg19.pathogenic.vus.txt',sep='\t',stringsAsFactors = F)
vus_diffexp <- vus[vus$Ref.Gene %in% diffexp$Genes,]
vus_diffexp <- vus_diffexp[vus_diffexp$gnomAD_genome_ALL < 0.001,]
vus_diffexp <- merge(vus_diffexp,diffexp,by.x="Ref.Gene",by.y="Genes")
vus <- vus_diffexp[vus_diffexp$Func.refGene %in% c("exonic","splicing"),]
ranks <- deframe(diffexp)
pathways.hallmark <- gmtPathways("~/research/resources/h.all.v7.0.symbols.gmt")
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

pathways.kegg <- gmtPathways("~/research/resources/c2.cp.kegg.v7.0.symbols.gmt")
fgseaKegg <- fgsea(pathways=pathways.kegg, stats=ranks, nperm=1000)


pathways.reactome <- gmtPathways("~/research/resources/c2.cp.reactome.v7.0.symbols.gmt")
fgseaReactome <- fgsea(pathways=pathways.reactome, stats=ranks, nperm=1000)

pathways.tftargs <- gmtPathways("~/research/resources/c3.tft.v7.0.symbols.gmt")
fgseaTFtargs <- fgsea(pathways=pathways.tftargs, stats=ranks, nperm=1000)

pathways.onco <- gmtPathways("~/research/resources/c6.all.v7.0.symbols.gmt")
fgseaOnco <- fgsea(pathways=pathways.onco, stats=ranks, nperm=1000)

pathways.go <- gmtPathways("~/research/resources/c5.bp.v7.0.symbols.gmt")
fgseaGo <- fgsea(pathways=pathways.go, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

