setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/DEG breast Primary_CNS")

# Enrichment analysis: data preparation
# add the symbol column
rm(list = ls())
options(stringsAsFactors = F)
# load("comprehensive_DEG.Rdata")
nrDEG = read.csv("/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/DEG_3521.csv")
rownames(nrDEG) <- nrDEG$X
nrDEG <- nrDEG[,-1]

library(dplyr)
deg = nrDEG
deg = mutate(deg,symbol = rownames(deg))
head(deg)

# add the change column: up and down regulation
logFC_t <- with(deg, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
change = ifelse(deg$logFC > logFC_t,'up',
                ifelse(deg$logFC < -logFC_t,'down','stable'))
deg <- mutate(deg, change)
head(deg)
table(deg$change)

# add the ENTREZID column: map gene symbols to ENTREZID
library(ggplot2)
library(R.utils)
library(clusterProfiler)
library(org.Hs.eg.db) # Load human gene annotation
R.utils::setOption("clusterProfiler.download.method","auto")
s2e <- bitr(unique(deg$symbol),fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
head(s2e)
head(deg)
deg <- inner_join(deg,s2e,by=c("symbol" = "SYMBOL"))

head(deg)
save(deg, file = "deg_with_label_3521.Rdata")


# Enrichment analysis #

# KEGG pathway analysis
rm(list = ls())
options(stringsAsFactors = F)
library(R.utils)
library(clusterProfiler)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method","auto")
load("deg_with_label_3521.Rdata")
gene_up = deg[deg$change == 'up','ENTREZID']
gene_down = deg[deg$change == 'down', 'ENTREZID']
gene_diff = c(gene_up, gene_down)
gene_all = deg[,'ENTREZID']
kk.up <- enrichKEGG(gene = gene_up,
                    organism = 'hsa',
                    universe = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
head(kk.up)[,1:6]
dim(kk.up)
kk.down <- enrichKEGG(gene = gene_down,
                      organism = 'hsa',
                      universe = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
head(kk.down)[,1:6]
dim(kk.down)
kk.diff <- enrichKEGG(gene = gene_diff,
                      organism = 'hsa',
                      pvalueCutoff = 0.05)
head(kk.diff)[,1:6]
class(kk.diff)
dim(kk.diff)
# Change entrez id to gene symbol
kk.up =DOSE::setReadable(kk.up, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk.down =DOSE::setReadable(kk.down, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
kk.diff =DOSE::setReadable(kk.diff, OrgDb='org.Hs.eg.db',keyType='ENTREZID')

# extract enrichment result data frame
kegg_diff_dt <- kk.diff@result
write.csv(kegg_diff_dt,"/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/kegg_3521.csv")



# Visualization
dotplot(kk.diff)
barplot(kk.diff)
# kegg_plot <- function(up_kegg, down_kegg){
#   dat=rbind(up_kegg,down_kegg)
#   colnames(dat)
#   dat$pvalue = -log10(dat$pvalue)
#   dat$pvalue = dat$pvalue*dat$group
#   
#   dat = dat[order(dat$pvalue, decreasing = F),]
#   
#   g_kegg <- ggplot(dat, aes(x=reorder(Description,order(pvalue,decreasing = F)),y = pvalue, fill = group)) +
#     geom_bar(stat="identity", width = 0.8) +
#     scale_fill_gradient(low="blue",high="red",guide = FALSE) +
#     scale_x_discrete(name = "Pathway names") +
#     scale_y_continuous(name = "log10P-value")+
#     coord_flip() + theme_bw() + theme(plot.title = element_text(hjust = 0.5))+
#     ggtitle("Pathway Enrichment")
# }
# 
# g_kegg <- kegg_plot(up_kegg, down_kegg)
# g_kegg

# GO database analysis #
rm(list = ls())
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) # load human gene annotation
options(stringsAsFactors = F)
load("deg_with_label.Rdata")

gene_up = deg[deg$change == 'up', 'ENTREZID']
gene_down = deg[deg$change == 'down', 'ENTREZID']
gene_diff = c(gene_up,gene_down)
head(deg)

# Cellular component
ego_CC <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)



# Biological process
ego_BP <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)


# Molecular function
ego_MF <- enrichGO(gene = gene_diff,
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)


# All
ego_ALL <- enrichGO(gene = gene_diff,
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.01,
                    readable = TRUE)

# Plot bar
barplot(ego_CC, showCategory=10, font.size = 8, title = "EnrichmentGO_CC")
barplot(ego_BP, showCategory=10, font.size = 8, title = "EnrichmentGO_BP")
barplot(ego_MF, showCategory=10, font.size = 8, title = "EnrichmentGO_MF")
barplot(ego_ALL, font.size = 8, title = "EnrichmentGO_ALL",split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free")
# Plot dot
dotplot(ego_CC, showCategory=10, font.size = 8, title="EnrichmentGO_CC")
dotplot(ego_BP, showCategory=10, font.size = 8, title="EnrichmentGO_BP")
dotplot(ego_MF, showCategory=10, font.size = 8, title="EnrichmentGO_MF")
dotplot(ego_ALL, showCategory=8, font.size = 8, title = "EnrichmentGO_ALL",split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY~.,scale="free")

# Get enriched genes and logFC in GO
library(dplyr)
library(tidyr)
GO = ego_ALL@result
GO_ID <- data.frame("gene_id"=GO$geneID)
GO_ID <- GO_ID %>% separate(gene_id, sep = "/", into = as.character(1:max(GO$Count)))
GO_ID["logFC"] = deg[match(GO_ID[,1],deg$symbol),"logFC"]
GO_ID["logFC"] = substr(GO_ID$logFC,1,4)
GO_ID = tidyr::unite(GO_ID,"symbol_logFC",as.character(1),logFC,sep = "_")

for (i in 2:max(GO$Count)){
  GO_ID["logFC"] = deg[match(GO_ID[,2],deg$symbol),"logFC"]
  GO_ID["logFC"] = substr(GO_ID$logFC,1,4)
  GO_ID = tidyr::unite(GO_ID,"temp",2,logFC,sep = "_",na.rm = T)
  GO_ID = tidyr::unite(GO_ID,"symbol_logFC",symbol_logFC,temp,sep = "/", na.rm = T)
}
GO_ID["symbol_logFC"] = gsub("/*$","",GO_ID$symbol_logFC)
GO["symbol_logFC"] = GO_ID$symbol_logFC

write.table(GO,"/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/GO_integrate.csv",row.names=FALSE,col.names=TRUE,sep=",")


# GSEA analysis

rm(list = ls())
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) # load human gene annotation
options(stringsAsFactors = F)
load("deg_with_label_3521.Rdata")
geneList <- deg$logFC
names(geneList) <- deg$ENTREZID
geneList <- sort(geneList, decreasing = T)
halmt<-read.gmt("/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/Primary_CNS/h.all.v2023.1.Hs.entrez.gmt") #读gmt文件
hallmark <-GSEA(geneList, TERM2GENE = halmt,eps = 0) #GSEA分析
print(hallmark@result)
dotplot(hallmark,font.size = 10) #出点图 
dotplot(hallmark,color="pvalue")  #按p值出点图 
write.table(hallmark@result,"/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/HALLMARK_3521.csv",row.names=FALSE,col.names=TRUE,sep=",")












