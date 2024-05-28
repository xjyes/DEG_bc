setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/DEG breast Primary_CNS")


library(sva)
library("dplyr")
load("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/Primary_CNS/exp.Rdata")
load("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/Primary_CNS/group_list.Rdata")

data3521$symbol = rownames(data3521)
data14682$symbol = rownames(data14682)
data14683$symbol = rownames(data14683)
merge_eset = inner_join(data3521,data14682, by = "symbol")
dim(merge_eset)
merge_eset = inner_join(merge_eset,data14683, by = "symbol")
dim(merge_eset)

rownames(merge_eset) = merge_eset$symbol
merge_eset = merge_eset[, c(-7)]


exp = as.matrix(merge_eset)
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
boxplot(data)

batchType = c(rep(1,6), rep(2,60), rep(3,20))
modType = c(gl3521, gl14682, gl14683)
mod = model.matrix(~as.factor(modType))
outTab = data.frame(ComBat(data, batchType, mod, par.prior = TRUE))
boxplot(outTab)

save(outTab, file = "aggreg_exp.Rdata")


rm(list = ls())
load("aggreg_exp.Rdata")
load("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/Primary_CNS/group_list.Rdata")
library(limma)
group_list = c(gl3521, gl14682, gl14683)
design <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exp)

contrast.matrix <- makeContrasts(paste0(c("treat","control"), collapse = "-"),levels = design)
contrast.matrix

fit <- lmFit(outTab, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput = topTable(fit2, coef = 1, n = Inf)
nrDEG = na.omit(tempOutput)
head(nrDEG)

save(nrDEG, file = "comprehensive_DEG.Rdata")
write.csv(nrDEG, file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/DEG_composed.csv")


# Volcano plot
rm(list = ls())
library(ggplot2)
load("comprehensive_DEG.Rdata")
DEG = nrDEG
logFC_cutoff <- with(DEG, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff, "up", "down"), "not"))
title <- ("TREAT vs CONTROL")
sfg = DEG[DEG$change != "not", ]
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(P.value)") +
  ggtitle(title) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('deepskyblue','black','red'), name = "p.value<0.05") +
  scale_size_manual(values = c(1.75, 0.5, 1.75)) +
  guides(size = FALSE)
g

write.csv(sfg, file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/results/DEG breast Primary_CNS/sfg_composed.csv")


# heatmap
rm(list = ls())
options(stringsAsFactors = F)
load("comprehensive_DEG.Rdata")
load("aggreg_exp.Rdata")
load("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/Primary_CNS/group_list.Rdata")
library(pheatmap)
choose_gene = rownames(nrDEG)
choose_matrix = outTab[choose_gene,]
matrix_scale = t(scale(t(choose_matrix)))
annotation_col = data.frame(sampletype = c(gl3521, gl14682, gl14683))
row.names(annotation_col) <- colnames(choose_matrix)
pheatmap(matrix_scale, scale = "row", annotation_col = annotation_col,
         cluster_rows = T, cluster_cols = F,show_rownames = F,show_colnames = F,border_color = NA)




























