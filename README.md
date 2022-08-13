# TCGA2
options("repos" = c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

cran_packages <- c('tidyr',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   'factoextra',
                   'FactoMineR',
                   'pheatmap',
                   "survival",
                   "survminer",
                   "patchwork",
                   "ggstatsplot",
                   "ggplotify",
                   "cowplot",
                   "glmnet",
                   "ROCR",
                   "caret",
                   "randomForest",
                   "survminer",
                   "Hmisc",
                   "e1071",
                   "deconstructSigs",
                   "AnnoProbe",
                   "timeROC",
                   "circlize",
                   "tinyarray"
) 
Biocductor_packages <- c("limma",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "SummarizedExperiment",
                         "DESeq2",
                         "edgeR",
                         "ggpubr",
                         "rtracklayer",
                         "genefilter",
                         "maftools",
                         "ComplexHeatmap",
                         "GDCRNATools"
)


for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


# use BiocManager to install
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}



for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}
proj = "TCGA-UCEC"
if(F){
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".htseq_counts.tsv.gz"),destfile = paste0(proj,".htseq_counts.tsv.gz"))
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".GDC_phenotype.tsv.gz"),destfile = paste0(proj,".GDC_phenotype.tsv.gz"))
  download.file(url = paste0("https://gdc.xenahubs.net/download/",proj, ".survival.tsv"),destfile = paste0(proj,".survival.tsv"))
}
clinical = read.delim(paste0(proj,".GDC_phenotype.tsv.gz"),fill = T,header = T,sep = "\t")
surv = read.delim(paste0(proj,".survival.tsv"),header = T)
clinical[1:4,1:4]
head(surv)
dat = read.table(paste0(proj,".htseq_counts.tsv.gz"),check.names = F,row.names = 1,header = T)
dat = as.matrix(2^dat - 1)
dat[1:4,1:4]
dat[97,9]
as.character(dat[97,9])
exp = apply(dat, 2, as.integer)
exp[1:4,1:4]
rownames(exp) = rownames(dat)
exp[1:4,1:4]


library(stringr)
head(rownames(exp))
library(AnnoProbe)
annoGene(rownames(exp),ID_type = "ENSEMBL")
rownames(exp) = str_split(rownames(exp),"\\.",simplify = T)[,1];head(rownames(exp))
library(tinyarray)
exp = trans_array(exp,ids = re,from = "ENSEMBL",to = "SYMBOL")
exp[1:4,1:4]
nrow(exp)
exp1 = exp[rowSums(exp)>0,]
nrow(exp1)
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
nrow(exp)
table(str_sub(colnames(exp),14,15))
Group = ifelse(as.numeric(str_sub(colnames(exp),14,15)) < 10,'tumor','normal')
Group = factor(Group,levels = c("normal","tumor"))
table(Group)
save(exp,Group,proj,clinical,surv,file = paste0(proj,".Rdata"))


rm(list = ls())
load("TCGA-UCEC.Rdata")
table(Group)
library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
if(!file.exists(paste0(proj,"_dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(proj,"_dd.Rdata"))
}
load(file = paste0(proj,"_dd.Rdata"))
class(dds)
library(DESeq2)
colData <- data.frame(row.names =colnames(exp), 
                      condition=Group)
if(!file.exists(paste0(proj,"_dd.Rdata"))){
  dds <- DESeqDataSetFromMatrix(
    countData = exp,
    colData = colData,
    design = ~ condition)
  dds <- DESeq(dds)
  save(dds,file = paste0(proj,"_dd.Rdata"))
}
load(file = paste0(proj,"_dd.Rdata"))
class(dds)
res <- results(dds, contrast = c("condition",rev(levels(Group))))
c("condition",rev(levels(Group)))
#> [1] "condition" "tumor"     "normal"
class(res)
DEG1 <- as.data.frame(res)
DEG1 <- DEG1[order(DEG1$pvalue),] 
DEG1 = na.omit(DEG1)
DEG1$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG1$change)
head(DEG1)
logFC_t = 2
pvalue_t = 0.05

k1 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange < -logFC_t);table(k1)

k2 = (DEG1$pvalue < pvalue_t)&(DEG1$log2FoldChange > logFC_t);table(k2)
head(DEG1)

library(edgeR)

dge <- DGEList(counts=exp,group=Group)
dge$samples$lib.size <- colSums(dge$counts)
dge <- calcNormFactors(dge) 

design <- model.matrix(~Group)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
fit <- glmLRT(fit) 

DEG2=topTags(fit, n=Inf)
class(DEG2)

DEG2=as.data.frame(DEG2)
head(DEG2)

k1 = (DEG2$PValue < pvalue_t)&(DEG2$logFC < -logFC_t);table(k1)

k2 = (DEG2$PValue < pvalue_t)&(DEG2$logFC > logFC_t);table(k2)

DEG2$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))

head(DEG2)

table(DEG2$change)

library(limma)
dge <- DGEList(counts=exp)
dge <- calcNormFactors(dge)
v <- voom(dge,design, normalize="quantile")

design <- model.matrix(~Group)
fit <- lmFit(v, design)
fit= eBayes(fit)

DEG3 = topTable(fit, coef=2, n=Inf)
DEG3 = na.omit(DEG3)

k1 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC < -logFC_t);table(k1)
k2 = (DEG3$P.Value < pvalue_t)&(DEG3$logFC > logFC_t);table(k2)
DEG3$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG3$change)
head(DEG3)
tj = data.frame(deseq2 = as.integer(table(DEG1$change)),
                edgeR = as.integer(table(DEG2$change)),
                limma_voom = as.integer(table(DEG3$change)),
                row.names = c("down","not","up")
);tj
save(DEG1,DEG2,DEG3,Group,tj,file = paste0(proj,"_DEG.Rdata"))

#可视化
library(ggplot2)
library(tinyarray)
exp[1:4,1:4]
dat = log2(cpm(exp)+1)
pca.plot = draw_pca(dat,Group);pca.plot
save(pca.plot,file = paste0(proj,"_pcaplot.Rdata"))
cg1 = rownames(DEG1)[DEG1$change !="NOT"]
cg2 = rownames(DEG2)[DEG2$change !="NOT"]
cg3 = rownames(DEG3)[DEG3$change !="NOT"]

h1 = draw_heatmap(dat[cg1,],Group,n_cutoff = 2)
h2 = draw_heatmap(dat[cg2,],Group,n_cutoff = 2)
h3 = draw_heatmap(dat[cg3,],Group,n_cutoff = 2)

v1 = draw_volcano(DEG1,pkg = 1,logFC_cutoff = logFC_t)
v2 = draw_volcano(DEG2,pkg = 2,logFC_cutoff = logFC_t)
v3 = draw_volcano(DEG3,pkg = 3,logFC_cutoff = logFC_t)

library(patchwork)
(h1 + h2 + h3) / (v1 + v2 + v3) +plot_layout(guides = 'collect') &theme(legend.position = "none")
ggsave(paste0(proj,"_heat_vo.png"),width = 15,height = 10)

#三大R包差异基因对比
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

up = intersect(intersect(UP(DEG1),UP(DEG2)),UP(DEG3))
down = intersect(intersect(DOWN(DEG1),DOWN(DEG2)),DOWN(DEG3))
dat = log2(cpm(exp)+1)
hp = draw_heatmap(dat[c(up,down),],Group,n_cutoff = 2)
#
install.packages('VennDiagram')
library('VennDiagram')
#上调、下调基因分别画维恩图
up_genes = list(Deseq2 = UP(DEG1),
                edgeR = UP(DEG2),
                limma = UP(DEG3))

down_genes = list(Deseq2 = DOWN(DEG1),
                  edgeR = DOWN(DEG2),
                  limma = DOWN(DEG3))

up.plot <- draw_venn(up_genes,"UPgene")
down.plot <- draw_venn(down_genes,"DOWNgene")
#维恩图拼图

library(patchwork)
#up.plot + down.plot
# 拼图
pca.plot + hp+up.plot +down.plot+ plot_layout(guides = "collect")
ggsave(paste0(proj,"_heat_ve_pca.png"),width = 15,height = 10)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#f87669","#2fa1dd")),
                       labels = c("tumor","normal"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

m = Heatmap(t(scale(t(exp[c(up,down),]))),name = " ",
            col = col_fun,
            top_annotation = top_annotation,
            column_split = Group,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL)
m



