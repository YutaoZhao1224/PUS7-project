library(DESeq2)
library(ggplot2)
library(dplyr)
library(FactoMineR)
library(tibble)
DEseq_custom <- function(RPKM, Count, CTR, TRE, FC = 1){
  cond <- colnames(Count) %>% strsplit('_')
  Condition <- c()
  for (i in cond){
    Condition <- append(Condition, i[1])
  }
  Condition <- factor(Condition)
  ColData <- data.frame(row.names=colnames(Count), Condition)
  dds <- DESeqDataSetFromMatrix(Count, colData = ColData, design = ~ Condition)
  dds$Condition <- relevel(dds$Condition,ref=as.character(Condition[1]))
  dds <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res <- results(dds, contrast=c("Condition",TRE,CTR))
  res_de <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res_de <- rownames_to_column(res_de,"Gene_name")
  RPKM <- RPKM %>% filter(rowMeans(.)> 0.4) %>% rownames_to_column("Gene_name")
  mutate(res_de, direction = if_else(abs(log2FoldChange) < FC | pvalue > 0.05, 'ns', 
                                     if_else(log2FoldChange >= FC, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>% 
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, pvalue, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE) %>% summarise(counts = n()) %>% print()
  
  output <- mutate(res_de, direction = if_else(abs(log2FoldChange) < FC | pvalue > 0.05, 'ns', 
                                               if_else(log2FoldChange >= FC, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>%
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, pvalue, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE)
  name <- paste(CTR,'_vs_',TRE,'.DEresults.txt', sep = '')
  write.table(output,name,quote = F, 
              sep = "\t", row.names = F, col.names = T)
  return(output)
}

DEseq_custom_padj <- function(RPKM, Count, CTR, TRE, FC = 1){
  cond <- colnames(Count) %>% strsplit('_')
  Condition <- c()
  for (i in cond){
    Condition <- append(Condition, i[1])
  }
  Condition <- factor(Condition)
  ColData <- data.frame(row.names=colnames(Count), Condition)
  dds <- DESeqDataSetFromMatrix(Count, colData = ColData, design = ~ Condition)
  dds$Condition <- relevel(dds$Condition,ref=as.character(Condition[1]))
  dds <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
  res <- results(dds, contrast=c("Condition",TRE,CTR))
  res_de <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res_de <- rownames_to_column(res_de,"Gene_name")
  RPKM <- RPKM %>% filter(rowMeans(.)> 0.4) %>% rownames_to_column("Gene_name")
  #print(res_de)
  mutate(res_de, direction = if_else(abs(log2FoldChange) < FC | padj > 0.05, 'ns', 
                                     if_else(log2FoldChange >= FC, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>% 
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, padj, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE) %>% summarise(counts = n()) %>% print()
  
  output <- mutate(res_de, direction = if_else(abs(log2FoldChange) < FC | padj > 0.05, 'ns', 
                                               if_else(log2FoldChange >= FC, 'up', 'down'))) %>% 
    left_join(RPKM,by = "Gene_name") %>% na.omit() %>%
    select(Gene_name, contains(c(CTR,TRE)),log2FoldChange, padj, direction) %>% 
    group_by(direction) %>% arrange(desc(abs(log2FoldChange)),.by_group = TRUE)
  name <- paste(CTR,'_vs_',TRE,'.DEresults.txt', sep = '')
  write.table(output,name,quote = F, 
              sep = "\t", row.names = F, col.names = T)
  return(output)
}

sicon_vs_si7sk_RPKM <- read.delim("~/workspace/RNA-seq/HCT116-PUS7/mapping/sicon_vs_si7sk.RPKM.txt",row.names = 1)
sicon_vs_si7sk_Counts <- read.delim("~/workspace/RNA-seq/HCT116-PUS7/mapping/sicon_vs_si7sk.Counts.txt",row.names = 1)
sicon_vs_sipus7_RPKM <- read.delim("~/workspace/RNA-seq/HCT116-PUS7/mapping/sicon_vs_sipus7.RPKM.txt",row.names = 1)
sicon_vs_sipus7_Counts <- read.delim("~/workspace/RNA-seq/HCT116-PUS7/mapping/sicon_vs_sipus7.Counts.txt",row.names = 1)
total_RPKM <- cbind(sicon_vs_si7sk_RPKM,sicon_vs_sipus7_RPKM[,4:6])
total_Counts <- cbind(sicon_vs_si7sk_Counts,sicon_vs_sipus7_Counts[,4:6])

###PCA
RPKM1 <- filter(total_RPKM, rowMeans(select(total_RPKM,1:9))>1)
gene.pca <- PCA(t(RPKM1), ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(gene.pca)

###DEseq2
DE_sicon_vs_si7sk <- DEseq_custom_padj(sicon_vs_si7sk_RPKM,sicon_vs_si7sk_Counts,"sicon","si7sk")
DE_sicon_vs_sipus7 <- DEseq_custom_padj(sicon_vs_sipus7_RPKM,sicon_vs_sipus7_Counts,"sicon","sipus7")
### overlap
si7sk_down <- filter(DE_sicon_vs_si7sk,direction=="down") %>% pull(Gene_name)
si7sk_up <- filter(DE_sicon_vs_si7sk,direction=="up") %>% pull(Gene_name)
sipus7_down <- filter(DE_sicon_vs_sipus7,direction=="down") %>% pull(Gene_name)
sipus7_up <- filter(DE_sicon_vs_sipus7,direction=="up") %>% pull(Gene_name)

intersect(si7sk_down,sipus7_down) %>% length()
write.table(intersect(si7sk_down,sipus7_down),file = "down.gene.txt",quote = F, row.names = F)
intersect(si7sk_up,sipus7_up) %>% length()
write.table(intersect(si7sk_up,sipus7_up),file = "up.gene.txt",quote = F, row.names = F)


intersect(si7sk_up,sipus7_down) %>% length()
intersect(si7sk_down,sipus7_up) %>% length()
### GO analysis
library(clusterProfiler)
#BiocManager::install("org.Hs.eg.db")
GO_analysis <- function(Gene_list){
  library(org.Hs.eg.db)
  Transformed <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                        keys = Gene_list,
                        column = "ENTREZID",
                        multiVals = "first",
                        keytype ="SYMBOL" ) %>% unname()
  ggo <- enrichGO(gene     = Transformed,
                  OrgDb    = org.Hs.eg.db,
                  ont      = "BP")
  return(ggo)
}
up_GO <- GO_analysis(sipus7_up)
head(up_GO,10)
down_GO <- GO_analysis(sipus7_down)
head(down_GO,10)

si7sk_up_GO <- GO_analysis(si7sk_up)
head(up_GO,10)
si7sk_down_GO <- GO_analysis(si7sk_down)
head(down_GO,10)

### KEGG-GSEA
GSEA_function <- function(Gene_list, Gene_list_FC){
  Transformed <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                                       keys = Gene_list,
                                       column = "ENTREZID",
                                       multiVals = "first",
                                       keytype ="SYMBOL" ) %>% unname()
  names_for_GSEA <- data.frame(value = Gene_list_FC, Gene_name=Transformed) %>% 
    na.omit() %>% pull(Gene_name)
  values_for_GSEA <- data.frame(value = Gene_list_FC, Gene_name=Transformed) %>% 
    na.omit() %>% pull(value)
  names(values_for_GSEA) <- names_for_GSEA
  values_for_GSEA <- sort(values_for_GSEA, decreasing = T)
  KEGG_kk_entrez <- gseKEGG(geneList     = values_for_GSEA,
                            organism     = "hsa", 
                            pvalueCutoff = 0.05)
  KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                               OrgDb=org.Hs.eg.db,
                               keyType='ENTREZID')
  kk_gse <- KEGG_kk
  
  kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.05 & abs(kk_gse$NES)>0.5]
  kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
  kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
  
  down_gsea <- kk_gse_cut_down[order(kk_gse_cut_down$NES,decreasing = F),]
  up_gsea <- kk_gse_cut_up[order(kk_gse_cut_up$NES,decreasing = T),]
  return(list(down=down_gsea,up=up_gsea))
}
a <- GSEA_function(DE_sicon_vs_sipus7$Gene_name,DE_sicon_vs_sipus7$log2FoldChange)

View(a$down)
View(a$up)




### GSEA for customed gene set-Hallmark
GSEA_function <- function(Gene_list, Gene_list_FC){
  Transformed <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                                       keys = Gene_list,
                                       column = "ENTREZID",
                                       multiVals = "first",
                                       keytype ="SYMBOL" ) %>% unname()
  names_for_GSEA <- data.frame(value = Gene_list_FC, Gene_name=Transformed) %>% 
    na.omit() %>% pull(Gene_name)
  values_for_GSEA <- data.frame(value = Gene_list_FC, Gene_name=Transformed) %>% 
    na.omit() %>% pull(value)
  names(values_for_GSEA) <- names_for_GSEA
  values_for_GSEA <- sort(values_for_GSEA, decreasing = T)
  gmt_set <- read.gmt("~/Genome/tools_custom/RNA-seq/GSEA_set/h.all.v2023.2.Hs.entrez.gmt")
  KEGG_kk_entrez <- GSEA(geneList     = values_for_GSEA,
                         TERM2GENE = gmt_set,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
  KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez, 
                               OrgDb=org.Hs.eg.db,
                               keyType='ENTREZID')
  kk_gse <- KEGG_kk
  
  kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.05 & abs(kk_gse$NES)>0.5]
  kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
  kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
  
  down_gsea <- kk_gse_cut_down[order(kk_gse_cut_down$NES,decreasing = F),]
  up_gsea <- kk_gse_cut_up[order(kk_gse_cut_up$NES,decreasing = T),]
  return(list(all=KEGG_kk, down=down_gsea,up=up_gsea))
}
a <- GSEA_function(DE_sicon_vs_sipus7$Gene_name,DE_sicon_vs_sipus7$log2FoldChange)
b <- GSEA_function(DE_sicon_vs_si7sk$Gene_name,DE_sicon_vs_si7sk$log2FoldChange)
siPUS7_P53_pathway <- enrichplot::gseaplot2(a$all,
                                            a$up$ID[6],#富集的ID编号
                                            title = a$up$ID[6],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T) #显示pvalue等信息
siP7SK_P53_pathway <- enrichplot::gseaplot2(b$all,
                                            b$up$ID[4],#富集的ID编号
                                            title = b$up$ID[4],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T)

siPUS7_TNF_pathway <- enrichplot::gseaplot2(a$all,
                                            a$up$ID[2],#富集的ID编号
                                            title = a$up$ID[2],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T) #显示pvalue等信息
siP7SK_TNF_pathway <- enrichplot::gseaplot2(b$all,
                                            b$up$ID[1],#富集的ID编号
                                            title = b$up$ID[1],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T)

siPUS7_E2F_pathway <- enrichplot::gseaplot2(a$all,
                                            a$down$ID[1],#富集的ID编号
                                            title = a$down$ID[1],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T) #显示pvalue等信息

si7sk_E2F_pathway <- enrichplot::gseaplot2(b$all,
                                            b$down$ID[1],#富集的ID编号
                                            title = b$down$ID[1],#标题
                                            color = "red", #GSEA线条颜色
                                            base_size = 20,#基础字体大小
                                            rel_heights = c(0.8, 0.1, 1),#副图的相对高度
                                            subplots = 1:2,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                                            ES_geom = "line", #enrichment score用线还是用点"dot"
                                            pvalue_table = T) #显示pvalue等信息




write.table(a$down,file = "GSEA_sipus7.down.gene.txt",quote = F, row.names = F,col.names = T,sep = "\t")
write.table(a$up,file = "GSEA_sipus7.up.gene.txt",quote = F, row.names = F,col.names = T,sep = "\t")
write.table(b$down,file = "GSEA_si7sk.down.gene.txt",quote = F, row.names = F,col.names = T,sep = "\t")
write.table(b$up,file = "GSEA_si7sk.up.gene.txt",quote = F, row.names = F,col.names = T,sep = "\t")
View(b$down)
View(a$down)
View(b$up)
View(a$up)