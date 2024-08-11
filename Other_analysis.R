#Some codes for dc atlas analyses.

#Xiaojing Chu (tsuki1128@outlook.com)

#2024-08-09

################################################
# 1. caluculate signature scores for lamp3+ dcs.
# 2. infer cellular origins of lamp3+ dcs by singleR.
# 3. show distribution of lamp3+ dc cellular origins across cancer types.
# 4. identify marker genes for lamp3+ dcs with different cellular origins
# 5. calculate enriched pathways for the marker genes.
# 6. visualization of IL15, CCL19 and CXCL9 in lamp3+ dcs of different origins.
# 7. Ternary plot showing cDC1, cDC2, LC-like signature scores for cDC1, cDC2, LC-like and LAMP3+ DCs.
# 8. Correlation analyses between cDC1-derived LAMP3+ DC, cDC2-derived LAMP3+ DC and LC-like-derived LAMP3+ DC proportions and CD8 T proportions.
# 9. monocle2 for lamp3+ dcs, cDC1, cDC2 and LC-like dcs.
# 10. Scenic for AS DC.
# 11. gsva for cdc2 subsets.
# 12. singleR for lc-like DC subset.
################################################


library(Seurat)
library(dplyr)
library(ggrepel)

options(ggrepel.max.overlaps = Inf)

dc<-readRDS("DC_raw_0611.rds")

dc<-NormalizeData(dc)

t<-dc

##### get marker genes of cDC1_CLEC9A, cDC2_CD1C and LC_CD207.

markers<-read.table("DC_celltype_markers.txt",header=T,row.names = 1,sep="\t")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

genes<-list(
  `cDC1_CLEC9A`=top10$gene[which(top10$cluster=="cDC1_CLEC9A")],
  `cDC2_CD1C`=top10$gene[which(top10$cluster=="cDC2_CD1C")],
  `LC_CD207`=top10$gene[which(top10$cluster=="LC_CD207")]
)

##### caluculate signature scores for lamp3+ dc.

t<-AddModuleScore(t,features = genes)

dat<-FetchData(subset(t,celltype %in% c("DC_LAMP3","cDC1_CLEC9A","cDC2_CD1C","LC_CD207")),vars = c("Cluster1","Cluster2","Cluster3","celltype"))

##### infer cellular origins of lamp3+ dc by singleR.

ref<-as.SingleCellExperiment(subset(dc,celltype %in% c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207") & Tissue %in% "Tumor"))

Query<-subset(dc,celltype=="DC_LAMP3" & Tissue %in% "Tumor")

pred.scRNA <- SingleR(test = GetAssayData(Query, assay = "RNA"), ref = ref,labels = ref$celltype, fine.tune = F,de.method = "wilcox",genes=genes,BPPARAM = BiocParallel:::MulticoreParam(workers=2))

col<-c("#1d78b8","#fb7d0c","#2d9a67","#ae3ef6","#8c5650",
       "#b6bc6a","#18bfd7","#acc7ed")

res1<-table(pred.scRNA$pruned.labels[which(apply(pred.scRNA$scores,1,max)>0.2)])
res1<-res1/sum(res1)
res2<-table(pred.scRNA$pruned.labels[which(apply(pred.scRNA$scores,1,max)>0.3)])
res2<-res2/sum(res2)
res3<-table(pred.scRNA$pruned.labels[which(apply(pred.scRNA$scores,1,max)>0.4)])
res3<-res3/sum(res3)

results<-as.data.frame(c(res1,res2,res3))

results<-transform(results,
                   CellType=c(names(res1),names(res2),names(res3)),
                   Threshold=c(rep("Score> 0.2",length(res1)),rep("Score> 0.3",length(res2)),rep("Score> 0.4",length(res3))))
names(results)<-c("Proportion","CellType","Threshold")

ggplot(results,aes(x=Threshold,y=Proportion,fill=CellType))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = col)+theme_classic()

##### use 0.3 as final threshold for downstream analyses.

res2<-pred.scRNA[which(apply(pred.scRNA$scores,1,max)>0.3),]
Query<-AddMetaData(Query,metadata = as.data.frame(res2))

get<-FetchData(Query,vars = c("labels","cancer"))


##### remove cancer types with less than 100 cells.

flag<-names(which(table(get$cancer)>100))

get<-get[which(get$cancer %in% flag),]

tab<-table(get$cancer,get$labels)
tab<-tab/apply(tab,1,sum)
tab<-na.omit(tab)

tab<-as.data.frame(tab)

tab<-tab[order(tab$Freq),]

tab$Var1<-factor(tab$Var1,levels = c(tab$Var1[tab$Var2=="cDC1_CLEC9A"]))

##### show distribution of lamp3+ dc cellular origins across cancer types.

ggplot(tab,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = col)+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))


##### identify marker genes for lamp3+ dc with different cellular origins

tmp<-subset(Query,labels %in% c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207"))

Idents(tmp)<-"labels"

markers <- FindAllMarkers(tmp, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)

markers <- markers[which(markers$p_val_adj<0.05),]

##### calculate enriched pathways for the marker genes.

library(clusterProfiler)
library(org.Hs.eg.db)

for (i in c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207")){
  
  geneList<-markers$gene[which(markers$cluster==i)]
  
  geneID_list <- bitr(geneList,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = 'org.Hs.eg.db')
  
  ego <- enrichGO(
    gene = geneID_list$ENTREZID,
    # universe = allgene
    OrgDb         = 'org.Hs.eg.db',
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable    = TRUE
  )
  
  write.table(ego,file=paste0("ego_LAMP3_",i,".txt"),quote=F,sep="\t")
  
}

##### visualization of IL15, CCL19 and CXCL9 in lamp3+ dc of different origins.

test<-FetchData(tmp,vars = c("IL15","CCL19","CXCL9","labels"))

my_comp1<-list(c("cDC1_CLEC9A","cDC2_CD1C"),c("cDC1_CLEC9A","LC_CD207"))

v1<-VlnPlot(tmp,features = "IL15",pt.size = 0)+stat_compare_means(method = "wilcox",comparisons = my_comp1)+
  scale_fill_manual(values=c("#b6bc6a","#18bfd7","#e178c5"))+NoLegend()+theme(axis.title.x = element_blank())
v2<-VlnPlot(tmp,features = "CCL19",pt.size = 0)+stat_compare_means(method = "wilcox",comparisons = my_comp1)+
  scale_fill_manual(values=c("#b6bc6a","#18bfd7","#e178c5"))+NoLegend()+theme(axis.title.x = element_blank())
v3<-VlnPlot(tmp,features = "CXCL9",pt.size = 0)+stat_compare_means(method = "wilcox",comparisons = my_comp1)+
  scale_fill_manual(values=c("#b6bc6a","#18bfd7","#e178c5"))+NoLegend()+theme(axis.title.x = element_blank())

v1|v2|v3

##### Ternary plot showing cDC1, cDC2, LC-like signature scores for cDC1, cDC2, LC-like and LAMP3+ DCs.

library("Ternary")

tmp<-apply(dat[,1:3],2,function(w){(w-min(w))/(max(w)-min(w))})

tmp2<-apply(tmp,1,function(w){w/sum(w)})
tmp2<-as.data.frame(t(tmp2))

TernaryPlot(axis.labels = seq(0, 1, by = 0.1),alab = "cDC1 score",blab = "cDC2 score",clab = "LC score",
            grid.minor.lines = 0, grid.lty = "dashed")

# Add points

bg<-tmp2[which(dat$celltype %in% c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207")),]
celltype<-data.frame(celltype=dat[which(dat$celltype %in% c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207")),"celltype"])

col<-rep(NA,length(celltype$celltype))

col[which(celltype$celltype=="cDC1_CLEC9A")]<-"#b6bc6a"
col[which(celltype$celltype=="cDC2_CD1C")]<-"#18bfd7"
col[which(celltype$celltype=="LC_CD207")]<-"#e178c5"

TernaryPoints(bg, col = col, pch = ".")

TernaryPoints(tmp2[which(dat$celltype %in% c("DC_LAMP3")),], col = "#2d9a67", pch = ".")

##### Correlation analyses between cDC1-derived LAMP3+ DC, cDC2-derived LAMP3+ DC and LC-like-derived LAMP3+ DC proportions and CD8 T proportions.

cd4<-read.csv("CD4_all_new_index_tissue_0617.csv",header=T,row.names = 1)
cd8<-read.csv("CD8_all_new_index_tissue_0617.csv",header=T,row.names = 1)

cd4<-cd4[cd4$Tissue%in% "Tumor",]
cd8<-cd8[cd8$Tissue%in% "Tumor",]

cd8$celltype<-"CD8T"

data<-rbind(cd4,cd8)

flag<-names(which(table(data$sample_tissue)>10))

data<-data[data$sample_tissue %in% flag,]

data<-split(data,data$sample_tissue)

x<-lapply(data,function(w){table(w$celltype)[2]/nrow(w)})

x<-as.data.frame(do.call(rbind,x))

names(x)<-"CD8_T"

x<-na.omit(x)

##
dcprp<-FetchData(Query,vars = c("sample_tissue","labels"))

dcprp<-na.omit(dcprp)

flag<-names(which(table(dcprp$sample_tissue)>10))

dcprp<-dcprp[which(dcprp$sample_tissue %in% flag),]

sam<-as.character(unique(dcprp$sample_tissue))

res<-matrix(NA,nrow=length(sam),ncol=1)
res<-as.data.frame(res)
rownames(res)<-sam

for (i in sam){
  if(length(table(dcprp$labels[which(dcprp$sample_tissue==i)]))!=3){res[i,1]<-NA}
  else{res[i,1]<-table(dcprp$labels[which(dcprp$sample_tissue==i)])[1]/nrow(dcprp[which(dcprp$sample_tissue==i),])}
}

res<-na.omit(res)
names(res)<-"cDC1_derived_LAMP3DC"

ol<-intersect(rownames(res),rownames(x))

plot<-data.frame(CD8T=x[ol,],cDC1_derived_LAMP3DC=res[ol,])

rownames(plot)<-ol

p1=ggplot(data=plot,aes(x=cDC1_derived_LAMP3DC,y=CD8T))+geom_point()+
  geom_smooth(method = "glm",se=F)+theme_classic()
#4.5,4.2

####cDC2_CD1C
res<-matrix(NA,nrow=length(sam),ncol=1)
res<-as.data.frame(res)
rownames(res)<-sam

for (i in sam){
  if(length(table(dcprp$labels[which(dcprp$sample_tissue==i)]))!=3){res[i,1]<-NA}
  else{res[i,1]<-table(dcprp$labels[which(dcprp$sample_tissue==i)])[2]/nrow(dcprp[which(dcprp$sample_tissue==i),])}
}

res<-na.omit(res)
names(res)<-"cDC2_derived_LAMP3DC"

ol<-intersect(rownames(res),rownames(x))

plot<-data.frame(CD8T=x[ol,],cDC2_derived_LAMP3DC=res[ol,])

rownames(plot)<-ol

p2=ggplot(data=plot,aes(x=cDC2_derived_LAMP3DC,y=CD8T))+geom_point()+
  geom_smooth(method = "glm",se=F)+theme_classic()

####LC
res<-matrix(NA,nrow=length(sam),ncol=1)
res<-as.data.frame(res)
rownames(res)<-sam

for (i in sam){
  if(length(table(dcprp$labels[which(dcprp$sample_tissue==i)]))!=3){res[i,1]<-NA}
  else{res[i,1]<-table(dcprp$labels[which(dcprp$sample_tissue==i)])[3]/nrow(dcprp[which(dcprp$sample_tissue==i),])}
}

res<-na.omit(res)
names(res)<-"LC_derived_LAMP3DC"

ol<-intersect(rownames(res),rownames(x))

plot<-data.frame(CD8T=x[ol,],LC_derived_LAMP3DC=res[ol,])

rownames(plot)<-ol

p3=ggplot(data=plot,aes(x=LC_derived_LAMP3DC,y=CD8T))+geom_point()+
  geom_smooth(method = "glm",se=F)+theme_classic()

p3+p4

#####
Sobj<-readRDS("DC_raw_0611.rds")

tmp<-subset(Sobj,Tissue %in% "Tumor")

tmp2<-FetchData(tmp,vars = c("celltype"))

cells<-tmp2 %>% mutate(cells=rownames(tmp2)) %>% 
  group_by(celltype) %>% sample_n(n()/5)

tmp<-subset(tmp, cells = cells$cells)


library(monocle)
expr_matrix <- as(as.matrix(tmp@assays$RNA$counts), 'sparseMatrix')

p_data <- tmp@meta.data

f_data <- data.frame(gene_short_name = row.names(tmp),row.names = row.names(tmp))

pd <- new('AnnotatedDataFrame', data = p_data)

fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix, phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)

cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.25)

markers<-read.table("DC_celltype_markers.txt",header=T,row.names = 1,sep="\t")

markers<-markers[which(markers$cluster %in% c("cDC1_CLEC9A","cDC2_CD1C","LC_CD207","DC_LAMP3")),]

cds<- setOrderingFilter(cds, markers$gene)

cds <- reduceDimension(cds,max_components = 2,reduction_method = 'DDRTree')

cds <- orderCells(cds)

plot_cell_trajectory(cds,color_by="celltype",size=1,show_backbone=TRUE,cell_size = 0.4)+
  scale_color_manual(values = c("#b6bc6a","#18bfd7","#2d9a67","#e178c5"))

###### Scenic analysis

##data preprocessing
library(Seurat)
library(dplyr)
library(patchwork)

Sobj<-readRDS("DC_raw_0611.rds")

Sobj<-subset(Sobj,Tissue %in% "Tumor")

set.seed(1)

tmp<-FetchData(Sobj,vars = c("celltype"))

cells1<-tmp %>% mutate(cells=rownames(tmp)) %>% 
  group_by(celltype) %>% filter(n()>1000) %>% sample_n(1000)

cells2<-tmp %>% mutate(cells=rownames(tmp)) %>%
  group_by(celltype) %>% filter(n()<=1000 & n()>=300)

cells<-rbind(cells1,cells2)

Sobj<-subset(Sobj, cells = cells$cells)

sceasy::convertFormat(Sobj,from="seurat",to="sce",outFile="dc_ace_1w.rds")

##scenic

library(SingleCellExperiment)
library(SCENIC)

ace<-readRDS("dc_ace_1w.rds")

exprMat <- counts(ace)

cellInfo <- colData(ace)

cellInfo<-data.frame(row.names=rownames(cellInfo),CellType=cellInfo$celltype)

saveRDS(cellInfo, file="int/cellInfo.Rds")

exprMat<-as.matrix(exprMat)

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

scenicOptions <- initializeScenic(org="hgnc", dbDir="scenic", nCores=10)

scenicOptions@inputDatasetInfo$cellInfo<-"int/cellInfo.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions)

exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

##### cdc2 GSVA

library(GSEABase)
library(GSVA)

raw<-readRDS("cDC2_raw_0610.rds")

raw<-NormalizeData(raw)

c2gmt <- getGmt("/data1/chuxj/ref/GOBP_Hs.gmt")

mat<-raw@assays$RNA@data

res.ssgsea <- gsva(as.matrix(mat), c2gmt, method = "ssgsea", kcdf = "Gaussian", min.sz = 10)

res.ssgsea<-as.data.frame(t(res.ssgsea))

group<-FetchData(raw,vars = "celltype")

ttest<-matrix(NA,nrow=ncol(res.ssgsea),ncol=4)
ttest<-as.data.frame(ttest)
rownames(ttest)<-names(res.ssgsea)
names(ttest)<-c("pvalue","mean_other","mean_cxcl9","tvalue")

ttest2<-matrix(NA,nrow=ncol(res.ssgsea),ncol=4)
ttest2<-as.data.frame(ttest2)
rownames(ttest2)<-names(res.ssgsea)
names(ttest2)<-c("pvalue","mean_other","mean_c1qc","tvalue")

ttest3<-matrix(NA,nrow=ncol(res.ssgsea),ncol=4)
ttest3<-as.data.frame(ttest3)
rownames(ttest3)<-names(res.ssgsea)
names(ttest3)<-c("pvalue","mean_other","mean_isg15","tvalue")

for (i in names(res.ssgsea)){
  
  message(paste0("working with ",i))
  
  t<-t.test(res.ssgsea[which(group$celltype!="C07_cDC2_CXCL9"),i],res.ssgsea[which(group$celltype=="C07_cDC2_CXCL9"),i]) 
  ttest[i,1]<-t$p.value
  ttest[i,2]<-t$estimate[1]
  ttest[i,3]<-t$estimate[2]
  ttest[i,4]<-t$statistic
  
  t<-t.test(res.ssgsea[which(group$celltype!="C02_cDC2_C1QC"),i],res.ssgsea[which(group$celltype=="C02_cDC2_C1QC"),i]) 
  ttest2[i,1]<-t$p.value
  ttest2[i,2]<-t$estimate[1]
  ttest2[i,3]<-t$estimate[2]
  ttest2[i,4]<-t$statistic
  
  t<-t.test(res.ssgsea[which(group$celltype!="C06_cDC2_ISG15"),i],res.ssgsea[which(group$celltype=="C06_cDC2_ISG15"),i]) 
  ttest3[i,1]<-t$p.value
  ttest3[i,2]<-t$estimate[1]
  ttest3[i,3]<-t$estimate[2]
  ttest3[i,4]<-t$statistic
}

ttest<-transform(ttest,padj=p.adjust(ttest$pvalue,method = "fdr"))
ttest2<-transform(ttest2,padj=p.adjust(ttest2$pvalue,method = "fdr"))
ttest3<-transform(ttest3,padj=p.adjust(ttest3$pvalue,method = "fdr"))

ttest<-ttest[ttest$padj<0.05,]
ttest2<-ttest2[ttest2$padj<0.05,]
ttest3<-ttest3[ttest3$padj<0.05,]

ttest<-ttest[rev(order(ttest$tvalue)),]
ttest2<-ttest2[rev(order(ttest2$tvalue)),]
ttest3<-ttest3[rev(order(ttest3$tvalue)),]

write.table(ttest,file="gsva_cxcl9_gobp.txt",quote=F,sep="\t")
write.table(ttest2,file="gsva_c1qc_gobp.txt",quote=F,sep="\t")
write.table(ttest3,file="gsva_isg15_gobp.txt",quote=F,sep="\t")

##### singleR for LC-like

library(Seurat)
library(SingleR)
library(tidyverse)
library(ggplot2)

dc<-readRDS("DC_raw_0611.rds")
dc<-NormalizeData(dc)

Idents(dc)<-"celltype"

markers<-read.table("DC_celltype_markers.txt",header=T,row.names = 1,sep="\t")

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

#######
genes<-list(`AS_DC_CD5`=top10$gene[which(top10$cluster=="AS_DC_CD5")],
            `cDC1_CLEC9A`=top10$gene[which(top10$cluster=="cDC1_CLEC9A")],
            `cDC2_CD1C`=top10$gene[which(top10$cluster=="cDC2_CD1C")],
            `DC_LAMP3`=top10$gene[which(top10$cluster=="DC_LAMP3")],
            `DC_PKIB`=top10$gene[which(top10$cluster=="DC_PKIB")],
            `DC_PRDM16`=top10$gene[which(top10$cluster=="DC_PRDM16")],
            `DC_ZEB2`=top10$gene[which(top10$cluster=="DC_ZEB2")],
            `pDC_LILRA4`=top10$gene[which(top10$cluster=="pDC_LILRA4")]
)

ref<-as.SingleCellExperiment(subset(dc,celltype %in% c("AS_DC_CD5","cDC1_CLEC9A","cDC2_CD1C",
                                                       "DC_LAMP3","DC_PKIB",
                                                       "DC_PRDM16","DC_ZEB2","pDC_LILRA4")))
Query<-subset(dc,celltype=="LC_CD207")
pred.scRNA <- SingleR(test = GetAssayData(Query, assay = "RNA"), ref = ref,labels = ref$celltype, fine.tune = F,de.method = "wilcox",genes=genes,BPPARAM = BiocParallel:::MulticoreParam(workers=2))

res<-table(pred.scRNA$pruned.labels[which(apply(pred.scRNA$scores,1,max)>0.3)])
res<-res/sum(res)
