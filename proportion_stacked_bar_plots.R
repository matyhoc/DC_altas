# plot cell numbers
library(patchwork)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(scales)
##### DC_proportion_by_cancer
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")

celltype_colors <- c("pDC_LILRA4"="#aec7e8",
"LC_CD207"="#e377c2",
"DC_LAMP3"="#279e68",
"cDC1_CLEC9A"="#b5bd61",
"AS_DC_CD5"="#1f77b4",
"DC_RUNX3"="#aa40fc",
"DC_GTF2IRD1"="#8c564b",
"DC_TRIM33"="#ff7f0e",
"DC_MKI67"="#d62728",
"cDC2_CD1C"="#17becf")
metadata_cancer <- metadata[metadata$Tissue=='Tumor',]
plotC <- table(metadata_cancer$cancer, metadata_cancer$celltype) %>% melt()
colnames(plotC) <- c("Sample", "CellType","Number")
plotC <- plotC %>%
  group_by(Sample) %>%
  mutate(Total = sum(Number),
         cDC2_CD1C_Proportion = Number[CellType == "cDC2_CD1C"] / Total) %>%
  ungroup() %>%
  arrange(desc(cDC2_CD1C_Proportion)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
PC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=celltype_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))
ggsave(PC2,file="DC_proportion_by_cancer.pdf",width = 12, height = 5)

##### cDC2_proportion_by_cancer
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/cDC2_raw_0610.csv")
celltype_colors <- c('cDC2_CD1C'='#9DC1D2',
'cDC2_C1QC'='#72BCC8',
'cDC2_CCR7'='#279e68',
'cDC2_FCN1'='#5E9A93',
'cDC2_HSP'='#aa40fc',
'cDC2_ISG15'='#8c564b',
'cDC2_CXCL9'='#e377c2'
)
metadata_cancer <- metadata[metadata$Tissue=='Tumor',]
plotC <- table(metadata_cancer$cancer, metadata_cancer$celltype) %>% melt()
colnames(plotC) <- c("Sample", "CellType","Number")
celltype_order <- c("cDC2_CD1C","cDC2_C1QC","cDC2_CCR7","cDC2_FCN1","cDC2_HSP","cDC2_ISG15","cDC2_CXCL9")

plotC <- plotC %>%
  group_by(Sample) %>%
  mutate(Total = sum(Number),
         cDC2_CD1C_Proportion = Number[CellType == "cDC2_CXCL9"] / Total) %>%
  ungroup() %>%
  arrange(desc(cDC2_CD1C_Proportion)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)),
         CellType = factor(CellType, levels = celltype_order))
PC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.8, aes(group = CellType), position = "fill") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "", y = "Cell proportion") +
  scale_y_continuous(labels = percent) +
  theme(axis.text = element_text(size = 12, colour = "black")) +
  theme(axis.title.y = element_text(size = 12, colour = "black")) +
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.6)) +
  theme(legend.text = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15))
ggsave(PC2,file="cDC2_proportion_by_cancer.pdf",width = 12, height = 5)
##### Tissue_type_by_cancer
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
tissue_type_colors <- c(
'Patient Metastatic LN'='#b5bd61',
'PBMC_Patient'='#d62728',
'Patient Normal LN'='#aec7e8',
'Tumor thrombus'='#17becf',
'PBMC_Healthy'='#e377c2',
'Tumor'='#1f77b4',
'Metastatic tissue'='#8c564b',
'Non-malignant lesion'='#aa40fc',
'Adjacent normal tissue'='#ff7f0e',
'Healthy tissue'='#279e68')
plotC <- table(metadata$cancer, metadata$Tissue) %>% melt()
colnames(plotC) <- c("Sample", "Tissue_type","Number")
plotC <- plotC %>%
  group_by(Sample) %>%
  mutate(Total = sum(Number),
         Tumor_Proportion = Number[Tissue_type == "Tumor"] / Total) %>%
  ungroup() %>%
  arrange(desc(Tumor_Proportion)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
PC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = Tissue_type)) +
  geom_bar(stat = "identity", width=0.8,aes(group=Tissue_type),position="fill")+  
  scale_fill_manual(values=tissue_type_colors) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))
ggsave(PC2,file="DC_proportion_by_tissue.png",width = 12, height = 5)
##### Immune_cell_proportion_by_dataset
DC <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
CD4T <- read.csv("/data2/maty/DC/script/Treg/file/CD4_all_new_index_tissue_0617.csv")
CD8T <- read.csv("/data2/maty/DC/script/Treg/file/CD8_all_new_index_tissue_0617.csv")
B <- read.csv("/data2/maty/DC/script/Treg/file/B_0619.csv")
Mye <- read.csv("/data2/maty/DC/script/Treg/file/Myeloid_cell_without_DC_0619.csv")
NK <- read.csv("/data2/maty/DC/script/Treg/file/NK_0618.csv")

DC['Celltype'] = 'Dendritic cells'
CD4T['Celltype'] = 'T cells'
CD8T['Celltype'] = 'T cells'
B['Celltype'] = 'B/Plasma cells'
Mye['Celltype'] = 'Other myeloid cells'
NK['Celltype'] = 'NK/ILC'

DC <- subset(DC, select = c(X, dataset_id, cancer, Celltype))
CD4T <- subset(CD4T, select = c(X, dataset_id, cancer, Celltype))
CD8T <- subset(CD8T, select = c(X, dataset_id, cancer, Celltype))
B <- subset(B, select = c(X, dataset_id, cancer, Celltype))
Mye <- subset(Mye, select = c(X, dataset_id, cancer, Celltype))
NK <- subset(NK, select = c(X, dataset_id, cancer, Celltype))

combined_df <- bind_rows(DC, CD4T, CD8T, B, Mye, NK)
combined_df <- combined_df[combined_df$cancer!='PCNSL',]
combined_df <- subset(combined_df, select = c(X, dataset_id, Celltype))

celltype_colors <- c(
'B/Plasma cells'="#279e68",
'NK/ILC'="#1f77b4",
'Other myeloid cells'="#ff7f0e",
"Dendritic cells"="#d62728",
'T cells'="#b5bd61")
plotC <- table(combined_df$dataset_id, combined_df$Celltype) %>% melt()
colnames(plotC) <- c("Sample", "CellType","Number")
desired_order <- c('T cells','B/Plasma cells','NK/ILC','Other myeloid cells','Dendritic cells')
plotC$CellType <- factor(plotC$CellType, levels = desired_order)
plotC <- plotC %>%
  group_by(Sample) %>%
  mutate(Total = sum(Number),
         DC_Proportion = Number[CellType == "Dendritic cells"] / Total) %>%
  ungroup() %>%
  arrange(desc(DC_Proportion)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
PC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=celltype_colors,breaks = c('T cells','B/Plasma cells','NK/ILC','Other myeloid cells','Dendritic cells'))+
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent,breaks = seq(0, 1, 0.1))+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.8, vjust = 0.6))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()
ggsave(PC2,file="DC_proportion_by_dataset_swap_x_y.pdf",width = 8, height = 20)
##### Immune_cell_proportion_by_cancer_type
DC <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
CD4T <- read.csv("/data2/maty/DC/script/Treg/file/CD4_all_new_index_tissue_0617.csv")
CD8T <- read.csv("/data2/maty/DC/script/Treg/file/CD8_all_new_index_tissue_0617.csv")
B <- read.csv("/data2/maty/DC/script/Treg/file/B_0619.csv")
Mye <- read.csv("/data2/maty/DC/script/Treg/file/Myeloid_cell_without_DC_0619.csv")
NK <- read.csv("/data2/maty/DC/script/Treg/file/NK_0618.csv")

DC['Celltype'] = 'Dendritic cells'
CD4T['Celltype'] = 'T cells'
CD8T['Celltype'] = 'T cells'
B['Celltype'] = 'B/Plasma cells'
Mye['Celltype'] = 'Other myeloid cells'
NK['Celltype'] = 'NK/ILC'

DC <- subset(DC, select = c(X, cancer, Tissue,  Celltype))
CD4T <- subset(CD4T, select = c(X, cancer, Tissue,  Celltype))
CD8T <- subset(CD8T, select = c(X, cancer, Tissue,  Celltype))
B <- subset(B, select = c(X, cancer, Tissue,  Celltype))
Mye <- subset(Mye, select = c(X, cancer, Tissue,  Celltype))
NK <- subset(NK, select = c(X, cancer, Tissue,  Celltype))

combined_df <- bind_rows(DC, CD4T, CD8T, B, Mye, NK)

celltype_colors <- c(
'B/Plasma cells'="#279e68",
'NK/ILC'="#1f77b4",
'Other myeloid cells'="#ff7f0e",
"Dendritic cells"="#d62728",
'T cells'="#b5bd61")

combined_df <- combined_df[combined_df$Tissue=='Tumor',]
combined_df <- combined_df[combined_df$cancer!='PCNSL',]

plotCC <- table(combined_df$cancer, combined_df$Celltype) %>% melt()
colnames(plotCC) <- c("Sample", "CellType","Number")

desired_order <- c('T cells','B/Plasma cells','NK/ILC','Other myeloid cells','Dendritic cells')
plotCC$CellType <- factor(plotCC$CellType, levels = desired_order)

plotCC <- plotCC %>%
  group_by(Sample) %>%
  mutate(Total = sum(Number),
         DC_Proportion = Number[CellType == "Dendritic cells"] / Total) %>%
  ungroup() %>%
  arrange(desc(DC_Proportion)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
PC2 <- ggplot(data = plotCC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=celltype_colors,breaks = c('T cells','B/Plasma cells','NK/ILC','Other myeloid cells','Dendritic cells')) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent,breaks = seq(0, 1, 0.1))+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.title = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 10))
ggsave(PC2,file="DC_proportion_by_CANCER_only_tumor.pdf",width = 20, height = 5)