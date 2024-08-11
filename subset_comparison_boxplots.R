library(plyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
rm(list=ls())
##### AS DC and LC-like across various cancer types
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
metadata_filt <- metadata[metadata$Tissue=='Tumor',]

col_names <- colnames(metadata_filt)
print(col_names)

cell_num <- ddply(metadata_filt, .(sample_tissue,cancer), nrow)
colnames(cell_num) <- c("sample","cancer","Total")

cell_num_cluster <- ddply(metadata_filt, .(sample_tissue,cancer,celltype), nrow)
colnames(cell_num_cluster) <- c("sample","cancer","celltype","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,sample+cancer+Total~celltype,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("sample","cancer","Total"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=10,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)
# Color panel -----------------
c54 <- c('BLCA' = 'dodgerblue2',"BRCA"='green4',"CESC"='#E31A1C',"CRC"='#6A3D9A',"CSCC"='#FF7F00',
         "CTCL"='#FB9A99',"ESCA"='#CAB2D6',"GIST"='khaki2',"LIVER"='deeppink1',"HNSC"='blue1',
         'KIDNEY'='steelblue4','LUNG'='green1','MEL'='yellow4','OS'='forestgreen',
         'OV-FTC'='orange','UCEC'='cornflowerblue', 'PDAC'='magenta','PRAD'='darkolivegreen4',
         'STAD'='indianred1','THCA'='tan4','GCTB'='mediumorchid1','TGCT'='firebrick4',
         'GEJ'='lightsalmon','NB'='tan1','ACP'='darkgray',
         'EPN'='wheat4','GLIOMA'='#DDAD4B','MB'='chartreuse','PCNSL'='seagreen1','PitNET'='moccasin',
         'VS'='mediumvioletred','MEN'='seagreen','UVM'='cadetblue1')
# AS DC
type <- "AS_DC_CD5"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 30, label.x = 2) + ggtitle(paste0("Proportion of ", type))+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )+xlab("")
p1
ggsave("/data2/maty/DC/script/DC_merge/figures/Part1/AS_DC_cancer_type.pdf", plot = p1, width = 10, height = 4)
# LC-like
type <- "LC_CD207"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 70, label.x = 2) + ggtitle(paste0("Proportion of ", type))+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )+xlab("")
p1
ggsave("/data2/maty/DC/script/DC_merge/figures/Part1/LC_cancer_type.pdf", plot = p1, width = 10, height = 4)
##### Subset comparison
# AS DC
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]

metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("CRC","CSCC",'ESCA','HNSC','STAD'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("AS_DC_CD5"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE, width=0.5)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, ncol = 5)
p <- p+stat_compare_means(size=3,label.y=17)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part1/AS_DC_compare.pdf", plot = p, width = 10.5, height = 3.8)
# LAMP3+ DC
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]

metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c('CRC','CSCC','ESCA','LIVER','HNSC','STAD'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("DC_LAMP3"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, scales = "free_y", ncol = 6)
p <- p+stat_compare_means(size=3)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part1/LAMP3_DC_compare.pdf", plot = p, width = 12, height = 4)
# LC-like
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c('CRC','KIDNEY','LUNG','UCEC'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("LC_CD207"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE,width=0.5)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, scales = "free_y", 
                    ncol = 2)
p <- p+stat_compare_means(size=3)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part1/LC_compare.pdf", plot = p, width = 5, height = 7)
#cDC2_CXCL9
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/cDC2_raw_0610.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c('CRC','CSCC','ESCA','HNSC','KIDNEY'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("cDC2_CXCL9"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE,width=0.5)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, scales = "free_y", 
                    ncol = 5)
p <- p+stat_compare_means(size=3)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part2-cDC2/cDC2_CXCL9_compare.pdf", plot = p, width = 10, height = 4)
#cDC2_ISG15
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/cDC2_raw_0610.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c('CSCC','ESCA','HNSC'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("cDC2_ISG15"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE,width=0.5)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, scales = "free_y", 
                    ncol = 4)
p <- p+stat_compare_means(size=3)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part2-cDC2/cDC2_ISG15_compare.pdf", plot = p, width = 6.5, height = 4)
#cDC2_C1QC
metadata <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/cDC2_raw_0610.csv")
metadata_filt_N_T <- metadata[metadata$Tissue %in% c("Tumor","Adjacent normal tissue"),]
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c('CSCC','ESCA','LUNG','MEN'),]
metadata_filt_N_T_cancer$Tissue <- factor(metadata_filt_N_T_cancer$Tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num) <- c("sample_tissue","cancer","Tissue","Total")

metadata_filt_LC <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$celltype %in% c("cDC2_C1QC"),]

cell_num_cluster <- ddply(metadata_filt_LC, .(sample_tissue,cancer,Tissue), nrow)
colnames(cell_num_cluster) <- c("sample_tissue","cancer","Tissue","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)

cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=10,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

compare_means(Proportion~Tissue, data=cell_num_cluster_summary, group.by = "cancer")
p <- ggboxplot(cell_num_cluster_summary, x="Tissue", y="Proportion", color = "Tissue", 
palette = "jco", add = "jitter", facet.by = "cancer", short.panel.labs = FALSE,width=0.5)+ ylab("Proportion")+ xlab("")
p <- p + facet_wrap(~ cancer, scales = "free_y", 
                    ncol = 4)
p <- p+stat_compare_means(size=3)
p <- p + theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
              strip.text = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              axis.title.y = element_text(size = 10),legend.position = "null")
p
ggsave("/data2/maty/DC/script/DC_merge/figures/Part2-cDC2/cDC2_C1QC_compare.pdf", plot = p, width = 8.5, height = 4)