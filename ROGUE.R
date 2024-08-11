# Calculate ROGUE index
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
data <- readRDS('/data2/maty/DC/script/DC_merge/R/DC_raw_celltype_tissue.rds')
sample_counts <- table(data@meta.data$sample_tissue)
selected_samples <- names(sample_counts[sample_counts >= 100])
seurat_obj <- subset(data, subset = sample_tissue %in% selected_samples)
data <- seurat_obj
data[['RNA']]@counts <- as.matrix(data[['RNA']]@counts)
ent.res <- SE_fun(data[['RNA']]@counts)
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.res <- rogue(data[['RNA']]@counts, labels = data@meta.data$celltype, samples = data@meta.data$sample_tissue, platform = "UMI", span = 0.9)
write.csv(rogue.res,'rogue.res.csv')
# Visulization
library(ggpubr)
library(patchwork)
library(ggsci)
library(reshape2)
library(plyr)
rogue.res <- read.csv("/data2/maty/DC/script/DC_merge/R/rogue.res.csv")
rogue.res_melted_data <- melt(rogue.res, id.vars = "X")
df = rogue.res_melted_data
df <- na.omit(df)
c54 <- c('cDC2_CD1C'='#b5bd61',
'cDC1_CLEC9A'='#d62728',
'DC_MKI67'='#aec7e8',
'DC_TRIM33'='#17becf',
'DC_GTF2IRD1'='#e377c2',
'pDC_LILRA4'='#1f77b4',
'DC_RUNX3'='#8c564b',
'AS_DC_CD5'='#aa40fc',
'LC_CD207'='#ff7f0e',
'DC_LAMP3'='#279e68')
median_table <- ddply(df,.(variable), function(x){median(x$value)})
median_table_sorted <- median_table[order(median_table$V1),]
df$variable <- factor(df$variable, levels=median_table_sorted$variable)
p1 <- ggboxplot(df, x = "variable", y = "value",
                color = "variable", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 1, label.x = 1) + ggtitle(paste0(''))+ theme(
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
      )+xlab("")+ ylab("ROGUE index")
p1
ggsave("/data2/maty/DC/script/DC_merge/figures/Part2-cDC2/ROGUE.pdf", plot = p1, width = 8, height = 5.5)