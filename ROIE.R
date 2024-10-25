*********************************************************************
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
library(pheatmap)
*********************************************************************
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

# cDC2
df <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/cDC2_raw_0610.csv")
all <- df[(df$Tissue == "Tumor" | df$Tissue == "Adjacent normal tissue"), ]
summary <- table(all[,c('celltype','Tissue')])
roe <- as.data.frame(ROIE(summary))
desired_col_order <- c('Adjacent normal tissue','Tumor')
roe_reordered <- roe[, match(desired_col_order, colnames(roe))]
p1 <- pheatmap(roe_reordered,  #要绘制热图的矩阵
               color = colorRampPalette(c('#B0C4DE',"#e9e9e9","brown3"))(201),
               breaks= c(0,0.0065,0.013,0.0195,0.026,0.0325,0.039,0.0455,0.052,0.0585,0.065,0.0715,0.078,0.0845,0.091,0.0975,0.104,0.1105,0.117,0.1235,0.13,0.1365,0.143,0.1495,0.156,0.1625,0.169,0.1755,0.182,0.1885,0.195,0.2015,0.208,0.2145,0.221,0.2275,0.234,0.2405,0.247,0.2535,0.26,0.2665,0.273,0.2795,0.286,0.2925,0.299,0.3055,0.312,0.3185,0.325,0.3315,0.338,0.3445,0.351,0.3575,0.364,0.3705,0.377,0.3835,0.39,0.3965,0.403,0.4095,0.416,0.4225,0.429,0.4355,0.442,0.4485,0.455,0.4615,0.468,0.4745,0.481,0.4875,0.494,0.5005,0.507,0.5135,0.52,0.5265,0.533,0.5395,0.546,0.5525,0.559,0.5655,0.572,0.5785,0.585,0.5915,0.598,0.6045,0.611,0.6175,0.624,0.6305,0.637,0.6435,0.65,0.6565,0.663,0.6695,0.676,0.6825,0.689,0.6955,0.702,0.7085,0.715,0.7215,0.728,0.7345,0.741,0.7475,0.754,0.7605,0.767,0.7735,0.78,0.7865,0.793,0.7995,0.806,0.8125,0.819,0.8255,0.832,0.8385,0.845,0.8515,0.858,0.8645,0.871,0.8775,0.884,0.8905,0.897,0.9035,0.91,0.9165,0.923,0.9295,0.936,0.9425,0.949,0.9555,0.962,0.9685,0.975,0.9815,0.988,0.9945,1.001,1.0075,1.014,1.0205,1.027,1.0335,1.04,1.0465,1.053,1.0595,1.066,1.0725,1.079,1.0855,1.092,1.0985,1.105,1.1115,1.118,1.1245,1.131,1.1375,1.144,1.1505,1.157,1.1635,1.17,1.1765,1.183,1.1895,1.196,1.2025,1.209,1.2155,1.222,1.2285,1.235,1.2415,1.248,1.2545,1.261,1.2675,1.274,1.2805,1.287,1.2935,1.3, Inf),
               border_color = "black",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               legend = TRUE,
               legend_breaks = c(0, 2),
               legend_labels = c("low","high"),
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 15, 
               fontsize_col = 15,
               display_numbers = T,
               fontsize_number = 10,
               number_color='black',
               cellwidth = 60,cellheight = 40,
               angle_col = 45,
               width = 100,
               height =50
       )
pdf("/data2/maty/DC/script/DC_merge/figures/Part2-cDC2/ROIE_cDC2_T_P.pdf", width = 12, height = 8)
print(p1)
dev.off()

# LC-like excluding cutaneous cancers
df <- read.csv("/data2/maty/DC/script/DC_merge/h5ad/DC_raw_0611.csv")
df <- df[!df$cancer %in% c('CSCC', 'MEL', 'CTCL'), ]
all <- df[(df$Tissue == "Tumor" | df$Tissue == "Adjacent normal tissue" | df$Tissue == "Healthy tissue" | df$Tissue == "Non-malignant lesion" | df$Tissue == "PBMC_Patient" | df$Tissue == "PBMC_Healthy"), ]
summary <- table(all[,c('celltype','Tissue')])
roe <- as.data.frame(ROIE(summary))
roe_LC_CD207 <- subset(roe, rownames(roe) == "LC_CD207")
roe <- roe_LC_CD207
desired_col_order <- c("Healthy PBMC", "Patient PBMC", "Healthy tissue","Non-malignant lesion", 'Adjacent normal tissue','Tumor')
roe_reordered <- roe[, match(desired_col_order, colnames(roe))]
roe_transposed <- t(roe_reordered)
p1 <- pheatmap(roe_transposed,
               color = colorRampPalette(c('white','orange','red'))(101),
               breaks= c(0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.315,0.33,0.345,0.36,0.375,0.39,0.405,0.42,0.435,0.45,0.465,0.48,0.495,0.51,0.525,0.54,0.555,0.57,0.585,0.6,0.615,0.63,0.645,0.66,0.675,0.69,0.705,0.72,0.735,0.75,0.765,0.78,0.795,0.81,0.825,0.84,0.855,0.87,0.885,0.9,0.915,0.93,0.945,0.96,0.975,0.99,1.005,1.02,1.035,1.05,1.065,1.08,1.095,1.11,1.125,1.14,1.155,1.17,1.185,1.2,1.215,1.23,1.245,1.26,1.275,1.29,1.305,1.32,1.335,1.35,1.365,1.38,1.395,1.41,1.425,1.44,1.455,1.47,1.485,1.5,Inf),
               border_color = "black",
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               legend = TRUE,
               legend_breaks = c(0, 1.5),
               legend_labels = c("low","high"),
               show_rownames = TRUE,
               show_colnames = TRUE,
               fontsize_row = 15, 
               fontsize_col = 15,
               display_numbers = T,
               fontsize_number = 10,
               number_color='black',
               cellwidth = 40,cellheight = 40,
               angle_col = 0,
               width = 10,
               height =12
       )
pdf("/data2/maty/DC/script/DC_merge/figures/Part1/ROIE_LC.pdf", width = 12, height = 8)
print(p1)
dev.off()

