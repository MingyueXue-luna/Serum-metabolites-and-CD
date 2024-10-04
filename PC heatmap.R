rm(list = ls())

setwd("C:/Users/74511/Desktop/Projects data/Serum_metabomic/Income new batch/PC_cor/")
data <- read.csv("Olink_Cor_heatmap.csv")
#data <- read.csv("heatmap_Cor.csv")

data <- as.matrix(data)


names <- read.csv("heatmap_rownames_olink.csv")
#names <- read.csv("heatmap_rownames.csv")

colnames(data) <- names[1:10,1]
rownames(data) <- names[11:20,1]



#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)



Heatmap(data)

# change color：
library(circlize)
col_fun <- colorRamp2(c(-0.45, -0.225,0.225, 0.45), c("#5296cc", "#cad5f9", "#fdedf6","#f064af"))

p_data <- read.csv("Olink_q_heatmap.csv")
#p_data <- read.csv("heatmap_FDR.csv")

p_data <- as.matrix(p_data)


# 再画：
Heatmap(data)




Heatmap(data, 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        # 设置颜色：
        col = col_fun,
        column_title = "Spearman correlation coefficient between metabolomics-derived and protiomics-derived PCs",
        # 调整热图格子的边框颜色和粗细：
        rect_gp = gpar(col = "white", lwd = 1),
        border_gp = gpar(col = "black", lwd = 1),
        # 行列名：
        row_names_gp = gpar(fontsize = 10, fontface = "bold.italic"), # 调整字体为斜体：
        column_names_gp = gpar(fontsize = 10,fontface = "bold"),
        column_names_rot = 45,
        row_names_side = "left",
        
        #row_dend_width = unit(2, "cm"),
        # 设置单元格中添加文字：
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", data[i, j]), x, y ,gp = gpar(fontsize = 10,fontface = "bold"))
          if (p_data[i,j] < 0.001) {
            grid.text(sprintf("            ***", data[i, j]), x, y, gp = gpar(fontsize = 10,fontface = "bold"))
          } else if (p_data[i,j]<0.01) {
            grid.text(sprintf("            **", data[i, j]), x, y, gp = gpar(fontsize = 10,fontface = "bold"))
          } else if (p_data[i,j]<0.05) {
            grid.text(sprintf("         *", data[i, j]), x, y, gp = gpar(fontsize = 10,fontface = "bold"))
          } else {
            grid.text(sprintf("", data[i, j]), x, y, gp = gpar(fontsize = 5))
          }
        },
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(col_fun = col_fun, 
                                    at = c(-0.5,-0.25, 0, 0.25, 0.5),
                                    # labels = c("low", "zero", "high"),
                                    title = " ",
                                    legend_height = unit(15.5, "cm"),
                                    border = TRUE,
                                    #title_position = "topcenter",
                                    #title_gp = gpar(fontsize = 5),
                                    #labels_gp = gpar(fontsize = 5),
                                    direction = "vertical"
                                    #grid_height = unit(20, "mm")
        ),
)


T_data <- read.csv("P_olink.csv")
T_data <- as.matrix(T_data)
T_data <- T_data < 0.05
