# monocle 3 ppseudotime analysis
# convert seurat object to dcs object
cds <- as.cell_data_set(eye)
# cluster cells
cds <- cluster_cells(cds, resolution=1e-3)
# learn graph
cds <- learn_graph(cds) 
# label roots and branch points here. chose AUnd as the root  

plot_cells(cds, color_cells_by = "cluster")

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 'AUnd']))  
# pseudotime plot
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph = T) 


# monocle 3 ppseudotime analysis of MF and PRs
# convert seurat object to dcs object
cds_prs <- as.cell_data_set(eye)
# cluster cells
cds_prs <- cluster_cells(cds_prs, resolution=1e-3)
# learn graph
cds_prs <- learn_graph(cds_prs) 
# label roots and branch points here. chose AUnd as the root  

plot_cells(cds_prs, color_cells_by = "cluster")

cds_prs <- order_cells(cds_prs, reduction_method = "UMAP", root_cells = colnames(cds_prs[, clusters(cds_prs) == 'MF']))  
# pseudotime plot
plot_cells(cds_prs, color_cells_by = "pseudotime",
           show_trajectory_graph = T) 
# snATACseq pseudotime

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)

DefaultAssay(eye_ATAC) <- "peaks"

DimPlot(eye_ATAC, label = T)

cds.ATAC <- as.cell_data_set(eye_ATAC)
cds.ATAC <- cluster_cells(cds = cds.ATAC, reduction_method = "UMAP")
cds.ATAC <- learn_graph(cds.ATAC, use_partition = TRUE)

cds.ATAC <- order_cells(cds.ATAC, reduction_method = "UMAP")

plot_cells(
  cds = cds.ATAC,
  color_cells_by = "pseudotime",
  show_trajectory_graph = T, cell_size = 1,trajectory_graph_segment_size = 2,trajectory_graph_color = 'black'
)
