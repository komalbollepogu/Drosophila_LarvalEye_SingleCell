library(Seurat)
library(ggplot2)
library(sctransform)
library(ggplot2)
library(scales)
library(monocle3)
library(SeuratWrappers)


# Load the data in R. two biological repeats larval1 and larval2
larval1 <- Read10X(data.dir = "E:/male1/outs/filtered_feature_bc_matrix")


wpp <- CreateSeuratObject(counts = larval1, project = "wpp")


larval2 <- Read10X(data.dir = "E:/male2/outs/filtered_feature_bc_matrix")


wpp2 <- CreateSeuratObject(counts = larval2, project = "wpp2")

#merge both biological repeats

wpp_merge <- merge(larval1, y = larval2, add.cell.ids = c("wpp1","wpp2"), project = "combined_wpp")


# store mitochondrial percentage in object meta data. Perform QC using the plots below

wpp_merge <- PercentageFeatureSet(wpp_merge, pattern = "mt:", col.name = "percent.mt")


VlnPlot(wpp_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(wpp_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(wpp_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


CombinePlots(plots = list(plot1))

CombinePlots(plots = list(plot2))

plot3 <- plot1 + ylim(0,20)
CombinePlots(plots = list(plot1,plot2))


#now subset based on mitochindrial gene expression and number of features. filter cells that have unique feature counts over 5500 or less than 200

wpp_subset <- subset(wpp_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 45)

wpp_subset

# run sctransform wheich performes normalization and scaling. regress out percent mito

wpp_subset <- SCTransform(wpp_subset,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE)

# Perform linear dimensional reduction using PCA
wpp_subset <- RunPCA(wpp_subset, verbose = TRUE)

# run non-linear dimensional reduction umap with the top 50 dimensions and then run clustering analysis with a resolution of 1.
wpp_subset <- RunUMAP(wpp_subset, dims = 1:50, verbose = FALSE)

wpp_subset <- FindNeighbors(wpp_subset, dims = 1:50, verbose = TRUE)

wpp_subset <- FindClusters(wpp_subset, verbose = TRUE, resolution=1)

# dimplot gives a umap with clusters. 
# split the clusters as biological repeats
DimPlot(wpp_subset, label = TRUE, split.by = "orig.ident", pt.size = 1.5) 

DimPlot(wpp_subset, label = TRUE, pt.size = 1.5) 


# subset only cells from eye disc. clusters 20,30,26,14,22,23,9 and 31 are from brain, glia and optic lobe
eye <- subset(wpp_10, idents =c(20,30,26,14,22,24,9,31) )
# run sctransform, PCA and UMAP
eye <- SCTransform(eye,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE)

eye <- RunPCA(eye, verbose = TRUE)

eye <- RunUMAP(eye, dims = 1:50, verbose = FALSE)

eye <- FindNeighbors(eye, dims = 1:50, verbose = TRUE)

eye <- FindClusters(eye, verbose = TRUE, resolution = 1)

DimPlot(eye, label = TRUE, pt.size = 1.5) 

# for supplemental figure 1 L, M
DimPlot(eye, label = TRUE, pt.size = 1.5, split.by = 'orig.ident') 

eye <- RenameIdents(
  
  object = eye,
  
  '0' = 'Cones',
  
  '1' = 'PUnd','2' = 'AUnd',
  '3'= 'R8',  '4'= 'R2/5', '5' = 'R3/4',
  '6'= 'R1/6', '7'='R7','8' = 'PPN', '9' = 'MF', '10' = 'Convergence',
  '11' = 'PC', '12'='PPD','13'='Oc','14'='LM')

# colors for figure 1
DimPlot(eye, label = F, pt.size = 1.5,
        cols = c('R3/4' = '#00B4f0','PPN'='#00BA38','MF'='#DE8C00','R8'='#FF64B0',
                 'R2/5'='#F564E3','AUnd'='#00BFC4','R1/6'='#619CFF','R7'='#00C08B',
                 'Cones'='#7CAE00','PUnd'='#F8766D','SMW'='#00C08B','PC'='#DE8C00','PPD'='#B79F00',
                 'Convergence'='#C77CFF','Oc'='#F564E3','LM'='#F564E3')) + NoAxes()+NoLegend() 



# R34 subcluster
R34 <- subset(eye, idents =c('R3/4') )
# run PCA and UMAP
R34 <- RunPCA(R34, verbose = TRUE)

R34 <- RunUMAP(R34, dims = 1:50, verbose = FALSE)

R34 <- FindNeighbors(R34, dims = 1:50, verbose = TRUE)

R34 <- FindClusters(R34, verbose = TRUE, resolution = 1)

DimPlot(R34, label = TRUE, pt.size = 3.5) 



