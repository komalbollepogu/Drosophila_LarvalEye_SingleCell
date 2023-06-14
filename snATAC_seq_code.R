# snATAC-seq data analysis code
# load packages
library(Cicero)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(IRanges)
library(hdf5r)
library(ensembldb)
library(GenomicRanges)
library(GenomicFeatures)
library(Biobase)
library(rtracklayer)

# read the counts into R
counts <- Read10X_h5(filename = "E:/outs/filtered_peak_bc_matrix.h5")



# read the metadata into R
metadata <- read.csv(  file = "E:/outs/singlecell.csv",
                       header = TRUE,
                       row.names = 1)

# create chromatin assay

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'dm6',
  fragments = 'E:/outs/fragments.tsv.gz'
  
)

# create suerat object using chromatin assay

larval_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# import GTF file into R

gtf <- rtracklayer::import('Drosophila_melanogaster.BDGP6.28.99.chr.gtf.gz')

# gene coordinates

gene.coords <- gtf[gtf$type == c('gene','exon')]

# change gene coordinates to ucsc format

seqlevelsStyle(gene.coords) <- 'UCSC'

gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

genome(gene.coords) <- "dm6"

# add the genomic coordinates annotation to the seurat object

Annotation(larval_ATAC) <- gene.coords

larval_ATAC[['peaks']]

# computing QC metircs
#  nucleosome signal score per cell

larval_ATAC <- NucleosomeSignal(object = larval_ATAC)
# TSS enrichment score per cell

larval_ATAC <- TSSEnrichment(object = larval_ATAC, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks. dm6 blacklist

larval_ATAC$pct_reads_in_peaks <- larval_ATAC$peak_region_fragments / larval_ATAC$passed_filters * 100
larval_ATAC$blacklist_ratio <- larval_ATAC$blacklist_region_fragments / larval_ATAC$peak_region_fragments

larval_ATAC$high.tss <- ifelse(larval_ATAC$TSS.enrichment > 2, 'High', 'Low')

TSSPlot(larval_ATAC, group.by = 'high.tss') + NoLegend()

larval_ATAC$nucleosome_group <- ifelse(larval_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

FragmentHistogram(object = larval_ATAC, group.by = 'nucleosome_group')

# plot the distribution of each QC metric 

VlnPlot(
  object = larval_ATAC,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

#remove cells that are outliers for these QC metrics

P0 <- subset(
  x = larval_ATAC,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

#Normalization and linear dimensional reduction
P0 <- RunTFIDF(P0)
P0 <- FindTopFeatures(P0, min.cutoff = 'q0')
P0 <- RunSVD(P0)

# correlation between depth and reduced dimension components
DepthCor(P0)

#Non-linear dimension reduction and clustering

P0 <- RunUMAP(object = P0, reduction = 'lsi', dims = 2:50)
P0 <- FindNeighbors(object = P0, reduction = 'lsi', dims = 2:50)
P0 <- FindClusters(object = P0, verbose = FALSE, algorithm = 3, resolution = 4)
DimPlot(object = P0, label = TRUE, pt.size = 1.5) + NoLegend()

# integration with scRNA seq


# quantify gene activity. P0 is snATAC and eye is scRNA
gene.activities <- GeneActivity(P0, features = VariableFeatures(eye))

# add gene activities as a new assay
P0[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(P0) <- "ACTIVITY"
P0 <- NormalizeData(P0)
P0 <- ScaleData(P0, features = rownames(P0))

# identify anchors between the two datasets

transfer.anchors <- FindTransferAnchors(reference = eye, query = P0, features = VariableFeatures(object = eye),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

# Annotate scATAC-seq cells using label transfer

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = eye$seurat_clusters,
                                     weight.reduction = P0[["lsi"]], dims = 2:50)

# add cell type predictons as metadata to snATAC-seq

P0 <- AddMetaData(P0, metadata = celltype.predictions)

#ground-truth annotation used for evaluation.
P0$annotation_correct <- P0$predicted.id == P0$seurat_clusters

p1 <- DimPlot(P0, group.by = "predicted.id", label = TRUE, pt.size = 2) + NoLegend() + ggtitle("Predicted annotation")

p2 <- DimPlot(P0, group.by = "seurat_clusters", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2 # annotations are extremely similar


predictions <- table(P0$seurat_clusters, P0$predicted.id)

predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

library(cowplot)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.75, hjust = 1))

correct <- length(which(P0$seurat_clusters == P0$predicted.id))

incorrect <- length(which(P0$seurat_clusters != P0$predicted.id))

data <- FetchData(P0, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                                                                                  labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2

#Co-embedding scRNA-seq and scATAC-seq datasets

genes.use <- VariableFeatures(eye)

refdata <- GetAssayData(eye, assay = "RNA", slot = "data")[genes.use, ]

# refdata has a scRNA-seq expression matrix for the scRNA-seq cells.  
# The output has an imputed scRNA-seq matrix for each of the ATAC cells

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = P0[["lsi"]],
                           dims = 2:50)
P0[["RNA"]] <- imputation

coembed <- merge(x = eye, y = P0)

# run PCA and UMAP on this object, to visualize the co-embedding of scRNA and snATAC

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)

coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)

coembed <- RunUMAP(coembed, dims = 1:50)

DimPlot(coembed, group.by = c("orig.ident", "seurat_clusters"))

# imputation of RNA on snATAC clusters

# predict gene expression values

rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(eye, assay = "RNA", slot = "data"),
  weight.reduction = P0[["lsi"]],
  dims = 2:50
)

# add predicted values as a new assay
P0[["predicted"]] <- rna

# subset only eye. remove brain10, glia17, optic stalk24

eye_ATAC <- subset(P0, idents = c(10,17,24,))

# run non-linear dimension reduction
eye_ATAC <- RunUMAP(object = eye_ATAC, reduction = 'lsi', dims = 2:50)
eye_ATAC <- FindNeighbors(object = eye_ATAC, reduction = 'lsi', dims = 2:50)
eye_ATAC <- FindClusters(object = eye_ATAC, verbose = FALSE, algorithm = 3)
DimPlot(object = eye_ATAC, label = TRUE, pt.size = 1.5) + NoLegend()



eye_ATAC <- RenameIdents(
  
  object = eye_ATAC,
  
  '0' = 'PUnd',
  
  '1' = 'PUnd',
  
  '2' = 'PUnd','3' = 'SMW',
  '4'= 'PUnd',
  '5'= 'AUnd',
  '6' = 'PUnd',
  '7'= 'AUnd',
  '8'='MF+PPN',
  '9' = 'Cones', '10' = 'MF+PPN',
  
  '11' = 'AUnd','12' = 'AUnd', '13'='Cones','14'='SMW','15'='R7',
  '16'='AUnd','17'='PUnd','18'='PUnd','19'='SMW','20'='Convergence','21'='PUnd','22'='R2/5','23'='R1/6',
  '24'='Convergence','25'='R3/4','26'='R7','28'='Convergence','29'='R8','30'='SMW')



# to make coverage plots fig 6E
CoveragePlot(
  object = eye_ATAC,
  region = "dac",  extend.upstream = 5000,
  extend.downstream =5000
)

# figure 6F
CoveragePlot(
  object = eye_ATAC,
  region = "CAP",
  extend.upstream = 3000,
  extend.downstream =3000, annotation = T)

# to make coverage plots fig 7 and for other genes, only the gene name was changed in the command below

CoveragePlot(
  object = eye_ATAC,
  region = c("sens",'Wnt2'),  extend.upstream = 5000,
  extend.downstream =5000
)





saveRDS(eye_ATAC, file = "eye_ATAC.rds")


saveRDS(eye, file = "eye_RNA.rds")


# for motif analyses

eye_ATAC <- readRDS("/eye_ATAC.rds")

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)
# Get a list of motif position frequency matrices for drosophila from the JASPAR database

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

library(BSgenome.Dmelanogaster.UCSC.dm6)

# add motif information

eye_ATAC <- AddMotifs(
  object = eye_ATAC,
  genome = BSgenome.Dmelanogaster.UCSC.dm6,
  pfm = pfm
)

# perform differential accessibility analyses
Cones_peaks <- FindMarkers(
  object = eye_ATAC,
  ident.1 = 'Cones',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.COnes.peak <- rownames(Cones_peaks[Cones_peaks$p_val < 0.005, ])


# find enriched motifs
cones.motifs <- FindMotifs(
  object = eye_ATAC,
  features = top.Cones.peak
)

#  visualize the different motif sequences
MotifPlot(
  object = eye_ATAC,
  motifs = head(rownames(cones.motifs))
)






