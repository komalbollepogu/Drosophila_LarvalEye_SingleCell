# testing different random seeds and hyperparameters for UMAP
set.seed(42)
larval <- eye

larval <- SCTransform(larval,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE, seed.use = 42)

larval <- RunPCA(larval, verbose = TRUE, seed.use = 42)

larval <- RunUMAP(larval, dims = 1:50, verbose = FALSE, seed.use = 42)

larval <- FindNeighbors(larval, dims = 1:50, verbose = TRUE)

larval <- FindClusters(larval, verbose = TRUE, resolution = 1)

DimPlot(larval,label = F, pt.size = 1.5)

# seed 123

set.seed(123)
larval <- eye

larval <- SCTransform(larval,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE, seed.use = 123)

larval <- RunPCA(larval, verbose = TRUE, seed.use = 123)

larval <- RunUMAP(larval, dims = 1:50, verbose = FALSE, seed.use = 123)

larval <- FindNeighbors(larval, dims = 1:50, verbose = TRUE)

larval <- FindClusters(larval, verbose = TRUE, resolution = 1)

DimPlot(larval,label = F, pt.size = 1.5)

# seed 1000. 
set.seed(1000)
larval <- eye

larval <- SCTransform(larval,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE, seed.use = 1000)

larval <- RunPCA(larval, verbose = TRUE, seed.use = 1000)

larval <- RunUMAP(larval, dims = 1:50, verbose = FALSE, seed.use = 1000)

larval <- FindNeighbors(larval, dims = 1:50, verbose = TRUE)

larval <- FindClusters(larval, verbose = TRUE, resolution = 1)

DimPlot(larval,label = F, pt.size = 1.5)

#seed 2000
set.seed(2000)
larval <- eye

larval <- SCTransform(larval,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE, seed.use = 2000)

larval <- RunPCA(larval, verbose = TRUE, seed.use = 2000)

larval <- RunUMAP(larval, dims = 1:50, verbose = FALSE, seed.use = 2000)

larval <- FindNeighbors(larval, dims = 1:50, verbose = TRUE)

larval <- FindClusters(larval, verbose = TRUE, resolution = 1)

DimPlot(larval,label = F, pt.size = 1.5)

# testing different dimensions

set.seed(42)
larval <- eye

larval <- SCTransform(larval,variable.features.n = 5000, vars.to.regress =c("percent.mt") ,verbose = TRUE, seed.use = 42)

larval <- RunPCA(larval, verbose = TRUE, seed.use = 42)

# ran this code with different dimension each time. Tested 10,20,30,40 and ,50 dimensions
larval <- RunUMAP(larval, dims = 1:50, verbose = FALSE, seed.use = 42)

larval <- FindNeighbors(larval, dims = 1:50, verbose = TRUE)

larval <- FindClusters(larval, verbose = TRUE, resolution = 1)

DimPlot(larval,label = F, pt.size = 1.5)
