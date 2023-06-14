
library(SCENIC)

library(Seurat)


# converting eye1 into matrix

scenic_eye <- as.matrix(GetAssayData(eye, slot = "counts"))

cellInfo <- data.frame(seuratCluster=Idents(eye))



org="dmel" 

dbDir="cisTarget_databases" # RcisTarget databases location

myDatasetTitle="SCENIC on larval eye disc" # choose a name for your analysis

data(defaultDbNames)

dbs <- defaultDbNames[[org]]

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle) 

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

# seurat cell info

cellInfo <- data.frame(seuratCluster=Idents(eye))

head(cellInfo)

cellInfo$nGene <- colSums(scenic_eye>0)


cellInfo <- data.frame(cellInfo)


cellTypeColumn <- "seuratCluster"

colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"

cbind(table(cellInfo$CellType))

saveRDS(cellInfo, file="int/cellInfo.Rds")

# assign cell colors

colVars <- list(CellType=c("PUnd"="forestgreen","AUnd"="yellow","SMW"="darkorange", "PPD"="magenta4","MF"="hotpink", 
                           "Cones"="red3", 
                           "PPN"="skyblue", 
                           "PC"="darkblue",
                           "Convergence"="moccasin","R16"="darkviolet","R34"="mediumaquamarine",
                           "R7"= "tan4","R8"="orange4","R25"="plum4"))





colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

#convert dgc matrix to matrix

genesKept <- geneFiltering(scenic_eye, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(scenic_eye),
                           minSamples=ncol(scenic_eye)*.01)

interestingGenes <- c("dpp", "h", "emc")

interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- scenic_eye[genesKept, ]

dim(exprMat_filtered)

dim(scenic_eye)

runCorrelation(exprMat_filtered, scenicOptions)


exprMat_filtered <- log2(exprMat_filtered+1)

runGenie3(exprMat_filtered, scenicOptions)

#reloading expr matrix

exprMat <- as.matrix(fr8.data)

logMat <- log2(exprMat+1)

dim(exprMat)


scenicOptions <- readRDS("int/scenicOptions.Rds")

scenicOptions@settings$verbose <- TRUE

scenicOptions@settings$nCores <- 4

scenicOptions@settings$seed <- 123


runSCENIC_1_coexNetwork2modules(scenicOptions)


runSCENIC_2_createRegulons(scenicOptions)

help("Defunct")
library(BiocParallel)

library(AUCell)

SnowParam()

runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
nPcs <- c(5) 
# nPcs <- c(5,15,50)

fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


#par(bg = "black")
par(mfrow=c(1,2))

tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_filtered, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("ato", "onecut")],], plots="Expression")

regulonNames <- c( "ato","onecut")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)

cellCol <- SCENIC::plotEmb_rgb(scenicOptions,regulonNames = regulonNames,aucType = "Binary",aucMaxContrast = 1,offColor = "lightgray",showLegend = F)
#zthis line works for me.

regulonNames <- list(red=c("ato", "onecut"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- SCENIC::plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

library(Cairo)
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()

library(KernSmooth)
library(RColorBrewer)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)


#par(bg = "black")
par(mfrow=c(1,2))

regulonNames <- c( "ato","onecut")
cellCol <- SCENIC::plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)

regulonNames <- list(red=c("ato", "onecut"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- SCENIC::plotEmb_rgb(scenicOptions, regulonNames, aucType="Binary")

regulons <- loadInt(scenicOptions, "regulons")
regulons[c("ato", "onecut")]

regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="onecut" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5)) 

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="onecut"]
viewMotifs(tableSubset) 


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

library(ComplexHeatmap)

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity",column_km = 2)
column_order(ht)

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)


minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))

topRegulators <- reshape2::melt(regulonActivity_byCellType_Binarized)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)

# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

library(plotly)
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)
fig<- plotly::ggplotly(rssPlot$plot)

for #pdf just run rssPlot

plotly::export(p = fig, #the graph to export
               file = "1.pdf") 

htmlwidgets::saveWidget(
  widget = fig, #the plotly object
  file = "figure.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)


plotRSS_oneSet(rss, setName = "Convergence")