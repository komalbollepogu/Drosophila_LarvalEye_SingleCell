
# convergence subcluster and PCA

Convergence <- subset(eye, idents =c('Convergence') )

Convergence <- RunPCA(Convergence, verbose = TRUE)

# eigen value for each PC for Convergence
pca1 = Convergence$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained1 = eigValues / sum(eigValues)
write.csv(varExplained1, file = "Convergence_eigen_values.csv")

# top 1000 Convergence PC1 genes with weights/loadings

Convergenceloadings <-Loadings(object = Convergence[["pca"]])[1:1000,1:2]

write.csv(Convergenceloadings, file = "Convergenceloadings.csv")

#R8 subcluster

R8 <- subset(eye, ident.use =c('R8') )

R8 <- RunPCA(R8, verbose = TRUE)

# eigen value for each PC for R8
pca1 = R8$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained2 = eigValues / sum(eigValues)
write.csv(varExplained2, file = "R8_eigen_values.csv")

# top 1000 R8 PC1 genes with weights/loadings

R8loadings <-Loadings(object = R8[["pca"]])[1:1000,1:2]

write.csv(R8loadings, file = "R8loadings.csv")

# R3/4 subcluster

R34 <- subset(eye, ident.use =c('R3/4') )

R34 <- RunPCA(R34, verbose = TRUE)

# eigen value for each PC for R34
pca1 = R34$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained3 = eigValues / sum(eigValues)
write.csv(varExplained3, file = "R34_eigen_values.csv")

# top 1000 R34 PC1 genes with weights/loadings

R34loadings <-Loadings(object = R34[["pca"]])[1:1000,1:2]

write.csv(R34loadings, file = "R34loadings.csv")

# subset R25

R25 <- subset(eye, ident.use =c('R2/5') )

R25 <- RunPCA(R25, verbose = TRUE)

# eigen value for each PC for R25
pca1 = R25$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained4 = eigValues / sum(eigValues)
write.csv(varExplained4, file = "R25_eigen_values.csv")

# top 1000 R25 PC1 genes with weights/loadings

R25loadings <-Loadings(object = R25[["pca"]])[1:1000,1:2]

write.csv(R25loadings, file = "R25loadings.csv")


# R1/6 subcluster

R16 <- subset(eye, ident.use =c('R1/6') )

R16 <- RunPCA(R16, verbose = TRUE)

# eigen value for each PC for R16
pca1 = R16$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained5 = eigValues / sum(eigValues)
write.csv(varExplained5, file = "R16_eigen_values.csv")

# top 1000 R16 PC1 genes with weights/loadings

R16loadings <-Loadings(object = R16[["pca"]])[1:1000,1:2]

write.csv(R16loadings, file = "R16loadings.csv")

# R7 subcluster 
R7 <- subset(eye, ident.use =c('R7') )

R7 <- RunPCA(R7, verbose = TRUE)
# eigen value for each PC for R7
pca1 = R7$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained6 = eigValues / sum(eigValues)
write.csv(varExplained6, file = "R7_eigen_values.csv")

# top 1000 R7 PC1 genes with weights/loadings

R7loadings <-Loadings(object = R7[["pca"]])[1:1000,1:2]

write.csv(R7loadings, file = "R7loadings.csv")

# Cones subcluster and PCA

Cones <- subset(eye, idents =c('Cones') )

Cones <- RunPCA(Cones, verbose = TRUE)

# eigen value for each PC for Cones
pca1 = Cones$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained7 = eigValues / sum(eigValues)
write.csv(varExplained7, file = "Cones_eigen_values.csv")

# top 1000 Cones PC1 genes with weights/loadings

Conesloadings <-Loadings(object = Cones[["pca"]])[1:1000,1:2]

write.csv(Conesloadings, file = "Conesloadings.csv")

# AUnd subcluster and PCA

AUnd <- subset(eye, idents =c('AUnd') )

AUnd <- RunPCA(AUnd, verbose = TRUE)

# eigen value for each PC for AUnd
pca1 = AUnd$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained8 = eigValues / sum(eigValues)
write.csv(varExplained8, file = "AUnd_eigen_values.csv")

# top 1000 AUnd PC1 genes with weights/loadings

AUndloadings <-Loadings(object = AUnd[["pca"]])[1:1000,1:2]

write.csv(AUndloadings, file = "AUndloadings.csv")


# PUnd subcluster and PCA

PUnd <- subset(eye, idents =c('PUnd') )

PUnd <- RunPCA(PUnd, verbose = TRUE)

# eigen value for each PC for PUnd
pca1 = PUnd$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained9 = eigValues / sum(eigValues)
write.csv(varExplained9, file = "PUnd_eigen_values.csv")

# top 1000 AUnd PC1 genes with weights/loadings

PUndloadings <-Loadings(object = PUnd[["pca"]])[1:1000,1:2]

write.csv(PUndloadings, file = "PUndloadings.csv")

# MF subcluster and PCA

MF <- subset(eye, idents =c('MF') )

MF <- RunPCA(MF, verbose = TRUE)

# eigen value for each PC for MF
pca1 = MF$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained10 = eigValues / sum(eigValues)
write.csv(varExplained10, file = "MF_eigen_values.csv")

# top 1000 MF PC1 genes with weights/loadings

MFloadings <-Loadings(object = MF[["pca"]])[1:1000,1:2]

write.csv(MFloadings, file = "MFloadings.csv")

# top 1000 AUnd PC1 genes with weights/loadings

PUndloadings <-Loadings(object = PUnd[["pca"]])[1:1000,1:2]

write.csv(PUndloadings, file = "PUndloadings.csv")

# PPN subcluster and PCA

PPN <- subset(eye, idents =c('PPN') )

PPN <- RunPCA(PPN, verbose = TRUE)

# eigen value for each PC for PPN
pca1 = PPN$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained11 = eigValues / sum(eigValues)
write.csv(varExplained11, file = "PPN_eigen_values.csv")

# top 1000 PPN PC1 genes with weights/loadings

PPNloadings <-Loadings(object = PPN[["pca"]])[1:1000,1:2]

write.csv(PPNloadings, file = "PPNloadings.csv")

# SMW subcluster and PCA

SMW <- subset(eye, idents =c('SMW') )

SMW <- RunPCA(SMW, verbose = TRUE)

# eigen value for each PC for SMW
pca1 = SMW$pca
eigValues = (pca1@stdev)^2  ## EigenValues
varExplained12 = eigValues / sum(eigValues)
write.csv(varExplained12, file = "SMW_eigen_values.csv")

# top 1000 SMW PC1 genes with weights/loadings

SMWloadings <-Loadings(object = SMW[["pca"]])[1:1000,1:2]

write.csv(SMWloadings, file = "SMWloadings.csv")



























































