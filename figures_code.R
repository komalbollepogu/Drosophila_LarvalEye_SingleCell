#figure 1D

DimPlot(eye, label = F, pt.size = 1.5)+ NoLegend()

#Figure 2 A-I

FeaturePlot(eye, "h", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "Optix", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "dpp", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "sens", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "ro", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "svp", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "B-H1", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "pros", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

FeaturePlot(eye, "lz", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#Figure 2J dotplot

features1 <- c('toy','hth','eya','so','ey','tsh','eyg','Optix','h','dac','E(spl)m8-HLH','E(spl)m4-BFM','dpp','rn','ato','sca','sens','boss','run','ro',
               'svp','salm','Claspin','PCNA','B-H1','pros','sv','ct','lz','sev')


DotPlot(eye,dot.scale = 8, features = features1, idents =c("AUnd", "PPN",
                                                           "MF","R8", 
                                                           "R2/5", "R3/4",'SMW',
                                                           "R1/6", "R7","Cones", 
                                                           "PUnd",'Convergence' ), 
        cluster.idents = F)  + theme(axis.text.x = element_text(angle = 45,  size = 15),axis.text.y = element_text(angle = 0, 
                                                                                                                   hjust = 1, size = 12, face = 'italic'),axis.title.y.right = element_text(size = 9, face = 'bold'),
                                     axis.title = element_text(size=15,face="bold"),legend.text=element_text(size=10),
                                     legend.title=element_text(size=10),axis.line = element_line(size=0.6))+coord_flip()

#figure 3 B-E


FeaturePlot(eye, "CG42458", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "liprin-gamma", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG34347", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "fipi", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#Fig 3G dot plot


features2 <- c('CG6414','AdamTS-A',"CG32447",'a','Nplp1','CG30460','Achl',
               'CG42458','Lim3','Liprin-gamma','scb','CG34347','CG34371',
               'CG11382','aurA','RnrS','CG1218','CG31221','CG15630','pyr','ths','Cpr49Ac',
               'CAP','Lrt','SoxN','Wnt2','NijC','wat','Gasp','CG11029','magu',
               'inv','Rbp6','Ace','brp')

DotPlot(eye,dot.scale = 8, features = features1, idents =c("AUnd", "PPN",
                                                           "MF","R8", 
                                                           "R2/5", "R3/4",'SMW',
                                                           "R1/6", "R7","Cones", 
                                                           "PUnd",'Convergence' ), 
        cluster.idents = F)  + theme(axis.text.x = element_text(angle = 45,  size = 15),axis.text.y = element_text(angle = 0, 
                                                                                                                   hjust = 1, size = 12, face = 'italic'),axis.title.y.right = element_text(size = 9, face = 'bold'),
                                     axis.title = element_text(size=15,face="bold"),legend.text=element_text(size=10),
                                     legend.title=element_text(size=10),axis.line = element_line(size=0.6))+coord_flip()



#Figure 4 A-D and H


DimPlot(R34, label = F, pt.size = 3.5)+ NoLegend()

FeaturePlot(R34, "svp", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(R34, "DIP-delta", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(R34, "CG4341", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

features4 <- c('pk','salm','svp','klg','Drl-2','Swip-1','sda','DIP-delta','DIP-theta','blanks','px','Magi',
               'CG4341','dpr1','ru','Hmx','beat-IIIc')

# fig 4H dot plot
DotPlot(R34,dot.scale = 10, features = features4, idents = c('Early R3/4','R3','R4') )  + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20, face = 'italic'),
                                                                                                axis.text.y = element_text(angle = 0, hjust = 1, size = 20),
                                                                                                axis.title.y.right = element_text(size = 20),
                                                                                                axis.title = element_text(size=20,face="bold"),
                                                                                                legend.text=element_text(size=10),
                                                                                                legend.title=element_text(size=10),
                                                                                                axis.line = element_line(size=1.25))



#Figure 5A-D, I,K,M

DimPlot(furrowPRs, label = F, pt.size = 5.5)+ NoLegend()

FeaturePlot(furrowPRs, "ato", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowPRs, "sens", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowPRs, "boss", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowPRs, "onecut", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowPRs, "chp", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowPRs, "qvr", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')


features3 <- c('dpp','sens','ro','svp','B-H1','pros','brp','Rbp6','Ace','jeb','chp','mspo','dpr9','nolo','Rim',
               'tau','dpr12','2mit',
               'qvr','Pde6',
               'Octbeta3R','mtt','Pde1c','Trpm','CG4168')

DotPlot(eye,dot.scale = 8, features = features2, idents =c("AUnd", "PPN",
                                                           "MF","R8", 
                                                           "R2/5", "R3/4",'SMW',
                                                           "R1/6", "R7","Cones", 
                                                           "PUnd",'Convergence',"PC",'PPD','LM','Oc' ), 
        cluster.idents = F)  + theme(axis.text.x = element_text(angle = 45,  size = 15),axis.text.y = element_text(angle = 0, 
                                                                                                                   hjust = 1, size = 12, face = 'italic'),axis.title.y.right = element_text(size = 9, face = 'bold'),
                                     axis.title = element_text(size=15,face="bold"),legend.text=element_text(size=10),
                                     legend.title=element_text(size=10),axis.line = element_line(size=0.6))+coord_flip()

#Figure6
#B
DimPlot(eye_ATAC, label = F, pt.size = 1.5)+ NoLegend()

#D
FeaturePlot(eye, "dac", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#E
CoveragePlot(
  object = eye_ATAC,
  region = "dac",  extend.upstream = 15000,
  extend.downstream =5000)

#F
CoveragePlot(
  object = eye_ATAC,
  region = "CAP",  extend.upstream = 5000,
  extend.downstream =5000
)

#G
FeaturePlot(eye, "CAP", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
#H
FeaturePlot(eye_ATAC, "CAP", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#Figure 7
#A
CoveragePlot(
  object = eye_ATAC,
  region = "sens",  extend.upstream = 5000,
  extend.downstream =5000
)

#B
FeaturePlot(eye, "sens", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#B'

FeaturePlot(eye_ATAC, "sens", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
#C
CoveragePlot(
  object = eye_ATAC,
  region = "Wnt2",  extend.upstream = 5000,
  extend.downstream =5000
)

#D
FeaturePlot(eye, "Wnt2", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#D'

FeaturePlot(eye_ATAC, "Wnt2", pt.size = 3.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

# figure 8 complex heatmap

ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))


# Supplemental figures
#suppl. Figure1 B-K

FeaturePlot(eye, "salm", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "mirr", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Ance", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "pnr", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "odd", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "sob", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "drm", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "oc", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "en", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "hh", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

#Figure 1 M and N
DimPlot(eye, label = F, pt.size = 1.5, split.by = 'orig.ident')+ NoLegend()


#Suppl.Figure2 A-s

FeaturePlot(eye, "ase", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG42313", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "side", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "scb", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Lim3", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG11382", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG34371", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG31221", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG45105", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Cpr49Ac", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "pyr", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "ths", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Lrt", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "SoxN", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "jv", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "dnd", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG45263", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "aay", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Cad96Cb", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')



#Supplemental Figure 4 A-M


FeaturePlot(eye, "CG6414", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG32447", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "arc", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Nplp1", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG30460", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Achl", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "aurA", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "RnrS", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG1218", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "inv", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Gasp", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG11029", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "magu", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')


#Supplemental Figure 5A-E
DimPlot(R34, label = F, pt.size = 3.5)+ NoLegend()

FeaturePlot(R34, "E(spl)mdelta-HLH", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "E(spl)malpha-BFM", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "E(spl)mbeta-HLH", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "E(spl)mgamma-HLH", pt.size = 1.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')


#Supplemental Figure 7 A-C

DimPlot(furrowprs, label = F, pt.size = 3.5)+ NoLegend()

FeaturePlot(furrowprs, "svp", pt.size =5.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(furrowprs, "B-H1", pt.size = 5.5, cols = c('gray72','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

# Supplemental Figure 8 A-k

FeaturePlot(eye_ATAC, "Optix", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "dpp", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "ato", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "sens", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "ro", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "svp", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "B-H1", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "pros", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "ct", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "lz", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')


DimPlot(eye_ATAC, label = F, pt.size = 1.5)+ NoLegend()


#Supplemental Figure 9 A-C, A'-E'

FeaturePlot(eye, "dac", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "ato", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "lz", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "sv", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "pros", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

CoveragePlot(
  object = eye_ATAC,
  region = "dac",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "ato",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "lz",  extend.upstream = 5000,
  extend.downstream =5000
)


CoveragePlot(
  object = eye_ATAC,
  region = "sv",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "pros",  extend.upstream = 5000,
  extend.downstream =5000
)

# Supplemental Figure 10 A-G

FeaturePlot(eye, "sv", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "sv", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "lz", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

CoveragePlot(
  object = eye_ATAC,
  region = "2mit",  extend.upstream = 5000,
  extend.downstream =5000
)

FeaturePlot(eye, "2mit", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye_ATAC, "2mit", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')



# Supplementary Figure 11 A-E and A'-E'



FeaturePlot(eye, "Lim3", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "lnRNA:CR45401", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "lnRNA:CR44344", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG13323", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG46301", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

CoveragePlot(
  object = eye_ATAC,
  region = "Lim3",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "lnRNA:CR45401",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "lnRNA:CR44344",  extend.upstream = 5000,
  extend.downstream =5000
)


CoveragePlot(
  object = eye_ATAC,
  region = "CG13323",  extend.upstream = 5000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "CG46301",  extend.upstream = 5000,
  extend.downstream =5000
)


# Supplementary Figure 12 A-D and A'-D'

FeaturePlot(eye, "fipi", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "pyr", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "Cpr49Ac", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "CG45105", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

CoveragePlot(
  object = eye_ATAC,
  region = "fipi",  extend.upstream = 25000,
  extend.downstream =5000
)

CoveragePlot(
  object = eye_ATAC,
  region = "pyr",  extend.upstream = 5000,
  extend.downstream =55000
)

CoveragePlot(
  object = eye_ATAC,
  region = "Cpr49Ac",  extend.upstream = 5000,
  extend.downstream =5000
)


CoveragePlot(
  object = eye_ATAC,
  region = "CG45105",  extend.upstream = 5000,
  extend.downstream =5000
)

# Supplementary Figure 13 A-B
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"])
rssPlot <- plotRSS(rss)
fig<- plotly::ggplotly(rssPlot$plot)

FeaturePlot(eye, "Rbp6", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

# Supplementary Figure 14 D-L

FeaturePlot(eye, "Optix", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Ariss, "Optix", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Gonzales, "Optix", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "dpp", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Ariss, "dpp", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Gonzales, "dpp", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(eye, "sens", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Ariss, "sens", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')
FeaturePlot(Gonzales, "sens", pt.size = 1.5, cols = c('gray','blue'), max.cutoff = 'q90',min.cutoff = 'q0')

















