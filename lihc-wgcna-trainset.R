rm(list = ls())
options(stringsAsFactors = F)
setwd('C:/Users/12sigma/Desktop/paper-revised-exp')

load('trainset_mrna_survival.RData')

library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)

datExpr = trainset_mrna[, 1:5001]
survival = subset(trainset_mrna, select = c('PATIENT_ID', 'OS_STATUS', 'OS_MONTHS'))
rownames(datExpr) = datExpr$PATIENT_ID
rownames(survival) = survival$PATIENT_ID
datExpr = datExpr[,-1]

## step 1 
if(T){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  png("step2-beta-value.png",width = 800,height = 600)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  abline(h=0.90,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}


## step2 
if(T){
  net = blockwiseModules(
    datExpr,
    power = sft$powerEstimate,
    maxBlockSize = 6000,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = F, 
    verbose = 3
  )
  table(net$colors) 
}


## step 3
if(T){
  mergedColors = labels2colors(net$colors)
  table(mergedColors)
  moduleColors=mergedColors
  
  png("step3-genes-modules.png",width = 800,height = 600)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}

if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  par(mar = c(0,5,2,0))
  plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1,cex.lab=1)

  sample_colors <- numbers2colors(as.numeric(factor(survival$OS_STATUS)), 
                                  colors = c("blue","red"),signed = FALSE)
  par(mar = c(1,4,3,1),cex=0.8)
  
  png("sample-survival-cluster.png",width = 800,height = 600)
  plotDendroAndColors(datExpr_tree, sample_colors,
                      groupLabels = colnames(sample),
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
}


## step 4
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  moduleColors <- labels2colors(net$colors)
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  write.csv(MEs, "trainset_eigengene.csv")
}


## step 5
if(T){
  datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")
  module = "turquoise"
  modNames = substring(names(MEs), 3)
  column = match(module, modNames)
  moduleGenes = moduleColors==module
  turquoise_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
  names(turquoise_module)="genename"
  turquoise_KME<-as.data.frame(datKME[moduleGenes,column]) 
  names(turquoise_KME)="KME"
  rownames(turquoise_KME)=turquoise_module$genename
  FilterGenes = abs(turquoise_KME$KME) > 0.8
  table(FilterGenes)

  turquoise_hub<-subset(turquoise_KME, abs(turquoise_KME$KME)>0.8)
  write.csv(turquoise_hub, "hubgene_KME_turquoise_0.8.csv")
  
  HubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
  write.csv (HubGenes,file = "TopHubGenes_of_each_module.csv",quote=F)
}


