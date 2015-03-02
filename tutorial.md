# WGCNA (26/02/2015)
# Contrução da rede de co-regulação com os dados de RNASeq
# Lembrar de update os pacotes e, se necessário, até o R
# Lembrar de salvar as tabelas no formato .csv (write.csv)
# Instalando o pacote WGCNA

source("http://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("WGCNA")

workingDir="."
setwd("/home/user/Documents/dir")
getwd()
options(stringsAsFactors = FALSE)

#########################################################
1. DATA INPUT AND CLEANING

rpkm=read.csv("file.csv",header=TRUE)
dim(rpkm)
names()

# Note that each row corresponds to a gene and column to a sample or auxiliary information.
# We now remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(rpkm[, -1]))
names(datExpr0) = rpkm$Genes                      #Usar o nome da primeira coluna após $
rownames(datExpr0) = names(rpkm)[-1]

# Checking data for excessive missing values and identification of outlier microarray samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# Se apareceu TRUE, não há outliers na amostra. Caso contrário, rodar...
if (!gsg$allOK)
    {
    if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
    }

# Clustering the samples...
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
# Add a linha cutoff para identificar e excluir os outliers
abline(h = 15, col = "red") 
# Salvar o plot no formato pdf
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, ght = 15, minSize = 10)              #Em 'ght' colocar o mesmo valor de h
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Save the objects and the file to a .RData format
save(rpkm,datExpr0, datExpr, gsg, clust, nGenes, nSamples,
keepSamples, sampleTree,file="dir/input.RData")

####################################################
2.STEP-BY-STEP NETWORK CONSTRUCTION

lnames=load("dir/input.RData")
lnames

# Choosing a soft-threshold to fit a scale-free topology to the network
powers=c(c(1:10)),seq(from=12,to=20,by=2))
sft=pickSoftThreshold(datExpr,dataIsExpr=TRUE,
    powerVector=powers,corFnc = cor,corOptions = list(use = 'p'),
    networkType = "unsigned")
# Plotting the results
sizeGrWindow(9,5)
# Save the plot
pdf(file="dir/scaleindependence.pdf",w=9,h=5)
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",
    ylab="Scale Free Topology Model Fit,signed R²",
    type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1],sft$fitIndices[,5],
    xlab="Soft Threshold (power)",
    ylab="Mean Connectivity",type="n",
    main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,
    cex=cex1,col="red")
dev.off()

# After chosing the power value, calculate the co-expression similarity and ajacency
softPower=10
adjacency=adjacency(datExpr,power=softPower,type="unsigned")
# Clustering using TOM
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM
geneTree=flashClust(as.dist(dissTOM),method="average")
sizeGrWindow(12,9)
pdf(file="dir/geneclustering.pdf")
plot(geneTree,xlab="",sub="",
    main="Gene clustering on TOM-based dissimilarity",
    labels=FALSE,hang=0.04)
minModuleSize=30
dynamicMods=cutreeDynamic(dendro=geneTree,method="tree", distM=dissTOM,
    deepSplit=2,pamRespectsDendro=FALSE,
    minClusterSize = minModuleSize)
table(dynamicMods)
write.csv(dynamicMods,file="dynamicMods.csv")
dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
write.csv(dynamicColors,file="dynamicColors.csv")
sizeGrWindow(8,6)
pdf(file="dir/dendogram_dynamictree.pdf",w=8,h=6)
plotDendroAndColors(geneTree,dynamicColors,
    "Dynamic Tree Cut",dendroLabels=FALSE,
    hang=0.03,addGuide=TRUE,guideHang=0.05,
    main="Gene dendrogram and module colors")
dev.off()

# Merging of modules whose expression profiles are very similar
MEList=moduleEigengenes(datExpr,colors=dynamicColors)
MEs=MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss=1-cor(MEs)
# Cluster module eigengenes
METree=flashClust(as.dist(MEDiss),method="average")
sizeGrWindow(7,6)
pdf(file="dir/clusteringeigengenes.pdf",w=7,h=6)
plot(METree,main="CLustering of module eigengenes",
    xlab="",sub="")
# A height cut of 0.25 correspond to correlation of 0.75    
MEDissThres=0.25
abline(h=MEDissThres,col="red")
dev.off()

merge=mergeCloseModules(datExpr,dynamicColors,
    cutHeight=MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs=merge$newMEs
sizeGrWindow(12,9)
pdf(file="dir/dendrogram_merged.pdf",w=12,h=9)
plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),
    c("Dynamic Tree Cut","Merged dynamic"),
    dendroLabels=FALSE,hang=0.03,
    addGuide=TRUE,guideHang=0.05)
dev.off()

moduleColors=mergedColors
colorOrder=c("grey",standardColors(50))
moduleLabels=match(moduleColors,colorOder)-1
MEs=mergedMEs
save(rpkm,datExpr0,datExpr,powers,sft,adjacency,TOM,dissTOM,
    geneTree,dynamicMods,dynamicColors,MEList,
    MEs,METree,merge,mergedColors,mergedMEs,
    moduleColors,colorOrder,moduleLabels,
    file="dir/networkconstruction_stepbystep.pdf")

# Networkheatmap construction
# Generating the heatmap plot for all genes take a substantial amount of time.So because this 
# it was necessary to restrict the number of genes
nSelect=1000
set.seed(10)
select=sample(nGenes,size=nSelect)
selectTOM=dissTOM[select,select]
selectTree=flashClust(as.dist(selectTOM),method="average")
selectColors=moduleColors[select]
sizeGrWindow(9,9)
pdf(file="dir/heatmap.pdf",w=9,h=9)
plotDIss=selectTOM^7
diag(plotDiss)=NA
TOMplot(plotDiss,selecTree,selectColors,
    main="Network heatmap plot, selected genes")
dev.off()

# Export the network into an edge list file VisANT can read
module="antiquewhite2"
probes=names(datExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modTOM=TOM[inModule,inModule]
dimnames(modTOM)=list(modProbes,modProbes)
vis=exportNetworkToVisANT(modTOM,file=paste("VisANTInput-",module,
".txt",sep=""),weighted=TRUE,threshold=0.2)
save(nSelect,select,selectTOM,selectTree,
    selectColors,plotDiss,file="heatmap_exportingdata.RData")
    





