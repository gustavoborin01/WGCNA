# WGCNA (26/02/2015)
# Contrução da rede de co-regulação com os dados de RNASeq
#Lembrar de update os pacotes e, se necessário, até o R
#Instalando o pacote WGCNA

source("http://bioconductor.org/biocLite.R")
biocLite("impute")
install.packages("WGCNA")

workingDir="."
setwd("/home/user/Documents/dir")
getwd()
options(stringsAsFactors = FALSE)

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

#Se apareceu TRUE, não há outliers na amostra. Caso contrário, rodar...
if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#Clustering the samples...
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
#Add a linha cutoff para identificar e excluir os outliers
abline(h = 15, col = "red") 
#Salvar o plot no formato pdf
dev.off()






