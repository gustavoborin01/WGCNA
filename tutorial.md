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


rpkm=read.csv("file.csv")





