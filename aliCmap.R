#Ali Cmap
library(gdata)
library(PharmacoGx)
require(xtable)
library(GSA)

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains")
turnUp = FALSE

if(turnUp == TRUE)
{
  genes = read.xls("ListTwo_Ensembl_Turn_Up.xls", sheet = 1)
  direc = rep(1, length(genes$Ensembl.Gene.ID))
  fileNameAli = "drugs to turn up list two genes.xls"
}
if(turnUp == FALSE)
{
  genes = read.xls("ListOne_Ensembl_Turn_Down.xls", sheet = 1)
  direc = rep(-1, length(genes$Ensembl.Gene.ID))
  fileNameAli = "drugs to turn down list one genes.xls"
}

ensembleVec = paste(as.character(genes$Ensembl.Gene.ID), "_at", sep = "")


#below two lines are time consuming, so just use result as its data independent
#data(CMAPsmall)
#drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna")
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
load("cmap_sig_rna.RData")
geneSigCmap = as.data.frame(NULL)
geneSigCmap[1:length(ensembleVec),1] = ensembleVec
rownames(geneSigCmap) = ensembleVec
geneSigCmap[,2] = direc
colnames(geneSigCmap) = c("feature", "direction")
#geneSigCmap is formatted like HDAC_genes which is used in the pharmaco documentation example
data(HDAC_genes)

#remove 1:2 from column part of drug.perurbation to do all the drugs
res <- apply(drug.perturbation[,,c("tstat", "fdr")],2, function(x, genesCmap)
{
  return(connectivityScore(x=x,y=genesCmap[,2,drop=FALSE],method="gsea", nperm=100))
}, genesCmap=geneSigCmap)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)
res <- res[order(res[,1], decreasing=TRUE),]
xtable(res,caption='Connectivity Score results for the gene signature.')

resFix = as.data.frame(res)
resFix[,2] = replace(resFix[,2], resFix[,2] == 0, 1/(100+1))

write.table(resFix, fileNameAli, sep="\t")


