
# Integrated genomic analyses of ovarian carcinoma. <-- TCGA Paper
# import limma package and lumi package that was mentioned in analysis report
# seems a time consuming step when run manually, fix to only be done if need be (aka dont update all always)
#.libPaths("D://rLib")
#source("https://bioconductor.org/biocLite.R");
#biocLite("limma", suppressUpdates=TRUE);
library(limma);
#biocLite("genefu")
library(genefu)
library(illuminaHumanv4.db)
#biocLite("TCGAbiolinks")
library(TCGAbiolinks)
library(matrixStats)
library(ROCR)
library(org.Hs.eg.db)
library(parallel)
#library(jetset)
library(Biobase)
library(gdata)
library(GSVA)
#MetaGx Packages
#biocLite("Biobase")
#biocLite("genefilter")
#biocLite("logging")
#biocLite("survcomp")

#query <- GDCquery(project = "TCGA-OV",
#                  data.category = "Transcriptome Profiling",
#                  data.type = "Gene Expression Quantification", 
#                  workflow.type = "HTSeq - Counts")

#previous report has sample names in figures so interchange barcode with sample names when comparing

#sample name
#Sample number
#Array Barcode 
sampNameMat = matrix(c(
  "72199_49e+", "A12", "9421742080_A",
  "70535_133-", "A32", "9421742080_B",
  "71423_49e+","A36", "9421742080_C",
  "71377_133-", "A26", "9421742080_D",
  "71029_49e+", "A21", "9421742080_E",
  "72143_133+","A4", "9421742080_F",
  "71853_49e+", "A30", "9421742080_G",
  "71377_133+", "A25", "9421742080_H",
  "72143_49e+", "A6", "9421742080_I",
  "72130_133+", "A16", "9421742080_J",
  "72130_49e+", "A18", "9421742080_K",
  "70924_133-", "A23", "9421742080_L",
  "65846_133+", "A1", "9421742096_A",
  "65846_133-", "A2", "9421742096_B",
  "68584_133+", "A13", "9421742096_C",
  "72199_133+", "A10", "9421742096_D",
  "72143_133-", "A5", "9421742096_E",
  "71029_133-", "A20", "9421742096_F",
  "72199_133-", "A11", "9421742096_G",
  "71029_133+", "A19", "9421742096_H",
  "70535_133+", "A31", "9421742096_I",
  "67794_133-", "A8", "9421742096_J",
  "65846_49e+", "A3", "9421742096_K",
  "67794_133+", "A7", "9421742096_L",
  "70535_49e+", "A33", "9421898014_A",
  "67794_49e+", "A9", "9421898014_B",
  "68584_133-", "A14", "9421898014_C",
  "70924_49e+", "A24", "9421898014_D",
  "68584_49e+", "A15", "9421898014_E",
  "71853_133-", "A29", "9421898014_F",
  "71423_133-", "A35", "9421898014_G",
  "71377_49e+", "A27", "9421898014_H",
  "72130_133-", "A17", "9421898014_I",
  "71423_133+", "A34", "9421898014_J",
  "70924_133+", "A22", "9421898014_K",
  "71853_133+", "A28", "9421898014_L"), 
  36, 3, byrow = TRUE);



sampNames = sampNameMat[,1]; 

#should use read.idat function to get data matrix from scratch later using natalie's raw idat files



#below has entrex gene ids and gene names for dataframe x used below in the analysis, cant get it to load with read.ilmn though
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Natalie Files");
rawData = read.delim("sample_probe_nonnorm_minusbackground_FinalReportNoHeaderNotepad.txt");
geneNameColInd = which(grepl("TargetID",colnames(rawData)))
entrezIDColInd = which(grepl("ENTREZ_GENE_ID",colnames(rawData)))
accessionInd = which(grepl("ACCESSION",colnames(rawData)))

#GeneSpring files seem to be working fine
#set directory below to one containing probe summary profiles from GeneSpring
setwd("C:\\Users\\micha\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data GeneSpring");
#47300 probes
x <- read.ilmn("Ali HT-12 expression_Sample_Probe_Profile_nonnormNoHeader.txt"); 
colnames(x$E) = sampNames;
x$geneNames = as.vector(rawData[, geneNameColInd]);
x$entrezIds = as.vector(rawData[, entrezIDColInd]);
x$probeIds = rownames(x$E);
x$acc = as.vector(rawData[, accessionInd]);

xOrig = x;

#sigVec = c()
#for(i in 1:dim(rawData)[2])
#{
#  print(grep("AVG_Signal", colnames(rawData)[i]))
#    if(length(grep("AVG_Signal", colnames(rawData)[i])) > 0)
#      sigVec = c(sigVec, i) 
#}
#rawExp = rawData[,sigVec]
#colnames(rawExp) = sampNames
#x$E = rawExp

#Normalization and log Step
x <- normalizeBetweenArrays(x,method="quantile")
#results no better at number genes = 20 using neqc, although the control probes were not supplied
#x <- neqc(x)
x$E = log2(x$E)

boxplot(log2(x$E),range=0,ylab="log2 intensity")
plotDensity(x$E)


plus49 = which(grepl("49", sampNames));
plus133 = which(grepl("133\\+", sampNames));
minus133 = which(grepl("133-", sampNames));

plus49Genes = x[, plus49];
plus133Genes = x[, plus133];
minus133Genes = x[, minus133];

# decide to focus on CAF relevant cell line from this point and onwards
#33411 probes ater this loop/step (if skip log(x$E) only? no because below says the same)
x = x[, plus49];

#note name botTwenty is a relic from the analysis report, value of 5 implies bottom 95 percentile
percentile = 5;
numPatients = 12;
numProbes = dim(x$E)[1];
samples = dim(x$E)[2];
numbCellLines = samples/numPatients
botTwentyPercentile = ceiling(numProbes*((100 - percentile)/100))
probeCount = matrix(0, numProbes, numbCellLines)
twentyPercOfGroup = ceiling(numPatients*(percentile/100));

for(i in 1:samples)
{
  samp = x$E[1:numProbes,i];
  sampSort = sort(samp, decreasing = TRUE);
  sortInd = sort(samp, decreasing = TRUE, index.return=TRUE)$ix;
  cellLineNumb = floor((i - 1)/numPatients) + 1;
  bottomTwentyTile = sortInd[botTwentyPercentile:numProbes];
  probeCount[bottomTwentyTile, cellLineNumb] = probeCount[bottomTwentyTile, cellLineNumb] + 1;
}

removeProbes = rowSums(probeCount >= twentyPercOfGroup);
removeProbes = which( removeProbes == numbCellLines, arr.ind=TRUE);
x$E = x$E[-removeProbes,];
x$entrezIds = x$entrezIds[-removeProbes];
x$probeIds = x$probeIds[-removeProbes];
x$geneNames = x$geneNames[-removeProbes];
x$other$Detection = x$other$Detection[-removeProbes, ]


#44000
missingEntrezIds = which(is.na(x$entrezIds));
x$E = x$E[-missingEntrezIds,];
x$entrezIds = x$entrezIds[-missingEntrezIds];
x$probeIds = x$probeIds[-missingEntrezIds];
x$geneNames = x$geneNames[-missingEntrezIds];
x$other$Detection = x$other$Detection[-missingEntrezIds, ]

xBeforeNorm = x;

boxplot(x$E,range=0,ylab="log2 intensity")

#Remove duplicate probes based on IQR

sortInd = sort(x$entrezIds, decreasing = FALSE, index.return=TRUE)$ix;
x$entrezIds = x$entrezIds[sortInd];
x$E = x$E[sortInd, ];
x$probeIds = x$probeIds[sortInd];
x$geneNames = x$geneNames[sortInd];

bestProbes = c()
i = 1;
while(i < dim(x$E)[1])
{
  entrezId = x$entrezIds[i];
  origInd = i;
  probesWithId = c();
  while(x$entrezIds[i] == entrezId)
  {
    probesWithId = c(probesWithId, i);
    i = i + 1;
  }
  if(length(probesWithId) > 1)
  {
    probes = x$E[probesWithId,]
    iqrs = rowIQRs(probes);
    keepProbe = which(iqrs == max(iqrs)) + (origInd - 1);
    bestProbes = c(bestProbes, keepProbe)
  }
  if(length(probesWithId) == 1)
  {
    bestProbes = c(bestProbes, origInd)
  }
}

x = x[bestProbes, ];
x$geneNames = x$geneNames[bestProbes];
x$entrezIds = x$entrezIds[bestProbes];
x$probeIds = x$probeIds[bestProbes];

xNormAndRemoved = x;
xGen = x;

dissimilarity <- 1 - cor(x$E);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Clustering of Probes", xlab="");
heatmap(x$E, Rowv=NA, Colv=as.dendrogram(hc))

clustOneOut = FALSE;

if(clustOneOut == TRUE)
{
  xFull = x;
  for(i in 1:dim(xFull$E)[2])
  {
    print(i)
    x = xFull
    x = x[, -i];
    dissimilarity <- 1 - cor(x$E);  #cor default is pearson
    distance <- as.dist(dissimilarity);
    #clusters <- hclust(distance, method = "average");
    #clusterCut <- cutree(clusters, 2)
    plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="");
    readline()
    plot(1, i)
    print("next")
  }
  
}

clusterPlot = FALSE;
if(clusterPlot == TRUE)
{
#install.packages("latticeExtra")
library(latticeExtra)
library(reshape2)
setwd("C:\\Users\\Michael\\Documents\\OVC Project");
groupingMat = as.matrix(read.csv("Leave One Out Clustering Diag.csv"));
rownames(groupingMat) = colnames(groupingMat)
gfMatMelt = melt(groupingMat)
cloud(value~Var1+Var2, gfMatMelt, panel.3d.cloud=panel.3dbars, col.facet='grey', xbase=0.4, ybase=0.4, scales=list(arrows=FALSE, col=1),par.settings = list(axis.line = list(col = "transparent")), xlab = "Samples Left Out of Clustering", ylab = "Sample", zlab = "Group")
}

fapBoxPlots = FALSE;
if(fapBoxPlots == TRUE)
{
  
  #Take 2
  #Group 1 Fap Lo = 72143, 72130, 70535, 71377, 71029 semi outlier 72199, 68584. col = 5,6,8,12,3,1,11 
  #Group 2 Fap Hi = 71423, 71853, 67794, 70924, 65846. col = 2,4,9,10,7
  
  #these are the fap hi percentages measured by ali
  aliFapLoCafGroupProt = c(18.3, 7.16, 14, 19.8, 0, 2.13, 51.1)
  aliFapHiCafGroupProt = c(29.6, 23, 41.2, 75.7, 15.6)
  boxplot(aliFapLoCafGroupProt, aliFapHiCafGroupProt, names = c("Fap Lo Group", "Fap Hi Group"), ylab = "Lab Measured Fap Hi %")
  
  fapRow = which(x$entrezIds == 2191);
  fapGeneExp = x$E[fapRow, ]
  aliFapHiGroupExp = c(fapGeneExp[2], fapGeneExp[4], fapGeneExp[7], fapGeneExp[9], fapGeneExp[10]);
  aliFapLoGroupExp = c(fapGeneExp[3], fapGeneExp[12], fapGeneExp[8], fapGeneExp[6], fapGeneExp[5], fapGeneExp[11], fapGeneExp[1]);
  #although order is switched presumably Fap Lo group is meant to be fap hi group
  boxplot(aliFapLoGroupExp, aliFapHiGroupExp, names = c("Fap Lo Group", "Fap Hi Group"), ylab = "FAP Gene Expression")
  
#Bar plot
fapMeasured <- c(29.6, 23, 15.6, 41.2, 75.7, 0, 19.8, 14, 7.16, 18.3, 51.1, 2.13)
colours = c("blue","blue","blue","blue","blue","red","red","red","red","red","red","red")
labelsPat = c("71423", "71853", "65846", "67794", "70924", "71029", "71377","70535", "72130", "72143", "68584", "72199")
barplot(fapMeasured, col=colours, xlab = "Samples (Blue is Fap Hi Red is Fap Lo)", ylab = "Ali Fap High % Measured",  names.arg=labelsPat)
}

c
diffExpress = TRUE

if(diffExpress == TRUE)
{
cafGroupOneGenes = x[,c(2, 4, 7, 9, 10)];
cafGroupTwoGenes = x[,c(1, 3, 5, 6, 8, 11, 12)];
#Ali's outliers
#cafGroupOtherGenes = x[,c(1, 11)];
cafPValsOneTwo = rep(0, dim(x$E)[1]);
cafGroupDiff = rep(0, dim(x$E)[1]);

for(i in 1:dim(x$E)[1])
{
  #wilcox is var1 - var2 so cafGroupOne - cafGroupTwo, so fapHi - fapLo according to Ali's measurements
  cafGroupOneGene = c(cafGroupOneGenes$E[i,]);
  cafGroupTwoGene = c(cafGroupTwoGenes$E[i,]);
  
  wilResult = wilcox.test(cafGroupOneGene, cafGroupTwoGene, paired=FALSE, conf.int = TRUE)
  cafPValsOneTwo[i] = wilResult[3];
  cafGroupDiff[i] = wilResult[9]
}

#adjusted p values all above 0.05??? non adjusted has 883 differentially expressed genes at p=0.05 or 246 at p=0.01
minP = 0.01;
adjustedCafPValsOneTwo = p.adjust(cafPValsOneTwo, method = "bonferroni", n = length(cafPValsOneTwo));
sigGeneIndsNoAdj = which(cafPValsOneTwo < minP)
sigGeneInds = which(adjustedCafPValsOneTwo < minP);
sigGeneIndsFive = which(cafPValsOneTwo < 0.05);
sigGeneIndsOne = which(cafPValsOneTwo < 0.01);
}

#define a gene signature
patientOutliers = c(1, 11);

#Begin regression analysis/gene signature creation. xSig is just diff expressed genes
#xGen is all genes after normalizing, removing duplicates and low percentile. unbiased approach uses xGen

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\R Code");
source("one-fold-CrossValidation.R");
source("oneFolCVRankings.R");
source("oneFoldValRoc.R");

xGenReg = t(xGen[sigGeneIndsNoAdj,]$E);

#fapHi values for patients in same order as col names in xGenReg and xSigReg, Ali should provid NA values soon
#fapHi = c(2.13, 70.4, 100, 77, 81.7, 92.8, 15.6, 86, 58.8, 24.2, 51.1, 19.8)
#fapGroup= c(1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1);
fapGroup = c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0);

#naVals = which(is.na(fapHi)==TRUE);
#xGenReg = xGenReg[-naVals, ];
#fapHi = fapHi[-naVals];

#using parallelizing takes around 15 minutes for 6

minGeneNumb = 1;
maxGeneNumb = 1000;
geneInc = 1;
corVsGenes = c();
numGenes = c();

regressionMod = FALSE;

if(regressionMod == TRUE)
{
geneNumVec = seq(from=minGeneNumb, to=maxGeneNumb, by=geneInc);
no_cores <- detectCores() - 2;
cl <- makeCluster(no_cores)
clusterExport(cl, list("onefoldcv", "myRankingOneFold"))
oneFoldResults = clusterApplyLB(cl, geneNumVec, function(x, y, z) onefoldcv(x, y, z), y = xGenReg, z = fapGroup)
stopCluster(cl);

for(i in 1:length(oneFoldResults))
{
  oneFoldResult = oneFoldResults[[i]][, 1];
  corVsGenes = c(corVsGenes, cor(oneFoldResult, fapGroup));
}
plot(geneNumVec, corVsGenes, xlab = "number of genes used in sig score", ylab = "sample's sig score and fapGroup correlation");
lines(geneNumVec, corVeGenes)
}

#FAP Hi group
cafGroupOneGenes = x[,c(2, 4, 7, 9, 10)]; #right cluster
#cafGroupTwoGenes = x[,c(3, 5, 6, 8, 12)];
#FAP Lo Group
cafGroupTwoGenes = x[,c(1, 3, 5, 6, 8, 11, 12)];
#Ali's outliers
#cafGroupOtherGenes = x[,c(1, 11)];
xBeforeRemove = x; 

#could make things quicker by checking which genes are never in top 500 once and removing them
#for(i in seq(minGeneNumb, maxGeneNumb, geneInc))
#{
#  x = xBeforeRemove;
#  scores = NULL

#get the diffeentially expressed gene ranks once for the case of having each sample removed
#then get the scores for the top n genes

xOneOutOrder = list()
cafGroupDiffOneOut = list()
jEnd = dim(x)[2]
  for(j in 1:jEnd)
  {
    x = xBeforeRemove;
    print(j)
    xSig = xBeforeRemove$E;
    testSamp = t(xSig[, j, drop = FALSE]);
    x = x[, -j]
    cafPValsOneTwo = rep(0, dim(x$E)[1]);
    cafGroupDiff = rep(0, dim(x$E)[1]);
    cafGroupLogFc = rep(0, dim(x$E)[1]);
    
    for(k in 1:dim(x$E)[1])
    {
      cafGroupOneGene = c(cafGroupOneGenes$E[k,]);
      cafGroupTwoGene = c(cafGroupTwoGenes$E[k,]);
      
      #group one is fap Hi, group 2 is fap lo, so this is fap hi - fap lo in test
      wilResult = wilcox.test(cafGroupOneGene, cafGroupTwoGene, paired=FALSE, conf.int = TRUE)
      cafPValsOneTwo[k] = wilResult[3];
      cafGroupDiff[k] = wilResult[9]
      cafGroupLogFc[k] = log(median(cafGroupOneGene)/median(cafGroupTwoGene))
    }
    
    cafPValsOneTwo = as.numeric(cafPValsOneTwo);
    cafGroupDiff = as.numeric(cafGroupDiff);
    cafGroupLogFc = as.numeric(cafGroupLogFc)
    
    sigGenes = which(cafPValsOneTwo <= 0.05)
    cafPValsOneTwoSig = cafPValsOneTwo[sigGenes]
    cafGroupDiffSig = cafGroupDiff[sigGenes]
    cafGroupLogFcSig = cafGroupLogFc[sigGenes]
    xOrd =  x[sigGenes, ]
    
    #sigRank = order(abs(cafGroupLogFcSig), decreasing = TRUE);
    sigRank = order(abs(cafGroupDiffSig), decreasing = TRUE);
    xOrd = xOrd[sigRank, ]
    cafGroupDiffSig = cafGroupDiffSig[sigRank];
    cafPValsOneTwoSig = cafPValsOneTwoSig[sigRank];
    cafGroupLogFcSig = cafGroupLogFcSig[sigRank]
    
    xOneOutOrder[[j]] = xOrd;
    cafGroupDiffOneOut[[j]] = cafGroupDiffSig;
  }
    
scoreList = list();
aucVec = c();
corVec = c();
x = xBeforeRemove;
posCont = FALSE;

for(i in seq(minGeneNumb, maxGeneNumb, geneInc))
{
  scores = NULL;
for(j in 1:dim(x)[2])
{
  xSig = xBeforeRemove$E;
  testSamp = t(xSig[, j, drop = FALSE]);
  
  topGenes = rownames(xOneOutOrder[[j]])[1:i]
  mysig <- cbind("probe"=topGenes, "EntrezGene.ID"=NA, "coefficient"=(cafGroupDiffOneOut[[j]][1:i]) / sum(abs(cafGroupDiffOneOut[[j]][1:i]), na.rm=TRUE))
  if(posCont == TRUE)
  {
    randGenes = sample(1:dim(x)[1], i, replace = T)
    topGenes = rownames(xOneOutOrder[[j]])[randGenes]
    mysig <- cbind("probe"=topGenes, "EntrezGene.ID"=NA, "coefficient"=(cafGroupDiffOneOut[[j]][randGenes]) / sum(abs(cafGroupDiffOneOut[[j]][randGenes]), na.rm=TRUE))
  }
  rownames(mysig) = topGenes;
  scores <- rbind(scores, cbind("score"=genefu::sig.score(x=mysig, testSamp, annot=geneNameFrame, do.mapping=FALSE, signed=TRUE)$score, "fold"=j))
}
  rownames(scores) = colnames(x);
  scoreList[[i - (minGeneNumb - 1)]] = scores;
  corVec = c(corVec, cor(as.vector(scores[, 1]), fapGroup))
 
 rocPred = prediction(as.vector(scores[, 1]), fapGroup)
 rocPerf <- performance(rocPred,"tpr","fpr")
 #plot(rocPerf,col="black")
 #abline(a=0, b= 1, col ="green")
 auc <- performance(rocPred,"auc")
 auc <- unlist(slot(auc, "y.values"))
 aucVec = c(aucVec, auc)

}

#tested this and saw expected result, as one goes farther into xOneOutOrder the means differ less
geneId = rownames(xOneOutOrder[[1]])[1]
probeRow = which(rownames(cafGroupOneGenes) == geneId)
geneExpGroupOne = cafGroupOneGenes$E[probeRow, ]
geneExpGroupTwo = cafGroupTwoGenes$E[probeRow, ]
boxplot(geneExpGroupOne, geneExpGroupTwo, names = c("Group One", "Group Two"), ylab = "Gene Expression")

#Get the gene signature, take top 50 diff expressed genes among all samples

cafPValsOneTwo = rep(0, dim(x$E)[1]);
cafGroupDiff = rep(0, dim(x$E)[1]);
cafGroupLogFc = rep(0, dim(x$E)[1]);
#k = 9 is a tie case
for(k in 1:dim(x$E)[1])
{
  cafGroupOneGene = c(cafGroupOneGenes$E[k,]);
  cafGroupTwoGene = c(cafGroupTwoGenes$E[k,]);
  
  #Fap Hi (Group 1) - FAP Lo (Group 2)
  wilResult = wilcox.test(cafGroupOneGene, cafGroupTwoGene, paired=FALSE, conf.int = TRUE)
  cafPValsOneTwo[k] = wilResult[3];
  cafGroupDiff[k] = wilResult[9]
  cafGroupLogFc[k] = log(median(cafGroupOneGene)/median(cafGroupTwoGene))
}

xPreGeneSig = x
cafPValsOneTwo = as.numeric(cafPValsOneTwo);
cafGroupDiff = as.numeric(cafGroupDiff);
cafGroupLogFc = as.numeric(cafGroupLogFc)

sigGenes = which(cafPValsOneTwo <= 0.05)
cafPValsOneTwoSig = cafPValsOneTwo[sigGenes]
cafGroupDiffSig = cafGroupDiff[sigGenes]
cafGroupLogFcSig = cafGroupLogFc[sigGenes]
xOrd =  x[sigGenes, ]
xOrd$geneNames = xOrd$geneNames[sigGenes]

#sigRank = order(abs(cafGroupLogFcSig), decreasing = TRUE);
sigRank = order(abs(cafGroupDiffSig), decreasing = TRUE)
cafGroupDiffSig = cafGroupDiffSig[sigRank];
cafPValsOneTwoSig = cafPValsOneTwoSig[sigRank];
cafGroupLogFcSig = cafGroupLogFcSig[sigRank]
xOrd$geneNames = xOrd$geneNames[sigRank]

xOrd = xOrd[sigRank, ]

#positive u-stat and fold change means it is expressed more in FAP Hi/Group 1
#aliTab = cbind(cafGroupDiffSig, cafGroupLogFcSig, cafPValsOneTwoSig)
#rownames(aliTab) = xOrd$geneNames
#colnames(aliTab) = c("u-stat", "log(FC)", "u-stat p value")
#orderTab = order(abs(aliTab[, 1]), decreasing = TRUE)
#aliTabSorted = aliTab[orderTab, ]

#setwd("C://Users/Michael//Documents//OVC Project")
#write.table(aliTabSorted, "fibroblast differentially expressed genes sorted by abs(u-stat).xls", sep="\t")

#variable data is loaded into is xPreGeneSig
#setwd("C://Users/Michael//Documents//OVC Project//R Code")
#load("cd133+sig.RData")
#load("fibroblastSig.RData")

#70/1000 with better p value/D index with 50 gene sig, try 30 to get < 50/1000
#log rank test p value is decent at 30 and bad (0.24) at 25, 35 seems to be the best trade off (better D index p value actually)
#35 seems best for fib sig
#10 seems best for cd133+sig
#something is messed up in the code!!! same numGenesInSig gives different results after a while
#pvalcomd = 0.00188

#xOrd$geneNames = xOrd$geneNames[sigRank]

#Top 40 genes
#[1] "SERPINB2" "TSPAN8"   "DHRS2"    "GJB2"     "EGFL6"    "MMP11"    "PEG3AS"   "C7"       "PEG3"     "MFAP5"   
#[11] "SFRP2"    "VCAN"     "COL11A1"  "D4S234E"  "SAA1"     "SLC12A8"  "TMEM158"  "FAP"      "MMP1"     "GREM1"   
#[21] "MATN2"    "EPYC"     "RASL11B"  "CXCL6"    "PLAU"     "PLCXD3"   "SEMA3D"   "DIO2"     "C18ORF51" "ABLIM1"  
#[31] "TDO2"     "GSTM5"    "KRT8"     "KIF5C"    "SPRR2F"   "MEGF10"   "TPR"      "ITM2A"    "ITGB2"    "C4BPB"   

#FAP Sig Genes actually in the data (29)
#[1] "SERPINB2" "TSPAN8"   "DHRS2"    "GJB2"      "MMP11"    "C7"       "PEG3"     "MFAP5"    "SFRP2"    "VCAN"    
#[11] "COL11A1"  "SAA1"     "FAP"      "MMP1"     "MATN2"    "EPYC"     "CXCL6"    "PLAU"     "DIO2"     "ABLIM1"  
#[21] "TDO2"     "GSTM5"    "KRT8"     "KIF5C"    "SPRR2F"   "TPR"      "ITM2A"    "ITGB2"    "C4BPB"  

#Caf results good at 40, good but worse at 50, decent but worse than 40 and 50 at 30, bad at 60
numGenesInSig = 40;
#xOrd$geneNames = x[sigRank]
#fapStat = cafGroupDiff[28];  #check gene names to see its in the signature, i.e top 35.
#results worse with 30 genes
geneSig = rownames(xOrd)[1:numGenesInSig]
geneSigCoef = cafGroupDiffSig[1:numGenesInSig];


posOnly = FALSE
negOnly = FALSE
test = FALSE
if(posOnly == TRUE)
{
  geneSigCoef = geneSigCoef[which(geneSigCoef > 0)]
  geneSig = geneSig[which(geneSigCoef > 0)]
}
if(negOnly == TRUE)
{
  geneSigCoef = geneSigCoef[which(geneSigCoef < 0)]
  geneSig = geneSig[which(geneSigCoef < 0)]
}
if(test == TRUE)
{
  geneSigCoef = rep(-1, numGenesInSig)
}


#Begin survival analysis using data from MetaGX
#should install each time because eset objects are very odd,
#unlike other r variables it seems you dont store a copy in a list but store a reference
#so changes to one variable effect another unexpectedly
#install.packages("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\MetaGxOvarian_0.99.0.tar.gz", repos = NULL, type="source")

#install.packages("MetaGxOvarian_0.99.0.tar.gz", type ="source", repos = NULL)
#biocLite("MetaGxOvarian")
#install.packages("Matrix")
#install.packages("lattice")
library("Matrix")
library("lattice")
library(MetaGxOvarian)
source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
source(system.file("extdata", "createEsetList.R", package="MetaGxOvarian"))
#esets likely contains all the data in one object, TCGA is just one dataset
#No missing entrez ids
threeGroups = FALSE
censoring = TRUE

getSurvInfo <- function(survObjData, bestProbes, geneSignature, geneSigCoef, dataValsOrig, dataInfo, dataTimeToDeath, dataVitalStat)
{
  #IMPORTANT: must sort the probes and then use bestprobes as best probes indices are obtained
  #after all the data is sorted
  dataInfoEntrezGene.ID = as.numeric(dataInfo$EntrezGene.ID)
  sortInd = sort(dataInfoEntrezGene.ID, decreasing = FALSE, index.return=TRUE)$ix;
  
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[sortInd];
  dataVals = dataValsOrig[sortInd, ];
  dataInfoProbeset = as.vector(dataInfo$probeset[sortInd])
  dataInfoGene = as.vector(dataInfo$gene[sortInd])
  
  dataVals = dataVals[bestProbes, ]
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[bestProbes]
  #print("hello")
  dataInfoProbeset = dataInfoProbeset[bestProbes]
  dataInfoGene = dataInfoGene[bestProbes]

  geneSigInds = c()
  geneSigEntrez = c()
  geneSigNames = c()
  dataGeneSigInds = c()
  genesFoundInd = c()
  for(i in 1:length(geneSignature))
  {
    geneSigInd = which(xOrig$probeIds == geneSignature[i]);
    geneSigInds = c(geneSigInds, geneSigInd)
    geneSigEntrez = c(geneSigEntrez, xOrig$entrezIds[geneSigInd])
    geneSigNames = c(geneSigNames, xOrig$geneNames[geneSigInd])
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
    if(length(which(dataInfoEntrezGene.ID == geneSigEntrez[i])) > 0)
      genesFoundInd = c(genesFoundInd, i)
    dataGeneSigInds = c(dataGeneSigInds, which(dataInfoEntrezGene.ID == geneSigEntrez[i]));
    #print(i)
    #print(which(dataInfoEntrezGene.ID == geneSigEntrez[i]))
  }
  #geneSigInds = which(xOrig$probeIds %in% geneSig);
  geneSigEntrez = xOrig$entrezIds[geneSigInds]
  geneSigNames = geneSigNames[genesFoundInd]
  #print(geneSigNames)
  #tcgaGeneSigInds = which(tcgaInfo$EntrezGene.ID %in% geneSigEntrez);
  
  scores = NULL
  topGenes = dataInfoProbeset[dataGeneSigInds]
# print(length(topGenes))
# print(length(geneSignature))
# print(length(topGenes) > 0)
  #if none of the randomly selected genes are in the data set then move on
  if(length(topGenes) < round(length(geneSignature)/2) | length(topGenes) < 1)
  {
   return(survObjData)
  }
  
  print(paste("number of genes in data set from signature found = ", length(topGenes)))
  
  #if even one gene expression value is NA than the score is NA, for now just remove NAs at end
  # may be better to check how many of the genes are NA and than remove NA genes while getting
  # the score with the remaining genes if enough genes remain
  
  rownames(dataVals) = dataInfoProbeset
  for(i in 1:dim(dataVals)[2])
  {
    testSamp = t(dataVals[, i, drop = FALSE]);
    mysig <- cbind("probe"=topGenes, "EntrezGene.ID"=NA, "coefficient"= (as.numeric(geneSigCoef[genesFoundInd])/sum(abs(as.numeric(geneSigCoef[genesFoundInd])))))
    rownames(mysig) = topGenes;
    scores <- rbind(scores, cbind("score"=genefu::sig.score(x=mysig, testSamp, annot=NULL, do.mapping=FALSE, signed=TRUE)$score, "fold"=j))
  }
  
  scoreVals = as.numeric(scores[,1])
  
  #probably should adjust to not assume NA in time to death is the same as NA in vital stat
  #dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
  #dataVitalStat = data@phenoData@data$vital_status;
  unknown = which(is.na(dataTimeToDeath) == TRUE)
  if(length(unknown) > 0)
  {
    dataTimeToDeath = dataTimeToDeath[-unknown];
    dataVitalStat = dataVitalStat[-unknown]
    scoreVals = scoreVals[-unknown]
  }
  
  missingGenes = which(is.na(scoreVals) == TRUE)
  if(length(missingGenes > 0))
  {
    dataTimeToDeath = dataTimeToDeath[-missingGenes];
    dataVitalStat = dataVitalStat[-missingGenes]
    scoreVals = scoreVals[-missingGenes] 
  }
  
  threeGroups = FALSE
  if(threeGroups == TRUE)
  {
    sortInd = sort(scoreVals, decreasing = TRUE, index.return=TRUE)$ix;
    thirdData = round(length(scoreVals)/3)
    dataGroups = matrix(0, length(scoreVals))
    dataGroups[sortInd[1:thirdData]] = 0;
    dataGroups[sortInd[thirdData:(2*thirdData)]] = 1;
    dataGroups[sortInd[(2*thirdData + 1):length(sortInd)]] = 2;
  }
  if(threeGroups == FALSE)
  {
    scoreMed = median(scoreVals)
    dataGroups = as.integer((scoreVals <= scoreMed));
  }
  #print(sum(is.na(dataGroups)))

  dataVitalStat = as.character(dataVitalStat)
  dataVitalStat[dataVitalStat == "living"] = 0
  dataVitalStat[dataVitalStat == "deceased"] = 1
  dataVitalStat = as.numeric(dataVitalStat)
  
  scoreValsOrig = scoreVals
  quantVals = quantile(scoreVals, c(.025, .975))
  oldLow = quantVals[1]
  oldHigh = quantVals[2]
  newLow = -1
  newHigh = 1
  scoreVals = newLow*(1 - (scoreVals - oldLow)/(oldHigh - oldLow)) + newHigh*((scoreVals - oldLow)/(oldHigh - oldLow))
  
  dataMatrix = matrix(c(dataTimeToDeath, dataVitalStat, dataGroups, scoreVals, scoreValsOrig), nrow = length(dataVitalStat), ncol = 5)
  survObjData[[1]] = rbind(survObjData[[1]], dataMatrix)
  
  dind <- D.index(x=scoreVals, surv.time=dataTimeToDeath, surv.event=dataVitalStat)
  
  if(length(survObjData) == 1)
  {
    survObjData[[2]] = c(dind$d.index)
    survObjData[[3]] = c(dind$se)
    survObjData[[4]] = c(dind$p.value)
  }
  if(length(survObjData) > 1)
  {
    survObjData[[2]] = c(survObjData[[2]], dind$d.index)
    survObjData[[3]] = c(survObjData[[3]], dind$se)
    survObjData[[4]] = c(survObjData[[4]] ,dind$p.value)
  }
  
  return(survObjData)
}

#MetaGx Data that is not missing vital status or days to death
ovDataListOrig = list()
ovDataListBestProbes = list()
ovDataListOrig[[1]] = TCGA
ovDataListOrig[[2]] = E.MTAB.386
ovDataListOrig[[3]] = GSE13876
ovDataListOrig[[4]] = GSE14764
ovDataListOrig[[5]] = GSE17260
ovDataListOrig[[6]] = GSE18520
#many patients (42) with almost no expression values for 90% of probes in just GSE19829
#Solution: remove patients with NA scores. In future may want to just remove probes and get score
ovDataListOrig[[7]] = GSE19829
ovDataListOrig[[8]] = GSE26193
ovDataListOrig[[9]] = GSE26712
ovDataListOrig[[10]] = GSE30009
ovDataListOrig[[11]] = GSE30161
ovDataListOrig[[12]] = GSE32062
ovDataListOrig[[13]] = GSE32063
ovDataListOrig[[14]] = GSE49997
ovDataListOrig[[15]] = GSE51088
ovDataListOrig[[16]] = GSE8842
ovDataListOrig[[17]] = GSE9891
ovDataListOrig[[18]] = PMID17290060
ovDataListOrig[[19]] = PMID19318476

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\R Code");
source("probeSelection.R");
source("survivalStatsAndPlot.R");

ovDataList = ovDataListOrig;
#below works the first time for some i but not the second, changing original arrays permenantly somehow? See eset weirdness note
#checked that running below multiple times doesnt efect curves (as should be the case)
#Fix, seeing multiple entrez ids in getsurvinfo function
for(i in 1:length(ovDataList))
{
  print(i)
  #data = ovDataListOrig[[i]];  still changing permenantly :s
  ovDataListBestProbes[[i]] = getBestProbes(ovDataList[[i]])
}

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Subtyping")
source("verhaakSubtypeFunction.R")
supplementary.data.sheet7 <- read.xls("JCI65833sd1.xls", sheet=7, skip=1,  perl = "C:\\Perl64\\bin\\perl.exe")
supplementary.data.sheet1 <- read.xls("JCI65833sd1.xls", skip=1,  perl = "C:\\Perl64\\bin\\perl.exe")

#Note that dataset 7 and 10 dont work, see subtyping function for more details
#suspect that genes arent present to generate scores for some signatures in those datasets
#ovDataList[[i]]$Verhaak.subtypes[1:5]
for(i in 1:length(ovDataList))
{
  print(i)
  if(i != 10 && i!= 7)
  ovDataList[[i]] = getVerhaakSubtypes(ovDataList[[i]], supplementary.data.sheet1, supplementary.data.sheet7)
}

subtyping = FALSE
threeGroups = FALSE
tcgaOnly = FALSE
#using actual gene signature

if(subtyping == FALSE)
{
  survInfo = list()
  survInfo[[1]] = as.data.frame(NULL)
  iEnd = length(ovDataList)
  if(tcgaOnly == TRUE)
    iEnd = 1  #TCGA data is first 
  for(i in 1:iEnd)
  {
    print(i)
    data = ovDataList[[i]];
    dataVals = data@assayData$exprs;
    dataInfo = data@featureData@data;
    dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
    dataVitalStat = data@phenoData@data$vital_status;
    
    survInfo = getSurvInfo(survInfo, ovDataListBestProbes[[i]], geneSig, geneSigCoef,  dataVals, dataInfo, dataTimeToDeath, dataVitalStat)
  }
  #D ind combined 1.145, 14/1000 random signatures better pval = 0.014
  survStats = getSurvivalStats(survInfo, TRUE, threeGroups, censoring)
  pvalCombD = 2*pnorm(-abs(log(survStats[1])/survStats[2]))
}
if(subtyping == TRUE)
{
  survInfoImr = list()
  survInfoImr[[1]] = as.data.frame(NULL)
  survInfoPro = list()
  survInfoPro[[1]] = as.data.frame(NULL)
  survInfoDif = list()
  survInfoDif[[1]] = as.data.frame(NULL)
  survInfoMes = list()
  survInfoMes[[1]] = as.data.frame(NULL)
  
  imr = list()
  imr[[1]] = as.data.frame(NULL)
  pro = list()
  pro[[1]] = as.data.frame(NULL)
  dif = list()
  dif[[1]] = as.data.frame(NULL)
  mes = list()
  mes[[1]] = as.data.frame(NULL)
  
  iEnd = length(ovDataList)
  if(tcgaOnly == TRUE)
    iEnd = 1  #TCGA data is first 
  for(i in 1:iEnd)
  {
    if(i != 10 && i!= 7)
    {
      subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes);
      
      data = ovDataList[[i]]
      dataVals = data@assayData$exprs;
      dataInfo = data@featureData@data;
      dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
      dataVitalStat = data@phenoData@data$vital_status;
      
      imrInds = which(subtypes == "IMR")
      proInds = which(subtypes == "PRO")
      difInds = which(subtypes == "DIF")
      mesInds = which(subtypes == "MES")
      
      survInfoImr = getSurvInfo(survInfoImr, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, imrInds], dataInfo, dataTimeToDeath[imrInds], dataVitalStat[imrInds])
      survInfoPro = getSurvInfo(survInfoPro, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, proInds], dataInfo, dataTimeToDeath[proInds], dataVitalStat[proInds])
      survInfoDif = getSurvInfo(survInfoDif, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, difInds], dataInfo, dataTimeToDeath[difInds], dataVitalStat[difInds])
      survInfoMes = getSurvInfo(survInfoMes, ovDataListBestProbes[[i]], geneSig, geneSigCoef, dataVals[, mesInds], dataInfo, dataTimeToDeath[mesInds], dataVitalStat[mesInds])
      
      imr[[i]] = survInfoImr[[1]][4][,1]
      pro[[i]] = survInfoPro[[1]][4][,1]
      dif[[i]] = survInfoDif[[1]][4][,1]
      mes[[i]] = survInfoMes[[1]][4][,1]
      
      print(i) 
    }
  }  
  survStatsImr = getSurvivalStats(survInfoImr, TRUE, threeGroups, censoring)
  survStatsPro = getSurvivalStats(survInfoPro, TRUE, threeGroups, censoring)
  survStatsDif = getSurvivalStats(survInfoDif, TRUE, threeGroups, censoring)
  survStatsMes = getSurvivalStats(survInfoMes, TRUE, threeGroups, censoring)
  
  boxplot(survInfoImr[[1]][4][,1], survInfoPro[[1]][4][,1], survInfoDif[[1]][4][,1], survInfoMes[[1]][4][,1], names = c("Immunoreactive", "Proliferative", "Differentiated", "Messenchymal"), ylab = "Sigscores")
  compSubtypes = TRUE
}



if(compSubtypes == TRUE)
{
  compList = list()
  compList[[1]] = as.data.frame(NULL)
  sizeList = list()
  sizeList[[1]] = as.data.frame(NULL)
  iEnd = length(ovDataList)
  
  for(i in 1:iEnd)
  {
    if(i != 10 && i!= 7)
    {
      classes = as.factor(c(rep("imr", length(imr[[i]])), rep("pro", length(pro[[i]])), rep("dif", length(dif[[i]])), rep("mes", length(mes[[i]]))))
      compMat = as.data.frame(NULL)
      sizeMat = as.data.frame(NULL)
      wilResGr = pairwise.wilcox.test(c(imr[[i]], pro[[i]], dif[[i]], mes[[i]]), classes, alternative = "greater", p.adjust.method = "fdr")
      wilResLe = pairwise.wilcox.test(c(imr[[i]], pro[[i]], dif[[i]], mes[[i]]), classes, alternative = "less", p.adjust.method = "fdr")
      
     compMat = rbind(compMat, c(wilcox.test(imr[[i]], pro[[i]], alternative = "less")[3], wilcox.test(imr[[i]], dif[[i]], alternative = "less")[3], wilcox.test(imr[[i]], mes[[i]], alternative = "less")[3]))  
     compMat = unname(compMat)
     compMat = rbind(compMat, c(wilcox.test(pro[[i]], imr[[i]], alternative = "less")[3], wilcox.test(pro[[i]], dif[[i]], alternative = "less")[3], wilcox.test(pro[[i]], mes[[i]], alternative = "less")[3])) 
     compMat = unname(compMat)
     compMat = rbind(compMat, c(wilcox.test(dif[[i]], imr[[i]], alternative = "less")[3], wilcox.test(dif[[i]], pro[[i]], alternative = "less")[3], wilcox.test(dif[[i]], mes[[i]], alternative = "less")[3]))  
     compMat = unname(compMat)
     compMat = rbind(compMat, c(wilcox.test(mes[[i]], imr[[i]], alternative = "less")[3], wilcox.test(mes[[i]], pro[[i]], alternative = "less")[3], wilcox.test(mes[[i]], dif[[i]], alternative = "less")[3]))  

     sizeMat = rbind(sizeMat, c(length(imr[[i]]) + length(pro[[i]]), length(imr[[i]]) + length(dif[[i]]), length(imr[[i]]) + length(mes[[i]])))
     sizeMat = rbind(sizeMat, c(length(pro[[i]]) + length(imr[[i]]), length(pro[[i]]) + length(dif[[i]]), length(pro[[i]]) + length(mes[[i]])))
     sizeMat = rbind(sizeMat, c(length(dif[[i]]) + length(imr[[i]]), length(dif[[i]]) + length(pro[[i]]), length(dif[[i]]) + length(mes[[i]])))
     sizeMat = rbind(sizeMat, c(length(mes[[i]]) + length(imr[[i]]), length(mes[[i]]) + length(pro[[i]]), length(mes[[i]]) + length(dif[[i]])))
     
     compList[[i]] = compMat
     sizeList[[i]] = sizeMat
    }
    #remove null elements since subtyping function doesnt work on the 2 datasets
    compList[[7]] = NULL
    compList[[9]] = NULL
    sizeList[[7]] = NULL
    sizeList[[9]] = NULL
  }
  finalPval = c(NULL)
  for(i in 1:dim(compList[[1]])[1])
  {
    for(j in 1:dim(compList[[1]])[2])
    {
      pval = c(NULL)
      weights = c(NULL)
      for(k in 1:length(compList))
      {
        pval = c(pval, compList[[k]][i,j])
        weights = c(weights, sizeList[[k]][i,j])
      }
      finalPval = c(finalPval, combine.test(p=pval, weight = weights, method="z.transform"))
    }
  }
  names(finalPval)= c("imr - pro", "imr - diff", "imr - mes", "pro - imr", "pro - diff", "pro - mess", "diff - imr", "dif - pro", "diff - mes", "mes - imr", "mes - pro", "mes - diff")
}
#check this works, weights appear not to factor in for fisher but do for z.transform, choice of method seems to make difference?


#testing random gene signatures
randStats = as.data.frame(NULL)
geneList = list()
#getSurvInfo <- function(survObjData, bestProbes, geneSignature, geneSigCoef, dataValsOrig, dataInfo, dataTimeToDeath, dataVitalStat)

numTests = 1000;
for(i in 1:numTests)
{
  print(i)
  randInds = sample(1:dim(xOrd)[1], numGenesInSig, replace = T)
  geneSigRand = rownames(xOrd)[randInds]
  #geneSig = rownames(xOrd)[1:numGenesInSig]
  geneSigCoef = cafGroupDiffSig[randInds];
  geneList[[i]] = geneSigRand
  survInfo = list()
  survInfo[[1]] = as.data.frame(NULL)
  
  iEnd = length(ovDataList)
  if(tcgaOnly == TRUE)
    iEnd = 1  #TCGA data is first 
  for(j in 1:iEnd)
  {
    data = ovDataList[[j]];
    dataVals = data@assayData$exprs;
    dataInfo = data@featureData@data;
    dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
    dataVitalStat = data@phenoData@data$vital_status;

    survInfo = getSurvInfo(survInfo, ovDataListBestProbes[[j]], geneSigRand, geneSigCoef, dataVals, dataInfo, dataTimeToDeath, dataVitalStat)
    #survInfo = getSurvInfo(survInfo, ovDataListBestProbes[[i]], geneSigRand, geneSigCoef, dataVals, dataInfo, dataTimeToDeath, dataVitalStat)
    
    #survInfo = getSurvInfo(survInfo, ovDataListBestProbes[[i]], geneSig, geneSigCoef,  dataVals, dataInfo, dataTimeToDeath, dataVitalStat)
    
  }
  
  randStats = rbind(randStats, getSurvivalStats(survInfo, FALSE, FALSE, censoring))
}
#combined D, combined D se, log rank P, all data D p value, all data D value 
#stats = c(dindComb$estimate, dindComb$se, (1 - pchisq(bb$chisq, 1)), dind$p.value, dind$d.index)
#pval = sum(randStats[,1] >= survStats[1])/length(randStats[,1])
#randstats 1 is exp(coef), getting p value of combined D estimate, stat centered about 0.
pvalCombD = sum(2*pnorm(-abs(log(randStats[,1])/randStats[,2])) <= 2*pnorm(-abs(log(survStats[1])/survStats[2])))/length(randStats[,1])
pvalAllDataNoCombD = sum(randStats[,4] < survStats[4])/length(randStats[,1])

#testing random gene signatures on subtypes
randStatsImr = as.data.frame(NULL)
randStatsPro = as.data.frame(NULL)
randStatsDif = as.data.frame(NULL)
randStatsMes = as.data.frame(NULL)
geneListSub = list()

numTests = 20;
for(i in 1:numTests)
{
  print(i)
  randInds = sample(1:dim(xOrd)[1], numGenesInSig, replace = T)
  geneSigRand = rownames(xOrd)[randInds]
  geneListSub[[i]] = geneSigRand
  survInfo = list()
  survInfo[[1]] = as.data.frame(NULL)
  iEnd = length(ovDataList)
  if(tcgaOnly == TRUE)
    iEnd = 1  #TCGA data is first 
  for(j in 1:iEnd)
  {
    if(j != 10 && j!= 7)
    {
      subtypes = as.vector(ovDataList[[i]]$Verhaak.subtypes);
      
      data = ovDataList[[i]]
      dataVals = data@assayData$exprs;
      dataInfo = data@featureData@data;
      dataTimeToDeath = data@phenoData@data$days_to_death/365.25;
      dataVitalStat = data@phenoData@data$vital_status;
      
      #if(censoring == TRUE)
      #{
      #  missingPat = which(is.na(dataTimeToDeath))
      #  if(length(missingPat) > 0)
      #  {
      #    dataTimeToDeath = dataTimeToDeath[-missingPat]
      #    dataVitalStat = dataVitalStat[-missingPat]
      #    dataVals = dataVals[, -missingPat]
      #    subtypes = subtypes[-missingPat]
      #  }
      #  censoredDat = censor.time(dataTimeToDeath, dataVitalStat, time.cens = 10)
      #  dataTimeToDeath = censoredDat$surv.time.cens
      #  dataVitalStat = censoredDat$surv.event.cens
      #  #patients followed up after time.cens get vital status of NA
      #  missingPat = which(is.na(dataVitalStat))
      #  if(length(missingPat) > 0)
      #  {
      #    dataTimeToDeath = dataTimeToDeath[-missingPat]
      #    dataVitalStat = dataVitalStat[-missingPat]
      #    dataVals = dataVals[, -missingPat]
      #    subtypes = subtypes[-missingPat]
      #  }
      #}
      imrInds = which(subtypes == "IMR")
      proInds = which(subtypes == "PRO")
      difInds = which(subtypes == "DIF")
      mesInds = which(subtypes == "MES")
      
      
      survInfoImr = getSurvInfo(survInfoImr, geneSigRand, dataVals[, imrInds], dataInfo, dataTimeToDeath[imrInds], dataVitalStat[imrInds])
      survInfoPro = getSurvInfo(survInfoPro, geneSigRand, dataVals[, proInds], dataInfo, dataTimeToDeath[proInds], dataVitalStat[proInds])
      survInfoDif = getSurvInfo(survInfoDif, geneSigRand, dataVals[, difInds], dataInfo, dataTimeToDeath[difInds], dataVitalStat[difInds])
      survInfoMes = getSurvInfo(survInfoMes, geneSigRand, dataVals[, mesInds], dataInfo, dataTimeToDeath[mesInds], dataVitalStat[mesInds])
    }
  }  
  randStatsImr = rbind(randStatsImr, getSurvivalStats(survInfoImr, FALSE, FALSE, censoring))
  randStatsPro = rbind(randStatsPro, getSurvivalStats(survInfoPro, FALSE, FALSE, censoring))
  randStatsDif = rbind(randStatsDif, getSurvivalStats(survInfoDif, FALSE, FALSE, censoring))
  randStatsMes = rbind(randStatsMes, getSurvivalStats(survInfoMes, FALSE, FALSE, censoring))
}


#done


if(regressionMod == TRUE)
{
  geneNumVec = seq(from=minGeneNumb, to=maxGeneNumb, by=geneInc);
  no_cores <- detectCores() - 2;
  cl <- makeCluster(no_cores)
  clusterExport(cl, list("onefoldcv", "myRankingOneFold"))
  oneFoldResults = clusterApplyLB(cl, geneNumVec, function(x, y, z) onefoldcv(x, y, z), y = xGenReg, z = fapGroup)
  stopCluster(cl);
  
  for(i in 1:length(oneFoldResults))
  {
    oneFoldResult = oneFoldResults[[i]][, 1];
    corVsGenes = c(corVsGenes, cor(oneFoldResult, fapGroup));
  }
  plot(geneNumVec, corVsGenes, xlab = "number of genes used in sig score", ylab = "sample's sig score and fapGroup correlation");
  lines(geneNumVec, corVeGenes)
}

getRandStats = function(geneSignature, ovarianDataList, randStatsClust)
{
  survInfo = as.data.frame(NULL)
  for(j in 1:length(ovarianDataList))
  {
    survInfo = getSurvInfo(ovarianDataList[[j]], survInfo, geneSignature)
  }
  survInfoOrig = survInfo
  #D ind = 1.139343 with se of 0.03926089
  #should look into what is causing occasional NA value
  dindCombRand = combine.est(survInfoOrig[[2]][!is.na(survInfoOrig[[2]])], survInfoOrig[[3]][!is.na(survInfoOrig[[3]])])
  survInfo = survInfo[[1]]
  colnames(survInfo) = c("timesToDeath", "vitalStats", "groups", "scores")
  survInfo$groups = as.integer((survInfo$scores <= median(survInfo$scores)));
  survObj <- survfit(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
  bb <- survdiff(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups,rho=0)
  dind <- D.index(x=survInfo$scores, surv.time=survInfo$timesToDeath, surv.event=survInfo$vitalStat)
  
  randStatsClust = rbind(randStatsClust, c(dindCombRand$estimate, dindCombRand$se, (1 - pchisq(bb$chisq, 1)), dind$p.value, dind$d.index))
  
  return(randStatsClust)
}

library(GSA)
library(piano)
setwd("C:\\Users\\Michael\\Documents\\OVC Project")
gSets <- GSA.read.gmt("c5.all.v5.1.entrez.gmt.txt")
dfgSets <- as.data.frame(cbind(unlist(gSets$genesets))) ## genes
dfgSNames <- as.data.frame(cbind(unlist(gSets$geneset.names)))
listS <- lapply(gSets$genesets,length)
gTogs <- data.frame(dfgSets$V1,rep(dfgSNames$V1,listS))
names(gTogs) <- c("V1","V2")
gTogs <- gTogs[gTogs$V1!="",]
a <- loadGSC(gTogs)

 
stats.str <- as.vector(cafGroupDiff) # here it will be your t-stat
names(stats.str) <- as.vector(x$entrezIds) # give names to each element
#gsares.str <- runGSA(geneLevelStats=stats.str,gsc=a,nPerm=1000,geneSetStat="gsea",adjMethod="none",gsSizeLim=c(4,Inf))
#analyze resulting GSEA object that venkatta sent as my computer ran the above too slowly

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\GSEA")
gsares.str = load("GSEA-IPA-Fib.RData")

#removed str from venkattas code inGSAsummaryTable(gsares.str)
gsasummary.str <- GSAsummaryTable(gsares)
gsasummary.str <- gsasummary.str[order(gsasummary.str[,4]),]
#fix the zeros, remove NA slash collapse, p.adjust
gsasummaryneat.str = gsasummary.str
gsasummaryneat.str[,4] = replace(gsasummaryneat.str[,4],gsasummaryneat.str[,4] == 0, 1/(1000+1));
gsasummaryneat.str[,6] = replace(gsasummaryneat.str[,6],gsasummaryneat.str[,6] == 0, 1/(1000+1));
gsasummaryneat.str[firstNA:length(gsasummary.str[,1]), 4] = gsasummaryneat.str[firstNA:length(gsasummaryneat.str[,1]), 6]
adjustedP = p.adjust(gsasummaryneat.str[,4], method = "fdr", n = length(gsasummaryneat.str[,4]));
gsasummaryneat.str[,5] = adjustedP
gsasummaryneat.str = gsasummaryneat.str[,1:5]
colnames(gsasummaryneat.str) = c("names", "Genes (Total)", "Stat (dist.dir)", "p value", "adjusted p value")
gsasummaryneat.str <- gsasummaryneat.str[order(gsasummaryneat.str[,5]),]

write.table(gsasummaryneat.str, "gseaTableNeat.xls", sep="\t")

#CMAP

cmapGsea = FALSE
cmapGwc = TRUE
if(cmapGsea == TRUE)
{
geneSigInds = c()
geneSigEntrez = c()
geneSigNames = c()
dataGeneSigInds = c()
for(i in 1:length(geneSig))
{
  geneSigInd = which(xOrig$probeIds == geneSig[i]);
  geneSigInds = c(geneSigInds, geneSigInd)
  geneSigEntrez = c(geneSigEntrez, xOrig$entrezIds[geneSigInd])
  geneSigNames = c(geneSigNames, xOrig$geneNames[geneSigInd])
  #print(which(dataInfo$EntrezGene.ID == geneSigEntrez[i]))
}
#geneSigInds = which(xOrig$probeIds %in% geneSig);
geneSigEntrez = xOrig$entrezIds[geneSigInds]

library(illuminaio)
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Annotation Files");
#Note: probeId in excel file is/x data is array_address_location variable of annot$gene and target id is the ilmn_gene var
#downloaded the file off of the internet, correct for the bead chip used during the experiments
annot = readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx");

excelTargetIds = annot$probes$ILMN_Gene;
excelProbeIds = annot$probes$Array_Address_Id;
excelEntrezIds = annot$probes$Entrez_Gene_ID;
excelGeneNames = annot$probes$Symbol;

y <- illuminaHumanv4ENSEMBL
# Get the manufacturer identifiers that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(y)
# Convert to a list
yy <- as.list(y[mapped_genes])

##Can do the mapping from array address to illumina ID using a revmap

z <- revmap(illuminaHumanv4ARRAYADDRESS)
mapped_probes <- mappedkeys(z)
# Convert to a list
zz <- as.list(z[mapped_probes])

# Convert the object to a list
ww <- as.list(illuminaHumanv4ALIAS2PROBE)


yy[1]
zz[1]

#get ensemble ids for cmap
ensembleVec = c();
genesPresent = c()
for (i in 1:length(geneSigEntrez))
{
  print(i)
  #multiple indexEnt for same entrez lead to same ensemble id
  indexEnt = which(excelEntrezIds == geneSigEntrez[i])
  probeId = excelProbeIds[indexEnt[1]]
  #if one uses the zz line instead of the ww line the results differ by 2 ensemble gene ids?
  #mappings should be equivalent from documentation, likely just multiple ensemble ids per gene issues
  #ilmnId = zz[which(names(zz) == probeId)][[1]]
  ilmnId = ww[which(names(ww) == excelGeneNames[indexEnt][1])][[1]][1]
  yyInd = which(names(yy) == ilmnId);
  if(length(yyInd) > 0)
  {
    #one instance had multiple ensemble Ids (2) for the same entrez gene id, take first one??
    ensembleId = yy[yyInd][[1]][1]
    print(ensembleId)
    ensembleId = paste0(ensembleId, "_at")
    ensembleVec = c(ensembleVec, ensembleId)
    genesPresent = c(genesPresent, i)
  }
}



library(PharmacoGx)
require(xtable)
#below two lines are time consuming, so just use result as its data independent
#data(CMAPsmall)
#drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna")
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
load("cmap_sig_rna.RData")
geneSigCmap = as.data.frame(NULL)
geneSigCmap[1:length(ensembleVec),1] = ensembleVec
rownames(geneSigCmap) = ensembleVec
geneSigCmap[,2] = cafGroupDiff[genesPresent]
colnames(geneSigCmap) = c("feature", "direction")
#flip the directions as we query for drugs that result in the reversing the signature
geneSigCmap[,2] = replace(geneSigCmap[,2], geneSigCmap[,2] > 0, 1)
geneSigCmap[,2] = replace(geneSigCmap[,2], geneSigCmap[,2] < 0, -1)
#geneSigCmap is formatted like HDAC_genes which is used in the pharmaco documentation example
data(HDAC_genes)

#remove 1:2 to do all the drugs
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
#adjustedP = p.adjust(resFix[,2], method = "fdr", n = length(resFix[,2]));
#resFix[,3] = NULL
#resFix[1:length(resFix[,2]),3] = adjustedP
#resFix <- resFix[order(resFix[,3]),]
#colnames(resFix) = c("Connectivity", "P Value", "FDR Adjusted P Value")

write.table(resFix, "ovarian fibroblast project cmap clean no inverse.xls", sep="\t")
}
if(cmapGwc == TRUE)
{
  library(illuminaio)
  library(PharmacoGx)
  setwd("C:\\Users\\Michael\\Documents\\PMH Research\\BenNeelCmap\\Data\\cmap pset")
  load("CMAP.RData")
  drugNotes = drugInfo(CMAP)
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
  load("cmap_sig_rna.RData")
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Annotation Files");
  #Note: probeId in excel file is/x data is array_address_location variable of annot$gene and target id is the ilmn_gene var
  #downloaded the file off of the internet, correct for the bead chip used during the experiments
  annot = readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx");
  
  excelTargetIds = annot$probes$ILMN_Gene;
  excelProbeIds = annot$probes$Array_Address_Id;
  excelEntrezIds = annot$probes$Entrez_Gene_ID;
  excelGeneNames = annot$probes$Symbol;
  
  y <- illuminaHumanv4ENSEMBL
  # Get the manufacturer identifiers that are mapped to an Ensembl ID
  mapped_genes <- mappedkeys(y)
  # Convert to a list
  yy <- as.list(y[mapped_genes])
  
  ##Can do the mapping from array address to illumina ID using a revmap
  
  z <- revmap(illuminaHumanv4ARRAYADDRESS)
  mapped_probes <- mappedkeys(z)
  # Convert to a list
  zz <- as.list(z[mapped_probes])
  
  # Convert the object to a list
  ww <- as.list(illuminaHumanv4ALIAS2PROBE)
  
  
  yy[1]
  zz[1]
  
  #get ensemble ids for cmap
  ensembleVec = c();
  genesPresent = c()
  tStatVec = c()
  pValVec = c()
  for (i in 1:dim(x$E)[1])
  {
    if(i%%1000 == 0)
      print(i)
    #multiple indexEnt for same entrez lead to same ensemble id
    indexEnt = which(excelEntrezIds == x$entrezIds[i])
    probeId = excelProbeIds[indexEnt[1]]
    #if one uses the zz line instead of the ww line the results differ by 2 ensemble gene ids?
    #mappings should be equivalent from documentation, likely just multiple ensemble ids per gene issues
    #ilmnId = zz[which(names(zz) == probeId)][[1]]
    wwInd = which(names(ww) == excelGeneNames[indexEnt][1])
    if(length(wwInd) > 0)
    {
      ilmnId = ww[wwInd][[1]][1]
      yyInd = which(names(yy) == ilmnId);
      if(length(yyInd) > 0)
      {
        #one instance had multiple ensemble Ids (2) for the same entrez gene id, take first one??
        ensembleId = yy[yyInd][[1]][1]
        ensembleId = paste0(ensembleId, "_at")
        ensembleVec = c(ensembleVec, ensembleId)
        tStatVec = c(tStatVec, cafGroupDiff[i])
        pValVec = c(pValVec, cafPValsOneTwo[i])
        genesPresent = c(genesPresent, i)
      } 
    }
  }
  
  #deal with the 500 or so repeat ensemble ids
  repeatEns = duplicated(ensembleVec)
  repeatEns = which(repeatEns == TRUE)
  
  #keep highest abs(t-stat) probes of repeat ensemble id probes
  rmList = c()
  for(i in 1:length(repeatEns))
  {
    repInd = repeatEns[i]
    ensId = ensembleVec[repInd]
    probes = which(ensembleVec == ensId)
    tVals = abs(cafGroupDiff[probes])
    keeper = which(tVals == max(tVals))
    probesRm = probes[-keeper]
    rmList = c(rmList, probesRm)
  }
  #above code checks more then once per repeat id, remove non unique
  rmList = unique(rmList)
  ensembleVec = ensembleVec[-rmList]
  tStatVec = tStatVec[-rmList]
  pValVec = pValVec[-rmList]
  
  
  genesCmap = as.data.frame(NULL)
  #IMPORTANT, clarify this makes sense, also note that it just changes the sign of the results
  tStatVec = -1*tStatVec
  genesCmap[1:length(ensembleVec),1] = tStatVec
  genesCmap[1:length(ensembleVec),2] = pValVec
  rownames(genesCmap) = ensembleVec
  colnames(genesCmap) = c("tstat", "pvalue")
  
  cmapTab = as.data.frame(NULL)
  numDrugs = dim(drug.perturbation)[2]
  numbPerms = 20000
  drugNames = drug.perturbation[1, ,1][1:numDrugs]
  drugNames = names(drugNames)
  for(i in 1:numDrugs)
  {

    #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
    #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
    #connectivity score is high when either
    # drug t-stat high (drug raised expression of gene) & data t-stat high 
    #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
    #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
    #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
    
    drugStats = drug.perturbation[,i ,c("tstat", "pvalue")]
    dataCor <- connectivityScore(genesCmap, drugStats,method="gwc", gwc.method="spearman", nperm = numbPerms)
    drugInd = which(rownames(drugNotes) == drugNames[i])
    drugInfoTab = drugNotes[drugInd, ]
    cmapTab = rbind(cmapTab, c(dataCor, as.character(drugInfoTab)), stringsAsFactors = FALSE)
    
  }
  colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
  rownames(cmapTab) = drugNames
  sort = order(cmapTab[,"score"], decreasing = TRUE)
  cmapTabSort = cmapTab[sort,]
  cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
  write.table(cmapTabSort, "GWC CMAP Results Ovarian Fibroblast Project.xls", col.names = NA, sep="\t")
  
}

randStatsCl = as.data.frame(NULL)
geneList = list()

numTests = 1;
for(k in 1:numTests)
{
  randInds = sample(1:dim(xOrd)[1], numGenesInSig, replace = T)
  geneSigRand = rownames(xOrd)[randInds]
  geneList[[k]] = geneSigRand
}
#crashed last time, reduce cores to prevent crash
#no_cores <- detectCores() - 2;
no_cores = 2
cl <- makeCluster(no_cores)
clusterExport(cl, list("getRandStats", "getSurvInfo","rowIQRs", "xOrig", "j"))
randStatsCl = clusterApplyLB(cl, geneList[[1]], function(x, y, z) getRandStats(x, y, z), y = ovDataList, z = randStatsCl)
stopCluster(cl);


#External Data Sets With Potential
#System-Wide Analysis Reveals a Complex Network of Tumor-Fibroblast Interactions Involved in Tumorigenicity
#Prognostic Gene Expression Signature of Carcinoma Associated Fibroblasts in Non-Small Cell Lung Cancer
#Molecular Characterization of Breast Carcinoma Associate fibroblasts
#Carcinoma-associated fibroblasts' transcriptomic program predicts clinical outcome in stage II/III colorectal cancer
#Vitamin D receptor expression and associated gene signature in tumor stromal fibroblasts predict clinical outcome in colorectal cancer
#Expression data from cultured human esophageal squamous cell carcinoma cell lines and cultured human fibroblasts
#remaining data sets from array express website have n < 25

TCGACopy = TCGA
ovDataList[[1]] = TCGACopy
dat = TCGA@assayData$exprs
dim(ovDataList[[1]])
dim(TCGA@assayData$exprs)
TCGA@assayData$exprs = TCGA@assayData$exprs[, -dim(TCGA@assayData$exprs)[2]]
dim(ovDataList[[1]])
dim(TCGA@assayData$exprs)

#Molecular Characterization of Breast Carcinoma Associate fibroblasts

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\External Fibroblast Data\\E-GEOD-29270.processed.1")
arrayFiles = list.files()
fibArrays = read.delim(arrayFiles[27], header = FALSE)
fibArrays = fibArrays[1:44500,]
for(i in 2:length(arrayFiles))
{
  fibArray = read.delim(arrayFiles[i], header = FALSE)
  print(i)
  fibArrays[, i+1] = fibArray[1:44500, 2]
}
fibMat = fibArrays
fibOrig = fibArrays

fibMat = as.matrix(fibMat)
rownames(fibMat) = fibMat[,1]
cornerInds = which(rownames(fibMat) == "DarkCorner")
fibMat = fibMat[-cornerInds, ]
fibMat = fibMat[, -1]
fibMat = fibMat[-1, ]
class(fibMat) = "numeric"
fibMat = normalizeBetweenArrays(fibMat,method="quantile")
dissimilarity <- 1 - cor(fibMat);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Clustering Analysis More Samples Different Data", xlab="");
#heatmap(fibMat, Rowv=NA, Colv=as.dendrogram(hc))

#To Do
#acquire results for both signatures and put them on the wiki
#start creating a meta gx function that does survival analysis on genes one sends in, eventualy ass sweave
#Cmap on Ben Neel's data
#analyze prognostic value of all the genes from the yale project individually on metagx data
#get data from kevin and begin making pharmao object


