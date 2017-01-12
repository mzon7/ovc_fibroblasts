
# import limma package and lumi package that was mentioned in analysis report
# seems a time consuming step when run manually, fix to only be done if need be (aka dont update all always)
source("https://bioconductor.org/biocLite.R");
biocLite("limma", suppressUpdates=TRUE);
library(limma);
biocLite("lumi", suppressUpdates=TRUE);
library(lumi);
source("https://bioconductor.org/biocLite.R")
biocLite("genefu")
biocLite("illuminaio", supressUpdates = TRUE);
library(illuminaio)
biocLite("org.Hs.eg.db",  supressUpdates = TRUE);
library(matrixStats)
library(org.Hs.eg.db)
library(parallel)



#https://www.biostars.org/p/109248/ 
#gotta do above o try and get entrez gene ids
biocLite("illuminaHumanv4.db", supressUpdates = TRUE);
library(illuminaHumanv4.db)

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


# use limma function to load in files generated from using GenomeStudio on Illumina data
# probe profile.txt is main probe summary profile
# control probe profile.txt is name of file contaiing profiles for the control probes
# can use read.ilmn.targets if one has multiple probe summary profiles and the samples are summarized in a target frame
# if one has the control probe profiles data can be favorably background corrected and normalized using neqc or nec func
# otherwis can use stuff used on other single channel platforms


# GeneSpring raw data folder works but GenomeStudio folder sample...finalreport file doesnt work
#removed header from genespring file and it still works, so error has to do with data
#Error in `rownames<-`(`*tmp*`, value = list(ProbeID = c(6450255L, 2570615L,  : 
#length of 'dimnames' [1] not equal to array extent

f = blockCode(){
  
  #set directory below to one containing probe summary profiles and control probe profiles
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data and QC GenomeStudio");
  
  useCtrlProbe = 1;
  
  if(useCtrlProbe == 1){
    #z <- read.ilmn("sample_probe_nonnorm_FinalReportNewHeader.txt", ctrlfiles="control_probe_nonnorm_FinalReport.txt");  
  } else{
    #z <- read.ilmn("sample_probe_nonnorm_FinalReportNewHeader.txt");  
  }
  
}

#should use read.idat function to get data matrix from scratch later using natalie's raw idat files



#below has entrex gene ids and gene names for dataframe x used below in the analysis, cant get it to load with read.ilmn though
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Natalie Files");
rawData = read.delim("sample_probe_nonnorm_minusbackground_FinalReportNoHeaderNotepad.txt");
geneNameColInd = which(grepl("TargetID",colnames(rawData)))
entrezIDColInd = which(grepl("ENTREZ_GENE_ID",colnames(rawData)))
#accessionInd = which(grepl("ACCESSION",colnames(rawData)))

#GeneSpring files seem to be working fine
#set directory below to one containing probe summary profiles from GeneSpring
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data GeneSpring");
#47300 probes
x <- read.ilmn("Ali HT-12 expression_Sample_Probe_Profile_nonnormNoHeader.txt"); 
colnames(x$E) = sampNames;
x$geneNames = as.vector(rawData[, geneNameColInd]);
x$entrezIds = as.vector(rawData[, entrezIDColInd]);
x$probeIds = rownames(x$E);
#x$acc = as.vector(rawData[, accessionInd]);

xOrig = x;
#44000
missingEntrezIds = which(is.na(x$entrezIds));
x$E = x$E[-missingEntrezIds,];
x$entrezIds = x$entrezIds[-missingEntrezIds];
x$probeIds = x$probeIds[-missingEntrezIds];
x$geneNames = x$geneNames[-missingEntrezIds];



setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Annotation Files");
#Note: probeId in excel file is/x data is array_address_location variable of annot$gene and target id is the ilmn_gene var
#downloaded the file off of the internet, correct for the bead chip used during the experiments
annot = readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx");
probeID = annot$probes$Probe_Id;
geneNames = data.frame(Gene=unlist(mget(x = probeID,envir = illuminaHumanv4SYMBOL)));

#here excel refers to the column names of file sample_probe_nonnorm_minusbackground_FinalReport.txt which have the probe Ids
#in the same order and has the same number of probes as the Ali....txt file used in the analysis, aka files are aligned
excelTargetIds = annot$probes$ILMN_Gene;
excelProbeIds = annot$probes$Array_Address_Id;
excelEntrezIds = annot$probes$Entrez_Gene_ID;

#ex first probe in excel file/
start = "^"
end = "$"
i = 1;
probeId = x$probeIds[i]
probeId = paste(start, probeId, end, sep = "")
annotIndex = which(grepl(probeId, excelProbeIds));
verifyGeneName = excelTargetIds[annotIndex];
entrezId = excelEntrezIvds[annotIndex];


z <- illuminaHumanv4ARRAYADDRESS
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(y)
# Convert to a list
zz <- as.list(y[mapped_probes])
zzVars = names(zz);

y <- illuminaHumanv4ENTREZID
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(y)
# Convert to a list
yy <- as.list(y[mapped_probes])
yyVars = names(yy);

probeIds = rownames(geneNames)
#convert gene names to entrez ids
start = "^"
end = "$"

#alternative way to get entrez id from gene name/ target id in the excel files

geneName = "HS.100261"
geneName = paste(start,geneName,end,sep="")
probeIdInds = which(grepl(geneName, geneNames$Gene));
probeIdName = probeIds[probeIdInds]
probeIdName = paste(start, probeIdName, end, sep="");

#use first one as they all have the same entrez gene id
yyInd = which(grepl(probeIdName[1], yyVars));
entrezId = yy[yyInd][[1]];


#CONCLUSION: cannot obtain certain entrez gene ids as they are missing from the annotation file provided by Illumina company

#x <- normalizeBetweenArrays(x,method="quantile")
#x = neqc(x)

# Quality Assessment (these plots are the same as the ones mentioned in the report)
# evident in report they just called plot function of lumi package to generate the various plots
boxplot(log2(x$E),range=0,ylab="log2 intensity")
plotDensity(x$E)


#should fix this to put each arrays histogram on the plot
#hist(log(x$E))

#report says analysis performed on log2 data, this transform gives the same boxplot as theirs
#x$E = log2(x$E);

#f2 = blockCode(){
#block this code as report says the data is already normalized and it does not mention background corrections

#neqc likely not appropriate for genespring data as control probe profile is not provided as is for genomestudio files
#thus the mean and variance of the nagative control probe used to normalize and bacground correct are inferred from p vals
#y <- neqc(x)  
# background correct data, i.e remove noise
# model sign as exponential fucntion and noise as normal distribution, presumably subtract the noise
#can also use background correction methods from lumi package that are tailers to illumina data
#y <- backgroundCorrect(x,method="normexp")
# normalize data to make measurements from separate array comparable
# make the largest intensity in all the arrays have the same vale, the value being the mean of all the arrays largest intensities
# repeat process for the 2nd largest, 3rd largest, and so on values in the arrays so they all have the same distribution
#y <- normalizeBetweenArrays(y,method="quantile")

#}

plus49 = which(grepl("49", sampNames));
plus133 = which(grepl("133\\+", sampNames));
minus133 = which(grepl("133-", sampNames));

plus49Genes = x[, plus49];
plus133Genes = x[, plus133];
minus133Genes = x[, minus133];

#one way anova between each cell line, not important for CAF analysis, can likely ignore

#genePVals = rep(0, dim(x$E)[1]);

#for(i in 1:dim(x$E)[1])
#{
#  plus49Gene = c(plus49Genes$E[i,]);
#  names(plus49Gene) = NULL;
#  plus133Gene = c(plus133Genes$E[i,]);
#  names(plus133Gene) = NULL;
#  minus133Gene = c(minus133Genes$E[i,]);
#  names(minus133Gene) = NULL;

#  geneData = c(plus49Gene, plus133Gene, minus133Gene);
#  n = rep(numPatients, numbCellLines);
#  group = rep(c("49+", "133+", "133-"), n);

#tmpfn = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x));
#tapply(geneData, group, tmpfn);

#  data = data.frame(geneData = geneData, group = factor(group));
#  fit = lm(geneData ~ group, data);
#  genePVals[i] = anova(fit)$P[1];

#}

#minP = 0.05;
#adjustedGenePVals = p.adjust(genePVals, method = "BH", n = length(genePVals));
#sigGeneInds = which(adjustedGenePVals < minP);

#could quantile normalize to make all the density plots look the same? report says genespring quantile normalized the data
#but appears to be false based on the density distribution plot without quantile normalization in r?

# in report, keep probes that are > 20%tile of intensity distribution in 100 - percentile% (10 or more samples) of the samples in the groups
# so 1. break data into groups (here just the 49 e group), 2. give each gene a zero vector,
# 3. for each sample organize probes from most to least intense and each time probe <= 20%tile intensity add 1 to zero vector group
# 4. remove all probes that have an element of 3 vec >= 3

#do we agree with the logic of doing the below probe removal? check results with and without doing this
#makes sense to remove very bottom outliers as could be noise/an error. 20 % seems to high
#results nealry match after doing this, without it they look very different --> check


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

#33411 probes left at bot 5 percentile

xOrigTwo = x;

#get original data and do quantile normalization followed by per probe median normalization

#Quite similar to previous groups dendogram, extra genes --> slight difference?
#pearson centerd correlation = 1 - pearscorcoff, prev analysis used this as distance metric for clustering with average linkage rules
dissimilarity <- 1 - cor(x$E);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="");
heatmap(x$E, Rowv=NA, Colv=as.dendrogram(hc))

y <- normalizeBetweenArrays(x,method="quantile")
dissimilarity <- 1 - cor(y$E);
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation (quantile normalized data)", xlab="");

y <- normalizeBetweenArrays(x,method="quantile")
y <- backgroundCorrect(y,method="normexp")
dissimilarity <- 1 - cor(y$E);
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation (quantile normalized background corrected data)", xlab="");


#drastically different, note effect of method choice!?
clustMat = t(x$E);
hc2 <- hclust(dist(clustMat), method = "average")
plot(hc2)

#my.dist <- function(x) dist(x, method="binary")
#my.hclust <- function(d) hclust(d, method="average")
#hm <- heatmap(t(x$E), distfun=my.dist, hclustfun=my.hclust)

#have the two CAF groups



#without adjusting one gets 5038 significantly varying probes, oneway.._filtereddata.txt has 5435
#but report says adjusted BH method gives 1469 varying probes, here one has 2387 significant after adjusting
#if one quantile normalizes you get 8961 without adjusting and 5108 with adjusting, if just normalize before loop get similar results

#Note: the clustering done after anova appears to have 5000 ish genes, not 1469 as mentioned in the report
#Evidence is in hierachicalCluster_onewayANOVAcorpt05.xls, which has 5000+ probes that were used in clustering
#get same dendrogram regardless of normlization and adjusting p values. Nearly the same as prior group's

#xSig = x[sigGeneInds,];
#dissimilarity <- 1 - cor(xSig$E);  #cor default is pearson
#distance <- as.dist(dissimilarity);
#hc <- hclust(distance, method = "average");
#plot(hclust(distance), main="Dissimilarity = 1 - Correlation, Significant Genes", xlab="");
#heatmap(xSig$E, Rowv=NA, Colv=as.dendrogram(hc));

#Permitted to skip pointless one way anova part of report on 133+/- cell lines
#test 113+ vs 133- to show not many diff expressed genes?
#On to creating gene signatures to distinguish between the two subtypes

#handle multiple probes per gene, keep probe with highest iqr (venkata, not averaging?) in 49e cells

#1.order x$E, x$probeID, x$geneName, x$entrexGeneIDS by decreasing to increasing entrezgeneID
#2.run through entrezGeneId, when numbers are the same get iqr for each probe and determine probe with highest iqr
#3.record row index of probe to keep 4.move to next new entrez id probe in x and repeats 



xOrigTwo = x;
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
xGen = x;

dissimilarity <- 1 - cor(x$E);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="");
heatmap(x$E, Rowv=NA, Colv=as.dendrogram(hc))

y=x;
y <- normalizeBetweenArrays(x,method="quantile")
dissimilarity <- 1 - cor(y$E);
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation (quantile normalized data)", xlab="");

y <- normalizeBetweenArrays(x,method="quantile")
y <- backgroundCorrect(y,method="normexp")
dissimilarity <- 1 - cor(y$E);
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation (quantile normalized background corrected data)", xlab="");

#now results look extremely like the prior groups, samples at exact same spots in dendrograms


#cafGroupOneGenes = x[,c(3, 7, 23, 26, 28)];
#cafGroupTwoGenes = x[,c(5, 9, 11, 25, 32)];
#cafGroupOtherGenes = x[,c(1, 29)];

cafGroupOneGenes = x[,c(2, 4, 7, 9, 10)];
cafGroupTwoGenes = x[,c(3, 5, 6, 8, 12)];
cafGroupOtherGenes = x[,c(1, 11)];


cafPValsOneTwo = rep(0, dim(x$E)[1]);

for(i in 1:dim(x$E)[1])
{
  cafGroupOneGene = c(cafGroupOneGenes$E[i,]);
  names(cafGroupOneGene) = NULL;
  cafGroupTwoGene = c(cafGroupTwoGenes$E[i,]);
  names(cafGroupTwoGene) = NULL;
  
  geneData = c(cafGroupOneGene, cafGroupTwoGene);
  n = rep(5, 2);
  group = rep(c("Caf One", "Caf Two"), n);
  
  #tmpfn = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x));
  #tapply(geneData, group, tmpfn);
  
  data = data.frame(geneData = geneData, group = factor(group));
  fit = lm(geneData ~ group, data);
  cafPValsOneTwo[i] = anova(fit)$P[1];
  
}

#adjusted p values all above 0.05??? non adjusted has 883 differentially expressed genes at p=0.05 or 246 at p=0.01
minP = 0.01;
adjustedCafPValsOneTwo = p.adjust(cafPValsOneTwo, method = "bonferroni", n = length(cafPValsOneTwo));
sigGeneInds = which(adjustedCafPValsOneTwo < minP);
sigGeneIndsFive = which(cafPValsOneTwo < 0.05);
sigGeneIndsOne = which(cafPValsOneTwo < 0.01);

#define a gene signature
patientOutliers = c(1, 11);

xSig = x[sigGeneIndsOne, ];
xSig$E = xSig$E[ , -patientOutliers];
xSigTrans = t(xSig$E)
dissimilarity <- 1 - cor(xSigTrans);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="");
#heatmap(xSigTrans, Rowv=NA, Colv=as.dendrogram(hc))

xSig = x[sigGeneIndsFive, ];
xSig$E = xSig$E[ , -patientOutliers];
xSigTrans = t(xSig$E)
dissimilarity <- 1 - cor(xSigTrans);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation", xlab="");
#heatmap(xSigTrans, Rowv=NA, Colv=as.dendrogram(hc))

# need to email Ali to get patient response variable for making signature via regression correlation
# why not just use all genes and create model based on survival data and differentially expresed genes?
#if we think that these differentially expressed genes are corelated with survival why create 2 groups and look
#at the survival of the groups when we could predict how long they will survive based on the expression of the genes
#and then creat 3 groups, low medium and high survival and then see if were right in survival analysis


#Begin regression analysis/gene signature creation. xSig is just diff expressed genes
#xGen is all genes after normalizing, removing duplicates and low percentile. unbiased approach uses xGen

setwd("C:\\Users\\Michael\\Documents\\OVC Project\\R Code");
source("one-fold-CrossValidation.R");
source("oneFolCVRankings.R");

#xGen <- normalizeBetweenArrays(xGen,method="quantile")
xGen$E = log2(xGen$E);
xGenReg = t(xGen$E);
xSigReg = t(xsig$E);

#fapHi values for patients in same order as col names in xGenReg and xSigReg, Ali should provid NA values soon
fapHi = c(NA, 70.4, 100, 77, 81.7, 92.8, NA, 86, 58.8, 24.2, NA, NA)

naVals = which(is.na(fapHi)==TRUE);
xGenReg = xGenReg[-naVals, ];
xSigReg = xSigReg[-naVals, ];
fapHi = fapHi[-naVals];

#about 2.5 hours to run 30 loops

minGeneNumb = 1;
maxGeneNumb = 30;
geneInc = 1;
corVsGenes = c();
numGenes = c();

#if quantile normalizing go from 0.85 to 0.3 correlation but smoother looking plot
for(i in seq(from=minGeneNumb, to=maxGeneNumb, by=geneInc))
{
  print(i)
  oneFoldResult = onefoldcv(xGenReg, fapHi, i)
  corVsGenes = c(corVsGenes, cor(oneFoldResult[,1], fapHi));
  numGenes = c(numGenes, i);
}

plot(numGenes, corVsGenes, xlab = "number of genes used in sig score", ylab = "sample's sig score and fapHi% correlation");

# next steps
# 1. get raw data from natalie and analyze it, and determine how she handles multiple probes/gene 
# 2. apply the same normalization as the report
# 3. perform remaining quality checks (likely just call lumi plot function)
# 4. Perform same one way anova as was done and get about 1469 significantly varying probes
# 4. cluster analysis after anova and get same dendrogram as from alis ppt/_corrected jpeg
# 5. complete remaining things in the report, just some more anova stuff (I think)
# 5. create gene signatures to differentiate 2 groups (talk to Haibe-Kaines/Venkata about this)


# ALI/PREVIOUS ANALYST QUESTIONS
# Hi Ali, can you get me in touch with the person who previously analyzed the projects data so I could ask the following questions
# 1. Do you have the raw Illumina data that can be loaded into GeneSpring
# 2. As the Gene Expression Analysis Report does not mention background correction (only normalization) should I assume no background correction was applied>

