
# import limma package and lumi package that was mentioned in analysis report
# seems a time consuming step when run manually, fix to only be done if need be (aka dont update all always)
source("https://bioconductor.org/biocLite.R");
biocLite("limma", suppressUpdates=TRUE);
library(limma);
biocLite("lumi", suppressUpdates=TRUE)
library(lumi);

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
rawData = read.delim("sample_probe_nonnorm_minusbackground_FinalReportNoHeader.txt");
geneNameColInd = which(grepl("TargetID",colnames(rawData)))
entrezIDColInd = which(grepl("ENTREZ_GENE_ID",colnames(rawData)))

#GeneSpring files seem to be working fine
#set directory below to one containing probe summary profiles from GeneSpring
setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data GeneSpring");
x <- read.ilmn("Ali HT-12 expression_Sample_Probe_Profile_nonnormNoHeader.txt"); 
x$geneNames = as.vector(rawData[, geneNameColInd]);
x$entrezIds = as.vector(rawData[, entrezIDColInd]);
x$probeIds = rownames(x$E);

#x <- normalizeBetweenArrays(x,method="quantile")
#x = neqc(x)

# Quality Assessment (these plots are the same as the ones mentioned in the report)
# evident in report they just called plot function of lumi package to generate the various plots
boxplot(log2(x$E),range=0,ylab="log2 intensity")
plotDensity(x$E)


#should fix this to put each arrays histogram on the plot
#hist(log(x$E))

#report says analysis performed on log2 data, this transform gives the same boxplot as theirs
x$E = log2(x$E);

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


#could quantile normalize to make all the density plots look the same? report says genespring quantile normalized the data
#but appears to be false based on the density distribution plot without quantile normalization in r?

# in report, keep probes that are > 20%tile of intensity distribution in 80% (10 or more samples) of the samples in at least 1 of the 3 groups
# so 1. break data into 3 groups (or note first 12 is group 1, 2nd 2, 3rd 3), 2. give each gene a 3 zero vector,
# 3. for each sample organize probes from most to least intense and each time probe <= 20%tile intensity add 1 to 3 zero vector group
# 4. remove all probes that have an element of 3 vec >= 3

#do we agree with the logic of doing the below probe removal? check results with and without doing this
#results nealry match after doing this, without it they look very different --> check
numPatients = 12;
numProbes = dim(x$E)[1];
samples = dim(x$E)[2];
numbCellLines = samples/numPatients
botTwentyPercentile = ceiling(numProbes*((100 - 20)/100))
probeCount = matrix(0, numProbes, numbCellLines)
twentyPercOfGroup = ceiling(numPatients*(20/100));

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

xOrig = x;

sampNames = sampNameMat[,1]; 

colnames(x$E) = sampNames;

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

#drastically different, note effect of method choice!
clustMat = t(x$E);
hc2 <- hclust(dist(clustMat), method = "average")
plot(hc2)

#my.dist <- function(x) dist(x, method="binary")
#my.hclust <- function(d) hclust(d, method="average")
#hm <- heatmap(t(x$E), distfun=my.dist, hclustfun=my.hclust)

#handle multiple probes per gene, keep probe with highest iqr (venkata) not averaging?

#1.order x$E, x$probeID, x$geneName, x$entrexGeneIDS by decreasing to increasing entrezgeneID
#2.run through entrezGeneId, when numbers are the same get iqr for each probe and determine probe with highest iqr
#3.remove probes/rows from x$E, x$probeID, x$geneName, x$enetrez that dont have highest iqr of the probes for that 1 gene



#have the two CAF groups

plus49 = which(grepl("49", sampNames));
plus133 = which(grepl("133\\+", sampNames));
minus133 = which(grepl("133-", sampNames));

plus49Genes = x[, plus49];
plus133Genes = x[, plus133];
minus133Genes = x[, minus133];

#probably want to run one way anova between the 3 groups on each individual probe

genePVals = rep(0, dim(x$E)[1]);

for(i in 1:dim(x$E)[1])
{
  plus49Gene = c(plus49Genes$E[i,]);
  names(plus49Gene) = NULL;
  plus133Gene = c(plus133Genes$E[i,]);
  names(plus133Gene) = NULL;
  minus133Gene = c(minus133Genes$E[i,]);
  names(minus133Gene) = NULL;
  
  geneData = c(plus49Gene, plus133Gene, minus133Gene);
  n = rep(numPatients, numbCellLines);
  group = rep(c("49+", "133+", "133-"), n);
  
  #tmpfn = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x));
  #tapply(geneData, group, tmpfn);
  
  data = data.frame(geneData = geneData, group = factor(group));
  fit = lm(geneData ~ group, data);
  genePVals[i] = anova(fit)$P[1];
  
}

minP = 0.05;
adjustedGenePVals = p.adjust(genePVals, method = "BH", n = length(genePVals));
sigGeneInds = which(adjustedGenePVals < minP);

#without adjusting one gets 5038 significantly varying probes, oneway.._filtereddata.txt has 5435
#but report says adjusted BH method gives 1469 varying probes, here one has 2387 significant after adjusting
#if one quantile normalizes you get 8961 without adjusting and 5108 with adjusting, if just normalize before loop get similar results

#Note: the clustering done after anova appears to have 5000 ish genes, not 1469 as mentioned in the report
#Evidence is in hierachicalCluster_onewayANOVAcorpt05.xls, which has 5000+ probes that were used in clustering
#get same dendrogram regardless of normlization and adjusting p values. Nearly the same as prior group's

xSig = x[sigGeneInds,];
dissimilarity <- 1 - cor(xSig$E);  #cor default is pearson
distance <- as.dist(dissimilarity);
hc <- hclust(distance, method = "average");
plot(hclust(distance), main="Dissimilarity = 1 - Correlation, Significant Genes", xlab="");
heatmap(xSig$E, Rowv=NA, Colv=as.dendrogram(hc));

#Permitted to skip pointless one way anova part of report on 133+/- cell lines
#test 113+ vs 133- to show not many diff expressed genes?
#On to creating gene signatures to distinguish between the two subtypes

cafGroupOneGenes = x[,c(3, 7, 23, 26, 28)];
cafGroupTwoGenes = x[,c(5, 9, 11, 25, 32)];
cafGroupOtherGenes = x[,c(1, 29)];

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

#adjusted p values all above 0.05??? non adjusted has 3647 differentially expressed genes at p=0.05 or 801 at p=0.01
minP = 0.01;
adjustedCafPValsOneTwo = p.adjust(cafPValsOneTwo, method = "bonferroni", n = length(cafPValsOneTwo));
sigGeneInds = which(adjustedCafPValsOneTwo < minP);
sigGeneInds = which(cafPValsOneTwo < 0.01);

#define a gene signature




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

