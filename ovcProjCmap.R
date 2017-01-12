#global variables needed
args = commandArgs(trailingOnly = TRUE)
print(length(args))
args

## n shows the number of up and down regulated genes. can be 30, 50 or 70
gseaGwcComp =  as.numeric(args[1])
numGenesGsea = as.numeric(args[2])
cmapGwc = as.numeric(args[3])
cmapGsea = as.numeric(args[4])
numbPerms = as.numeric(args[5])
pCutOff = as.numeric(args[6])
mord = as.numeric(args[7])

if(is.na(mord))
  mord = FALSE

library(limma);
#biocLite("genefu")
library(genefu)
library(illuminaHumanv4.db)
#biocLite("TCGAbiolinks")
#library(TCGAbiolinks)
library(matrixStats)
library(ROCR)
library(org.Hs.eg.db)
library(parallel)
#library(jetset)
library(Biobase)
library(gdata)
library(GSVA)
library(illuminaio)

#Files needed

if(mord == TRUE)
{
  #source("gwcCmapFunction.R")
  #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Natalie Files");
  rawData = read.delim("sample_probe_nonnorm_minusbackground_FinalReportNoHeaderNotePad.txt");
  
  #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data GeneSpring");
  x <- read.ilmn("Ali HT-12 expression_Sample_Probe_Profile_nonnormNoHeader.txt"); 
  
  .libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
  library(PharmacoGx)
  library(parallel)
  require(xtable)
  
  #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
  load("cmap_sig_rna.RData")
  
  #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\BenNeelCmap\\Data\\cmap pset")
  load("CMAP.RData")
  
  #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Annotation Files");
  #Note: probeId in excel file is/x data is array_address_location variable of annot$gene and target id is the ilmn_gene var
  #downloaded the file off of the internet, correct for the bead chip used during the experiments
  annot = readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx");
  
  numDrugs = dim(drug.perturbation)[2]
  write(paste("pCutOff is ", pCutOff), "test.txt")
}
if(mord == FALSE)
{
  #setwd("C:\\Users\\Michael\\Documents\\PMH Research\\R Code General")
  #source("gwcCmapFunction.R")
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Natalie Files");
  rawData = read.delim("sample_probe_nonnorm_minusbackground_FinalReportNoHeaderNotePad.txt");
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Ali_Benjamin Haibe-Kains\\RAW Data\\Raw data GeneSpring");
  x <- read.ilmn("Ali HT-12 expression_Sample_Probe_Profile_nonnormNoHeader.txt"); 
  
  .libPaths(c("/mnt/work1/users/bhklab/Rlib", .libPaths()))
  library(PharmacoGx)
  library(parallel)
  require(xtable)
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
  load("cmap_sig_rna.RData")
  
  setwd("C:\\Users\\Michael\\Documents\\PMH Research\\BenNeelCmap\\Data\\cmap pset")
  load("CMAP.RData")
  
  setwd("C:\\Users\\Michael\\Documents\\OVC Project\\Annotation Files");
  #Note: probeId in excel file is/x data is array_address_location variable of annot$gene and target id is the ilmn_gene var
  #downloaded the file off of the internet, correct for the bead chip used during the experiments
  annot = readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx");
  
  numDrugs = dim(drug.perturbation)[2]
  
  gseaGwcComp = FALSE
  # the number of genes one picks from the top and bottom of the list
  numGenesGsea = 300
  cmapGwc = FALSE
  cmapGsea = TRUE
  numbPerms = 20000
  #pCutOff is for gwc
  pCutOff = 0.05
}
numbPerms = 1000

print(pCutOff)
#gseaGwcComp = FALSE
# the number of genes one picks from the top and bottom of the list
#numGenesGsea = 100
#cmapGwc = TRUE
#cmapGsea = FALSE
#numbPerms = 20000
#numDrugs = dim(drug.perturbation)[2]
#pCutOff is for gwc
#pCutOff = 0.05

#found this in natalies pdf, so hardcoded it
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
geneNameColInd = which(grepl("TargetID",colnames(rawData)))
entrezIDColInd = which(grepl("ENTREZ_GENE_ID",colnames(rawData)))
accessionInd = which(grepl("ACCESSION",colnames(rawData)))

#GeneSpring files seem to be working fine
#set directory below to one containing probe summary profiles from GeneSpring

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

#boxplot(log2(x$E),range=0,ylab="log2 intensity")
#plotDensity(x$E)


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


#define a gene signature
patientOutliers = c(1, 11);

#Begin regression analysis/gene signature creation. xSig is just diff expressed genes
#xGen is all genes after normalizing, removing duplicates and low percentile. unbiased approach uses xGen

#xGenReg = t(xGen[sigGeneIndsNoAdj,]$E);

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

if(gseaGwcComp == TRUE)
{
  sigGenes = which(cafPValsOneTwo < 0.01)
}

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

if(gseaGwcComp == FALSE)
{
numGenesTopOrBot = numGenesGsea/2;
diffSort = order(cafGroupDiffSig, decreasing = TRUE)
cafGroupDiffSort = cafGroupDiffSig[diffSort]
xOrd = xOrd[diffSort, ]
#xOrd$geneNames = x[sigRank]
#fapStat = cafGroupDiff[28];  #check gene names to see its in the signature, i.e top 35.
#results worse with 30 genes
geneSig = rownames(xOrd)[1:numGenesTopOrBot]
geneSig = c(geneSig, rownames(xOrd)[length(rownames(xOrd)):(length(rownames(xOrd)) - numGenesTopOrBot + 1)])
geneSigCoef = cafGroupDiffSort[1:numGenesTopOrBot];
geneSigCoef = c(geneSigCoef, cafGroupDiffSort[length(rownames(xOrd)):(length(rownames(xOrd)) - numGenesTopOrBot + 1)])
}
if(gseaGwcComp == TRUE)
{
  numGenesGsea = length(cafGroupDiffSig);
  #xOrd$geneNames = x[sigRank]
  #fapStat = cafGroupDiff[28];  #check gene names to see its in the signature, i.e top 35.
  #results worse with 30 genes
  geneSig = rownames(xOrd)[1:numGenesGsea]
  geneSigCoef = cafGroupDiffSig[1:numGenesGsea];
}

corWeighted <- 
  function (x, y, w, method=c("pearson", "spearman"), alternative=c("two.sided", "greater", "less"), nperm=0, nthread=1, setseed, na.rm=FALSE) {
    
    ######################
    wcor <- function (d, w, na.rm=TRUE) {
      ### NOTE::: THIS FORMULA CAN SUFFER CATASTROPHIC CANCELATION AND SHOULD BE FIXED!!!
      #     s <- sum(w, na.rm=na.rm)
      #     m1 <- sum(d[ , 1L] * w, na.rm=na.rm) / s
      #     m2 <- sum(d[ , 2L] * w, na.rm=na.rm) / s
      #     res <- (sum(d[ , 1L] * d[ , 2L] * w, na.rm=na.rm) / s - m1 * m2) / sqrt((sum(d[ , 1L]^2 * w, na.rm=na.rm) / s - m1^2) * (sum(d[ , 2L]^2 * w, na.rm=na.rm) / s - m2^2))
      CovM <- cov.wt(d, wt=w)[["cov"]]
      res <- CovM[1,2]/sqrt(CovM[1,1]*CovM[2,2])
      return (res)
    }
    
    ######################
    
    if (missing(w)) { w <- rep(1, length(x)) / length(x) }
    if (length(x) != length(y) || length(x) != length(w)) { stop("x, y, and w must have the same length") }
    method <- match.arg(method)
    if (method == "spearman") {
      x <- rank(x)
      y <- rank(y)
    }
    alternative <- match.arg(alternative)
    
    res <- c("rho"=NA, "p"=NA)
    
    ## remove missing values
    ccix <- complete.cases(x, y, w)
    if(!all(ccix) && !na.rm) { warning("Missing values are present") }
    if(sum(ccix) < 3) {
      return(res)
    }
    x <- x[ccix]
    y <- y[ccix]
    w <- w[ccix]
    
    wc <- wcor(d=cbind(x, y), w=w)
    res["rho"] <- wc
    if (nperm > 1) {
      if (!missing(setseed)) { set.seed(setseed) }
      splitix <- parallel::splitIndices(nx=nperm, ncl=nthread)
      if (!is.list(splitix)) { splitix <- list(splitix) }
      splitix <- splitix[sapply(splitix, length) > 0]
      mcres <- parallel::mclapply(splitix, function(x, xx, yy, ww) {
        pres <- sapply(x, function(x, xx, yy, ww) {
          ## permute the data and the weights
          d2 <- cbind(xx[sample(1:length(xx))], yy[sample(1:length(yy))])
          w2 <- ww[sample(1:length(ww))]
          return(wcor(d=d2, w=w2))
        }, xx=xx, yy=yy, ww=ww)
        return(pres)
      }, xx=x, yy=y, ww=w)
      perms <- do.call(c, mcres)
      
      switch (alternative,
              "two.sided" = { 
                if (res["rho"] < 0) { p <- sum(perms <= res, na.rm=TRUE) } else { p <- sum(perms >= res, na.rm=TRUE) }
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
                p <- p * 2
              },
              "greater" = {
                p <- sum(perms >= res, na.rm=TRUE) 
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
              },
              "less" = {
                p <- sum(perms <= res, na.rm=TRUE) 
                if (p == 0) { p <- 1 / (nperm + 1) } else { p <- p / nperm }
              })
      res["p"] <- p
    }
    return(res)
  }

combineTest <-
  function(p, weight, method=c("fisher", "z.transform", "logit"), hetero=FALSE, na.rm=FALSE) {
    if(hetero) { stop("function to deal with heterogeneity is not implemented yet!") }
    method <- match.arg(method)
    na.ix <- is.na(p)
    if(any(na.ix) && !na.rm) { stop("missing values are present!") }
    if(all(na.ix)) { return(NA) } ## all p-values are missing
    p <- p[!na.ix]
    k <- length(p)
    if(k == 1) { return(p) }
    if(missing(weight)) { weight <- rep(1, k); }
    switch(method,  
           "fisher"={
             cp <- pchisq(-2 * sum(log(p)), df=2*k, lower.tail=FALSE)
           }, 
           "z.transform"={
             z <- qnorm(p, lower.tail=FALSE)
             cp <- pnorm(sum(weight * z) / sqrt(sum(weight^2)), lower.tail=FALSE)
           }, 
           "logit"={
             tt <- (- sum(log(p / (1 - p)))) / sqrt(k * pi^2 * (5 * k + 2) / (3 * (5 * k + 4)))
             cp <- pt(tt,df=5*k+4, lower.tail=FALSE)
           })
    return(cp)
  }

#Begin GWC Code

if(cmapGwc == TRUE)
{
  drugNotes = drugInfo(CMAP)
  
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
        #ensembleId = paste0(ensembleId, "_at")
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
  
  fileNameCmap = "ovcProjectGwcCmap"
  ensembleVec = paste0(ensembleVec, "_at")
  genesCmapAll = as.data.frame(NULL)
  #Important to flip direction as files are disease - normal
  tStatVec = -1*tStatVec
  genesCmapAll[1:length(ensembleVec),1] = tStatVec
  #fdr no good for dep file as p min = 0.37
  genesCmapAll[1:length(ensembleVec),2] = pValVec
  rownames(genesCmapAll) = ensembleVec
  colnames(genesCmapAll) = c("tstat", "pvalue")
  
  sigGenes = which(genesCmapAll[, 2] < pCutOff)
  genesCmapAll = genesCmapAll[sigGenes,]
  #at 0.01 still have 1900 genes, test it out, but maybe too many and need lofFC threshold also
  orderingAll = order(abs(genesCmapAll[,1]), decreasing = TRUE)
  
  numSigIts = 15
  genesSigEnd = length(sigGenes)
  #if(length(sigRows) > 2000)
  #  genesSigEnd = 2000
  genesSigStart = 400
  genesSigInc = round((genesSigEnd - genesSigStart)/numSigIts) + 1
  
  #run below code twice basically, once on the significant genes and once on all the genes
  #choose drugs that are high on the list in both cases for a best of both worlds scenario
  
  #beginning and end included
  quickTest = FALSE
  
  cmapList = list()
  totRuns = round((genesSigEnd - genesSigStart)/genesSigInc)
  for(g in 1:totRuns)
  {
    genesUse = genesSigStart + (g-1)*genesSigInc
    sigRows = orderingAll[1:genesUse]
    genesCmap = genesCmapAll[sigRows, ]
    
    cmapTab = as.data.frame(NULL)
    numDrugs = dim(drug.perturbation)[2]
    drugNames = drug.perturbation[1, ,1][1:numDrugs]
    drugNames = names(drugNames)
    noCores <- detectCores() - 2
    # Initiate cluster
    cl <- makeCluster(noCores)
    pco = library(PharmacoGx)
    #clusterEvalQ(cl, library(PharmacoGx))
    
    if(quickTest == TRUE)
    {
      numDrugs = 100
      numbPerms = 100 
    }
    
    drugsList = list()
    for(i in 1:numDrugs)
    {
      drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
    }
    #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
    #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
    #connectivity score is high when either
    # drug t-stat high (drug raised expression of gene) & data t-stat high 
    #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
    #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
    #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
    clusterExport(cl, list("drugsList", "genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    
    dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = genesCmap, y = z,method="gwc", gwc.method="spearman", nperm = numbPerms))
    stopCluster(cl)
    
    for(i in 1:numDrugs)
    {
      drugInd = which(rownames(drugNotes) == drugNames[i])
      drugInfoTab = drugNotes[drugInd, ]
      cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
    }
    
    colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
    rownames(cmapTab) = drugNames[1:numDrugs]
    sorting = order(cmapTab[,"score"], decreasing = TRUE)
    cmapTabSort = cmapTab[sorting,]
    cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
    
    cmapTab = cmapTab[-which(is.na(cmapTab[,1])), ]
    cmapList[[g]] = cmapTab
    #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
    print(g)
  }
  
  rankFrame = as.data.frame(NULL)
  #below frame finds drugs that enhance phenotype the most
  rankFrameRev = as.data.frame(NULL)
  rankChangeFrame = as.data.frame(NULL)
  colnameVec = c()
  for(g in 1:totRuns)
  {
    ranking = c()
    rankingRev = c()
    cmapRes = cmapList[[g]]["score"]
    #added -1 so now highest rank --> most +ve conectivity score, rank = M --> higher score than M - 1 drugs
    ordering = sort(-1*as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    orderingRev = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    #ordering = order(as.numeric(cmapRes[,1]), decreasing = TRUE)
    for(i in 1:length(ordering))
    {
      ranking[ordering[i]] = i
      rankingRev[orderingRev[i]] = i
    }
    if(g == 1)
    {
      rankFrame = as.data.frame(ranking) 
      rankFrameRev = as.data.frame(rankingRev)
    }
    if(g > 1)
    {
      rankFrame = cbind(rankFrame, ranking) 
      rankFrameRev = cbind(rankFrameRev, rankingRev)
    }
    if(g == 2)
      rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
    if(g > 2)
      rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
    colnameVec = c(colnameVec, as.character(genesSigStart + (g-1)*genesSigInc))
    oldRanking = ranking
  }
  rownames(rankFrame) = rownames(cmapRes)
  colnames(rankFrame) = colnameVec
  
  rownames(rankFrameRev) = rownames(cmapRes)
  colnames(rankFrameRev) = colnameVec
  
  rownames(rankChangeFrame) = rownames(cmapRes)
  colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]
  
  rankChangeFrameOrig = rankChangeFrame
  for(g in 1:dim(rankChangeFrame)[2])
  {
    changes = rankChangeFrame[,g]
    #maybe use average instead? differnce from rank 1 to 2 still the same, 1 is just higher, min likely better
    changeRanks = rank(-1*changes, ties.method = "min")
    changeRanksOrig = rank(changes, ties.method = "min")
    rankChangeFrame[,g] = changeRanks
    rankChangeFrameOrig[,g] = changeRanksOrig
  }
  
  
  drugRankMeansPos = rowMeans(rankFrame)
  drugRankMeansNeg = rowMeans(rankFrameRev)
  drugRankSigns = c()
  drugRankMeans = c()
  for(g in 1:length(drugRankMeansPos))
  {
    if(drugRankMeansPos[g] > drugRankMeansNeg[g])
    {
      drugRankMeans[g] = drugRankMeansPos[g]
      drugRankSigns[g] = 1
    }
    if(drugRankMeansNeg[g] > drugRankMeansPos[g])
    {
      drugRankMeans[g] = drugRankMeansNeg[g]
      drugRankSigns[g] = -1
    }
  }
  #drugRankMeans = drugRankMeans
  drugRankChangeMeans = rowMeans(rankChangeFrame)
  drugRankChangeMeansOrig = rowMeans(rankChangeFrameOrig)
  sigFinalScores = drugRankSigns*(drugRankMeans + drugRankChangeMeans)/2
  
  #maybe some poor naming but pretty sure that drugRankMeansNeg = 1 means on average always had highest +score and best drug to reverse phenotype rank wise
  cmapResults = cbind(sigFinalScores, drugRankMeansNeg, drugRankChangeMeansOrig, cmapTab)
  sorting = order(cmapResults[,"sigFinalScores"], decreasing = TRUE)
  cmapResults = cmapResults[sorting,]
  
  write.table(cmapResults, paste("sig genes queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"in file", fileNameCmap,".xls"), col.names = NA, sep="\t")
  
  sigFinalScoresPos = (drugRankMeansPos + drugRankChangeMeans)/2
  sigFinalScoresNeg = -1*(drugRankMeansNeg + drugRankChangeMeans)/2
  pdf("ranking 1 final scores vs ranking 2 final scores")
  plot(sigFinalScoresPos, sigFinalScoresNeg)
  dev.off()
  
  topDrug = rownames(cmapResults)[1]
  rankFrameInd = which(rownames(rankFrame) == topDrug)
  rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
  rankVsSize = rankFrame[rankFrameInd, ]
  changeVsSize = rankChangeFrame[rankChangeInd,]
  
  sizeVec = c()
  for(g in 1:totRuns)
    sizeVec = c(sizeVec, genesSigStart + (g-1)*genesSigInc)
  
  pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"SigGenesRankVsQuerySizeTopDrug.pdf"))
  plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size", ylab = "Rank")
  dev.off()
  
  pdf(paste("sig gene queried", genesSigStart,"to",genesSigEnd,"by",genesSigInc,"SigGenesRankChangeFromQueryToQueryTopDrug.pdf"))
  plot(1:(totRuns-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank Change)")
  dev.off()
  
  save(rankFrame, file = paste("SigGenesRankFrame.RData", fileNameCmap))
  save(rankChangeFrame, file = paste("SigGenesRankChangeFrame.Rdata", fileNameCmap))
  
  
  #Run it again on the whole genome
  remainingGenome = FALSE
  if(remainingGenome == TRUE)
  {
  
  numGenomeIts = 20
  genesEnd = dim(genesCmapAll)[1]
  genesStart = genesSigEnd
  genesInc = round((genesEnd-genesStart)/numGenomeIts)
  
  
  cmapList = list()
  totRuns = round((genesEnd - genesStart)/genesInc) + 1
  for(g in 1:totRuns)
  {
    genesUse = genesStart + (g-1)*genesInc
    sigRows = orderingAll[1:genesUse]
    genesCmap = genesCmapAll[sigRows, ]
    
    cmapTab = as.data.frame(NULL)
    numDrugs = dim(drug.perturbation)[2]
    drugNames = drug.perturbation[1, ,1][1:numDrugs]
    drugNames = names(drugNames)
    noCores <- detectCores() - 2
    # Initiate cluster
    cl <- makeCluster(noCores)
    pco = library(PharmacoGx)
    #clusterEvalQ(cl, library(PharmacoGx))
    
    if(quickTest == TRUE)
    {
      numDrugs = 300
      numbPerms = 100 
    }
    
    drugsList = list()
    for(i in 1:numDrugs)
    {
      drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
    }
    #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
    #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
    #connectivity score is high when either
    # drug t-stat high (drug raised expression of gene) & data t-stat high 
    #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
    #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
    #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
    clusterExport(cl, list("drugsList", "genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
    
    dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = genesCmap, y = z,method="gwc", gwc.method="spearman", nperm = numbPerms))
    stopCluster(cl)
    
    for(i in 1:numDrugs)
    {
      drugInd = which(rownames(drugNotes) == drugNames[i])
      drugInfoTab = drugNotes[drugInd, ]
      cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
    }
    
    colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
    rownames(cmapTab) = drugNames[1:numDrugs]
    sorting = order(cmapTab[,"score"], decreasing = TRUE)
    cmapTabSort = cmapTab[sorting,]
    cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
    
    cmapTab = cmapTab[-which(is.na(cmapTab[,1])), ]
    cmapList[[g]] = cmapTab
    #write.table(cmapTabSort, paste("genes used is", genesUse, fileNameCmap), col.names = NA, sep="\t")
    print(g)
  }
  
  rankFrame = as.data.frame(NULL)
  rankChangeFrame = as.data.frame(NULL)
  colnameVec = c()
  for(g in 1:totRuns)
  {
    ranking = c()
    cmapRes = cmapList[[g]]["score"]
    ordering = sort(as.numeric(cmapRes[,1]), index.return=TRUE,decreasing = TRUE)$ix
    #ordering = order(as.numeric(cmapRes[,1]), decreasing = TRUE)
    for(i in 1:length(ordering))
      ranking[ordering[i]] = i
    if(g == 1)
      rankFrame = as.data.frame(ranking)
    if(g > 1)
      rankFrame = cbind(rankFrame, ranking)
    if(g == 2)
      rankChangeFrame = as.data.frame(abs(ranking - oldRanking))
    if(g > 2)
      rankChangeFrame = cbind(rankChangeFrame, abs(ranking - oldRanking))
    colnameVec = c(colnameVec, as.character(genesStart + (g-1)*genesInc))
    oldRanking = ranking
  }
  rownames(rankFrame) = rownames(cmapRes)
  colnames(rankFrame) = colnameVec
  
  rownames(rankChangeFrame) = rownames(cmapRes)
  colnames(rankChangeFrame) = colnameVec[2:length(colnameVec)]
  
  for(g in 1:dim(rankChangeFrame)[2])
  {
    changes = rankChangeFrame[,g]
    changeRanks = rank(changes, ties.method = "min")
    rankChangeFrame[,g] = changeRanks
  }
  
  
  drugRankMeans = rowMeans(rankFrame)
  drugRankChangeMeans = rowMeans(rankChangeFrame)
  genomeFinalScores = (drugRankMeans + drugRankChangeMeans)/2
  
  avgFinalScores = (genomeFinalScores + sigFinalScores)/2
  
  cmapResults = cbind(avgFinalScores, genomeFinalScores, sigFinalScores, drugRankMeans, drugRankChangeMeans, cmapTab)
  sorting = order(cmapResults[,"genomeFinalScores"], decreasing = FALSE)
  cmapResults = cmapResults[sorting,]
  
  write.table(cmapResults, paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"in file", fileNameCmap, ".xls"), col.names = NA, sep="\t")
  
  topDrug = rownames(cmapResults)[1]
  rankFrameInd = which(rownames(rankFrame) == topDrug)
  rankChangeInd = which(rownames(rankChangeFrame) == topDrug)
  rankVsSize = rankFrame[rankFrameInd, ]
  changeVsSize = rankChangeFrame[rankChangeInd,]
  
  sizeVec = c()
  for(g in 1:totRuns)
    sizeVec = c(sizeVec, genesStart + (g-1)*genesInc)
  
  pdf(paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"RankVsQuerySizeTopDrug.pdf"))
  plot(sizeVec, rankVsSize, main = paste("Rank as a Function of Query Size for the Top Drug", topDrug), xlab = "Query Size", ylab = "Rank")
  dev.off()
  
  pdf(paste("genome queried", genesStart,"to",genesEnd,"by",genesInc,"RankChangeFromQueryToQueryTopDrug.pdf"))
  plot(1:(totRuns-1), changeVsSize, main = paste("Change in Rank from Query to Query for the Top Drug", topDrug), xlab = "Query Number", ylab = "abs(Rank Change)")
  dev.off()
  
  save(rankFrame, file = paste(fileNameCmap, "GenomeGenesRankFrame.RData"))
  save(rankChangeFrame, file = paste(fileNameCmap, "GenomeGenesRankChangeFrame.Rdata"))
  }
  
  #gwcCmapAnalysis = function(drug.perturbation, drugNotes, fileNameCmap, ensembleVec, estimateVec, pVec, numSigIts, numGenomeIts, numbPerms = 1000, pCutOff = 0.05, quickTest = FALSE)
    
  #gwcCmapAnalysis(drug.perturbation, drugNotes, "ovcCmapAli", ensembleVec, -1*tStatVec, pValVec, numSigIts = 10, numGenomeIts = 20, quickTest = TRUE)
  oldGwc = FALSE
  if(oldGwc == TRUE)
  {
  genesCmap = as.data.frame(NULL)
  #IMPORTANT, clarify this makes sense
  tStatVec = -1*tStatVec
  genesCmap[1:length(ensembleVec),1] = tStatVec
  genesCmap[1:length(ensembleVec),2] = pValVec
  rownames(genesCmap) = ensembleVec
  colnames(genesCmap) = c("tstat", "pvalue")
  
  sigRows = which(genesCmap[, 2] < pCutOff)
  genesCmap = genesCmap[sigRows, ]
  
  cmapTab = as.data.frame(NULL)
  drugNames = drug.perturbation[1, ,1][1:numDrugs]
  drugNames = names(drugNames)
  
  #numDrugs = 5
  # Calculate the number of cores
  noCores <- detectCores() - 1
  # Initiate cluster
  cl <- makeCluster(noCores)
  #clusterExport(cl, list("genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "corWeighted"), envir = environment())
  pco = library(PharmacoGx)
  clusterExport(cl, list("genesCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
  
  drugsList = list()
  for(i in 1:numDrugs)
  {
    drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
  }
  #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
  #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
  #connectivity score is high when either
  # drug t-stat high (drug raised expression of gene) & data t-stat high 
  #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
  #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
  #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
  
  dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = genesCmap, y = z,method="gwc", gwc.method="spearman", nperm = numbPerms))
  stopCluster(cl)
  
  for(i in 1:numDrugs)
  {
    drugInd = which(rownames(drugNotes) == drugNames[i])
    drugInfoTab = drugNotes[drugInd, ]
    cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
  }
  
  #old serial way
  #for(i in 1:numDrugs)
  #{
  #  drugStats = drug.perturbation[,i ,c("tstat", "pvalue")]
  #  dataCor <- connectivityScore(genesCmap, drugStats,method="gwc", gwc.method="spearman", nperm = numbPerms)
  #  drugInd = which(rownames(drugNotes) == drugNames[i])
  #  drugInfoTab = drugNotes[drugInd, ]
  #  cmapTab2 = rbind(cmapTab2, c(dataCor, as.character(drugInfoTab)), stringsAsFactors = FALSE)
  #}
  colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
  rownames(cmapTab) = drugNames
  sort = order(cmapTab[,"score"], decreasing = TRUE)
  cmapTabSort = cmapTab[sort,]
  cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
  
  #setwd("C:\\Users\\Michael\\Documents\\OVC Project\\CMAP")
  write.table(cmapTabSort, paste("GWC CMAP Results Ovarian Fibroblast Project Using a P Threshold of =", pCutOff,".xls"), col.names = NA, sep="\t")
  }
  
}

#Begin GSEA Code

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
    #print(i)
    #multiple indexEnt for same entrez lead to same ensemble id
    indexEnt = which(excelEntrezIds == geneSigEntrez[i])
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
      print(ensembleId)
      if(length(grep(ensembleId, ensembleVec)) == 0)
      {
        ensembleId = paste0(ensembleId, "_at")
        ensembleVec = c(ensembleVec, ensembleId)
        genesPresent = c(genesPresent, i) 
      }
    }
    }
  }
  
  #below two lines are time consuming, so just use result as its data independent
  #data(CMAPsmall)
  #drug.perturbation <- drugPerturbationSig(CMAPsmall, mDataType="rna")

  geneSigCmap = as.data.frame(NULL)
  geneSigCmap[1:length(ensembleVec),1] = ensembleVec
  rownames(geneSigCmap) = ensembleVec
  geneSigCmap[,2] = cafGroupDiff[genesPresent]
  colnames(geneSigCmap) = c("feature", "direction")
  #flip the directions as we query for drugs that result in reversing the signature
  geneSigCmap[,2] = replace(geneSigCmap[,2], geneSigCmap[,2] > 0, 1)
  geneSigCmap[,2] = replace(geneSigCmap[,2], geneSigCmap[,2] < 0, -1)
  #geneSigCmap is formatted like HDAC_genes which is used in the pharmaco documentation example
  data(HDAC_genes)
  
  drugNames = drug.perturbation[1, ,1][1:numDrugs]
  drugNames = names(drugNames)
  
  # Calculate the number of cores
  noCores <- detectCores() - 3
  # Initiate cluster
  # Initiate cluster
  cl <- makeCluster(noCores)
  pco = library(PharmacoGx)
  #clusterEvalQ(cl, library(PharmacoGx))
  numDrugs = 5
  drugsList = list()
  for(i in 1:numDrugs)
  {
    drugsList[[i]] = drug.perturbation[,i ,c("tstat", "pvalue")]
  }
  #t-stat high in drug sig --> raised expression of gene (after drug - before drug)
  #and t-stat was flipped to query the reverse signature, i.e fap hi expression --> similar to fap lo expression
  #connectivity score is high when either
  # drug t-stat high (drug raised expression of gene) & data t-stat high 
  #i.e gene expression low in fap hi/disease group and drug raises the genes expression (fap hi low exp gene --> high exp gene)
  #drug t-stat low (drug lowered the expression of the gene) & data t-stat low
  #i.e gene expression high in fap Hi group and drug lowers the genes expression (fap hi high exp --> low gene exp)
  clusterExport(cl, list("drugsList", "geneSigCmap", "numbPerms", "connectivityScore", "gwc", "intersectList", "pco", "corWeighted", "combineTest"), envir = environment())
  
  dataCor = clusterApplyLB(cl, drugsList, function(z) connectivityScore(x = z, y = geneSigCmap[,2,drop=FALSE],method="gsea", nperm = numbPerms))
  stopCluster(cl)
  
  cmapTab = as.data.frame(NULL)
  for(i in 1:numDrugs)
  {
    drugInd = which(rownames(drugNotes) == drugNames[i])
    drugInfoTab = drugNotes[drugInd, ]
    cmapTab = rbind(cmapTab, c(dataCor[[i]], as.character(drugInfoTab)), stringsAsFactors = FALSE)
  }
  
  #old serial way
  #for(i in 1:numDrugs)
  #{
  #  drugStats = drug.perturbation[,i ,c("tstat", "pvalue")]
  #  dataCor <- connectivityScore(genesCmap, drugStats,method="gwc", gwc.method="spearman", nperm = numbPerms)
  #  drugInd = which(rownames(drugNotes) == drugNames[i])
  #  drugInfoTab = drugNotes[drugInd, ]
  #  cmapTab2 = rbind(cmapTab2, c(dataCor, as.character(drugInfoTab)), stringsAsFactors = FALSE)
  #}
  colnames(cmapTab) = c("score", "fdr adjusted pvalue", colnames(drugNotes[1,]))
  rownames(cmapTab) = drugNames[1:numDrugs]
  sort = order(cmapTab[,"score"], decreasing = TRUE)
  cmapTabSort = cmapTab[sort,]
  cmapTabSort[,2] = as.numeric(cmapTabSort[,2])
  cmapTabSort[,2] = replace(cmapTabSort[,2], cmapTabSort[,2] == 0, 1/(numbPerms+1))
  cmapTabSort[,2] = p.adjust(cmapTabSort[,2], method = "fdr")
  
  write.table(cmapTabSort, paste("GSEA CMAP Results Ovarian Fibroblast Project at Signature Size =", numGenesGsea,".xls"), col.names = NA, sep="\t")
  
  #old serial code
  #remove 1:5 to do all the drugs
  #res <- apply(drug.perturbation[,1:5,c("tstat", "fdr")],2, function(x, genesCmap)
  #{
  #  return(connectivityScore(x=x,y=genesCmap[,2,drop=FALSE],method="gsea", nperm=numbPerms))
  #}, genesCmap=geneSigCmap)
  #rownames(res) <- c("Connectivity", "P Value")
  #res <- t(res)
  #res <- res[order(res[,1], decreasing=TRUE),]
  #xtable(res,caption='Connectivity Score results for the gene signature.')
  
  #resFix = as.data.frame(res)
  #resFix[,2] = replace(resFix[,2], resFix[,2] == 0, 1/(numbPerms+1))
  #adjustedP = p.adjust(resFix[,2], method = "fdr", n = length(resFix[,2]));
  #resFix[,3] = NULL
  #resFix[1:length(resFix[,2]),3] = adjustedP
  #resFix <- resFix[order(resFix[,3]),]
  #colnames(resFix) = c("Connectivity", "P Value", "FDR Adjusted P Value")

}




