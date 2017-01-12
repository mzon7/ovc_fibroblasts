#get best probes for the data sets based on highest IQR

#not really using jetset so ignore that part, keeping for later just in case, doesn't work when TRUE
getBestProbes = function(data)
{
  dataValsOrig = data@assayData$exprs;
  #Result is much worse if the below line is included
  #dataVals <- normalizeBetweenArrays(dataVals,method="quantile")
  dataInfo = data@featureData@data;
  dataInfoEntrezGene.ID = as.numeric(dataInfo$EntrezGene.ID)
  dataInfoProbeset = as.vector(dataInfo$probeset)
  dataInfoGene = as.vector(dataInfo$gene)
  sortInd = sort(dataInfoEntrezGene.ID, decreasing = FALSE, index.return=TRUE)$ix;
  #definitly sorts least to greatest
  dataInfoEntrezGene.ID = dataInfoEntrezGene.ID[sortInd];
  dataVals = dataValsOrig[sortInd, ];
  dataInfoProbeset = dataInfoProbeset[sortInd];
  dataInfoGene = dataInfoGene[sortInd];
  
  #Dr. Haibe-Kains adapted code for selecting best probes on hgu133a chip
  useJetSet = FALSE
  
  if(useJetSet == TRUE)
  {
    
    js <- jscores(chip="hgu133a", probeset=dataInfoProbeset)
    js <- js[dataInfoProbeset, , drop=FALSE]
    #took out stripWhiteSpace as I couldnt find it
    geneid1 <- as.character(js[ ,"EntrezID"])
    names(geneid1) <- rownames(js)
    geneid2 <- sort(unique(geneid1))
    names(geneid2) <- paste("geneid", geneid2, sep=".")
    gix1 <- !is.na(geneid1)
    gix2 <- !is.na(geneid2)
    geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
    ## probes corresponding to common gene ids
    gg <- names(geneid1)[is.element(geneid1, geneid.common)]
    gid <- geneid1[is.element(geneid1, geneid.common)] ## duplicated gene ids
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
    ## unique gene ids
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
    ## which are the best probe for each gene
    js <- data.frame(js, "best"=FALSE, stringsAsFactors=FALSE)
    js[gg.uniq, "best"] <- TRUE
    ## data for duplicated gene ids
    if(length(gid.dupl) > 0) {    
      ## use jetset overall score to select the best probesets
      myscore <- js[gg.dupl,"overall"]
      myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
      myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
      js[myscore[ ,"probe"], "best"] <- TRUE
    }
    
    bestProbes = js["best"]
    bestProbeInds = which(bestProbes == TRUE)
    #dataVals = dataVals[bestProbeInds, ]
    #dataInfo = dataInfo[bestProbeInds, ]
    
  }
  
  if(useJetSet == FALSE)
  {
    
    bestProbes = c()
    i = 1;
    while(i < dim(dataVals)[1])
    {
      entrezId = dataInfoEntrezGene.ID[i];
      origInd = i;
      probesWithId = c();
      while(dataInfoEntrezGene.ID[i] == entrezId)
      {
        probesWithId = c(probesWithId, i);
        i = i + 1;
        if(i == dim(dataVals)[1])
        {
          if(dataInfoEntrezGene.ID[i] == entrezId)
            probesWithId = c(probesWithId, i)
          if(dataInfoEntrezGene.ID[i] != entrezId)
            bestProbes = c(bestProbes, i)
          
          break
        }
      }
      if(length(probesWithId) > 1)
      {
        probes = dataVals[probesWithId,]
        iqrs = rowIQRs(probes);
        keepProbe = which(iqrs == max(iqrs)) + (origInd - 1);
        bestProbes = c(bestProbes, keepProbe)
      }
      if(length(probesWithId) == 1)
      {
        bestProbes = c(bestProbes, origInd)
      }
    }
    
    #dataVals = dataVals[bestProbes, ]
    #dataInfo = dataInfo[bestProbes, ]
    
    #data@assayData$exprs = dataVals;
    #data@featureData@data = dataInfo;
    #data@assayData$exprsCopy = dataVals;
    return(bestProbes)
  }
  
}

getBestProbesOriginal = function(data)
{
  dataVals = data@assayData$exprs;
  #Result is much worse if the below line is included
  #dataVals <- normalizeBetweenArrays(dataVals,method="quantile")
  dataInfo = data@featureData@data;
  dataInfo$EntrezGene.ID = as.numeric(dataInfo$EntrezGene.ID)
  dataInfo$probeset = as.vector(dataInfo$probeset)
  dataInfo$gene = as.vector(dataInfo$gene)
  sortInd = sort(dataInfo$EntrezGene.ID, decreasing = FALSE, index.return=TRUE)$ix;
  dataInfo$EntrezGene.ID = dataInfo$EntrezGene.ID[sortInd];
  dataVals = dataVals[sortInd, ];
  dataInfo$probeset = dataInfo$probeset[sortInd];
  dataInfo$gene = dataInfo$gene[sortInd];
  
  #Dr. Haibe-Kains adapted code for selecting best probes on hgu133a chip
  useJetSet = FALSE
  
  if(useJetSet == TRUE)
  {
    
    js <- jscores(chip="hgu133a", probeset=dataInfo$probeset)
    js <- js[dataInfo$probeset, , drop=FALSE]
    #took out stripWhiteSpace as I couldnt find it
    geneid1 <- as.character(js[ ,"EntrezID"])
    names(geneid1) <- rownames(js)
    geneid2 <- sort(unique(geneid1))
    names(geneid2) <- paste("geneid", geneid2, sep=".")
    gix1 <- !is.na(geneid1)
    gix2 <- !is.na(geneid2)
    geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
    ## probes corresponding to common gene ids
    gg <- names(geneid1)[is.element(geneid1, geneid.common)]
    gid <- geneid1[is.element(geneid1, geneid.common)] ## duplicated gene ids
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
    ## unique gene ids
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
    ## which are the best probe for each gene
    js <- data.frame(js, "best"=FALSE, stringsAsFactors=FALSE)
    js[gg.uniq, "best"] <- TRUE
    ## data for duplicated gene ids
    if(length(gid.dupl) > 0) {    
      ## use jetset overall score to select the best probesets
      myscore <- js[gg.dupl,"overall"]
      myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
      myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
      js[myscore[ ,"probe"], "best"] <- TRUE
    }
    
    bestProbes = js["best"]
    bestProbeInds = which(bestProbes == TRUE)
    dataVals = dataVals[bestProbeInds, ]
    dataInfo = dataInfo[bestProbeInds, ]
    
  }
  
  if(useJetSet == FALSE)
  {
    
    bestProbes = c()
    i = 1;
    while(i < dim(dataVals)[1])
    {
      entrezId = dataInfo$EntrezGene.ID[i];
      origInd = i;
      probesWithId = c();
      while(dataInfo$EntrezGene.ID[i] == entrezId)
      {
        probesWithId = c(probesWithId, i);
        i = i + 1;
        if(i == dim(dataVals)[1])
        {
          break
        }
      }
      if(length(probesWithId) > 1)
      {
        probes = dataVals[probesWithId,]
        iqrs = rowIQRs(probes);
        keepProbe = which(iqrs == max(iqrs)) + (origInd - 1);
        bestProbes = c(bestProbes, keepProbe)
      }
      if(length(probesWithId) == 1)
      {
        bestProbes = c(bestProbes, origInd)
      }
    }
    
    dataVals = dataVals[bestProbes, ]
    dataInfo = dataInfo[bestProbes, ]
    
    data@assayData$exprs = dataVals;
    data@featureData@data = dataInfo;
    return(data)
  }
}
  