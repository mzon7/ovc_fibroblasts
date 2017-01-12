#routine to get survival stats and plot if need be

getSurvivalStats = function(survInfo, makePlot, threeGroups, censoring)
{
  
  #unknown2 is contained in unknown, but should write code to combine the two and then remove pateints
  #as some unnwon values are not in unknown2
  #unknown2 = which(is.na(tcgaVitalStat) == TRUE)
  if(threeGroups == FALSE)
  {
    survInfoOrig = survInfo
    #D ind = 1.139343 with se of 0.03926089
    dindComb = combine.est(survInfoOrig[[2]][!is.na(survInfoOrig[[2]])], survInfoOrig[[3]][!is.na(survInfoOrig[[3]])])
    survInfo = survInfo[[1]]
    colnames(survInfo) = c("timesToDeath", "vitalStats", "groups", "scores")
    #remove below line to get result of doing groupings by getting median of just that data set
    # and separating patients from data set into groups based on the median of their data set
    if(censoring == TRUE)
    {
      censoredDat = censor.time(survInfo$timesToDeath, survInfo$vitalStats, time.cens = 10)
      survInfo$timesToDeath = censoredDat$surv.time.cens
      survInfo$vitalStats = censoredDat$surv.event.cens 
      #patients followed up after time.cens get vital status of NA
      #missingPat = which(is.na(dataVitalStat))
      #if(length(missingPat) > 0)
      #{
      #  dataTimeToDeath = dataTimeToDeath[-missingPat]
      #  dataVitalStat = dataVitalStat[-missingPat] 
      #  dataVals = dataVals[, -missingPat] 
      #}
      
    }
   
    survInfo$groups = as.integer((survInfo$scores <= median(survInfo$scores)));
    survInfo$groups = survInfo[,3];
    survObj <- survfit(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
    bb <- survdiff(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups,rho=0)
    dind <- D.index(x=survInfo$scores, surv.time=survInfo$timesToDeath, surv.event=survInfo$vitalStat)
    
  }
  if(threeGroups == TRUE)
  {
    survInfoOrig = survInfo
    #D ind = 1.139343 with se of 0.03926089
    dindComb = combine.est(survInfoOrig[[2]][!is.na(survInfoOrig[[2]])], survInfoOrig[[3]][!is.na(survInfoOrig[[3]])])
    survInfo = survInfo[[1]]
    colnames(survInfo) = c("timesToDeath", "vitalStats", "groups", "scores")
    
    sortInd = sort(survInfo$scores, decreasing = TRUE, index.return=TRUE)$ix;
    thirdData = round(length(survInfo$scores)/3)
    survInfo$groups = matrix(0, length(survInfo$scores))
    survInfo$groups[sortInd[1:thirdData]] = 0;
    survInfo$groups[sortInd[thirdData:(2*thirdData)]] = 1;
    survInfo$groups[sortInd[(2*thirdData + 1):length(sortInd)]] = 2;
    
    if(censoring == TRUE)
    {
      censoredDat = censor.time(survInfo$timesToDeath, survInfo$vitalStats, time.cens = 10)
      survInfo$timesToDeath = censoredDat$surv.time.cens
      survInfo$vitalStats = censoredDat$surv.event.cens 
      #patients followed up after time.cens get vital status of NA
      #missingPat = which(is.na(dataVitalStat))
      #if(length(missingPat) > 0)
      #{
      #  dataTimeToDeath = dataTimeToDeath[-missingPat]
      #  dataVitalStat = dataVitalStat[-missingPat] 
      #  dataVals = dataVals[, -missingPat] 
      #}
      
    }
    
    survObj <- survfit(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
    bb <- survdiff(Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups,rho=0)
    dind <- D.index(x=survInfo$scores, surv.time=survInfo$timesToDeath, surv.event=survInfo$vitalStat)
    
  }
  datFrame = data.frame("surv.time" = survInfo$timesToDeath, "surv.event"=survInfo$vitalStats, "strat"=survInfo$groups)
  ddweights <- array(1, dim=nrow(datFrame))
  #km.coxph.plot(formula.s = Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups)
  km.coxph.plot(formula.s= Surv(survInfo$timesToDeath, survInfo$vitalStats) ~ survInfo$groups, data.s=datFrame,
                weight.s=ddweights, x.label="Time (years)", y.label="Percent Survival",
                main.title="", leg.text=paste(c("High Score Group", "Low Score Group"), "   ", sep=""),
                leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkred"),
               .lty=c(1,1), show.n.risk=TRUE, n.risk.step=2, n.risk.cex=0.85, verbose=FALSE, leg.bty = "n")
  #survObj <- survfit(Surv(tcgaTimeToDeath, tcgaVitalStat) ~ tcgaGroups)
  #bb <- survdiff(Surv(tcgaTimeToDeath, tcgaVitalStat) ~ tcgaGroups,rho=0)
  #dind <- D.index(x=scoreVals, surv.time=tcgaTimeToDeath, surv.event=tcgaVitalStat)
  if(makePlot == TRUE)
  {
    temp = NULL
    if(threeGroups == FALSE)
    {
      plot(survObj, lty = 2:3,col=3:4, main=paste("Kaplan-Meier Survival Curves with the Number of genes in the signature =", numGenesInSig), xlab="time (years)", ylab="survival function")
      #high sigscore class = 0, low sigscore class = 1
      legend("top",legend=c(paste("Group 1 - High Sigscore, n =", sum(survInfo$groups == 0)), paste("Group 2 - Low Sigscore, n =", sum(survInfo$groups == 1))), col=3:4, lty=2:3, horiz=FALSE, bty='n')
    }
    if(threeGroups == TRUE)
    {
      plot(survObj, lty = 2:4,col=3:5, main=paste("Kaplan-Meier estimate with 95% confidence bounds with numGenesInSig =", numGenesInSig), xlab="time (years)", ylab="survival function")
      legend("topright",legend=c(paste("Group 2 - High Sigscore, n =", sum(survInfo$groups == 0)), paste("Group 2 - Medium Sigscore, n =", sum(survInfo$groups == 1)), paste("Group 2 - Low Sigscore, n =", sum(survInfo$groups == 2))), col=3:5, lty=2:4, horiz=FALSE, bty='n')
    }
    
    temp$x = 4000.5
    temp$y = 0.5
    #text(temp, cat("D.index p value", dind$p.value))
    Corner_text <- function(text, location="topright"){
      legend(location,legend=text, bty ="n", pch=NA) 
    }
    textVal = paste("Combined D.index ", sprintf("%.4g", dindComb$estimate));
    Corner_text(text= textVal, location = "bottomleft")
    
    logRankP = 1 - pchisq(bb$chisq, 1)
    temp$y = 0.6
    textVal = paste("log rank test p value", sprintf("%.4g", logRankP));
    Corner_text(text= textVal, location = "bottomright")
    
    #get gene names for significant genes
    #sigGeneNames = xOrig$geneNames[geneSigInds]
    #xOrig$geneNames[which(xOrig$probeIds %in% geneListSecond500[[472]])]
    
    #for survivak curve, strat is group number, survival time obvious,x is my sigscores,
    #log rank test is enough for now, D.index later maybe
  }
  stats = c(dindComb$estimate, dindComb$se, (1 - pchisq(bb$chisq, 1)), dind$p.value, dind$d.index)
  return(stats)
}
