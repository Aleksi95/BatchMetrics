
#source("../plotFunctions.R")
library(pspearman)
library(reshape2)

meanRankCorr = function(resdata, rank = 1:11){
  corr2 = apply(t(resdata), 2, function(x) spearman.test(x, rank)$estimate)
  return(mean(corr2))
}

rankCorrelations = function(resdata, type = "spearman", score = "ratio", levels = seq(0,1, by=0.1), metrics = NULL){
  if(is.null(metrics)){
    metrics = names(resdata[[score]])
  }
  meanCorrs = c()
  pvalues = c()
  data = resdata[[score]][metrics]
  for(d in data){
    #data = data.frame(method_class, t(d))
    vec = unlist(c(t(d)))
    lvl = rep(levels, nrow(d))
    if(type == "spearman"){
  
      corrs = spearman.test(vec, lvl)$estimate
      pvals = spearman.test(vec, lvl)$p.value
    }
    if(type == "pearson"){
      corrs = cor.test(vec, lvl)$estimate
      pvals = cor.test(vec, lvl)$p.value
    }
    meanCorrs = c(meanCorrs, corrs)
    pvalues = c(pvalues, pvals)
  }
  names(meanCorrs) = names(data)
  names(pvalues) = names(data)
  return(list(estimates = meanCorrs, p.values = pvalues))
}


#Ftest2 = function(data){
#  data = cbind(group = as.factor(1:nrow(data)), as.data.frame(data))
  
  #is.na(data) = sapply(data, is.infinite)
  
#  data = na.omit(data)
  
#  melt_df = melt(data, id.vars = "group")
  
#  fit = aov(value~group, data = melt_df)
  
#  return(summary(fit)[[1]][1,"F value"])
#}

Ftest2 = function(data, test = "AOV", alternative = "increasing"){
  
  data = as.matrix(data)
  if(!is.null(colnames(data))){
    colnames(data) = NULL
  }
  data = cbind(group = as.factor(1:nrow(data)), as.data.frame(data))
  
  #is.na(data) = sapply(data, is.infinite)
  
  data = na.omit(data)
  
  melt_df = melt(data, id.vars = "group")
  
  if(test == "AOV"){
    fit = aov(value~group, data = melt_df)
  
    return(summary(fit)[[1]][1,"F value"])
  } 
  if(test == "JT"){
    fit = DescTools::JonckheereTerpstraTest(melt_df$value, g = ordered(melt_df$group), alternative = alternative)
    
    return(c(JT = fit$statistic, pvalue = fit$p.value))
  }
}



write_results = function(Results, filename, score = "ratio", levels = seq(0,1, by=0.1), metrics = NULL, latex = FALSE){
  for(i in 1:length(Results)){
  res_data = Results[[i]]
  if(is.null(metrics)){
    metrics = names(res_data$results[[score]])
  }
  scoresdf = data.frame(matrix(ncol = length(metrics), nrow = 0))
  #colnames(scoresdf) = metrics
  
  #fscores1 = sapply(metrics, function(i){
    #data = t(res_data$results[[score]][[i]])
    #data[is.infinite(data)] = NA
    #Ftest2(data)
  #})
  fscores1 = sapply(metrics, function(i) Ftest2(t(res_data$results[[score]][[i]])))
  #jtScores = sapply(metrics, function(i) Ftest2(t(res_data$results[[score]][[i]]), 
                                               # test = "JT",
                                               # alternative = c("increasing", "decreasing")[(score == "batch.signal") + 1])[1])
  
  pears = rankCorrelations(res_data$results, type = "pearson", score = score, levels = levels, metrics = metrics)$estimates
  spearman = rankCorrelations(res_data$results, score = score, levels = levels, metrics = metrics)$estimates
  scoresdf = rbind(scoresdf, pears)
  scoresdf = rbind(scoresdf, spearman)
  scoresdf = rbind(scoresdf, fscores1)
  scoresdf = rbind(scoresdf, log(fscores1))
  logFxRho = spearman * log(fscores1)
  scoresdf = rbind(scoresdf, logFxRho)
  #scoresdf = rbind(scoresdf, jtScores)
  #scoresdf = rbind(scoresdf, log(jtScores)*spearman)
  sign = c(-1,1)[(score == "batch.signal") + 1]
  scoresdf = rbind(scoresdf, rank(sign*logFxRho))
  #scoresdf = rbind(scoresdf, rank(-(log(jtScores)*spearman)))
  colnames(scoresdf) = metrics
  rownames(scoresdf) = c("Pearson correlation",  
                         "Rank correlation", 
                         "ANOVA F-statistic", 
                         "log(F)", 
                         "\\rho \\times log(F)", 
                         "Ranking")
  
 if(latex == TRUE){
   writeLatexTable(round(scoresdf, digits = 3), filename = paste(names(Results)[i], filename, sep = "_"))
 } else {
    write.table(round(scoresdf, digits = 3), file = paste(names(Results)[i], filename, sep = "_"))
 }
  }
}

writeLatexTable = function(table, filename){
  fileConn = file(filename)
  lines = c()
  lines = c(lines, paste("\\begin{tabular}{", paste0(rep("|l", ncol(table)+1), collapse = ""), "|}", sep = ""))
  lines = c(lines, paste0("\\hline"))
  lines = c(lines, paste0(paste("&",paste(colnames(table), collapse = " & "), "\\\\",  "\\hline ")))
  for(i in 1:nrow(table)){
    newline = paste0(paste(rownames(table)[i], "&", paste(table[i,], collapse = " & "),  "\\\\", "\\hline"))
    lines = c(lines, newline)
  }
  lines = c(lines, "\\end{tabular}")
  writeLines(lines, con = fileConn)
  close(fileConn)
}



