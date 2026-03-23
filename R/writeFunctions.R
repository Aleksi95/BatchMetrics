
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

rankCorrs2 = function(resdata, metrics = NULL, types = NULL){
  if(is.null(metrics)) metrics = names(resdata[[1]][[1]])
  if(is.null(types)) types = names(resdata[[1]])
  
  sapply(metrics, function(m) {
    sapply(types, function(t){
      mscores = as.vector(sapply(resdata, function(l) l[[t]][[m]]))
      ranks = rep(1:11, length(resdata))
      rcor = cor.test(mscores, ranks, method = "spearman")$estimate * ifelse(t == "batch", -1, 1)
    })
  })
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




