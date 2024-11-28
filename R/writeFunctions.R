
#source("../plotFunctions.R")
library(pspearman)
library(reshape2)

rankCorrelations = function(resdata, rank, type = "ratio"){
  meanCorrs = c()
  pvalues = c()
  for(d in resdata[[type]]){
    #data = data.frame(method_class, t(d))
    
    corrs = apply(t(d), 2, function(x) spearman.test(x, rank)$estimate)
    pvals = apply(t(d), 2, function(x) spearman.test(x, rank)$p.value)
    meanCorrs = c(meanCorrs, mean(corrs))
    pvalues = c(pvalues, mean(pvals))
  }
  names(meanCorrs) = names(resdata[[type]])
  names(pvalues) = names(resdata[[type]])
  return(list(estimates = meanCorrs, p.values = pvalues))
}


Ftest2 = function(data){
  data = as.matrix(data)
  if(!is.null(colnames(data))){
    colnames(data) = NULL
  }
  data = cbind(group = as.factor(1:nrow(data)), as.data.frame(data))
  
  #is.na(data) = sapply(data, is.infinite)
  
  data = na.omit(data)
  
  melt_df = melt(data, id.vars = "group")
  
  fit = aov(value~group, data = melt_df)
  
  return(summary(fit)[[1]][1,"F value"])
}


write_results = function(res_data, filename, type = "ratio"){
  metrics = names(res_data$results[[type]])
  scoresdf = data.frame(matrix(ncol = length(metrics), nrow = 0))
  #colnames(scoresdf) = metrics
  
  fscores1 = sapply(metrics, function(i) Ftest2(t(res_data$results[[type]][[i]])))
  scoresdf = rbind(scoresdf, rankCorrelations(res_data$results, 1:11)$estimates)
  print(fscores1)
  scoresdf = rbind(scoresdf, log(fscores1))
  scoresdf = rbind(scoresdf, rankCorrelations(res_data$results, 1:11)$estimates * log(fscores1))
  colnames(scoresdf) = metrics
  rownames(scoresdf) = c("rho", "logF", "rhoXlogF")
  
  
  write.table(round(scoresdf, digits = 3), file = filename)
  
}

