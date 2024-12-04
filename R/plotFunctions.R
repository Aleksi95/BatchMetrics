library(ggplot2)
library(reshape2)
library(gridExtra)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

require('ggplot2')
require('reshape2')
require('patchwork')

linePlot = function(dat, title = "title"){
  
  meanI = seq(1, ncol(dat), 2)
  means = dat[,meanI]
  mean_df = data.frame(corr_level = seq_along(means[, 1]),
                       means)
  mean_df = melt(mean_df, id.vars = 'corr_level', variable.name = 'method', value.name = 'mean')
  
  sd_df = melt(dat[,-meanI], variable.name = 'method1', value.name = 'sd')
  df = cbind(mean_df, sd_df)
  
  p = ggplot(data = df, aes(x = corr_level, y = mean, color = method)) + ggtitle(label = title)
  
  p = p + geom_line() + geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd, color = method1),width=0.1) + theme(legend.position = "bottom")
  
  return(p)
}

std = function(d) apply(d, 1, function(x) (x-min(x))/(max(x)-min(x)))

linePlot2 = function(result_data, filename = "linePlot.png", metrics =NULL, order = "leftToRight", standardize = TRUE){
  
  if(standardize == TRUE){
    result_data = lapply(result_data, function(d) lapply(d, function(X) t(std(X))))
  }
  if(is.null(metrics)){
    metrics = names(result_data[[1]])
  }
  
  plotlist = list()
  for(m in metrics){
    metr_data = lapply(result_data, function(d) d[[m]])
    if(order == "leftToRight"){
      levels = ncol(metr_data[[1]]):1
    } else if(order == "rightToLeft"){
      levels = 1:ncol(metr_data[[1]])
    } else {
      stop("provide ordering of x-axes!")
    }
    
    #medians = sapply(metr_data, colMedians)
    #maxs = sapply(metr_data, colMaxs)
    #mins = sapply(metr_data, colMins)
    
    
    colMedians = function(x) apply(x, 2, median)
    colMaxs = function(x) apply(x, 2, max)
    colMins = function(x) apply(x, 2, min)
    
    mediansdf = data.frame(noise_level= levels, sapply(metr_data, colMedians))
    maxsdf = data.frame(noise_level = levels,sapply(metr_data, colMaxs))
    minsdf = data.frame(noise_level = levels, sapply(metr_data, colMins))
   #q75s = data.frame(corr_level = levels, std(sapply(metr_data, function(d) colQuantiles(d, probs = 0.75))))
    #q25s = data.frame(corr_level = levels, std(sapply(metr_data, function(d) colQuantiles(d, probs = 0.25))))
    
    df = rbind(cbind(quantile = "median", melt(mediansdf, id.vars = 'noise_level', variable.name = 'signal')),
          cbind(quantile = "max", melt(maxsdf, id.vars = 'noise_level', variable.name = 'signal')),
          cbind(quantile = "min", melt(minsdf, id.vars = 'noise_level', variable.name = 'signal')))
          #cbind(quantile = "q75", melt(q75s, id.vars = 'corr_level', variable.name = 'signal')),
          #cbind(quantile = "q25", melt(q25s, id.vars = 'corr_level', variable.name = 'signal')))
    
  #sd_df = melt(dat[,-meanI], variable.name = 'method1', value.name = 'sd')
  
  #df = cbind(mean_df, max = as.vector(maxs), min = as.vector(mins), q75 = as.vector(q75), q25 = as.vector(q25))
  
      
    p = ggplot(data = df, aes(x = noise_level, y = value, linetype = quantile, color = signal)) + ggtitle(label = m) + theme(plot.title = element_text(size=20))
    
    p = p + geom_line() + scale_linetype_manual(name = "",
                                                labels = c('max', 'median', 'min'),
                                                values = c(2,1,2))
    
    if(m != metrics[length(metrics)]){
      p = p + theme(legend.position = "none")#+ geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd, color = method1),width=0.1)
    
    } else {
      legendp = p + theme(legend.position = "bottom", legend.key.size = unit(2, 'cm'), legend.text = element_text(size=15),
                          legend.title = element_blank())
      p = p + theme(legend.position = "none")
    } 
    plotlist[[m]] = p
  
  }
  legend = extract_legend(legendp)
  ln1 = ceiling(length(metrics)/2)
  ln2 = length(metrics) - ln1
  png(filename = filename, width = 1400, height = 900)
  grid.arrange(arrangeGrob(grobs = plotlist, ncol = 3, nrow = 2),
               legend, nrow = 2, heights = c(ln2*10, ln1))
  #ln1 = ceiling(length(metrics)/2)
  #ln2 = length(metrics) - ln1
  #par(mfrow = c(ln1,ln2))
  #multiplot(plotlist = plotlist, cols = 3)
  
  dev.off()
  return(plotlist)
}

plotResults = function(res_data, metrics = c("F-scores", "Davies-Bouldin index", "Chi-square test", "kBET rejection rate")){
  
  dat = res_data$means
  n = ncol(dat)
  
  title.ind = 1
  
  #p = ggplot() + theme_void()
  plotlist = list()
  for(i in seq(1, n, 6)){
    p = linePlot(dat[,i:(i+5)], title = metrics[title.ind])
    plotlist[[title.ind]] = p
    title.ind = title.ind + 1
  }
  return(plotlist)
}

plotMetrics = function(object, filename){
  dat = object$means
  metrics = names(object$results[[1]])
  n = ncol(dat)
  for(m in metrics){
    #largestX = max(sapply(c("bio", "batch", "ratio"), function(t) c(max(dat[[paste(m, "mean", t)]]) - min(dat[[paste(m, "mean", t)]]))))
    for(t in c("bio", "batch", "ratio")){
      x = dat[,paste(m, "mean", t)]
      dat[,paste(m, "mean", t)] = (x-min(x))/(max(x)-min(x))
    }
  }
  png(filename = filename, width = 1400, height = 700)
  ln1 = ceiling(length(metrics)/2)
  ln2 = length(metrics) - ln1
  par(mfrow = c(ln1,ln2))
  plotlist = plotResults(list(means = dat), metrics = metrics)
  multiplot(plotlist = plotlist, cols = 3)
  dev.off()
}


plotDim = function(data, class, method, title = NULL){
  
  if(is.null(dim(data))){
    stop("data must have at least 2 dimensions")
  }
  
  if(method == "MDS"){
    mds = limma::plotMDS(t(data), plot = FALSE)
    x = mds$x
    y = mds$y
  }
  if(method == "tsne"){
    tsne = Rtsne::Rtsne(data, perplexity = (nrow(data) - 1)/3)
    x = tsne$Y[,1]
    y = tsne$Y[,2]
  }
  if(method == "umap"){
    um = umap::umap(data)
    x = um$layout[,1]
    y = um$layout[,2]
  }
  if(method == "PCA"){
    pcs = prcomp(t(data))
    x = pcs$rotation[,1]
    y = pcs$rotation[,2]
  }
  d = data.frame(x = x, y=y, class = class, names = rownames(data))
  ggplot(d, aes(x,y)) + geom_text(aes(label = names, color = class)) + ggtitle(paste(title, ",", method, "plot 2 dimensions"))
}


