###BOOTSTRAP FROR RUNNING QUALITY CONTROL METRICS ON BATCH CORRECTED DATASETS


#'@import vegan
#'@import limma
#'@import kBET
#'@import FNN
#'@import DESeq2
#'@import sva
#'
#'@importFrom stats setNames
#'@importFrom stats sd
#'@importFrom utils write.csv2
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 geom_boxplot
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 aes
#'@importFrom reshape2 melt
#'@importFrom grDevices pdf
#'@importFrom grDevices dev.off


## Bootstrapping the genes in the data
### Wrapper for many datasets
#requireNamespace("kBET")
requireNamespace("FNN")
requireNamespace("limma")
requireNamespace("DESeq2")
requireNamespace("sva")


runMetrics = function(data, sample_types, batch, y = NULL,distMatrix=NULL, metrics, zeroRows = FALSE){

  tmp = list()
  for(m in metrics){
    if(m %in% c("F-score", "Davies-Bouldin", "kBET", "Silhouette", "kNN", "mindist")){
      dist = TRUE
    } else {
      dist = FALSE
    }
    if(dist){
      if(m == "F-score"){
        metric = getFscore
      } else if(m == "Davies-Bouldin"){
        metric = DaviesBouldinIndex
      } else if(m == "kBET"){
        metric = getkBET
      } else if(m == "Silhouette"){
        metric = getSilhouette
      } else if(m == "kNN"){
        metric = kNN_proportions
      } else if(m == "mindist"){
        metric = avgMinDist
      }
      if(is.null(distMatrix)){
        distMatrix = GetDistMatrix(data)
      }
      metrBio = metric(distMatrix, sample_types)
      metrBatch = metric(distMatrix, batch)
    } else {
      samplefactor = as.factor(as.numeric(as.factor(sample_types)))
      bfactor = as.factor(as.numeric(as.factor(batch)))

      if(m == "avedist"){
        metric = avedist
      } else if(m == "kldist"){
        metric = kldist
      } else if(m == "sepscore"){
        metric = sepscore
      } else if(m == "skewdiv"){
        metric = skewdiv
      } else if(m == "pvca"){
        metric = pvcam
      }



      if(m == "pvca"){
        metrBio = metric(t(data), samplefactor, y)
        metrBatch = metric(t(data), bfactor, y)
      } else {
        metrBio = metric(t(data), samplefactor)
        metrBatch = metric(t(data), bfactor)
      }

      if(m == "gPCA"){
        metrBio = gPCA_percentage(data, samplefactor)
        metrBatch = gPCA_percentage(data, bfactor)
      }

    }
    metrBio = as.numeric(metrBio)
    metrBatch = as.numeric(metrBatch)
    tmp[[m]] = c(bio = as.numeric(metrBio), batch = as.numeric(metrBatch), ratio = metrBio/metrBatch)
  }
  names(tmp) = metrics
  return(tmp)
}

QC_bootstrap = function(data_list, biol.groups, batches, method_names = NULL, metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN", "KL-divergence", "Silhouette", "mindist"), dist_method = "pearson", iters = 50, savefile = FALSE, filename = "evaluations.csv", plot = FALSE, zeroRows = FALSE, y = NULL){

  if(!is.list(data_list)){
    data_list = list(data_list)
  }



  ### the datasets should have EQUAL NUMBER OF ROWS and contain no NA values
  x_row = nrow(data_list[[1]])
  for(d in data_list){
    if(nrow(d) != x_row){
      stop("Error: datasets should have equal number of rows")
    }
    if(any(is.na(d))){
      stop("Error: datasets should not have NA values")
    }
    if(any(is.infinite(d))){
      stop("Error: datasets should contain only finite values")
    }
  }

  n = length(data_list)

  if(is.null(method_names)){
    if(!is.null(names(data_list))){
      method_names = names(data_list)
    } else {
      method_names = paste("method", as.character(1:n))
    }
  } else {
    if(length(data_list) != length(method_names)){
      stop("The list of method names has to be of the same length as the data list")
    }
  }

  r = cor(as.numeric(as.factor(biol.groups)), as.numeric(as.factor(batches)))
  if(abs(r) > 0.2){
    warning("Biological groups may be imbalanced across batches, results may be skewed")
  }

  tmp_bio = list()
  tmp_batch = list()
  tmp_ratio = list()


  for(k in 1:iters){

    ##sample the rows of the data
    #rowI = sample.int(x_row, replace = TRUE)
    newdata_list = list()
    dists = list()

    #metricsD = list()
    ##generate new datasets
    for(d in 1:n){

      if(zeroRows == FALSE){
        nZeroI = which(rowSds(data_list[[d]]) > 1e-6)
      } else {
        nZeroI = 1:nrow(newdata)
      }
      rowI = sample(nZeroI,size = length(nZeroI), replace = TRUE)
      newdata = data_list[[d]][rowI,]
      newdata_list[[d]] = newdata
      if(any(is.na(newdata))){
        print("data include NAs")
        newdata = na.omit(newdata)
      }

      dists[[d]] = GetDistMatrix(newdata, dist_method = dist_method)


    }

    for(m in metrics){
      if(m == "F-score"){

        ###Compute F-scores for biol. and batch signal
        fbatch = getFscores(dists, batches)
        fbio = getFscores(dists, biol.groups)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], fbio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], fbatch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], fbio/fbatch)
      }

      if(m == "Davies-Bouldin"){

        ###Compute Davies-Bouldin indices

        db_batches = DaviesBouldinScores(dists, batches)
        dbs = DaviesBouldinScores(dists, biol.groups)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], dbs)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], db_batches)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dbs/db_batches)
      }

      if(m == "Chi-square"){
        ###Compute chi-square p-values for clustered data

        chisq_bio = getChisq(dists, biol.groups)
        chisq_batch = getChisq(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], chisq_bio$pvals)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], chisq_batch$pvals)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], chisq_bio$pvals/chisq_batch$pvals)

      }

      if(m == "kBET"){
        ###Compute kBET rejection rates

        batch.estimates = getkBET(dists, batches)
        bio.estimates = getkBET(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], bio.estimates$rej.rates)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], batch.estimates$rej.rates)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], bio.estimates$rej.rates/batch.estimates$rej.rates)
      }

      if(m == "KL-divergence"){
        ###Compute Kullback-Leibler divergence of between and within cluster densities

        kld_bio = getKLdistances(dists, biol.groups)
        kld_batch = getKLdistances(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], kld_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], kld_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], kld_bio/kld_batch)
      }

      if(m == "Silhouette"){
        silh_bio = getSilhouettes(dists, biol.groups)
        silh_batch = getSilhouettes(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], silh_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], silh_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], silh_bio/silh_batch)
      }

      if(m == "kNN"){
        props_batch = kNN_proportions(dists, batches)

        props_bio = kNN_proportions(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], props_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], props_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], props_bio/props_batch)
      }


      if(m == "mindist"){

        dist_bio = getMinDists(dists, biol.groups)
        dist_batch = getMinDists(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], dist_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], dist_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dist_bio/dist_batch)

      }

      if(m %in% c("avedist", "kldist","sepscore","skewdiv","pvca")){
        bFactor = as.factor(as.numeric(as.factor(batches)))
        groupFactor = as.factor(as.numeric(as.factor(biol.groups)))

        if(m == "avedist"){
          aveDists_bio = sapply(newdata_list, function(d) avedist(t(d), groupFactor))
          aveDists_batch = sapply(newdata_list, function(d) avedist(t(d), bFactor))


          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

        if(m == "kldist"){
          print("computing kldist")
          t = proc.time()

          klDists_bio = sapply(newdata_list, function(d){
            kldist(t(d), groupFactor)
          })
          klDists_batch = sapply(newdata_list, function(d) kldist(t(d), bFactor))
          print((proc.time() - t)[3])

          tmp_bio[[m]] = rbind(tmp_bio[[m]], klDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], klDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], klDists_bio/klDists_batch)
        }

        if(m == "sepscore"){
          print("computing sepscore")
          t = proc.time()
          sepDists_bio = sapply(newdata_list, function(d) sepscore(t(d), groupFactor))
          sepDists_batch = sapply(newdata_list, function(d) sepscore(t(d), bFactor))
          print((proc.time() - t)[3])

          tmp_bio[[m]] = rbind(tmp_bio[[m]], sepDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], sepDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], sepDists_bio/sepDists_batch)
        }

        if(m == "skewdiv"){
          print("computing skewdiv")
          t = proc.time()
          aveDists_bio = sapply(newdata_list, function(d) skewdiv(t(d), groupFactor))
          aveDists_batch = sapply(newdata_list, function(d) skewdiv(t(d), bFactor))
          print((proc.time() - t)[3])

          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

        if(m == "pvca"){
          aveDists_bio = sapply(newdata_list, function(d) pvcam(t(d), groupFactor, y = as.factor(as.numeric(as.factor(y)))))
          aveDists_batch = sapply(newdata_list, function(d) pvcam(t(d), bFactor, y = as.factor(as.numeric(as.factor(y)))))


          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

      }


      ###to-do: add some other metrics (principal variance components, etc.)

    }



    ###compute the evaluation metric
  print(paste(k, "/", iters))
  }

  #tmp_bio = lapply(tmp_bio, setNames, nm = method_names)
  #tmp_batch = lapply(tmp_batch, setNames, nm = method_names)
  #tmp_ratio = lapply(tmp_ratio, setNames, nm = method_names)
  #names(tmp_bio) = metrics
  #names(tmp_batch) = metrics
  #names(tmp_ratio) = metrics


  qc_evaluations = data.frame(matrix(nrow = n, ncol = 0))
  types = c("bio", "batch", "ratio")

  for(i in 1:length(metrics)){
    for(tmp in list(tmp_bio, tmp_batch, tmp_ratio)){
      qc_evaluations = cbind(qc_evaluations, mean = colMeans(tmp[[i]], na.rm = TRUE))
      qc_evaluations = cbind(qc_evaluations, sd = apply(tmp[[i]], 2, sd))
    }
  }


  #str = paste(as.vector(sapply(metrics, function(x) rep(x, 4))), as.vector(sapply(c("batch", "bio"), function(x) rep(x,2))))
  #str = paste(str, c("mean", "sd"))

  row.names(qc_evaluations) = method_names
  colnames(qc_evaluations) = paste(rep(metrics, each=6), paste(colnames(qc_evaluations), rep(types, each=2)))

  #qc_evaluations = cbind(method_class, qc_evaluations)

  #reslist = list(FScores = FScores, Davies_Bouldin = dbScores, Chisq_test = chisq, Chisq_p = pvalue, kBET_p = kbet_pval, rej.rate = kbet_rej)

  reslist = list(biol.signal = tmp_bio, batch.signal = tmp_batch, ratio = tmp_ratio)


  #for(tmp in reslist){
    #tmp = lapply(tmp, setNames, nm = method_names)
  #}

  if(savefile){
    write.csv2(qc_evaluations, filename)
  }

  if(plot == TRUE){
    boxPlotPDF(reslist, "Batch metric evaluations boxplot.pdf")
  }


  return(list(means = qc_evaluations, results = reslist, method_class = method_names))

}


###functions for plotting results
requireNamespace("reshape2")
requireNamespace("ggplot2")


boxPlot = function(res_data, title = "metric", method_class){

  n = ncol(res_data)
  corr_method = as.factor(1:n)
  df = data.frame(corr_method, method_class, t(res_data))

  ##sort by method_class, raw data first

  df = rbind(df[1,], df[-1,][order(df[-1,]$method_class),])

  df = melt(df, id.vars = c("corr_method", "method_class"), variable.name = 'run')
  p = ggplot(df, aes(x = corr_method, y = value, fill = method_class)) + geom_boxplot() + ggtitle(label = title)


  return(p)
}

boxPlots = function(object){

  reslist = object$results
  method_class = object$method_class
  if(is.null(method_class)){
    method_class = object$means$method_class
  }
  metrics = names(reslist[[1]])

  title.ind = 1

  type.ind = 1
  types = c("for biological signal", "for batch signal", "ratio")

  for(tmp in reslist){
    title.ind = 1
    for(l in tmp){
      p = boxPlot(l, title = paste(metrics[title.ind], types[type.ind]), method_class)
      print(p)
      title.ind = title.ind + 1
    }
    type.ind = type.ind + 1
  }
}

boxPlotPDF = function(object, filename){
  pdf(file = filename, paper="a4r", height=9, width=11.5)
  boxPlots(object)
  dev.off()
}

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



getRanks = function(metrics_list, type = c("bio", "batch", "ratio")){


  metrics = names(metrics_list$results[[1]])
  ret = list()
  i = 1

  mean_results = metrics_list$means

  for(t in type){

    ranks = data.frame(row.names = rownames(mean_results))

    for(j in 1:length(metrics)){
      ord = order(mean_results[,paste(metrics, "mean", t)[j]], decreasing = !(t == "batch"))

      rank = 1:nrow(mean_results)
      rank[ord] <- 1:nrow(mean_results)
      ranks = cbind(ranks, rank)
    }

    colnames(ranks) = metrics

    ret[[i]] = ranks
    i = i+1
  }

  names(ret) = c(paste(type, "ranks"))
  return(ret)

}

diluteSeries = function(Data, batch, groups = NULL, corrData = NULL, corrData2 = NULL, batch2 = NULL, perc = seq(0,1, by = 0.1), corr_method = "ComBat"){

  if(is.null(corrData)){
    if(corr_method == "ComBat"){
    if(!is.null(groups)){
    	mod = GenerateDesignMatrices(groups)
    } else {
    	mod = NULL
    }
    corrData = ComBat(Data, batch = batch, mod = mod)
    } else {
    	stop("no other methods implemented yet!")
    }
  }
  if(is.null(batch2)){
    corrections1 = GetMeanCorrections(Data, corrData, batch)
    lapply(perc, function(elem) ApplyBatchCorrections(Data, batch, corrections1, p = elem))
  } else {
    corrections1 = GetMeanCorrections(Data, corrData, batch)
    corrections2 = GetMeanCorrections(corrData, corrData2, batch2)
    lapply(sqrt(perc), function(elem){

      ApplyBatchCorrections(ApplyBatchCorrections(Data,
                                                  batch,
                                                  corrections1, p = elem),
                                                            batch2,
                                                            corrections2,
                                                            p = elem)
      })
  }
}

ADS_A = function(data, samples, p = 1){
  sampleInds = lapply(unique(samples), function(i) which(samples == i))
  #means = sapply(sampleInds, function(I) data[,I] - rowMeans(data[,I]))
  for(I in sampleInds){
    data[,I] = (data[,I] - p*rowMeans(data[,I]))
  }
  #colnames(means) = unique(samples)
  return(data)
}

diluteSeries2 = function(data, samples, perc = seq(0,1, by = 0.1)){
  
  lapply(perc, function(p) ADS_A(data, samples, p = p))
}



#Generate datasets with variation, with deseq-normalization and batch-correction

QC_resample = function(CountData,coldata = NULL, batches, groups, metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist", "kldist"),
                       iters = 100, corrMethod = "combat", mod = NULL, Ctrl_sample, y = NULL, levels = seq(0,1, by = 0.1), dilute_samples = FALSE){

  tmp_bio = list()
  tmp_batch = list()
  tmp_ratio = list()

  n = length(levels)

  if(is.null(coldata)){
    coldata = data.frame(batch = batches, group = groups)

    rownames(coldata) = colnames(CountData)
  }

  for(i in 1:iters){
    newCounts = data.frame(matrix(nrow = nrow(CountData), ncol = 0))


    ##Sample the new count data
    for(j in 1:ncol(CountData)){
      probs = (gtools::rdirichlet(1, CountData[,j]))[1,]
      #probs = CountData[,j]/(sum(CountData[,j]))
      newCounts = cbind(newCounts, rmultinom(1, size = sum(CountData[,j]), probs))
    }

    rownames(newCounts) = rownames(CountData)
    colnames(newCounts) = colnames(CountData)

    DESeq_new = DESeq2_Wrapper(newCounts,
                               coldata, "group",
                               Design="batch + group", CtrlSample = Ctrl_sample)

    #print("error happens here?")
    vst_new = assay(vst(DESeq_new))
    #print("1")

    if(corrMethod == "combat"){
      correction = ComBat(vst_new, batch = batches, mod = mod)
    } #else if (corrMethod == "clustercenter"){
      #correction = ClusterCenter(vst_new, batches)
    #}
    #print("2")
    if(dilute_samples == TRUE){
      series = diluteSeries2(correction, groups, perc = levels)
    } else {

      series = diluteSeries(vst_new, correction, batches, perc = levels)
    }
    print("3")
    series = lapply(series, function(d) d[which(rowSds(d) > 1e-6),])
    dists = lapply(series, function(d) GetDistMatrix(d, dist_method = "pearson"))
    print("4")

    #metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN", "KL-divergence", "Silhouette", "kNN", "mindist")



    biol.groups = groups
    scaledF = FALSE
    print("computing metrics")
    for(m in metrics){
      if(m == "F-score" | m == "scaled F-score"){

        if(m == "scaled F-score"){
          scaledF = TRUE
        }

        ###Compute F-scores for biol. and batch signal
        fbatch = getFscores(dists, batches)
        fbio = getFscores(dists, biol.groups)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], fbio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], fbatch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], fbio/fbatch)
      }

      if(m == "Davies-Bouldin"){

        ###Compute Davies-Bouldin indices

        db_batches = DaviesBouldinScores(dists, batches)
        dbs = DaviesBouldinScores(dists, biol.groups)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], dbs)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], db_batches)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dbs/db_batches)
      }

      if(m == "Chi-square"){
        ###Compute chi-square p-values for clustered data

        chisq_bio = getChisq(dists, biol.groups)
        chisq_batch = getChisq(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], chisq_bio$pvals)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], chisq_batch$pvals)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], chisq_bio$pvals/chisq_batch$pvals)

      }

      if(m == "kBET"){
        ###Compute kBET rejection rates

        batch.estimates = getkBET(dists, batches)
        bio.estimates = getkBET(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], bio.estimates$rej.rates)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], batch.estimates$rej.rates)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], bio.estimates$rej.rates/batch.estimates$rej.rates)
      }

      if(m == "KL-divergence"){
        ###Compute Kullback-Leibler divergence of between and within cluster densities

        kld_bio = getKLdistances(dists, biol.groups)
        kld_batch = getKLdistances(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], kld_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], kld_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], kld_bio/kld_batch)
      }

      if(m == "Silhouette"){
        silh_bio = getSilhouettes(dists, biol.groups)
        silh_batch = getSilhouettes(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], silh_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], silh_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], silh_bio/silh_batch)
      }

      if(m == "kNN"){
        props_batch = kNN_proportions(dists, batches)

        props_bio = kNN_proportions(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], props_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], props_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], props_bio/props_batch)
      }


      if(m == "mindist"){

        dist_bio = getMinDists(dists, biol.groups)
        dist_batch = getMinDists(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], dist_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], dist_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dist_bio/dist_batch)

      }

      if(m %in% c("avedist", "kldist","sepscore","skewdiv","pvca", "gPCA")){
        bFactor = as.factor(as.numeric(as.factor(batches)))
        groupFactor = as.factor(as.numeric(as.factor(biol.groups)))
        newdata_list = series

        if(m == "avedist"){
          aveDists_bio = sapply(newdata_list, function(d) avedist(t(d), groupFactor))
          aveDists_batch = sapply(newdata_list, function(d) avedist(t(d), bFactor))


          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

        if(m == "kldist"){

          klDists_bio = sapply(newdata_list, function(d){
            kldist(t(d), groupFactor)
          })
          klDists_batch = sapply(newdata_list, function(d) kldist(t(d), bFactor))

          tmp_bio[[m]] = rbind(tmp_bio[[m]], klDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], klDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], klDists_bio/klDists_batch)
        }

        if(m == "sepscore"){
          sepDists_bio = sapply(newdata_list, function(d) sepscore(t(d), groupFactor))
          sepDists_batch = sapply(newdata_list, function(d) sepscore(t(d), bFactor))

          tmp_bio[[m]] = rbind(tmp_bio[[m]], sepDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], sepDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], sepDists_bio/sepDists_batch)
        }

        if(m == "skewdiv"){
          aveDists_bio = sapply(newdata_list, function(d) skewdiv(t(d), groupFactor))
          aveDists_batch = sapply(newdata_list, function(d) skewdiv(t(d), bFactor))


          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

        if(m == "pvca"){
          aveDists_bio = sapply(newdata_list, function(d) pvcam(t(d), groupFactor, y = as.factor(as.numeric(as.factor(y)))))
          aveDists_batch = sapply(newdata_list, function(d) pvcam(t(d), bFactor, y = as.factor(as.numeric(as.factor(y)))))


          tmp_bio[[m]] = rbind(tmp_bio[[m]], aveDists_bio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], aveDists_batch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], aveDists_bio/aveDists_batch)
        }

        if(m == "gPCA"){
          gpcaBio = sapply(newdata_list, function(d) gPCA_percentage(d, biol.groups, nperm = 250))
          gpcaBatch = sapply(newdata_list, function(d) gPCA_percentage(d, batches, nperm = 250))

          tmp_bio[[m]] = rbind(tmp_bio[[m]], gpcaBio)
          tmp_batch[[m]] = rbind(tmp_batch[[m]], gpcaBatch)
          tmp_ratio[[m]] = rbind(tmp_ratio[[m]], gpcaBio/gpcaBatch)
        }

      }


      ###to-do: add some other metrics (principal variance components, etc.)

    }

    print(paste(i, "/", iters, "done"))

  }

  method_names = paste(levels, corrMethod)

  tmp_bio = lapply(tmp_bio, setNames, nm = method_names)
  tmp_batch = lapply(tmp_batch, setNames, nm = method_names)
  tmp_ratio = lapply(tmp_ratio, setNames, nm = method_names)

  names(tmp_batch) = names(tmp_bio) = names(tmp_ratio) = metrics

  qc_evaluations = data.frame(matrix(nrow = n, ncol = 0))
  types = c("bio", "batch", "ratio")

  for(i in 1:length(metrics)){
    for(tmp in list(tmp_bio, tmp_batch, tmp_ratio)){
      qc_evaluations = cbind(qc_evaluations, mean = colMeans(tmp[[i]], na.rm = TRUE))
      qc_evaluations = cbind(qc_evaluations, sd = apply(tmp[[i]], 2, sd))
    }
  }

  row.names(qc_evaluations) = method_names
  colnames(qc_evaluations) = paste(rep(metrics, each=6), paste(colnames(qc_evaluations), rep(types, each=2)))

  reslist = list(biol.signal = tmp_bio, batch.signal = tmp_batch, ratio = tmp_ratio)
  return(list(means = qc_evaluations, results = reslist, method_class = method_names))

}

#wrapper function for running metric evaluations on generated variations of data

QC_wrapper = function(CountData, batch, group, y=NULL, metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist", "kldist", "gPCA") , iters = 100, perc = seq(0,1, by = 0.1), corrData = NULL, deseqData = NULL, parallel = FALSE, var_measure = c("bootstrap", "resampling", "bootstrapSamples", "resampleSamples"), cores = 4){
  
  if(is.null(y)){
    y = group
  }
  
  if(!is.character(group)){
    group = as.character(group)
  }
  
  sampleData = data.frame(batch = batch, group = group)
  
  rownames(sampleData) = colnames(CountData)
  
  if(is.null(deseqData)){
  
    deseqData = DESeq2_Wrapper(CountData, sampleData, "group", "batch + group", "control")
  
  }
  
  Data = assay(vst(deseqData))
  
  if(is.null(corrData)){
  
    corrData = ComBatAL(Data, batch, mod = model.matrix(~group))$data
  
  }

  series = diluteSeries(Data, corrData, batch, perc = perc)
  
  print("ADS 2 corrected")
  series2 = diluteSeries2(corrData, group, perc = perc)
  
  if(parallel == TRUE){
    print("using parallel computation:")
    library(parallel)
    reslist = mclapply(var_measure, function(c){
      if(c == "bootstrap"){
        res = QC_bootstrap(series, group, batch, metrics = metrics, iters = iters, y = y)
      }
      if(c == "resampling"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = "normal", mod =GenerateDesignMatrices(group), y = y, levels = perc)
      }
      if(c == "bootstrapSamples"){
        res = QC_bootstrap(series2, group, batch, metrics = metrics, iters = iters, y = y) 
      }
      if(c == "resampleSamples"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = "normal", mod =GenerateDesignMatrices(group), dilute_samples = TRUE, y = y, levels = perc)
      }
      
      write_results(res, filename = paste(c, "txt", sep = "."))
      return(res)
    }, mc.cores = cores)
    names(reslist) = var_measure
    return(reslist)
  } else {

    reslist = lapply(var_measure, function(c){
      if(c == "bootstrap"){
        res = QC_bootstrap(series, group, batch, metrics = metrics, iters = iters, y = y)
      }
      if(c == "resampling"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = "normal", mod =GenerateDesignMatrices(group), y = y, levels = perc)
      }
      if(c == "bootstrapSamples"){
        res = QC_bootstrap(series2, group, batch, metrics = metrics, iters = iters, y = y) 
      }
      if(c == "resampleSamples"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = "normal", mod =GenerateDesignMatrices(group), dilute_samples = TRUE, y = y, levels = perc)
      }
      
      write_results(res, filename = paste(c, "txt", sep = "."))
      return(res)
    })
    names(reslist) = var_measure
    return(reslist)
  }
}
