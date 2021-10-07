###BOOTSTRAP FROR RUNNING QUALITY CONTROL METRICS ON BATCH CORRECTED DATASETS


#'@import vegan
#'@import limma
#'@import kBET
#'@import FNN
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


QC_bootstrap = function(data_list, biol.groups, batches, method_names = NULL, metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN", "KL-divergence", "Silhouette", "kNN", "mindist"), Fscore_method = "scaled", iters = 50, savefile = FALSE, filename = "evaluations.csv", plot = FALSE){

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


  tmp_bio = list()
  tmp_batch = list()
  tmp_ratio = list()

  for(i in 1:length(metrics)){
    df = data.frame(matrix(nrow = 0, ncol = n))
    tmp_bio[[i]] = df
    tmp_batch[[i]] = df
    tmp_ratio[[i]] = df
  }


  for(k in 1:iters){

    ##sample the rows of the data
    rowI = sample.int(x_row, n=x_row, replace = TRUE)
    newdata_list = list()
    dists = list()


    ##generate new datasets
    for(d in 1:n){
      newdata = data_list[[d]][rowI,]
      newdata_list[[d]] = newdata
      dists[[d]] = GetDistMatrix(newdata, dist_method = "pearson")
      #print(any(is.na(dists[[d]])))
    }



    ###compute the evaluation metrics

    for(m in 1:length(metrics)){
      if(metrics[m] == "F-score"){

        ###Compute F-scores for biol. and batch signal
        fbatch = getFscores(dists, batches, method = Fscore_method)
        fbio = getFscores(dists, biol.groups, method = Fscore_method)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], fbio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], fbatch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], fbio/fbatch)
      }

      if(metrics[m] == "Davies-Bouldin"){

        ###Compute Davies-Bouldin indices

        db_batches = DaviesBouldinScores(dists, batches)
        dbs = DaviesBouldinScores(dists, biol.groups)
        tmp_bio[[m]] = rbind(tmp_bio[[m]], dbs)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], db_batches)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dbs/db_batches)
      }

      if(metrics[m] == "Chi-square"){
        ###Compute chi-square p-values for clustered data

        chisq_bio = getChisq(dists, biol.groups)
        chisq_batch = getChisq(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], chisq_bio$pvals)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], chisq_batch$pvals)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], chisq_bio$pvals/chisq_batch$pvals)

      }

      if(metrics[m] == "kBET"){
        ###Compute kBET rejection rates

        batch.estimates = getkBET(dists, batches)
        bio.estimates = getkBET(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], bio.estimates$rej.rates)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], batch.estimates$rej.rates)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], bio.estimates$rej.rates/batch.estimates$rej.rates)
      }

      if(metrics[m] == "KL-divergence"){
        ###Compute Kullback-Leibler divergence of between and within cluster densities

        kld_bio = getKLdistances(dists, biol.groups)
        kld_batch = getKLdistances(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], kld_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], kld_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], kld_bio/kld_batch)
      }

      if(metrics[m] == "Silhouette"){
        silh_bio = getSilhouettes(dists, biol.groups)
        silh_batch = getSilhouettes(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], silh_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], silh_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], silh_bio/silh_batch)
      }

      if(metrics[m] == "kNN"){
        props_batch = kNN_proportions(dists, batches)

        props_bio = kNN_proportions(dists, biol.groups)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], props_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], props_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], props_bio/props_batch)
      }


      if(metrics[m] == "mindist"){

        dist_bio = getMinDists(dists, biol.groups)
        dist_batch = getMinDists(dists, batches)

        tmp_bio[[m]] = rbind(tmp_bio[[m]], dist_bio)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], dist_batch)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], dist_bio/dist_batch)

      }


      ###to-do: add some other metrics (principal variance components, etc.)

    }

  }

  tmp_bio = lapply(tmp_bio, setNames, nm = method_names)
  tmp_batch = lapply(tmp_batch, setNames, nm = method_names)
  tmp_ratio = lapply(tmp_ratio, setNames, nm = method_names)
  names(tmp_bio) = metrics
  names(tmp_batch) = metrics
  names(tmp_ratio) = metrics


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


  for(tmp in reslist){
    tmp = lapply(tmp, setNames, nm = method_names)
  }

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



getRanks = function(metrics_list, type = c("bio", "batch", "ratio")){


  metrics = names(metrics_list$results[[1]])
  ret = list()
  i = 1

  mean_results = metrics_list$means

  for(t in type){
    #mean_results = data.frame()
    #for(l in metrics_list){
    #  mean_results = rbind(mean_results, l$means[nrow(l$means),])
    #}

    colnames(mean_results) = unlist(strsplit(colnames(mean_results), paste(" mean", t)))
    mean_results = mean_results[,metrics]

    ranks = data.frame(row.names = rownames(mean_results))

    for(j in 1:ncol(mean_results)){
      ord = order(mean_results[,j], decreasing = !(t == "batch"))

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

diluteSeries = function(Data, batch1, batch2 = NULL, corrections1, corrections2 = NULL, perc = seq(0,1, by = 0.1)){
  if(is.null(batch2)){
    lapply(perc, function(elem) ApplyBatchCorrections(Data, batch1, corrections1, p = elem))
  } else {
    lapply(sqrt(perc), function(elem) ApplyBatchCorrections(ApplyBatchCorrections(Data,
                                                                                  batch1,
                                                                                  corrections1, p = elem),
                                                            batch2,
                                                            corrections2,
                                                            p = elem))
  }
}

