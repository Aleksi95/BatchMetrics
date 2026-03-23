###BOOTSTRAP FROR RUNNING QUALITY CONTROL METRICS ON BATCH CORRECTED DATASETS


#'@import limma
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
requireNamespace("kBET")
requireNamespace("FNN")
requireNamespace("limma")
requireNamespace("DESeq2")
requireNamespace("sva")
requireNamespace("parallel")
requireNamespace("SummarizedExperiment")
requireNamespace("MatrixGenerics")
requireNamespace("SingleCellExperiment")
requireNamespace("scuttle")

evalBatchEffect = function(data, sample_types, batch, metric = "Davies-Bouldin",  y = NULL,distMatrix=NULL, zeroRows = FALSE){
  m = metric
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
    return(c(bio = as.numeric(metrBio), batch = as.numeric(metrBatch), ratio = metrBio/metrBatch))
}


runMetrics = function(data, sample_types, batch, y = NULL,distMatrix=NULL, metrics, zeroRows = FALSE){

  tmp = list()
  for(m in metrics){
   
    tmp[[m]] = evalBatchEffect(data, sample_types, batch, y = NULL,distMatrix=NULL, metric = m, zeroRows = FALSE)
  }
  names(tmp) = metrics
  return(tmp)
}

QC_bootstrap = function(data_list, biol.groups, batches, method_names = NULL, metrics = c("F-score", "Davies-Bouldin", "kNN", "KL-divergence", "Silhouette", "mindist"), dist_method = "pearson", iters = 50, savefile = FALSE, filename = "evaluations.csv", plot = FALSE, zeroRows = FALSE, y = NULL){

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

         if(m == "kBET"){
    ###Compute kBET rejection rates

        batch.estimates = sapply(newdata_list, function(d) kBET(t(d), batches)$summary[1,2])
        bio.estimates = sapply(newdata_list, function(d) kBET(t(d), biol.groups)$summary[1,2])

        tmp_bio[[m]] = rbind(tmp_bio[[m]], bio.estimates)
        tmp_batch[[m]] = rbind(tmp_batch[[m]], batch.estimates)
        tmp_ratio[[m]] = rbind(tmp_ratio[[m]], bio.estimates/batch.estimates)
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


  return(reslist)

}

QC_bootstrap_par = function(data_list, biol.groups, batches, method_names = NULL, metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN", "KL-divergence", "Silhouette", "mindist"), dist_method = "pearson", scaledF = FALSE, iters = 50, nCores=16, savefile = FALSE, filename = "evaluations.csv", plot = FALSE, zeroRows = FALSE, y = NULL){

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


  #registerDoParallel(numCores)
  #foreach(k=1:iters)%dopar%
   result = mclapply(1:iters, function(i) {

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
      #newdata = data_list[[d]][rowI,]
      newdata_list[[d]] = data_list[[d]][rowI,]
      if(any(is.na(newdata_list[[d]]))){
        print("data include NAs")
        newdata_list[[d]] = na.omit(newdata_list[[d]])
      }

      dists[[d]] = GetDistMatrix(newdata_list[[d]], dist_method = dist_method)


    }

     tmp_bio = list()
    tmp_batch = list()
    tmp_ratio = list()
    
    #biol.groups = groups
    scaledF = FALSE
    series = newdata_list
    print("computing metrics")

    for(m in metrics){
      if(m == "F-score" | m == "scaled F-score"){
        
        if(m == "scaled F-score"){
          scaledF = TRUE
        }
        
        ###Compute F-scores for biol. and batch signal
        fbatch = getFscores(dists, batches, scaled = scaledF)
        fbio = getFscores(dists, biol.groups, scaled = scaledF)
        tmp_bio[[m]] = fbio
        tmp_batch[[m]] = fbatch
        tmp_ratio[[m]] = fbio/fbatch
      }
      
      if(m == "Davies-Bouldin"){
        
        ###Compute Davies-Bouldin indices
        
        db_batches = DaviesBouldinScores(dists, batches)
        dbs = DaviesBouldinScores(dists, biol.groups)
        tmp_bio[[m]] = dbs
        tmp_batch[[m]] = db_batches
        tmp_ratio[[m]] = dbs/db_batches
      }
      
      if(m == "Chi-square"){
        ###Compute chi-square p-values for clustered data
        
        chisq_bio = getChisq(dists, biol.groups)
        chisq_batch = getChisq(dists, batches)
        
        tmp_bio[[m]] = chisq_bio$pvals
        tmp_batch[[m]] =  chisq_batch$pvals
        tmp_ratio[[m]] = chisq_bio$pvals/chisq_batch$pvals
        
      }
      
      if(m == "kBET"){
        ###Compute kBET rejection rates
        
        batch.estimates = sapply(series, function(d) kBET(t(d), batches)$summary[1,2])
        bio.estimates = sapply(series, function(d) kBET(t(d), biol.groups)$summary[1,2])
        
        tmp_bio[[m]] = bio.estimates
        tmp_batch[[m]] = batch.estimates
        tmp_ratio[[m]] =  bio.estimates/batch.estimates
      }
      
      if(m == "KL-divergence"){
        ###Compute Kullback-Leibler divergence of between and within cluster densities
        
        kld_bio = getKLdistances(dists, biol.groups)
        kld_batch = getKLdistances(dists, batches)
        
        tmp_bio[[m]] = kld_bio
        tmp_batch[[m]] =  kld_batch
        tmp_ratio[[m]] = kld_bio/kld_batch
      }
      
      if(m == "Silhouette"){
        silh_bio = getSilhouettes(dists, biol.groups)
        silh_batch = getSilhouettes(dists, batches)
        
        tmp_bio[[m]] = silh_bio
        tmp_batch[[m]] = silh_batch
        tmp_ratio[[m]] = silh_bio/silh_batch
      }
      
      if(m == "kNN"){
        props_batch = kNN_proportions(dists, batches)
        
        props_bio = kNN_proportions(dists, biol.groups)
        
        tmp_bio[[m]] =  props_bio
        tmp_batch[[m]] = props_batch
        tmp_ratio[[m]] = props_bio/props_batch
      }
      
      
      if(m == "mindist"){
        
        dist_bio = getMinDists(dists, biol.groups)
        dist_batch = getMinDists(dists, batches)
        
        tmp_bio[[m]] =  dist_bio
        tmp_batch[[m]] =  dist_batch
        tmp_ratio[[m]] = dist_bio/dist_batch
        
      }
      
      if(m %in% c("avedist", "kldist","sepscore","skewdiv","pvca", "gPCA")){
        bFactor = as.factor(as.numeric(as.factor(batches)))
        groupFactor = as.factor(as.numeric(as.factor(biol.groups)))
        
        
        if(m == "avedist"){
          aveDists_bio = sapply(series, function(d) avedist(t(d), groupFactor))
          aveDists_batch = sapply(series, function(d) avedist(t(d), bFactor))
          
          
          tmp_bio[[m]] =  aveDists_bio
          tmp_batch[[m]] =  aveDists_batch
          tmp_ratio[[m]] =  aveDists_bio/aveDists_batch
        }
        
        if(m == "kldist"){
          
          klDists_bio = sapply(series, function(d){
            kldist(t(d), groupFactor)
          })
          klDists_batch = sapply(series, function(d) kldist(t(d), bFactor))
          
          tmp_bio[[m]] =  klDists_bio
          tmp_batch[[m]] =  klDists_batch
          tmp_ratio[[m]] =  klDists_bio/klDists_batch
        }
        
        if(m == "sepscore"){
          sepDists_bio = sapply(series, function(d) sepscore(t(d), groupFactor))
          sepDists_batch = sapply(series, function(d) sepscore(t(d), bFactor))
          
          tmp_bio[[m]] =  sepDists_bio
          tmp_batch[[m]] = sepDists_batch
          tmp_ratio[[m]] =  sepDists_bio/sepDists_batch
        }
        
        if(m == "skewdiv"){
          aveDists_bio = sapply(series, function(d) skewdiv(t(d), groupFactor))
          aveDists_batch = sapply(series, function(d) skewdiv(t(d), bFactor))
          
          
          tmp_bio[[m]] =  aveDists_bio
          tmp_batch[[m]] =  aveDists_batch
          tmp_ratio[[m]] = aveDists_bio/aveDists_batch
        }
        
        if(m == "pvca"){
          aveDists_bio = sapply(series, function(d) pvcam(t(d), groupFactor, y = as.factor(as.numeric(as.factor(y)))))
          aveDists_batch = sapply(series, function(d) pvcam(t(d), bFactor, y = as.factor(as.numeric(as.factor(y)))))
          
          
          tmp_bio[[m]] = aveDists_bio
          tmp_batch[[m]] = aveDists_batch
          tmp_ratio[[m]] =  aveDists_bio/aveDists_batch
        }
        
        if(m == "gPCA"){
          gpcaBio = sapply(series, function(d) gPCA_percentage(d, biol.groups, nperm = 250))
          gpcaBatch = sapply(series, function(d) gPCA_percentage(d, batches, nperm = 250))
          
          tmp_bio[[m]] =  gpcaBio
          tmp_batch[[m]] =  gpcaBatch
          tmp_ratio[[m]] =  gpcaBio/gpcaBatch
        }
        
      }
      
      
      ###to-do: add some other metrics (principal variance components, etc.)
      
    }
    
    list(bio = tmp_bio,
         batch = tmp_batch,
         ratio = tmp_ratio)
    
  }, mc.cores = nCores)
  
  #stopCluster(cl)
  print(result[[1]])
  #result is a list with 100 elements, each element has 3 lists with results for m metric each
  bio_results = lapply(result, function(r) r[["bio"]])
  batch_results = lapply(result, function(r) r[["batch"]])
  ratio_results = lapply(result, function(r) r[["ratio"]])
  
  reslist = lapply(list(biol.signal = bio_results, batch.signal = batch_results, ratio = ratio_results), function(b){
    metrics = names(b[[1]])
    res = lapply(metrics, function(m){
      restable = NULL
      #print(length(b))
      for(i in 1:length(b)){
        #print(b[[i]])
        if(!is.null(b[[i]])){
          #print(b[[i]][[m]])
          restable = rbind(restable, b[[i]][[m]])
        }
      }
      restable
    })
    names(res) = metrics
    res
  })
  
  names(reslist) = c("biol.signal", "batch.signal", "ratio")
  
  return(reslist)

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


QC_resample = function(CountData, coldata = NULL, batches, groups, 
                           metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist", "kldist", "skewdiv"),
                           iters = 100, corrMethod = "ComBat", mod = NULL, Ctrl_sample, y = NULL, 
                           levels = seq(0,1, by = 0.1), dilute_samples = FALSE, 
                           design = NULL, normalization = "deseq", parallel = TRUE, nCores = 16, 
                           usePCA = FALSE, nPCs = 50, var.model = "multinom"){
  
  if(is.null(coldata)){
    coldata = data.frame(batch = batches, group = groups)
    rownames(coldata) = colnames(CountData)
  }
  
  # Ensure CountData is a matrix for faster access
  CountData <- as.matrix(CountData)
  
  if(var.model == "nbinom"){
    if(is.null(design)){
      design = as.formula("~batch + group")
    }
    #estimate parameters using DESEq2
    suppressMessages({
      DESeq_d = DESeqDataSetFromMatrix(CountData, coldata, design=design)
      DESeq_d = estimateSizeFactors(DESeq_d)
      DESeq_d = estimateDispersions(DESeq_d)
      
    })
  }
  # --- OPTIMIZATION 1: Pre-calculate Dirichlet parameters if possible or optimize loop ---
  # Depending on library, doing this inside the parallel worker might be safer, 
  # but we define the worker function cleanly below.
  
  run_iteration <- function(i) {
    
    
    # --- OPTIMIZATION 2: Vectorized Resampling (Removing the cbind loop) ---
    # Instead of growing a dataframe, generate all columns and combine once.
    # checking if we can use matrix operations. 
    
    n_genes <- nrow(CountData)
    n_samples <- ncol(CountData)
    
    # Pre-allocate matrix for speed
    newCounts <- matrix(0, nrow = n_genes, ncol = n_samples)
    colnames(newCounts) <- colnames(CountData)
    rownames(newCounts) <- rownames(CountData)
    
    l_theta = function(x, alpha) (x + alpha)/(sum(x) + alpha*length(x))
    
    for(j in 1:n_samples){
      # gtools::rdirichlet returns a matrix, we take the first row
      ##add Laplace smoothing?
      #theta = l_theta(CountData[,j], alpha = 1)
      if(var.model == "multinom"){
      probs <- gtools::rdirichlet(1, CountData[,j])[1,] 
      newCounts[,j] <- rmultinom(1, size = sum(CountData[,j]), prob = probs)
      } else if(var.model == "nbinom"){
        
        #estimate parameters using DESEq2
        newCounts[,j] = rnbinom(n = nrow(DESeq_d), mu = assays(DESeq_d)$mu[,j], size = dispersions(DESeq_d))
      }
    }
    
    if(anyNA(newCounts)) newCounts[is.na(newCounts)] = 0
    
    # --- Normalization ---
    if(normalization == "deseq"){
      if(is.null(design)){
        design = as.formula("~batch + group")
      }
      # Suppress messages to keep parallel logs clean
      suppressMessages({
        DESeq_new = DESeqDataSetFromMatrix(newCounts, coldata, design=design)
        DESeq_new = estimateSizeFactors(DESeq_new)
        DESeq_new = estimateDispersions(DESeq_new)
        vst_new = assay(vst(DESeq_new))
      })
    } else if(normalization == "scuttle"){
      se <- SingleCellExperiment(assays=list(counts = newCounts), colData=DataFrame(coldata))
      se = logNormCounts(se)
      vst_new = logcounts(se) # Result is a matrix
    }
    
    # 1. Keep genes with variance > 0
    keep_genes <- which(rowSds(vst_new) > 1e-6)
    vst_subset <- vst_new[keep_genes, ]
    
    # --- Batch Correction ---
    if(corrMethod == "ComBat"){
      # Use tryCatch to handle ComBat failures gracefully
      correction_subset <- tryCatch({
        sva::ComBat(dat = vst_subset, batch = batches, mod = mod)
      }, error = function(e) {
        return(vst_subset) # Fallback: No correction if ComBat fails
      })
    } else if (corrMethod == "limma"){
      correction_subset = limma::removeBatchEffect(vst_subset, batch = batches, group = factor(groups))
    } else if (corrMethod == "harmony"){
      usePCA = FALSE
      pca_res <- prcomp(t(vst_subset), rank. = nPCs, scale. = TRUE, center = TRUE)
      vst_subset <- t(pca_res$x) 
      correction_subset = t(harmony::RunHarmony(vst_subset, coldata, 'batch', verbose = FALSE))
    }
    
    if(all(abs(correction_subset - vst_subset) < 1e-10)) { warning("ComBat did not modify data") }
    
    
    vst_new <- vst_subset
    correction <- correction_subset
    
    # we process one level 'p', calculate metrics, and discard the data.
    
    biol.groups = groups
    bFactor = as.factor(as.numeric(as.factor(batches)))
    groupFactor = as.factor(as.numeric(as.factor(biol.groups)))
    
    # Initialize storage for this iteration
    iter_bio = list()
    iter_batch = list()
    iter_ratio = list()
    
    dilute_all = function(data, correct, p) (p)*correct + (1-p)*data
    
    # Loop over ADS levels sequentially
    for(ix in seq_along(levels)){
      
      p <- levels[ix]
      
      
      d_p = dilute_all(vst_new, correction, p)
      
      if(dilute_samples == TRUE){
        d_p <- ADS_A(correction, groups, p = p) 
      } 
      
      # Filter low variance rows (memory efficient subsetting)
      #if(corrMethod != "harmony"){ 
       # d_p <- d_p[which(rowSds(d_p) > 1e-6), , drop=FALSE]
      #}
      # 2. Dimensionality Reduction (Optional but HIGHLY recommended for SC)
      #print(usePCA)
      if(usePCA && nrow(d_p) > nPCs){
        # Calculate PCA and use it for distances/kBET
        pca_res <- prcomp(t(d_p), rank. = nPCs, scale. = TRUE, center = TRUE)
        data_for_dist <- t(pca_res$x) # Transpose because GetDistMatrix expects Genes x Samples
        # Note: kBET usually expects Samples x PCs. 
        # You might need to adjust kBET calls if using PCA to pass t(data_for_dist)
      } else {
        data_for_dist <- d_p
      }
      #print(dim(d_p))
      # 3. Calculate Distance Matrix (Once per level)
      # Force use of fast correlation if available in your Metrics.R logic
      dist_p <- GetDistMatrix(data_for_dist, dist_method = "pearson", fastCor = TRUE)
     
      if(corrMethod != "harmony"){ 
        doPCA = FALSE
      } else {
        doPCA = usePCA
      }
      
      # 4. Compute Metrics for this level
      # (Logic copied and adapted to work on single item instead of list)
      print("computing metrics")
      for(m in metrics){
        res_bio <- NA
        res_batch <- NA
        
        # Distance-based metrics
        if(m %in% c("F-score", "Davies-Bouldin", "Chi-square", "KL-divergence", "Silhouette", "kNN", "mindist")){
          if(m == "F-score") {
            res_batch = getFscore(dist_p, batches, scaled = FALSE)
            res_bio = getFscore(dist_p, biol.groups, scaled = FALSE)
          } else if(m == "Davies-Bouldin") {
            res_batch = DaviesBouldinIndex(dist_p, batches)
            res_bio = DaviesBouldinIndex(dist_p, biol.groups)
          } else if(m == "kNN") {
            res_batch = getKNNprop(dist_p, batches) # Assuming non-list wrapper exists or using sapply logic
            res_bio = getKNNprop(dist_p, biol.groups)
          } else if(m == "mindist") {
            res_batch = avgMinDist(dist_p, batches)
            res_bio = avgMinDist(dist_p, biol.groups)
          }
          # Add other distance metrics here...
        }
        
        # Data-based metrics (kBET, CMS, etc)
        # Note: kBET is slow. If using PCA, pass t(data_for_dist) which is Samples x PCs
        if(m == "kBET"){
          # kBET expects Samples x Features
          k_data <- t(data_for_dist) 
          
      
          res_batch = getkBET(k_data, batches, do.pca = doPCA)
          res_bio = getkBET(k_data, biol.groups, do.pca = doPCA)
        }
        
        if(m == "CMS"){
          # CMS usually runs on SCE object
          res_batch = getCMS(d_p, batches, doPCA = doPCA)
          res_bio = getCMS(d_p, biol.groups, doPCA = doPCA)
        }
        
        if(m == "kldist"){
          
   
          res_bio = kldist(t(data_for_dist), groupFactor)
          res_batch = kldist(t(data_for_dist), bFactor)
         
        }
        
     
        
        if(m == "gPCA"){
          
          res_bio = gPCA_percentage(d_p, biol.groups, nperm = 250)
          res_batch = gPCA_percentage(d_p, batches, nperm = 250)
          
        }
        
        # Store results
        if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
        if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
        if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
        
        # Find index of p
        #idx <- which(levels == p)
        iter_bio[[m]][ix] <- res_bio
        iter_batch[[m]][ix] <- res_batch
        iter_ratio[[m]][ix] <- res_bio / res_batch
      }
      
      # Clean up heavy objects immediately
      #print(ix)
      rm(d_p, dist_p, data_for_dist)
      gc(verbose=FALSE) 
    }
    print(paste(i, "/", iters))
    return(list(bio = iter_bio, batch = iter_batch, ratio = iter_ratio))
  }

	if(parallel){
  # Run in parallel
  	result = mclapply(1:iters, run_iteration, mc.cores = nCores)
  } else {
  	result = lapply(1:iters, run_iteration)
  }
  return(result)
}


DESeq2_Wrapper = function(CountData, coldata, designCols, ctrl_sample){
  
  
  design = as.formula(paste("~",paste0(designCols, collapse = "+"), sep = ""))
  
  deseq_mat_1 = DESeqDataSetFromMatrix(countData = CountData, colData = coldata, design = design) #create the DESeq-dataset
  deseqData1 = estimateSizeFactors(deseq_mat_1, quiet=TRUE) #estimate size factors for training data
  deseqData1 = estimateDispersions(deseqData1, quiet=TRUE) 
  return(deseqData1)
}

#wrapper function for running metric evaluations on generated variations of data

QC_wrapper = function(CountData, batch, group, y=NULL, metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist", "kldist", "gPCA") , iters = 100, perc = seq(0,1, by = 0.1), corrData = NULL, deseqData = NULL, parallel = FALSE, var_measure = c("bootstrap", "resampling", "bootstrapSamples", "resampleSamples"), cores = 4, design = NULL,  Ctrl_sample = "normal", rnaSeq = "bulk"){
  
  if(is.null(y)){
    y = group
  }
  
  if(!is.character(group)){
    group = as.character(group)
  }
  
  sampleData = data.frame(batch = batch, group = group)
  
  rownames(sampleData) = colnames(CountData)
  
  
  if(rnaSeq == "bulk"){
  
  if(is.null(deseqData)){
  
    deseqData = DESeq2_Wrapper(CountData, sampleData, "group", "batch + group", "normal")
  
  }
  
  Data = assay(vst(deseqData))
  
  normalization = "deseq"
  
  } else if (rnaSeq == "sc"){
    print("normalizing single cell data")
    se <- SingleCellExperiment(assays=list(counts = CountData),
                               colData=DataFrame(sampleData))
    
    se = logNormCounts(se)
    Data = logcounts(se)
    normalization = "scuttle"
  }
  
  if(is.null(corrData)){
  
    corrData = ComBatAL(Data, batch, mod = model.matrix(~group))$data
  
  }

  series = diluteSeries(Data, corrData, batch, perc = perc)
  
  print("ADS 2 corrected")
  series2 = diluteSeries2(corrData, group, perc = perc)
  
  if(parallel != TRUE){
   
    
    reslist = lapply(var_measure, function(c){
      if(c == "bootstrap"){
        res = QC_bootstrap(series, group, batch, metrics = metrics, iters = iters, y = y)
      }
      if(c == "resampling"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample =Ctrl_sample, mod =GenerateDesignMatrices(group), y = y, levels = perc, design = design, normalization = normalization)
      }
      if(c == "bootstrapSamples"){
        res = QC_bootstrap(series2, group, batch, metrics = metrics, iters = iters, y = y) 
      }
      if(c == "resampleSamples"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = Ctrl_sample, mod =GenerateDesignMatrices(group), dilute_samples = TRUE, y = y, levels = perc, design = design,normalization = normalization)
      }
      
      #write_results(res, filename = paste(c, "txt", sep = "."))
      return(res)
    })
    names(reslist) = var_measure
    return(reslist)
  } else {
 print("using parallel computation:")
 library(parallel)
    reslist = lapply(var_measure, function(c){
         if(c == "bootstrap"){
        res = QC_bootstrap_par(series, group, batch, metrics = metrics, iters = iters, y = y, nCores = cores)
      }
      if(c == "resampling"){
        print("parallel")
        res = QC_resample_par(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample =Ctrl_sample, mod =GenerateDesignMatrices(group), y = y, levels = perc, design = design, normalization = normalization)
      }
      if(c == "bootstrapSamples"){
        res = QC_bootstrap_par(series2, group, batch, metrics = metrics, iters = iters, y = y, nCores = cores) 
      }
      if(c == "resampleSamples"){
        res = QC_resample_par(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample = Ctrl_sample, mod =GenerateDesignMatrices(group), dilute_samples = TRUE, y = y, levels = perc, design = design, normalization = normalization)
      }
      
      
      #write_results(res, filename = paste(c, "txt", sep = "."))
      return(res)
    })
    names(reslist) = var_measure
    return(reslist)
  }
}
