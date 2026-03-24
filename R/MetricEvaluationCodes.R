###BOOTSTRAP FROR RUNNING QUALITY CONTROL METRICS ON BATCH CORRECTED DATASETS


#'@import limma
#'@import FNN
#'@import DESeq2
#'@import sva
#'
#'@importFrom stats setNames
#'@importFrom stats sd
#'@importFrom stats rnbinom
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

QC_bootstrap = function(DataNorm, biol.groups, batches, method_names = NULL, 
                            metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN"), levels = seq(0,1, by = 0.1), 
                            dist_method = "pearson", scaledF = FALSE, iters = 50, corrMethod = "ComBat", dilute_samples = FALSE, parallel = TRUE,
                            nCores=16, zeroRows = FALSE, usePCA = FALSE, nPCs=50){
  
  #if(!is.list(data_list)) data_list = list(data_list)
  coldata = data.frame(batch = batches, group = biol.groups)
  # Define worker function
  run_boot <- function(i) {
    
    iter_bio = list()
    iter_batch = list()
    iter_ratio = list()
    
    keep_genes <- which(rowSds(DataNorm) > 1e-6)
    DataNorm <- DataNorm[keep_genes, ]
    
    
    if(corrMethod == "ComBat"){
      # ComBatAL should be loaded/defined in your environment
      mod = model.matrix(~biol.groups)
      
      correction_subset = ComBat(DataNorm, batch = batches, mod = mod)
    } else if (corrMethod == "limma"){
      correction = limma::removeBatchEffect(DataNorm, batch = batches, group = factor(biol.groups))
    } else if (corrMethod == "harmony"){
      usePCA = FALSE
      pca_res <- prcomp(t(DataNorm), rank. = nPCs, scale. = TRUE, center = TRUE)
      DataNorm <- t(pca_res$x) 
      correction_subset = t(harmony::RunHarmony(DataNorm, coldata, 'batch', verbose = FALSE))
    }
    
    #vst_new <- vst_subset
    correction <- correction_subset
    
    # --- OPTIMIZATION 3: Sequential Processing of Levels ---
    # Instead of generating ALL diluted series and ALL dist matrices (Memory Hog),
    # we process one level 'p', calculate metrics, and discard the data.
    

    bFactor = as.factor(as.numeric(as.factor(batches)))
    groupFactor = as.factor(as.numeric(as.factor(biol.groups)))
    
    dilute_all = function(data, correct, p) (p)*correct + (1-p)*data
    
    # Get correction vectors once (outside the p loop)
    # Assuming GetMeanCorrections and ApplyBatchCorrections are available
    #if(dilute_samples == TRUE){
      # Logic for dilute samples (custom implementation based on provided code)
      # Pre-calculation for diluteSeries2 usually happens internally, 
      # but here we just call the logic per p.
    #} else {
     # corrections1 = GetMeanCorrections(DataNorm, correction, batches)
    #}
    
    
    # Loop over datasets (methods/levels)
    for(d_idx in seq_along(levels)){
      p = levels[d_idx]
      if(dilute_samples == TRUE){
        d_p <- ADS_A(correction, biol.groups, p = p) # Assuming ADS_A is the helper for diluteSeries2
      } else {
        d_p <- dilute_all(DataNorm, correction, p)
      }
      
      data_orig <- d_p
      
      # 1. Bootstrap Rows (Genes)
      if(zeroRows == FALSE){
        nZeroI = which(rowSds(data_orig) > 1e-6)
      } else {
        nZeroI = 1:nrow(data_orig)
      }
      
      rowI = sample(nZeroI, size = length(nZeroI), replace = TRUE)
      newdata = data_orig[rowI, ]
      
      if(any(is.na(newdata))) newdata = na.omit(newdata)
      
      # 2. PCA (Optional)
      if(usePCA){
        pca <- prcomp(t(newdata), rank. = nPCs, scale. = TRUE)
        data_for_calc <- t(pca$x)
        dist_m <- GetDistMatrix(data_for_calc, dist_method = dist_method, fastCor=TRUE)
      } else {
        data_for_calc <- newdata
        dist_m <- GetDistMatrix(newdata, dist_method = dist_method, fastCor=TRUE)
      }
      
      if(corrMethod != "harmony"){ 
        doPCA = FALSE
      } else {
        doPCA = usePCA
      }
      
      
      # 3. Calculate metrics immediately
      for(m in metrics){
        # ... [Insert single-metric calculation logic here] ...
        # e.g.:
        if(m == "F-score"){
          fb = getFscore(dist_m, batches, scaled=scaledF)
          fbi = getFscore(dist_m, biol.groups, scaled=scaledF)
          
          # Append result for this dataset index
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- fbi
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- fb
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- fbi/fb
          
          
        } else if(m == "Davies-Bouldin") {
          res_batch = DaviesBouldinIndex(dist_m, batches)
          res_bio = DaviesBouldinIndex(dist_m, biol.groups)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
        } else if(m == "kNN") {
          res_batch = getKNNprop(dist_m, batches) # Assuming non-list wrapper exists or using sapply logic
          res_bio = getKNNprop(dist_m, biol.groups)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
        } else if(m == "mindist") {
          res_batch = avgMinDist(dist_m, batches)
          res_bio = avgMinDist(dist_m, biol.groups)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
          # ... handle other metrics ...
        }
        
        if(m == "kBET"){
          # kBET expects Samples x Features
          k_data <- t(data_for_calc) 
    
          
          
          res_batch = getkBET(k_data, batches, do.pca = doPCA)
          res_bio = getkBET(k_data, biol.groups, do.pca = doPCA)

          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
        }
        
        if(m == "CMS"){
          # CMS usually runs on SCE object
          res_batch = getCMS(data_for_calc, batches, doPCA)
          res_bio = getCMS(data_for_calc, biol.groups, doPCA)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
        }
        
        if(m == "kldist"){
          
          
          res_bio = kldist(t(data_for_calc), groupFactor)
          res_batch = kldist(t(data_for_calc), bFactor)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <- numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
          
        }
        
        
        
        if(m == "gPCA"){
          
          res_bio = gPCA_percentage(data_for_calc, biol.groups, nperm = 250)
          res_batch = gPCA_percentage(data_for_calc, batches, nperm = 250)
          
          if(is.null(iter_bio[[m]])) iter_bio[[m]] <- numeric(length(levels))
          iter_bio[[m]][d_idx] <- res_bio
          
          if(is.null(iter_batch[[m]])) iter_batch[[m]] <-numeric(length(levels))
          iter_batch[[m]][d_idx] <- res_batch
          
          if(is.null(iter_ratio[[m]])) iter_ratio[[m]] <- numeric(length(levels))
          iter_ratio[[m]][d_idx] <- res_bio/res_batch
          
          
        }
      }
      
      # Explicitly clear memory
      rm(newdata, dist_m, data_for_calc)
    }
    print(paste(i, "/", iters))
    
    return(list(bio = iter_bio, batch = iter_batch, ratio = iter_ratio))
  }

  if(parallel){
  	result = mclapply(1:iters, run_boot, mc.cores = nCores)
  } else {
  	result = lapply(1:iters, run_boot, mc.cores = nCores)
  }
  # ... [Result Aggregation Logic] ...
  return(result)
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
      } else {
		d_p = dilute_all(vst_new, correction, p)

	  }
      
      # Filter low variance rows (memory efficient subsetting)
   
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


