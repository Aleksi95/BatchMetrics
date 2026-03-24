
QC_wrapper = function(CountData, batch, group, y=NULL, metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist", "kldist", "gPCA") , iters = 100, perc = seq(0,1, by = 0.1), corrData = NULL, corrMethod = "ComBat", deseqData = NULL, parallel = TRUE, var_measure = c("bootstrap", "resampling", "bootstrapSamples", "resampleSamples"), cores = 4, design = NULL,  Ctrl_sample = "normal", rnaSeq = "bulk", usePCA = TRUE,
                      nPCs = 50, var.model = "multinom"){
  
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
  
    deseqData = DESeq2_Wrapper(CountData, sampleData, "group", "batch + group", Ctrl_sample)
  
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
  

 
    reslist = lapply(var_measure, function(c){
         if(c == "bootstrap"){
           print("bootstrap")
        res = QC_bootstrap(Data, group, batch, metrics = metrics, iters = iters,parallel = parallel, nCores = cores, usePCA = usePCA, nPCs = nPCs, corrMethod = corrMethod)
      }
      if(c == "resampling"){
        print("resampling")
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample =Ctrl_sample, mod =GenerateDesignMatrices(group), y = y, levels = perc, design = design, normalization = normalization, parallel = parallel, nCores = cores, corrMethod = corrMethod, usePCA = usePCA, nPCs = nPCs, var.model = var.model)
      }
      if(c == "bootstrapSamples"){
        res = QC_bootstrap(Data, group, batch, metrics = metrics, iters = iters,parallel = parallel, nCores = cores, usePCA = usePCA, nPCs = nPCs, corrMethod = corrMethod, dilute_samples = TRUE)
      }
      if(c == "resampleSamples"){
        res = QC_resample(CountData, coldata = sampleData, batch, group, metrics = metrics, iters = iters, Ctrl_sample =Ctrl_sample, mod =GenerateDesignMatrices(group), y = y, levels = perc, design = design, normalization = normalization, parallel = parallel, nCores = cores, corrMethod = corrMethod, usePCA = usePCA, nPCs = nPCs, var.model = var.model, dilute_samples = TRUE)
      }
      
      gc()
      #write_results(res, filename = paste(c, "txt", sep = "."))
      return(res)
    })
    names(reslist) = var_measure
    return(reslist)
  
}
