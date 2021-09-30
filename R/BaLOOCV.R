#'@importFrom sva ComBat
#'@importFrom Harman harman
#'@importFrom Harman reconstructData
#'@importFrom limma removeBatchEffect

#'@param dataIn
#'@param batch1
#'@param batch2
#'@param sample_types
#'@param subsetInd1
#'@param subsetInd2
#'@param method
#'@param mse
#'@param true_params

#'@return a list of results

#source("BatchCorrectionRelatedCodes.R")
#library(sva)
#library(limma)
#library(Harman)


BaLOOCV = function(dataIn, batch1, batch2 = NULL, sample_types = NULL, subsetInd1 = NULL, subsetInd2 = NULL, method = "combat", mse = FALSE, true_params = NULL){


  test_data = NULL
  mean_mindists = c()
  mean_mses_mu = c()
  mean_mses_delta = c()
  if(!is.null(batch2)){
    batch1 = as.character(batch1)
    batch2 = as.character(batch2)
    batches = paste(batch1, batch2)
  } else {
    batches = as.character(batch1)
  }


  for(j in 1:ncol(dataIn)){
    X = dataIn
    Xtrain = X[,-j]
    Xtest = X[,j]
    batch = batches[-j]
    batchfirst = batch1[-j]
    if(!is.null(batch2)){
      batchsecond = batch2[-j]
    }
    if(!is.null(subsetInd1)){
      subset1 = subsetInd1[-j]
    }
    if(!is.null(subsetInd2)){
      subset2 = subsetInd2[-j]
    }

    #batch2 = JoinedSampleData[-j,22]
    #batch1 = JoinedSampleData[-j,16]
    if(!is.null(sample_types)){
      samples = sample_types[-j]
    } else {
      samples = NULL
    }
    #samples = sample_types[-j]
    #subset2 = c(TechReplInds[[1]],TechReplInds[[2]])[-j]

    if(method == "combat"){
      model = ComBat(Xtrain, batch)
    }
    if(method == "combat2x"){
      if(!is.null(batch2)){
        model = BatchCorr2(Xtrain, batchfirst, batchsecond, samples)
      } else {
        stop("Error: need second batch!")
      }
    }
    if(method == "samplecenter"){
      model = SamplesCenter(Xtrain, batch, samples)
    }
    if(method == "clustercenter"){
      model = ClusterCenter(Xtrain, batch)
    }
    if(method == "subsetcombat"){
      if(!is.null(batch2)){
        model = BatchCorr2(Xtrain, batchfirst, batchsecond, samples, subset1, subset2)
      } else {
        model = BatchCorr_WithSubset(Xtrain, batch, subset1, SampleTypes = samples)
      }
    }

    if(method == "harman"){
      if(!is.null(batch2)){
        model = BatchCorr2(Xtrain, batchfirst, batchsecond, samples, method = "harman")
      } else {
        model = harman(Xtrain, expt = samples, batch)
        model = reconstructData(model)
      }
    }
    if(method == "subsetharman"){
      if(!is.null(batch2)){
        model = BatchCorr2(Xtrain, batchfirst, batchsecond, samples, subset1, subset2, method = "harman")
      } else {
        model = BatchCorr_WithSubset(Xtrain, batch, subset1, SampleTypes = samples, CorrMethod = "harman")
      }
    }

    if(method == "limma"){
      model = removeBatchEffect(Xtrain, batch1, batch2, design = GenerateDesignMatrices(samples))
    }

    corrections = GetMeanCorrections(Xtrain, model, batch, DoSD = TRUE)
    if(method == "harman"){
      corrections = removeZeroSD(corrections, batch)
    }

    if(mse == TRUE & !is.null(true_params)){
      mu_mses = c()

      delta_mses = c()

      for(i in 1:length(unique(batch))){
        scale_estimate = CreateScalingFactor(corrections[[2]][[i]][,1],
                                             corrections[[2]][[i]][,2])
        scale_true = CreateScalingFactor(true_params$delta[[i]][,1],
                                         true_params$delta[[i]][,2])
        mu_mse = mean((rowMeans(corrections$mean.vals[[i]])- true_params$mu[[i]])^2)
        delta_mse = mean((scale_estimate - scale_true)^2)
        mu_mses = c(mu_mses, mu_mse)
        delta_mses = c(delta_mses, delta_mse)
      }
      mean_mses_mu = c(mean_mses_mu, mean(mu_mses))
      mean_mses_delta = c(mean_mses_delta, mean(delta_mses))

    }
    batchInd = which(unique(batches) == batches[j])
    batchInds = which(batch == batches[j])
    avg = rowMeans(corrections$mean.vals[[batchInd]])
    mean = rowMeans(Xtrain[,batchInds])
    ScaleVector = CreateScalingFactor(corrections[[2]][[batchInd]][,1],
                                      corrections[[2]][[batchInd]][,2])

    test_corrected = Xtest - mean
    test_corrected = test_corrected*ScaleVector
    test_corrected = test_corrected + mean + avg

    X[,j] = test_corrected

    dist = GetDistMatrix(X, dist_method = "pearson")
    mindists = c()
    for(I in which(unique(batches) != batches[j])){
      batchI = which(batches == unique(batches)[I])
      mindists = c(mindists, mean(dist[j, batchI], na.rm = TRUE))
    }
    mean_mindists = c(mean_mindists, mean(mindists))
    test_data = cbind(test_data, test_corrected)
  }

  colnames(test_data) = colnames(X)
  return(list(data = test_data, dists = mean_mindists, mu_mse = mean_mses_mu, delta_mse = mean_mses_delta))

}

#BaLOOCV(assay(VST_data), batches = JoinedSampleData[,22], method = "combat")
