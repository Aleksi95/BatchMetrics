#############################################
# I collect here codes that are related to
# Batch correction.
#
#############################################
#'@importFrom sva ComBat
#'@importFrom Harman harman
#'@importFrom Harman reconstructData
#'@importFrom stats kmeans
#'@importFrom stats sd
#'@importFrom stats var
#'@importFrom limma plotMDS
#'@importFrom grDevices pdf



BatchCorr_WithSubset <- function(DataIn, batch, SubsetIndex, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE){
   # Run ComBat on a subset of the data (use mean.only, by default)
   # Apply the generated correction on the rest of the data

   if(is.null(SampleTypes) | length(unique(SampleTypes[SubsetIndex]))== 1){
      ModelToCombat=NULL
   } else{
      ModelToCombat=GenerateDesignMatrices(SampleTypes[SubsetIndex])
   }
   CorrMethod = tolower(CorrMethod)
   if(CorrMethod =='combat'){
      ComBatOut <- ComBat(DataIn[,SubsetIndex], batch=batch[SubsetIndex],
                              mod=ModelToCombat, mean.only=mean.only)
      Corrections <- GetMeanCorrections(DataIn[,SubsetIndex], ComBatOut,
                                        batch=batch[SubsetIndex], DoSD=!mean.only )
   } else {
      if(CorrMethod =='harman'){
         if(is.null(SampleTypes)){
           stop("Error: harman requires sample types")
         }
         HarmanModel <- harman(DataIn[,SubsetIndex], batch=batch[SubsetIndex],
                               expt=SampleTypes[SubsetIndex])
	       HarmanOut <- reconstructData(HarmanModel)
	       mean.only=FALSE
	       Corrections <- GetMeanCorrections(DataIn[,SubsetIndex], HarmanOut,
                                           batch=batch[SubsetIndex], DoSD=!mean.only )
	       Corrections = removeZeroSD(Corrections, batches = batch[SubsetIndex])
      }


     if(CorrMethod == "samplecenter"){
       mean.only = TRUE
       CenterModel = SamplesCenter(DataIn[,SubsetIndex], samples = SampleTypes[SubsetIndex], batch = batch[SubsetIndex])
       Corrections = GetMeanCorrections(DataIn[,SubsetIndex], CenterModel, batch=batch[SubsetIndex], DoSD=!mean.only)
     }

     if(CorrMethod == "clustercenter"){

       CenterModel = ClusterCenter(DataIn[,SubsetIndex], batch = batch[SubsetIndex])
       Corrections = GetMeanCorrections(DataIn[,SubsetIndex], CenterModel, batch=batch[SubsetIndex], DoSD=!mean.only)
     }
   }
   #pcorrs_means = lapply(Corrections[[1]], function(d) p*d)
   #pcorrs_sds = lapply(Corrections[[2]], function(d) (p*(d-1) + 1))

   #Corrections = list(mean.vals = pcorrs_means, sd.vals = pcorrs_sds)


   CorrectedData <- ApplyBatchCorrections(DataIn, batch, Corrections, p, mean.only = mean.only)

   return(CorrectedData)

}





GetMeanCorrections <- function(OldData, NewData, batch, DoSD=TRUE){
   # Extract average shifts, created
   # in batch correction.
   # DoSD = Run Standard Deviation extraction on rows

   UniqueBatches = unique(batch)

   output = list(c())
   output2 = list(c())
   for(k in 1:length(UniqueBatches)){
     tmp_ind = which(batch == UniqueBatches[k])
     output[[k]] = NewData[,tmp_ind] - OldData[,tmp_ind]
     if(DoSD){
        output2[[k]] = cbind( apply( NewData[,tmp_ind], 1, sd),
                              apply( OldData[,tmp_ind], 1, sd))
     }
   }
   names(output) = UniqueBatches
   if(DoSD){
     names(output2) = UniqueBatches
   }
   return(list(mean.vals=output, sd.vals=output2))
}

###this function is to deal with exceptions with harman

removeZeroSD = function(corr_data, batches){
  uniq = unique(batches)
  sds = list()
  for(i in 1:length(uniq)){
    sd_vals = corr_data$sd.vals[[i]]
    sd_vals = apply(sd_vals, 1, function(r){
      if(r[1] != 0 & r[2] == 0 | r[1] == 0 & r[2] != 0){
        return(rep(0, 2))
      } else {
        return(r)
      }
    })
    sds[[i]] = (t(matrix(unlist(sd_vals), nrow = 2)))

  }
  names(sds) = uniq
  return(list(mean.vals = corr_data[[1]], sd.vals = sds))
}



FindMatchIndex <- function(vector1, vector2){
   # Which cases of vector1 are found in vector2
   # Return an index vector of these positions
   output = c()
   for(k in 1:length(vector1)){
      if(any(vector2 == vector1[k])){
        output = c(output, k)
      }
   }
   return(output)
}

CreateScalingFactor <- function(vector1, vector2){
   # Create a scaling factor.
   # Calculus is X=vector1/vector2
   # Created a separate function to handle all the exceptions
   #
   # 1st Exception is 0/0 cases. I set these to be 1
   # 2nd Exception: All X/0 cases are set as NA (X > 0)
   # 2nd Exception: These we should not have

   if(length(vector1) != length(vector2)){
     stop('Error in CreateScalingFactor. Inputs do not match')
   }
   output = rep(0, length(vector1))
   NA.ind = which(vector1 > 0 & vector2 == 0)
   one.ind = which(vector1 == 0 & vector2 == 0)
   output[NA.ind] = NA
   output[one.ind] = 1
   exclude = c(NA.ind, one.ind)
   output[-exclude] = vector1[-exclude]/vector2[-exclude]
   return(output)
}

ApplyBatchCorrections <- function(DataIn, Batches, corrections, p=1, mean.only = FALSE){
   # Code for creating batch corrections
   # This uses correction created on separate set
   # as a model
   # corrections = list of to lists.
   # 1st list is shift correction. 2nd list is scale correction

   # Added feature: p= percentage of batch correction (Aleksi)

   CtrlSubset = NULL

   Output = DataIn
   dat_sz = dim(DataIn)
   BatchNames = names(corrections[[1]])
   for(k in 1:length(corrections[[1]])){
      #Average profile
      tmp_aver = rowMeans(corrections[[1]][[k]]*p)
      tmp_ind = which(Batches == BatchNames[k])
      if(!is.null(CtrlSubset)){
        tmp_ind2 = FindMatchIndex(tmp_ind, CtrlSubset)
	      row.mean.v = rowMeans(DataIn[,tmp_ind[tmp_ind2]])
      } else{
        tmp_ind2 = NULL

        if(ncol(DataIn[,tmp_ind]) > 1){
	        row.mean.v = rowMeans(DataIn[,tmp_ind])
        } else {
          row.mean.v = DataIn[,tmp_ind]
        }
      }
      ShiftVector = (row.mean.v + tmp_aver)

      if(mean.only){
        ScaleVector = rep(1, length(corrections[[1]][[1]][,1]))
      } else {
        ScaleVector = CreateScalingFactor(corrections[[2]][[k]][,1]*p,
                                          corrections[[2]][[k]][,2]*p)
      }


      Output[,tmp_ind] = ReAdjustData(DataIn[,tmp_ind],
                                      ScaleVector, ShiftVector,
				      CtrlSubset=tmp_ind2)
   }
   return(Output)
}

ReAdjustData <- function(DataIn, ScaleVector, ShiftVector, p=1, CtrlSubset=NULL){
   # Rescale and shift the processed data
   # Code is for batch correction
   # ScaleVector is ratio between current and desired stand.dev.
   # ShiftVector contains the desired mean values

   dat_sz = dim(DataIn)
   if(length(ScaleVector) != dat_sz[1] | length(ShiftVector) != dat_sz[1]){
      stop('Inputs for RescaleData do not match in size')
   }
   row.means <- rowMeans(DataIn)
   if(is.null(CtrlSubset)){
     NewData <- DataIn - p*matrix(rowMeans(DataIn),dat_sz[1],1) %*% rep(1,dat_sz[2])
   } else{
     NewData <- DataIn - p*matrix(rowMeans(DataIn[,CtrlSubset]),dat_sz[1],1) %*% rep(1,dat_sz[2])
   }
   NewData <- NewData * ((p*(matrix(ScaleVector,dat_sz[1],1)-1) + 1) %*% rep(1,dat_sz[2]))
   NewData <- NewData + ( p*matrix(ShiftVector,dat_sz[1],1) %*% rep(1,dat_sz[2]) )
   return(NewData)
}

##############################################
# This code aims to generate design matrices
# from sample classifications
# This is for removeBatchEffect
##############################################

GenerateDesignMatrices <- function(SampleClassifications){


  if(is.vector(SampleClassifications)){
    DesignOut <- model.matrix(~SampleClassifications)[,-1,drop=FALSE]
    return(DesignOut)
  } else {
    dat_sz <- dim(SampleClassifications)
  }
  DesignOut <- matrix(0, dat_sz[1], 0)
  for(k in 1:dat_sz[2]){
     tmp = model.matrix(~SampleClassifications[,k])[,-1,drop=FALSE]
     DesignOut <- cbind(DesignOut, tmp)
  }
  return(DesignOut)
}




####
##This function uses means in each biological group of interest to center data and perform batch correction
## variance correction does not yet work
## Gives OK results


SamplesCenter = function(data, batch, samples, p = 1){

  mean.only = TRUE

  ###VARIANCE CORRECTION DOES NOT CURRENTLY WORK

  testcorr_data = data
  batchInds = lapply(unique(batch), function(j) which(batch == j))
  bioInds = lapply(unique(samples), function(j) which(samples == j))

  ###Center each biological group
  for(I in bioInds){
    for(j in batchInds){
      Ind = intersect(I, j)
      #print(Ind)
      if(length(Ind)>1){
        testcorr_data[,Ind] = (data[,Ind] - p*rowMeans(data[,Ind]))
        if(!mean.only){
          rowsds = apply(data[,Ind], 1, sd)
          testcorr_data[,Ind] = testcorr_data[,Ind]/rowsds
        }
      } else {
        testcorr_data[,Ind] = data[,Ind] - p*data[,Ind]
      }
    }
    testcorr_data[,I] = testcorr_data[,I] + p*rowMeans(data[,I])

    if(!mean.only){
      rowsds = apply(data[,I], 1, sd)
      testcorr_data[,I] = testcorr_data[,I]*rowsds
    }
  }

  return(testcorr_data)
}


##this correction method uses k-means clustering to define batch cluster centers, which are then used to normalize each batch
##gives a reasonable batch correction

ClusterCenter = function(data, batch, mean.only = FALSE, p = 1){

  testcorr_data = data
  batchInds = lapply(unique(batch), function(j) which(batch == j))


  init_centers = NULL
  batch.vars = NULL
  for(j in batchInds){
    if(length(j)> 1){
      init_centers = cbind(init_centers, rowMeans(data[,j]))
      if(!mean.only){
        rowsds = apply(data[,j], 1, sd)
        batch.vars = cbind(batch.vars, rowsds)
      }
    } else{
      init_centers = cbind(init_centers, data[,j])
      if(!mean.only){
        batch.vars = cbind(batch.vars, rep(1, length(data[,j])))
      }
    }
  }

  batch.vars[batch.vars == 0] <- 1

  clusts = kmeans(t(data), centers = t(init_centers))

  centers = t(clusts$centers)

  i = 1
  for(j in batchInds){
    if(length(j) > 1){
      if(!mean.only){
        testcorr_data[,j] = (p*(rowMeans(batch.vars) - 1) + 1)*((testcorr_data[,j] - p*rowMeans(data[,j]))/(p*(batch.vars[,i] - 1) +1)) + p*rowMeans(centers)

        i = i+1
      } else {
        testcorr_data[,j] = (testcorr_data[,j] - p*rowMeans(data[,j])) + p*rowMeans(centers)
      }
    } else {
      testcorr_data[,j] = (testcorr_data[,j] - p*(data[,j])) + p*rowMeans(centers)
    }
  }

  return(testcorr_data)

}

BatchCorr2 = function(dataIn, batch1, batch2, samples = NULL, subset1 = NULL, subset2 = NULL, p = 1, method = "combat"){
  method = tolower(method)
  if(!is.null(subset1)){
    if(!(method == "harman" | method == "combat")){
      stop("Error: subset correction defined only for combat and harman")
    }
    corr1 = BatchCorr_WithSubset(dataIn, batch1, subset1, p, SampleTypes = samples, CorrMethod = method)
    corr2 = BatchCorr_WithSubset(corr1, batch2, subset2, p, SampleTypes = samples, CorrMethod = method)
  } else {
    if(!is.null(samples)){
      mod = GenerateDesignMatrices(samples)
      if(method == "combat"){
        corr1 = ComBat(dataIn, batch1, mod)
        corr2 = ComBat(corr1, batch2, mod)
      }
      if(method == "harman"){
        harman1a = harman(dataIn, expt = samples, batch = batch1)
        harman2 = harman(reconstructData(harman1a), expt = samples, batch = batch2)
        corr2 = reconstructData(harman2)
      }
      if(method == "samplescenter"){
        corr1 = SamplesCenter(dataIn, batch1, samples, p)
        corr2 = SamplesCenter(corr1, batch2, samples, p)
      }

    } else {
      if(method == "harman"){
        stop("Error: harman requires sample types")
      }

      if(method == "samplescenter"){
        stop("Error: samplescenter requires sample types")
      }

      if(method == "combat"){
        corr1 = ComBat(dataIn, batch1)
        corr2 = ComBat(corr1, batch2)
      }

      if(method == "clustercenter"){
        corr1 = ClusterCenter(dataIn, batch1, mean.only = TRUE, p)
        corr2 = ClusterCenter(corr1, batch2, mean.only = TRUE, p)
      }


    }
  }
  return(corr2)
}

###Functions for plotting data

plotMDS_andLegend <- function(data, sample_types, LegendChoice="bottomleft",
                              SelectScores = NULL, MaxNum = 500, UseMax=TRUE, ...){
  # Create MDS plot
  # Define colors for sample types
  # Add legend to plot
  # LegendChoice selects where the legend goes
  # SelectScores = Scores used to select the used genes
  #

  #lls()
  if(!is.null(SelectScores)){
    if(length(SelectScores) != nrow(data)){
      stop('data and SelectScores do not match each other')
    }
    if(UseMax){
      data = data[order(SelectScores, decreasing=TRUE)[1:MaxNum],]
    } else{
      data = data[order(SelectScores)[1:MaxNum],]
    }
  }


  UniqNames = unique(sample_types)
  Matches = match(sample_types, UniqNames)
  NumbersTmp = list(Numbering=Matches, Keys=UniqNames)

  ColorNames = c('red','green','blue','black','gray','yellow',
                 'magenta','cyan','orange','brown','blueviolet',
                 'purple','pink','maroon','olivedrab','lightsteelblue',
                 'deeppink','palegreen','red4','green4','blue4')
  plotMDS(data, col=ColorNames[NumbersTmp[[1]]], ...)

  legend(LegendChoice, legend=NumbersTmp[[2]], fill=ColorNames[1:length(NumbersTmp[[2]])])

}

DoManyMDS_plots <- function(data, SampleTypeTable, ColVector, SelectScores = NULL, toPDF = FALSE, filename = NULL, ...){

  # This is wrapper for printing many MDS plots
  # in same document

  if(is.null(SelectScores)){
    SelectScores  = apply(data, 1, var)
  }
  if(toPDF){
    filename = paste(filename, "MDSplot.pdf", sep = "_")
    pdf(file=filename,paper="a4r", height=9, width=11.5)
  }
  for(k in 1:length(ColVector)){
    plotMDS_andLegend(data, SampleTypeTable[,ColVector[k]], SelectScores = SelectScores, ...)
  }
  if(toPDF){
    dev.off()
  }

}





