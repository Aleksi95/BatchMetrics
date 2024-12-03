
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#'@importFrom kBET kBET
#'@importFrom FNN get.knn
#'@importFrom stats model.matrix
#'@importFrom stats dist
#'@importFrom stats cor
#'@importFrom stats as.dist
#'@importFrom stats hclust
#'@importFrom stats cutree
#'@importFrom stats density
#'@importFrom stats weighted.mean
#'@importFrom FNN KL.dist
#'@importFrom FNN KL.divergence
#'@importFrom limma plotMDS

###BatchMetrics package 2021
#
##Authors: Petri TÃ¶rÃ¶nen, Aleksi Laiho
#
#
#Functions to measure group separation in batch correction and metrics to measure batch correction

#convert lower triangle matrices to symmetric
toSymmetric = function(dist, lower = TRUE){
  if(lower){
    dist[upper.tri(dist)] <- t(dist)[upper.tri(t(dist))]
  } else {
    dist[lower.tri(dist)] <- t(dist)[lower.tri(t(dist))]
  }
  return(dist)

}

requireNamespace("limma")

GetDistMatrix <- function(Data, dist_method = "pearson"){

  # Function for calculating pair-wise distance matrix
  # of samples
  # USAGE:
  # dist_matrix = GetDistMatrix(expr.data, dist.method)
  #
  # Where
  # expr.data is expression data
  # dist.method is distance method. Current alternatives below
  # ("euclidean","manhattan","pearson","spearman")

  #data should NOT contain NAs, because of the methods used
  if(any(is.na(Data))){
    stop("Data should not contain NAs!")
  }


  # Create DistMatrix
  if(dist_method == "euclidean" | dist_method == "manhattan"){
    #print(paste(dist_method,"distance used"))
    DistMatrix <- as.matrix(dist(t(Data), method = dist_method))
  }
  else if(dist_method == "pearson" | dist_method == "spearman"){
    #print(paste(dist_method,"distance used"))
    DistMatrix <- (1-cor(Data, method = dist_method))
  } else if(dist_method == "MDS"){
    #select 500 top genes with highest log-fold changes to compute dist matrix
    #counter curse of dimensionality
    DistMatrix = toSymmetric(plotMDS(Data, plot = FALSE)$distance.matrix)

  }
  # Do we want include other methods?
  #
  # See Jaskowiak, Campello, Costa (2014) BMC Bioinformatics
  # for alternatives
  return(DistMatrix)

}


SeparateDistTypes <- function(SampleTypes){
  # Create a matrix that separates
  # distances within the cluster and
  # between the clusters
  uniqSamples <- unique(SampleTypes)
  type_matrix <- matrix(0,length(SampleTypes),length(SampleTypes))
  counter = 1
  for( k in 1:length(uniqSamples)){
    ind_vect <- which(SampleTypes == uniqSamples[k])
    type_matrix[ind_vect,ind_vect] = counter
    counter = counter + 1
  }
  diag(type_matrix) = -1  # Remove diagonal from selection
  return(type_matrix)

}




DaviesBouldinIndex <- function(dist, sample_types){

  # Calculate Davies-Bouldin index
  # Davies-Bouldin compares distances in cluster to the nearest cluster
  #
  # Diagonal is removed using trick with na.rm.
  # Therefore data should not include NAs


  DistMatrix = dist
  SampleTypes = sample_types

  dist.size = dim( DistMatrix)
  UniqSamples = unique(SampleTypes)

  # GEnerate average distance between clusters
  # Generate also average of distances inside the cluster
  ClustDists = matrix(0, length(UniqSamples), length(UniqSamples))
  for( k in 1:length(UniqSamples)){
    others = 1:length(UniqSamples)
    others = others[others!=k]
    for( l in others){
      tmp = DistMatrix[ SampleTypes == UniqSamples[k],
                        SampleTypes == UniqSamples[l]]
      ClustDists[k,l] = mean(tmp)
    }
    tmp = DistMatrix[SampleTypes == UniqSamples[k],
                     SampleTypes == UniqSamples[k]]
    diag(tmp) = NA
    ClustDists[k,k] = mean(tmp, na.rm = TRUE)
  }
  score1 = rep(0, length(UniqSamples))
  score2 = score1
  for( k in 1:length(UniqSamples)){
    others = 1:length(UniqSamples)
    others = others[others!=k]
    tmp = -1
    tmp3 = -1
    for(l in others){
      #score 1
      tmp2 = (ClustDists[k,k] + ClustDists[l,l])/
        ClustDists[k,l]
      if(ClustDists[k,l] != 0) {
        if( tmp2 > tmp){
          tmp= tmp2
        }
        tmp2 = ClustDists[k,k]/ClustDists[k,l]
        if( tmp2 > tmp3){
          tmp3= tmp2
        }

      }

    }
    score1[k] = tmp
    score2[k] = tmp3
  }
  # BElow I vary how the mean is taken from DB-score
  output = list( DB_score=1/mean(score1), DB_score2=1/mean(score2))
  return(output)
}


getSilhouette <- function(dist, sample_types){

  # Calculate Davies-Bouldin index
  # Davies-Bouldin compares distances in cluster to the nearest cluster
  #
  # Diagonal is removed using trick with na.rm.
  # Therefore data should not include NAs
  DistTypeMatrix = SeparateDistTypes(sample_types)
  DistMatrix = dist

  dist.size = dim( DistMatrix)
  tmp = rep(0, dist.size[1])
  for( k in 1:dist.size[1] ){
    tmp2 = mean( DistMatrix[k,DistTypeMatrix[k,] > 0] )
    tmp3 = min( DistMatrix[k,DistTypeMatrix[k,] == 0] )
    tmp[k] = (tmp3 - tmp2)/max(tmp3,tmp2)
  }


  Output = list(MeanSilh=mean(tmp), SilhouetteVector=tmp)
  return(Output)

}

##Wrapper

getSilhouettes = function(dist, sample_types){
  if(!is.list(dist)){
    dists = list(dist)
  } else {
    dists = dist
  }

  return(sapply(dists, function(d) getSilhouette(d, sample_types)$MeanSilh))
}


###

### Compute the f-score by scaling the pairwise distances by the median or quartile of the pairwise distances of each datapoint

###INPUT: distance matrix, sample type vector containing group labels of interest, quantile of the pairwise distances to scale (numeric between 0 and 1, default 0.5)

###OUTPUT: scaled f-score of cluster separation (numeric)

##Compute f-score by computing the product of variance ratios of each cluster (biological group and batch effect) separately
##INPUT: Distance matrix, sample type vector containing labels of each biological group of interest
##OUTPUT: the  f-score for cluster separation (Numeric)

getFscore = function(dist, sample_types, scaled = TRUE, quantile = 0.5){

  if(scaled == TRUE){
    dist = apply(dist, 1, function(x) x/(quantile(x, quantile)))
  }

  DistTypeMatrix = SeparateDistTypes(sample_types)


  BetweenClust = sum(dist[DistTypeMatrix == 0])
  WithinClust = sum(dist[DistTypeMatrix > 0])

  return(BetweenClust/WithinClust)


}



###Wrapper for multiple datasets
##INPUT: Distance matrix or a list of distance matrices of datasets of interest, sample type vector or a list of vectors containing biological variables of interest

##OUTPUT: vector of f-scores for each dataset

getFscores = function(dists, sample_types, scaled = TRUE, quantile = 0.5){
  if(!is.list(dists)){
    dists = list(dists)
  }
  f_scores = c()
  for (k in 1:length(dists)){
    if(is.list(sample_types)){
      samples = sample_types[[k]]
    } else {
      samples = sample_types
    }

    if(length(quantile) > 1){
      quant = quantile[k]
    } else {
      quant = quantile
    }

    fscore = getFscore(dists[[k]], samples, scaled = scaled, quantile = quant)

    f_scores = c(f_scores, fscore)
  }
  return(f_scores)
}

###Clustering and chi-square test to measure batch effects. Performs hierarchical clustering and compares the clustering results by chi-square contingemcy table test to the actual clusters of interest (batches or biological groups)
##lower p-value means better cluster separation (worse batch correction)

##INPUT: distance matrix or a list of distance matrices, sample type vector containing the biological or batch labels of interest
##OUTPUT: a list of chi-square test statistics and p-values

getChisq = function(dists, sample_types){
  if(!is.list(dists)){
    dists = list(dists)
  }
  k_clusts = length(unique(sample_types))
  statistics = c()
  pvals = c()
  for(i in 1:length(dists)){
    dist = dists[[i]]
    clust = hclust(as.dist(dist), method = 'complete')
    chisq_test = chisq.test(table(cutree(clust, k = k_clusts), sample_types))
    statistics = c(statistics, chisq_test$statistic)
    pvals = c(pvals, chisq_test$p.value)
  }
  return(list(statistics = statistics, pvals =pvals))
}

###run kBET algorithm (BÃ¼ttner et al, 2019), a kNN-based metric for cluster separation, neighbourhood size is determined as the mean batch size
##INPUT: distance matrix or a list of distance matrices, sample type vector containing the biological or batch labels of interest

##OUTPUT: kBET rejection rate (see BÃ¼ttner et al 2019) and average p-value

getkBET = function(dists, sample_types){
  if(!is.list(dists)){
    dists = list(dists)
  }
  batch = as.numeric(as.factor(sample_types))
  initk =floor(mean(table(batch)))
  avg.pvals = c()
  rej.rates = c()
  for(i in 1:length(dists)){
    data = dists[[i]]
    knn = get.knn(data, k=initk, algorithm = 'cover_tree')
    batch.estimate = kBET(data, batch, plot = FALSE, k0 = initk, knn = knn)
    avg.pvals = c(avg.pvals, batch.estimate$average.pval)
    rej.rates = c(rej.rates, batch.estimate$summary$kBET.observed[1])
  }

  return(list(avg.pvals = avg.pvals, rej.rates = rej.rates))
}


###Compute Davies-Bouldin indices for multiple datasets

##INPUT: distance matrix or a list of distance matrices, sample type vector containing the biological or batch labels of interest
##OUTPUT: the davies-bouldin index of cluster separation (numeric)


DaviesBouldinScores = function(dists, sample_types){
  if(!is.list(dists)){
    dists = list(dists)
  }
  scores = sapply(dists, function(d) DaviesBouldinIndex(d, sample_types)$DB_score)
  return(scores)
}


### Compute Kullback-Leibler divergence of estimates (KDE) of within-cluster density and between-cluster density.
##INPUT: distance matrix or a list of distance matrices, sample type vector containing the biological or batch labels of interest
##OUTPUT: the Kullnack-leibler divergence of within-cluster density and between-cluster density


getKLD = function(dists, sample_types){

  if(!is.list(dists)){
    dists = list(dists)
  }

  disttype = SeparateDistTypes(sample_types)

  kld = sapply(dists, function(d){
    between = d[disttype == 0]
    within = d[disttype > 0]
    kldiv = mean(KL.divergence(density(between)$y, density(within)$y))
    return(kldiv)
  })

  return(kld)
}

### Compute Kullback-Leibler distance of estimates (KDE) of within-cluster density and between-cluster density.
##INPUT: distance matrix, sample type vector containing the biological or batch labels of interest
##OUTPUT: the weighted mean of kl-distances of  within-cluster density and between-cluster density

getKLdist = function(dists, sample_types, k = 5){

  batches = sample_types
  dist = dists
  enum = 1:length(unique(batches))
  disttype = SeparateDistTypes(batches)

  klds = c()
  for(i in enum){
    disti = dist[disttype == i]
    batchInd = which(batches == unique(batches)[i])
    for(j in enum[-i]){
      distj = dist[disttype == j]
      batchInd2 = which(batches == unique(batches)[j])
      distjj = dist[batchInd, batchInd2]
      kldiv1 = mean(KL.dist(density(disti)$y, density(distjj)$y, k = k))
      kldiv2 = mean(KL.dist(density(distj)$y, density(distjj)$y, k = k))
      kld = weighted.mean(c(kldiv1, kldiv2), c(length(batchInd), length(batchInd2)))
      klds = c(klds, kld)
    }
  }
  return(mean(klds))
}

###Wrapper for multiple datasets

getKLdistances = function(dists, sample_types, k = 5){

  if(!is.list(dists)){
    dists = list(dists)
  }

  ret = sapply(dists, function(d) getKLdist(d, sample_types, k))
  return(ret)
}


getKNNprop = function(dist, sample_types){

  k0=floor(mean(table(sample_types)))

  knn = get.knn(dist, k = k0)

  props = c()

  for(b in unique(sample_types)){
    neighbors = knn$nn.index[sample_types == b,]

    batchInd = which(sample_types == b)

    prop = mean(apply(neighbors, 1, function(x) mean(x %in% batchInd)))

    props = c(props, prop)

  }

  return(mean(props))
}

kNN_proportions = function(dist, sample_types){
  if(!is.list(dist)){
    dists = list(dist)
  } else {
    dists = dist
  }

  ret = sapply(dists, function(d) getKNNprop(d, sample_types))
  return(ret)
}


###Compute the average minimum distance between clusters
##INPUT: distance matrix, sample type vector containing the biological or batch labels of interest
##OUTPUT: the mean of minimum distances between the batches or biological clusters

avgMinDist = function(dists, batch){
  dist = dists
  min_dist = c()

  for(t in 1:length(unique(batch))){
    batchI = which(batch == unique(batch)[t])
    for(k in 1:length(unique(batch))){
      if (t < k){
        batchI2 = which(batch == unique(batch)[k])
        min_dist = c(min_dist, min(dist[batchI, batchI2]))
      }
    }

  }
  return(mean(min_dist))

}

###Wrapper for multiple datasets

getMinDists = function(dists, batch){
  if(!is.list(dists)){
    dists = list(dists)
  }

  sapply(dists, function(d) avgMinDist(d, batch))
}


EffectSizeEstimate = function(data, sample_types, dim.reduct.method = "MDS"){
  if(dim.reduct.method == "MDS"){
    MDS = limma::plotMDS(data, plot = FALSE)
    y = (MDS$x + MDS$y)/2
  } else if (dim.reduct.method == "PCA") {
    pca = prcomp(data)
    y = (pca$rotation[,1] + pca$rotation[,2])/2
  }
  model = lm(y~sample_types)
  sum = summary(model)
  return(list(coef = abs(sum$coefficients[2,1]), pval = sum$coefficients[2,4], R2 = sum$r.squared))
}

gPCA.batchdetect <-
  function(x,batch,filt=NULL,nperm=1000,center=FALSE,scaleY=FALSE,seed=NULL){

    # x : n x p matrix of genomic data
    # batch : length n vector indicating batch for each sample
    # filt : number of features to retain after filtering
    # nperm : number of permutations to perform during the permutation test
    # x2.imp : optional centered data matrix with imputed missing values
    # center : logical, is the data centered?
    # scaleY : should Y be scaled based on number of samples in each batch?
    # seed : ability to set the seed for random sampling
    if(!is.null(seed)){set.seed(seed)}

    # Permute batch:
    permute<-matrix(NA,ncol=length(batch),nrow=50000)
    for (j in 1:50000) {permute[j,]<-sample(batch,replace=FALSE)}
    samp<-sample(1:dim(permute)[1],nperm,replace=FALSE)
    permute.samp<-permute[samp,]

    # Center data:
    if(center==FALSE) {
      x2<-scale(x,center=T,scale=F)
    } else {
      x2<-x
    }

    # Test for missing values, impute if missing:
    # (save imputed x2 in working directory so don't have to run it again)
    if(sum(is.na(x))>0) {
      missing<-readline(prompt="Missing values detected. Continue with mean value imputation? (Note this may take a very long time, but it will automatically save in your working dir so you don't have to ever run it again.) [y/n] ")
      if (substr(missing,1,1)=="n") {
        stop("The PC cannot be calculated with missing values.")
      } else {
        x2.imp<-ifelse(is.na(x2),rowMeans(x2,na.rm=TRUE),x2)
        save(x2.imp,"x2.imputed.RData")
      }
    } else {
      x2.imp<-x2
    }

    # Filter data:
    if(is.null(filt)){
      data.imp<-x2.imp
    } else {
      sd<-apply(x2.imp,2,sd)
      rank<-rank(sd)
      keep<-(1:length(sd))[rank %in% (length(rank)-filt+1):length(rank)]
      data.imp<-x2.imp[,keep]
    }
    n<-dim(data.imp)[1]
    p<-dim(data.imp)[2]
    b<-length(unique(batch))
    n ; p ; b

    ## Test for dimensionality:
    if(length(batch)!=n) {stop("Matrices do not conform: length(batch)!=n")}

    # Establish Y matrix indicating batch:
    y<-matrix(nrow=length(batch),ncol=length(unique(batch)))
    for ( j in 1:length(unique(batch)) ){
      y[,j]<-ifelse(batch==j,1,0)
    }
    if (scaleY==FALSE){
      y2<-scale(y,center=T,scale=F) #y2.bat
    } else {
      ys<-matrix(nrow=length(batch),ncol=length(unique(batch)))
      nk<-apply(y,2,sum)
      for ( j in 1:length(unique(batch)) ){
        ys[,j]<-ifelse(batch==j,1/nk[j],0)
      }
      y2<-scale(ys,center=F,scale=F) #y2.bat
    }

    # Unguided SVD:
    svd.x<-svd(data.imp)

    ### Variance of unguided PCs
    PC.u<-data.imp%*%svd.x$v
    var.x<-var(PC.u)
    varPCu1<-diag(var.x)[1]/sum(diag(var.x))
    cumulative.var.u<-numeric()
    for( i in 1:dim(var.x)[1] ){
      cumulative.var.u[i]<-sum(diag(var.x)[1:i])/sum(diag(var.x))
    }

    # Guided SVD:
    svd.bat<-svd(t(y2)%*%data.imp)

    ### Variance of guided PCs
    PC.g<-data.imp%*% svd.bat$v
    var.bat<-var(PC.g)
    varPCg1<-diag(var.bat)[1]/sum(diag(var.bat))
    cumulative.var.g<-numeric()
    for( i in 1:dim(var.bat)[1] ){
      cumulative.var.g[i]<-sum(diag(var.bat)[1:i])/sum(diag(var.bat))
    }

    # Calculate test statistic delta:
    delta<-diag(var.bat)[1]/diag(var.x)[1]


    ##########################################################
    ### Begin loop for random sample of batch permutations ###
    delta.p<-numeric()
    for ( i in 1:nperm ){

      batch.p<-permute.samp[i,]

      y<-ys<-matrix(nrow=length(batch.p),ncol=length(unique(batch.p)))
      for ( j in 1:length(unique(batch.p)) ){
        y[,j]<-ifelse(batch.p==j,1,0)
      }
      if (scaleY==FALSE){
        y2<-scale(y,center=T,scale=F) #y2.bat
      } else {
        nk<-apply(y,2,sum)
        for ( j in 1:length(unique(batch.p)) ){
          ys[,j]<-ifelse(batch.p==j,1/nk[j],0)
        }
        y2<-scale(ys,center=F,scale=F) #y2.bat
      }



      # Perform gPCA
      svd.bat.p<-svd(t(y2)%*%data.imp)

      var.bat.p<-var(data.imp%*% svd.bat.p$v)
      PC.g.p<-diag(var.bat.p)[1]/sum(diag(var.bat.p))

      delta.p[i]<- diag(var.bat.p)[1]/diag(var.x)[1]  #Alternative test statistic


    } # end of permutation loop

    p.val<-sum(delta<delta.p)/length(delta.p)
    p.val
    p.val<-ifelse(p.val==0,"<0.001",round(p.val,3))

    out<-list(delta=delta,p.val=p.val,delta.p=delta.p,batch=batch,filt=filt,n=n,p=p,b=b,
              PCg=PC.g,PCu=PC.u,varPCu1=varPCu1,varPCg1=varPCg1,nperm=nperm,
              cumulative.var.u=cumulative.var.u,
              cumulative.var.g=cumulative.var.g)

  }


#calculates the percentage of variation explained by batch/biological groups
gPCA_percentage = function(data, sample_types, nperm = 1000){
  batch = as.numeric(as.factor(sample_types))
  out = gPCA.batchdetect(t(data), batch = batch, nperm = nperm)
  return(((out$varPCg1-out$varPCu1)/out$varPCg1))
}


kldist <-
  function(xba, batch) {
    
    xba <- scale(xba)
    
    allpairs <- combn(length(levels(batch)), 2)
    
    sum(apply(allpairs, 2, function(y) (sum(batch==y[1]) + sum(batch==y[2]))*kldistTwo(xba[batch==y[1],], xba[batch==y[2],]))/(nrow(xba)*(length(levels(batch))-1)))
    
  }

kldistTwo <-
  function(xb1, xb2) {
    
    ##require("FNN")
    
    distmat1 <- as.matrix(dist(xb1))
    dists1 <- distmat1[lower.tri(distmat1)]
    
    distmat2 <- as.matrix(dist(xb2))
    dists2 <- distmat2[lower.tri(distmat2)]
    
    distmat12 <- apply(xb1, 1, function(y) sqrt(rowSums(sweep(xb2, 2, y, "-")^2)))
    dists12 <- as.vector(distmat12)
    
    (nrow(xb1)*FNN::KL.dist(dists1, dists12, k=5)[5] + nrow(xb2)*
        FNN::KL.dist(dists2, dists12, k=5)[5])/(nrow(xb1) + nrow(xb2))
    
  }






