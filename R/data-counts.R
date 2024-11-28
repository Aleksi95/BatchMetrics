#' An example of raw RNA-Seq count data
#' 
#' This datasets contains 58884 rows (genes) and 107 samples. Both gene and sample names have been anonymized. 
#' 
#' The data contains the raw RNA-Seq count data matrix, and the DESeq-normalized, variance stabilized and batch corrected (ComBat) dataset.
#' 
#' This file also contains the batch variable (factor with 4 levels) and the sample group variable (5 sample types, one being the control samples)
#'
#' @author Aleksi Laiho \email{aleksi.laiho@@helsinki.fi}

#'Count data
#'@format matrix
"counts"

#'Normalized batch corrected data
#'@format matrix
"data"

#'batch variable (factor)
#'@format factor
"batch"

#'Sample group variable (factor)
#'#'@format factor
"group"
