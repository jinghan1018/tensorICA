#' @title 
#' Estimate dim and dJ by RMT (Random Matrix Theory)
#' 
#' @aliases EstDim
#'  
#' @description 
#' Auxiliary function to estimate dimensionality parameters
#' 
#' @param Data
#' The array or data-tensor format of input data. In this case, and in our specific applications, the first mode defines the tissue, cell or data-type, the second mode defines the samples and the third mode the features (e.g. CpGs or genes). 
#' 
#' @return dim
#' a vector which contains the number of significant components of each data matrix to search for. 
#' 
#' @return dJ
#' the number of significant components of joint variation across data/tissue types. 
#' 
#' @references 
#' Teschendorff AE, Han J, Paul D, Virta J, Nordhausen K.
#' \emph{Tensorial Blind Source Separation for Improved Analysis of Multi-Omic Data.}
#' Genome Biology (2018) 19:76.
#' 
#' 
#' @references
#' Plerou, V., Gopikrishnan, P., Rosenow, B., Amaral, L.A., Guhr, T., Stanley, H.E. 
#' #' \emph{ Random matrix approach tocross correlations in financial data.} 
#' Phys Rev E Stat Nonlin Soft Matter Phys (2002) 65(6), 066126 
#' 
#' Teschendorff, A.E., Zhuang, J., Widschwendter, M.
#' \emph{Independent surrogate variable analysis to deconvolve confounding factors in large-scale microarray profiling studies.}
#' Bioinformatics (2011) 27(11), 1496-1505
#' 
#' @examples 
#' data(buccalbloodtensor);
#' Dim.l <- EstDim(buccalbloodtensor$data);
#' dim <- Dim.l$dim;
#' dJ <- Dim.l$dJ;
#' 
#' @export
#'    
#'     
#' @import isva
EstDim <- function(Data){
  nt <- dim(Data)[1];
  data.l <- list();
  for(i in seq_len(nt)){
    data.l[[i]] <- t(Data[i,,]);
  }
  ### estimate joint variation
  data.m <- data.l[[1]];
  for(i in seq_len(nt-1)){
    data.m <- rbind(data.m,data.l[[i+1]]);
  }
  sd.v <- sqrt(apply(data.m,1,var));
  d <- EstDimRMT((data.m - rowMeans(data.m))/sd.v, plot=FALSE)$dim;
  dim <- vector();
  for(j in seq_len(nt)){
    dim[j] <- EstDimRMT(data.l[[j]] - rowMeans(data.l[[j]]), plot = FALSE)$dim;
  }
  return(list(dJ=d, dim=dim));
}



