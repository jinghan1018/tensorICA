#' @title 
#' Implement tensorial TICA 
#' 
#' @aliases DoTICA
#' 
#' @description
#' Implement tensorial ICA(tFOBI and tJADE) to decompose the tensor-valued multi-omic data 
#' 
#' @param Data
#' The data object, which contains the multi-way data in two different formats. The data gives the array or data-tensor format. In the former case, and in our specific applications, the first mode defines the tissue, cell or data-type, the second mode defines the samples and the third mode the features (e.g. CpGs or genes). In the latter case, each list entry corresponds to the cell/tissue or data-type and consists of the data-matrix which rows representing features and columns representing samples.  
#'
#' @param dim
#' A vector which contains the number of significant components of each data matrix to search for, and is typically obtained by applying RMT to each separate data/tissue-type matrix (i.e. to the individual entries of data$L above).
#'
#' @param method
#' Choose which TICA method,  tFOBI or tJADE.
#' 
#' @return TICA.S
#' The source array result from TICA, which is of the same size as input data array containing the principal components.
#' 
#' @return TICA.projS
#'  List of projected TICA.S array to each sliced matrix with transforming back to the first mode.
#'
#' @return TICA.W
#' 	List containing all the unmixing matrices of TICA result.
#' 
#' @return TPCA.S
#' The source array result from TPCA, which is of the same size as input data array containing the principal components.
#' 
#' @return TPCA.U 
#' 	List containing the rotation matrices of TPCA result.
#' 	
#' @references 
#' Teschendorff AE, Han J, Paul D, Virta J, Nordhausen K.
#' \emph{Tensorial Blind Source Separation for Improved Analysis of Multi-Omic Data.}
#' Genome Biology (2018) 19:76.
#' 
#' 
#' @references
#' Virta, J., Li, B., Nordhausen, K. and Oja, H.
#' \emph{Independent component analysis for tensorvalued data}
#' Journal of Multivariate Analysis (2017). 
#' 
#' 
#' @references
#' Virta, J., Li, B., Nordhausen, K. and Oja, H.
#' \emph{JADE for Tensor-Valued Observation}
#' Journal of Computational and Graphical Statistics. preprint available on ArXiv 
#' 
#' @examples
#' data(buccalbloodtensor);
#' Dim.l <- EstDim(buccalbloodtensor$data);
#' dim <- Dim.l$dim;
#' tica.o <- DoTICA(Data = buccalbloodtensor$data, dim = dim, method = "FOBI");
#' 
#' @export
#' 
#' @import tensorBSS
DoTICA <- function(Data, dim, method = c("FOBI","JADE")){
  data.a <- Data;
  cdata.a <- tensorCentering(data.a); ### center the data
  tpca.o <- tPCA(cdata.a,d = c(dim(data.a)[1],max(dim)));
  pdata.a <- tensorTransform(cdata.a,t(tpca.o$U[[2]]),2); ### whiten the data
  
  if(method=="FOBI"){    
    tica.o <- tFOBI(pdata.a);
  }else if(method=="JADE"){
    tica.o <- tJADE(pdata.a);
  }else{
    stop("Please input the correct 'method' name.")
  }
  Mix1.m <- solve(tica.o$W[[1]])
  
  projS.a <- tensorTransform(tica.o$S,Mix1.m,1)
  tica.projS.lm <- list();
  for(t in seq_len(dim(data.a)[1])){    
    tica.projS.lm[[t]] <- projS.a[t,,];
  }
  
  return(list(TICA.projS = tica.projS.lm, TPCA.S = tpca.o$S, 
              TPCA.U = tpca.o$U, TICA.S = tica.o$S, TICA.W = tica.o$W));
}
