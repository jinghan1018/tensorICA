#' @title 
#' Feature selection
#' 
#' @aliases feature_selection
#' 
#' @description
#' Feature selection.
#' 
#' @param tica.o
#' The TICA result object, which can be obtained by function *DoTICA*.
#'
#' @param component 
#' The component ID.
#' 
#' @param topN
#' The number of top-ranked features to select from each inferred component.
#' 
#' @param CLkurt
#' The confidence level estimated by top-ranked features' kurtosis.
#' 
#' @return pred.idx
#' The list of predicted feature index.
#' 
#' @return k.p
#' Visualizing the kurtosis of predicted feature.
#' 
#' 
#' @references 
#' Teschendorff AE, Han J, Paul D, Virta J, Nordhausen K.
#' \emph{Tensorial Blind Source Separation for Improved Analysis of Multi-Omic Data.}
#' Genome Biology (2018) 19:76.
#' 
#' @examples
#' data(buccalbloodtensor);
#' Dim.l <- EstDim(buccalbloodtensor$data);
#' dim <- Dim.l$dim;
#' tica.o <- DoTICA(Data = buccalbloodtensor$data, dim = dim, method = "FOBI");
#' feature.n <- feature_selection(tica.o = tica.o, component = 12, topN = 62);
#' buccalbloodtensor$testDMCs[feature.n$pred.idx[[1]]];
#' buccalbloodtensor$testDMCs[feature.n$pred.idx[[2]]];
#' feature.k <- feature_selection(tica.o = tica.o, component = 12, CLkurt = 0.95);
#' buccalbloodtensor$testDMCs[feature.k$pred.idx[[1]]];
#' buccalbloodtensor$testDMCs[feature.k$pred.idx[[2]]];
#' feature.k$k.p[[1]];
#' feature.k$k.p[[2]];
#'  
#' @export
#' 
#' @import dplyr tidyr ggplot2
#' @importFrom stats var
feature_selection <- function(tica.o, component, topN = 0, CLkurt = 0) {
  if(topN!=0 & CLkurt!= 0) {  
    stop("Please input only one way to do feature selection.");
  }
  pred.idx <- list();
  if(topN!=0){
    for(t in seq_len(dim(tica.o$TICA.S)[1])){
      tmp.s <- sort(abs(tica.o$TICA.S[t,component,]), decreasing = TRUE, index.return = TRUE);
      pred.idx[[t]] <- tmp.s$ix[seq_len(topN)];
    }
    return(list(pred.idx = pred.idx));
  }
  
  if(CLkurt!=0){
    k.p <- list();
    kurt <- function(v){
      z <- (v-mean(v,na.rm = TRUE))/sqrt(var(v, na.rm = TRUE));
      n <- length(which(is.na(v) == FALSE));
      k <- n*(n+1)*sum(z^4)/((n-1)*(n-2)*(n-3)) - (3*(n-1)^2)/((n-2)*(n-3));
      return(k);
    }
    k.v <- c();
    for(t in seq_len(dim(tica.o$TICA.S)[1])){
      tmp.s <- sort(abs(tica.o$TICA.S[t,component,]), decreasing = TRUE);
      for (i in seq_len(dim(tica.o$TICA.S)[3]-10)){
        k.v[i] <- kurt(tmp.s);
        tmp.s <- tmp.s[-1];
      }
      topN <- sum(k.v >= k.v[1]*(1-CLkurt));
      
      tmp.s <- sort(abs(tica.o$TICA.S[t,component,]), decreasing = TRUE, index.return = TRUE);
      pred.idx[[t]] <- tmp.s$ix[seq_len(topN)];
      tmp.df <- data.frame(y = k.v, x = seq_len(length(k.v)), color = rep(1,length(k.v)));
      tmp.df$color[seq_len(topN)] <- 2;
      k.p[[t]] <- ggplot(tmp.df, aes(x = x, y = y, color = as.factor(color))) + geom_point(pch = 21, size = 3) +
        scale_color_manual(values = c("black","red"),labels=c("Non-selected","Selected"), name = "") + 
        theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), 
                           legend.text = element_text(size = 14)) +
        labs( x = paste0("Sorted S[",t,",",component,",] ID"), y = "Kurtosis")
    }
    return(list(pred.idx = pred.idx, k.p = k.p));
  }
}
