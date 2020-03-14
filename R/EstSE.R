#' @title 
#' Estimate sensitivity
#' 
#' @aliases EstSE
#' 
#' @description
#' Estimate sensitivity.
#' 
#' @param tica.o
#' The TICA result object, which can be obtained by function *DoTICA*.
#'
#' @param tp
#' An index vector labeling the true positive features associated with a factor of interest, which we know drives joint variation in the data, and which therefore the algorithm should capture.
#'
#' @param topN
#' The number of top-ranked features to select from each inferred component. It must be specified and by default it equals the number of true positives.
#' 
#' 
#' @return se
#' The sensitivity list of all inferred components.
#' 
#' @return pv
#' The p value list of all inferred components.
#' 
#' @return se.p
#' A figure of 'se' list.
#' 
#' @return pv.p
#' A figure of 'pv' list.
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
#' SE <- EstSE(tica.o = tica.o, tp = seq(501,562));
#' SE$se.p;
#' SE$pv.p;
#'  
#' @export
#' 
#' @import dplyr tidyr ggplot2
#' @importFrom stats pbinom
EstSE <- function(tica.o, tp, topN = length(tp)){
  
  se.lv <- list(); pv.lv <- list();
  for(t in seq_len(dim(tica.o$TICA.S)[1])){
    se.lv[[t]] <- vector(); pv.lv[[t]] <- vector();
    for(cp in seq_len(dim(tica.o$TICA.S)[2])){
      
      tmp.s <- sort(abs(tica.o$TICA.S[t,cp,]), decreasing = TRUE, index.return = TRUE);
      ### select the topN features and declare these to be the positives for joint variation
      pred.li <- tmp.s$ix[seq_len(topN)];
      se.lv[[t]][cp] <- length(intersect(unique(pred.li),tp))/length(tp);
      n <- round(se.lv[[t]][cp]*length(tp));
      pv <- pbinom(n,size=topN,prob=length(tp)/dim(tica.o$TICA.S)[3], lower.tail = FALSE);
      pv.lv[[t]][cp] <- pv;
    }
  }
  tmp.df <- data.frame(se.lv);
  colnames(tmp.df) <- seq_len(length(se.lv));
  tmp.df <- gather(tmp.df,value = "y", key = "type");
  tmp.df$x <- rep(seq_len(length(se.lv[[1]])),length(se.lv))
  se.p <- ggplot(tmp.df, aes(x = x, y = y, fill = type)) + geom_bar(position = "dodge", stat = "identity") +
    scale_fill_discrete(name = "", labels = paste0("S[",seq_len(length(se.lv)),",cp,]"))  +
    scale_x_continuous(breaks = seq_len(length(se.lv[[1]]))) + 
    theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), legend.text = element_text(size = 14)) +
    labs(x = "Component", y = "Sensitivity")
  
  tmp.df <- data.frame(pv.lv);
  colnames(tmp.df) <- seq_len(length(pv.lv));
  tmp.df <- gather(tmp.df,value = "y", key = "type");
  tmp.df$x <- rep(seq_len(length(pv.lv[[1]])),length(pv.lv))
  tmp.df$y <- -log10(tmp.df$y)
  pv.p <- ggplot(tmp.df, aes(x = x, y = y, fill = type)) + geom_bar(position = "dodge", stat = "identity") +
    scale_fill_discrete(name = "", labels = paste0("S[",seq_len(length(pv.lv)),",cp,]"))  +
    scale_x_continuous(breaks = seq_len(length(pv.lv[[1]]))) +
    theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), legend.text = element_text(size = 14)) +
    labs( x = "Component", y = "P value (-log10)")
  names(se.lv) <- paste0("S[",seq_len(length(se.lv)),",cp,]");
  names(pv.lv) <- paste0("S[",seq_len(length(se.lv)),",cp,]");
  return(list(se=se.lv, pv=pv.lv, se.p = se.p, pv.p = pv.p));
}
