#' @title 
#' Visualize the weights of components
#' 
#' @aliases plot_component
#' 
#' @description
#' Visualize the weights of components to help with the feature selection.
#' 
#' @param tica.o
#' The TICA result object, which can be obtained by function *DoTICA*.
#'
#' @param component
#' The component ID to be plotted out. 
#' 
#' @param nt.idx
#' The first mode ID to be plotted out.
#' 
#' @param tp
#' An index vector labeling the true positive features associated with a factor of interest, 
#' which we know drives joint variation in the data, and which therefore the algorithm should capture. 
#' These indices must be in the order of the entries in the 3rd mode of the data$A object, 
#' or alternatively in the same order as the rows of the individual data matrices in data$L.
#' 
#' @return weight.p
#' 
#'
#' @return weight_abs.p
#' 
#' 
#' @return weight_prime.p
#' 
#' 
#' @return weight_prime_abs.p
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
#' cp13 <- plot_component(tica.o = tica.o, component = 12, nt.idx = c(1,2), tp = seq(501,562));
#' cp13$weight.p;
#' cp13$weight_abs.p;
#' cp13$weight_prime.p;
#' cp13$weight_prime_abs.p;
#'  
#' @export
#' 
#' @import ggplot2 dplyr tidyr
#' @importFrom cowplot plot_grid
#' @importFrom stats wilcox.test
plot_component  <- function(tica.o, component, nt.idx, tp = 0){
  if(length(nt.idx) != 2){
    stop("'nt.idx' should be length of 2.")
  }
  if(component > dim(tica.o$TICA.S)[2]){
    stop("Please input correct 'component' ID.")
  }
  if(nt.idx[1] > dim(tica.o$TICA.S)[1] | nt.idx[2] > dim(tica.o$TICA.S)[1]){
    stop("Please input correct 'nt.idx' ID.")
  }
  if(tp != 0){ 
    x <- tica.o$TICA.S[nt.idx[1],component,]
    y <- tica.o$TICA.S[nt.idx[2],component,]
    col.idx <- rep(1,dim(tica.o$TICA.S)[3])
    col.idx[tp] <- 2
    tmp.df <- data.frame(x = x, y = y, color = as.factor(col.idx))
    weight.p <- ggplot(tmp.df,aes(x = x, y = y,  colour = color)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = c("black","red"), name="", labels=c("Controls","Items of Interest"))  +
      theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), 
                         legend.text = element_text(size = 14),legend.position="top") + 
      labs(x = paste0("S[",nt.idx[1],",",component,",*]"), y = paste0("S[",nt.idx[2],",",component,",*]"), title="")
    abs.p <- ggplot(tmp.df,aes(x = abs(x), y = abs(y),  colour = color)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = c("black","red"), name="", labels=c("Controls","Items of Interest"))  +
      theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), 
                         legend.text = element_text(size = 14),legend.position="top") + 
      labs(x = paste0("|S[",nt.idx[1],",",component,",*]|"), y = paste0("|S[",nt.idx[2],",",component,",*]|"), title="")
    box_abs1.p <- ggplot(tmp.df, aes(x = color, y = abs(x), fill = color, colour = color)) + geom_boxplot(width = 0.5) + 
      scale_color_manual(values = c("black","red"), guide = FALSE) + scale_fill_manual(values = c("grey","pink"), guide = FALSE) +
      scale_x_discrete(labels = c("Controls","Items of Interest")) +
      annotate("text", x = 1.5, y = 0.75*max(tmp.df$y), colour = "darkred",
               label = paste0("Wilcox Test\nP = ", format(wilcox.test(abs(filter(tmp.df,color == 2)$x), 
                                                                      abs(filter(tmp.df,color == 1)$x), 
                                                                      alternative = "g")$p.value, digits = 3))) +
      theme_bw() + theme(axis.text=element_text(size=10), axis.title=element_text(size=10), 
                         legend.text = element_text(size = 10)) + 
      labs(y = paste0("|S[",nt.idx[1],",",component,",*]|"), x = "", title="")
    box_abs2.p <- ggplot(tmp.df, aes(x = color, y = abs(y), fill = color, colour = color)) + geom_boxplot(width = 0.5) + 
      scale_color_manual(values = c("black","red"), guide = FALSE) + 
      scale_fill_manual(values = c("grey","pink"), guide = FALSE) +
      scale_x_discrete(labels = c("Controls","Items of Interest")) +
      annotate("text", x = 1.5, y = 0.75*max(tmp.df$y), colour = "darkred",
               label = paste0("Wilcox Test\nP = ", format(wilcox.test(abs(filter(tmp.df,color == 2)$y), 
                                                                      abs(filter(tmp.df,color == 1)$y), 
                                                                      alternative = "g")$p.value, digits = 3))) +
      theme_bw() + theme(axis.text=element_text(size=10), axis.title=element_text(size=10), 
                         legend.text = element_text(size = 10)) + 
      labs(y = paste0("|S[",nt.idx[2],",",component,",*]|"), x = "", title="")
    box_abs_p <- plot_grid(box_abs1.p, box_abs2.p, nrow = 2)
    weight_abs.p <- plot_grid(abs.p, box_abs_p, ncol = 2, rel_widths =c(2,1))
    
    x_p <- tica.o$TICA.projS[[nt.idx[1]]][component,]
    y_p <- tica.o$TICA.projS[[nt.idx[2]]][component,]
    tmp.df <- data.frame(x = x_p, y = y_p, color = as.factor(col.idx))
    weight_prime.p <- ggplot(tmp.df,aes(x = x, y = y,  colour = color)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = c("black","red"), name="", labels=c("Controls","Items of Interest"))  +
      theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14),
                         legend.text = element_text(size = 14),legend.position="top") + 
      labs(x = paste0("S'[",nt.idx[1],",",component,",*]"), y = paste0("S'[",nt.idx[2],",",component,",*]"), title="")
    abs_prime.p <- ggplot(tmp.df,aes(x = abs(x), y = abs(y),  colour = color)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = c("black","red"), name="", labels=c("Controls","Items of Interest"))  +
      theme_bw() + theme(axis.text=element_text(size=14), axis.title=element_text(size=14), 
                         legend.text = element_text(size = 14),legend.position="top") + 
      labs(x = paste0("|S'[",nt.idx[1],",",component,",*]|"), y = paste0("|S'[",nt.idx[2],",",component,",*]|"), title="")
    box_abs1_prime.p <- ggplot(tmp.df, aes(x = color, y = abs(x), fill = color, colour = color)) + geom_boxplot(width = 0.5) + 
      scale_color_manual(values = c("black","red"), guide = FALSE) + 
      scale_fill_manual(values = c("grey","pink"), guide = FALSE) +
      scale_x_discrete(labels = c("Controls","Items of Interest")) +
      annotate("text", x = 1.5, y = 0.75*max(tmp.df$x), colour = "darkred",
               label = paste0("Wilcox Test\nP = ", format(wilcox.test(abs(filter(tmp.df,color == 2)$x), 
                                                                      abs(filter(tmp.df,color == 1)$x), 
                                                                      alternative = "g")$p.value, digits = 3))) +
      theme_bw() + theme(axis.text = element_text(size=10), axis.title = element_text(size=10), 
                         legend.text = element_text(size = 10)) + 
      labs(y = paste0("|S'[",nt.idx[1],",",component,",*]|"), x = "", title="")
    box_abs2_prime.p <- ggplot(tmp.df, aes(x = color, y = abs(y), fill = color, colour = color)) + geom_boxplot(width = 0.5) + 
      scale_color_manual(values = c("black","red"), guide = FALSE) + scale_fill_manual(values = c("grey","pink"), guide = FALSE) +
      scale_x_discrete(labels = c("Controls","Items of Interest")) + 
      annotate("text", x = 1.5, y = 0.75*max(tmp.df$y), colour = "darkred",
               label = paste0("Wilcox Test\nP = ", format(wilcox.test(abs(filter(tmp.df,color == 2)$y), 
                                                                      abs(filter(tmp.df,color == 1)$y), 
                                                                      alternative = "g")$p.value, digits = 3))) +
      theme_bw() + theme(axis.text = element_text(size=10), axis.title = element_text(size=10), 
                         legend.text = element_text(size = 10)) + 
      labs(y = paste0("|S'[",nt.idx[2],",",component,",*]|"), x = "", title="")
    box_abs_prime.p <- plot_grid(box_abs1_prime.p, box_abs2_prime.p, nrow = 2)
    weight_prime_abs.p <- plot_grid(abs_prime.p, box_abs_prime.p, ncol = 2, rel_widths =c(2,1))  
    
    return(list(weight.p = weight.p, weight_abs.p = weight_abs.p, 
                weight_prime.p = weight_prime.p, weight_prime_abs.p = weight_prime_abs.p)) 
  }
  if(tp == 0){
    x <- tica.o$TICA.S[nt.idx[1],component,]
    y <- tica.o$TICA.S[nt.idx[2],component,]
    tmp.df <- data.frame(x = x, y = y)
    weight.p <- ggplot(tmp.df,aes(x = x, y = y)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = "black", name = "", labels = "")  +
      theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                         legend.text = element_text(size = 14)) + 
      labs(x = paste0("S[",nt.idx[1],",",component,",*]"), y = paste0("S[",nt.idx[2],",",component,",*]"), title = "")
    abs.p <- ggplot(tmp.df,aes(x = abs(x), y = abs(y))) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = "black", name = "", labels = "")  +
      theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14), 
                         legend.text = element_text(size = 14),legend.position = "top") + 
      labs(x = paste0("|S[",nt.idx[1],",",component,",*]|"), y = paste0("|S[",nt.idx[2],",",component,",*]|"), title = "")
    
    x_p <- tica.o$TICA.projS[[nt.idx[1]]][component,]
    y_p <- tica.o$TICA.projS[[nt.idx[2]]][component,]
    tmp.df <- data.frame(x = x_p, y = y_p)
    weight_prime.p <- ggplot(tmp.df,aes(x = x, y = y)) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = "black", name = "", labels = "")  +
      theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
                         legend.text = element_text(size = 14)) + 
      labs(x = paste0("S'[",nt.idx[1],",",component,",*]"), y = paste0("S'[",nt.idx[2],",",component,",*]"), title="")
    abs_prime.p <- ggplot(tmp.df,aes(x = abs(x), y = abs(y))) + geom_point(pch = 20, size = 3) + 
      scale_color_manual(values = "black", name = "", labels = "")  +
      theme_bw() + theme(axis.text = element_text(size = 14), axis.title = element_text(size=14), 
                         legend.text = element_text(size = 14),legend.position = "top") + 
      labs(x = paste0("|S'[",nt.idx[1],",",component,",*]|"), y = paste0("|S'[",nt.idx[2],",",component,",*]|"), title="")
    
    return(list(weight.p = weight.p, weight_abs.p = abs.p, weight_prime.p = weight_prime.p, weight_prime_abs.p = abs_prime.p)) 
  }
}
