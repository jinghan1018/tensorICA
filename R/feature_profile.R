#' @title 
#' Visualize the feature profile
#' 
#' @aliases feature_profile
#' 
#' @description
#' Visualize the feature profile against phenotype .
#' 
#' 
#' @param Data
#' The data object used in TICA.
#' 
#' @param phenotype.v
#' A vector comtains phenotypes inputed by users to plot against, it should be in the same length with the number of sample.
#'
#' @param phenotype.is.categorical
#' A Boolean logical variable states whether the phenotype is categorical. 
#' 
#' @param feature.index
#' The index of features in data matrix used in TICA. It can be obtained by feature selection.
#' 
#' @param tissue.index
#' The index of tissue in data matrix used in TICA.
#' 
#' @return profile.p
#' All the profile figures.
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
#' feature.k <- feature_selection(tica.o = tica.o, component = 12, CLkurt = 0.95);
#'  (profile.p <- feature_profile(Data = buccalbloodtensor$data, 
#' phenotype.v = buccalbloodtensor$pheno.l$SmokingStatus, phenotype.is.categorical = TRUE,
#' feature.index = feature.k$pred.idx[[1]][seq_len(10)], tissue.index = 1))
#' @export
#' 
#' @import ggplot2 dplyr tidyr 
#' @importFrom cowplot plot_grid
#' @importFrom stats lm na.omit
#' 
feature_profile <- function(Data, phenotype.v, phenotype.is.categorical, feature.index, tissue.index){
  if(phenotype.is.categorical){
    profile.pl <- lapply(seq_len(length(feature.index)),function(j){
      tmp.df <- data.frame(x = as.factor(phenotype.v), y = Data[tissue.index,,feature.index[j]])
      tmp.df <- na.omit(tmp.df)
      profile.p <- ggplot(tmp.df, aes(x = x, y = y, fill = x)) + geom_boxplot(width = 0.5) +
        scale_fill_brewer(type = "div", palette = 5, aesthetics = "fill", guide = FALSE) +
        theme_bw() + theme(text = element_text(size = 8)) + 
        labs(x = "Phenotype", y = paste0("D[",tissue.index,",*,",feature.index[j],"]"), title="")
      return(profile.p)})
    profile.p <- plot_grid(plotlist = profile.pl, ncol=5)
  }
  else{
    profile.pl <- lapply(seq_len(length(feature.index)),function(j){
      tmp.df <- data.frame(x = phenotype.v, y = Data[tissue.index,,feature.index[j]])
      tmp.df <- na.omit(tmp.df)
      profile.p <- ggplot(tmp.df, aes(x = x, y = y)) + 
        geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) + 
        geom_point(pch = 20) +
        theme_bw() + theme(text = element_text(size = 8)) +
        labs(x = "Phenotype", y = paste0("D[",tissue.index,",*,",feature.index[j],"]"), title="")
      return(profile.p)})
    profile.p <- plot_grid(plotlist = profile.pl, ncol=5)
  }
  return(profile.p = profile.p)
}
