#' @title 
#' Correlate inferred components against phenotype
#' 
#' @aliases cor_phenotype
#' 
#' @description
#' Correlate inferred components against phenotype.
#' 
#' 
#' @param tica.o
#' the TICA result object, which can be obtained by function *DoTICA*.
#'
#' @param phenotype
#' A list comtains phenotype vectors inputed by users to correlate against, each should be in the same length with the number of sample.
#'
#' @param phenotype.is.categorical
#' A Boolean logical vector states whether the phenotype is categorical. It is of the same length of phenotype list.
#'
#' @param component
#' Choose the component significantly correlated to specific phenotype. It is an optional parameter for boxplot or scatter plot of chosen component against phenotype.
#' 
#' @return A_star_matrix.p
#' A heatmap of the second mode mixing matrix.
#'
#' @return pv.p
#' A heatmap of p value between phenotypes and the second mode mixing matrix components using linear regression model.
#' 
#' @return compheno.pl
#' A list of boxplots or scatter plots of chosen component against all phenotype.
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
#' phenotype.p <- cor_phenotype(tica.o = tica.o, phenotype = buccalbloodtensor$pheno.l,
#'  phenotype.is.categorical = buccalbloodtensor$pheno.i);
#' phenotype.p$A_star_matrix.p;
#' phenotype.p$pv.p;
#'  
#' @export
#' 
#' @import ggplot2 dplyr tidyr 
#' @importFrom cowplot theme_cowplot
#' @importFrom stats lm na.omit
cor_phenotype <- function(tica.o, phenotype, phenotype.is.categorical, component = 0){
  if(length(phenotype.is.categorical)!=length(phenotype)){
    stop("The 'phenotype' list is not matched to 'phenotype.is.categorical'.");
  }
  
  A_star.m <- t(solve(tica.o$TICA.W[[2]]))%*%t(tica.o$TPCA.U[[2]]);
  
  tmp.df <- gather(data.frame(A_star.m,id=seq_len(nrow(A_star.m))), key = "x", value = "weight", seq_len(ncol(A_star.m)))
  tmp.df$x <- as.numeric(as.factor(tmp.df$x))
  
  colnames(tmp.df) <- c("Components","Sample_ID","Weight")
  A_star_matrix.p <- ggplot(tmp.df , aes( x =  Sample_ID, y = Components , fill = Weight)) + 
    geom_tile() + theme_bw() + labs(title = "Mixing Matrix of Sample Mode") +
    scale_fill_gradient2(low = "blue", high = "red",na.value = "grey") +
    scale_y_continuous(limits = c(0,nrow(A_star.m)+1),breaks = seq_len(nrow(A_star.m)))
  
  pv.df <- data.frame()
  for(i in seq_len(dim(tica.o$TICA.S)[2])){
    pv.v <- unlist(lapply(seq_len(length(phenotype)),function(j){
      if(phenotype.is.categorical[j]){
        lm.o <- lm(A_star.m[i,] ~ factor(phenotype[[j]]))
      }else{
        lm.o <- lm(A_star.m[i,] ~ phenotype[[j]])
      }
      lm.s <- summary(lm.o)
      return(lm.s$coefficients[2,4])
    }))
    tmp.df <- data.frame(p = pv.v, components = i, phenotype = names(phenotype))
    pv.df <- rbind(pv.df,tmp.df)
  }
  pv.df$sig <- factor(as.numeric(cut(pv.df$p, breaks = c(0,1e-10,1e-5,0.001,0.05,1))), levels = seq_len(5))
  
  sig_labels.v <- c("1e-10","1e-5","1e-3")
  sig_labels.v[[4]] <- "0.05"
  sig_labels.v[[5]] <- "N.S"
  pv.p <- ggplot(pv.df, aes(x = components, y = phenotype, fill = sig)) + 
    geom_tile() + scale_fill_manual(values = c("darkred","red","tomato","pink","white"),
                                    limits=c("1","2","3","4","5"),
                                    name="P value", labels=sig_labels.v) +
    labs( y = "Phenotype", title="", x = "Component") +
    theme_cowplot(font_size = 10) +
    scale_x_continuous(name="Component", breaks=seq_len(dim(tica.o$TICA.S)[2]), expand = c(0,0)) + 
    theme(axis.line = element_blank(),
          legend.key=element_rect(color="black", linetype = 1, size= 0.5), 
          legend.key.size = unit(10, "pt"))
  if(component != 0){
    compheno.pl <- lapply(seq_len(length(phenotype)),function(j){
      if(phenotype.is.categorical[j]){
        tmp.df <- data.frame(x = as.factor(phenotype[[j]]), y = A_star.m[component,])
        tmp.df <- na.omit(tmp.df)
        sample_size <- tmp.df %>% group_by(x) %>% summarise(n())
        compheno.p <- ggplot(tmp.df, aes(x = x, y = y, fill = x)) + geom_boxplot(width = 0.5) +
          scale_fill_brewer(type = "div", palette = 5, aesthetics = "fill", name = "Sample Size", labels = paste0("n = ", sample_size$`n()`)) +
          theme_bw() + theme(text = element_text(size = 14)) + 
          labs(x = names(phenotype)[j], y = "Weight", title=paste0("Component ", component))
      }else{
        tmp.df <- data.frame(x = phenotype[[j]], y = A_star.m[component,])
        tmp.df <- na.omit(tmp.df)
        compheno.p <- ggplot(tmp.df, aes(x = x, y = y)) + geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) + geom_point(pch = 20, size = 3) +
          theme_bw() + theme(text = element_text(size = 14)) +
          labs(x = names(phenotype)[j], y = "Weight", title=paste0("Component ",component))
      }
      return(compheno.p)
    })
    names(compheno.pl) <- names(phenotype)
    return(list(A_star_matrix.p = A_star_matrix.p, pv.p = pv.p, compheno.pl = compheno.pl))
  }
  return(list(A_star_matrix.p = A_star_matrix.p, pv.p = pv.p))
}
