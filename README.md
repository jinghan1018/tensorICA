---
title: "tensorICA - Tensorial Independent Analysis decomposition methods on multi-omic tensor-valued data - R package"
author: "Andrew E. Teschendorff, Han Jing"
date: "`r Sys.Date()`"
package: "`r pkg_ver('tensorICA')`"
output:
  BiocStyle::html_document
bibliography: tensorICA.bib
vignette: >
  %\VignetteIndexEntry{Tensorial Independent Analysis decomposition methods on multi-omic tensor-valued data - R package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

The **tensorICA** package contains several utility functions for decomposing multi-omic tensor-valued data implementing two tensorial ICA methods, JADE and FOBI, and functions for visualizing the components against phenotype to help with feature selection.


# How to use **tensorICA** package

Here we use a small Illumina HumanMethylation450 BeadChip tensor-valued dataset as an example for tensor decomposition and feature selection. This tensor dataset consist of two matched blood-buccal data matrices defined over the same CpG sites (columns) and the same 152 women samples (rows). Because there are two distinct samples (blood and buccal) per individual, most of the variation is genetic. Hence, to reduce this background genetic variation, the tested CpG sites are obtained by combining 1000 non-smoking associated CpGs with the 62 smoking associated DMCs which are reported from a meta-analysis of several EWAS in blood associated with smoking exposure. We expect tensorICA can capture components corresponding smoking with large weight on the 62 smoking associated CpGs.

The input data *buccalbloodtensor* is stored in the package. 

First, We load **tensorICA** package and test dataset.
```{r load, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
library(tensorICA)
data(buccalbloodtensor)
```

Before the decomposition, we use RMT(Random Matrix Theory), function *EstDim*, to estimate the number of components.
```{r infer, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
require(isva)
Dim.l <- EstDim(buccalbloodtensor)
dim <- Dim.l$dim
dJ <- Dim.l$dJ

```
Here, *dim* is a vector containing components number of each data matrix and *dJ* is the number of significant components of joint variation across data.

Then, we decompose the tensor dataset using tensorial ICA, function *DoTICA*, to infer the sources of data variation. In this case, tICA uses tensorial PCA (TPCA) as a preliminary step. Here we choose tFOBI method in *DoTICA* function for example.
```{r check, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
tica.o <- DoTICA(data = buccalbloodtensor, dim = dim, method = "FOBI")
```

The *tica.o* contain respective rotation matrices *U*, unmixing matrices *W* and source matrices *S* which are of the same size as input matrices containing the principal components.

Then we use these components to help with the feature selection. The phenotype *phenotype.smk.v* we use, which is stored in the package, is the average number of smoking pack per year for each individual sample. We can use function *plot_phenotype* to visualize the components against phenotype. Here we set the parameter *dim.p* top be 2 to show the second mode components.
```{r check2, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
phenotype.p <- plot_phenotype(tica.o = tica.o, phenotype = phenotype.smk.v, dim.p  = 2)
phenotype.p$A_star_matrix.p;
phenotype.p$pv.p;
```
We observe from the p value heatmap that component 13 and component 14 are highly assiciated with smoking. So we use function *plot_component* to look deep into these two components.
```{r check3, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
cp13 <- plot_component(tica.o = tica.o, component = 13, nt.idx = c(1,2), tp);
cp13$weight.p;
cp13$weight_abs.p;
cp13$weight_prime.p;
cp13$weight_prime_abs.p;
cp14 <- plot_component(tica.o = tica.o, component = 14, nt.idx = c(1,2), tp);
cp14$weight.p;
cp14$weight_abs.p;
cp14$weight_prime.p;
cp14$weight_prime_abs.p;
```
To check the efficiency and accuracy of the methods, we can use function *EstSE* to estimate the sensitivity if the true positive set is known.
```{r check4, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
SE <- EstSE(tica.o = tica.o, tp = 1001:1062);
SE$se.p;
SE$pv.p;
```
Finally, we can do feature selection using function *feature_selection*. Here, we can directly input the number of predicted features. Or we can estimate the number of predicted features by inputing the percentage of nongaussianity which predicted feature can explain, so the function *feature_selection* can estimate it by kutosis.
```{r check5, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
feature.n <- feature_selection(tica.o = tica.o, component = 13, topN = 62);
testDMCs[feature.n$pred.idx[[1]]];
testDMCs[feature.n$pred.idx[[2]]];
feature.k <- feature_selection(tica.o = tica.o, component = 13, percentage = 95);
testDMCs[feature.k$pred.idx[[1]]];
testDMCs[feature.k$pred.idx[[2]]];
feature.k$k.p[[1]];
feature.k$k.p[[2]];
```
# Sessioninfo

```{r sessionInfo, echo=FALSE}
sessionInfo()
```

# References




