fFDR: Functional False Discovery Rate
============================

<!-- badges: start -->
[![R build status](https://github.com/StoreyLab/fFDR/workflows/R-CMD-check/badge.svg)](https://github.com/StoreyLab/fFDR/actions)
<!-- badges: end -->

This fFDR package implements the functional false discovery rate (FDR) methodology and estimates FDR from a collection of p-values and observations from an informative variable. The informative variable may affect the likelihood of a true null hypothesis or the power of a statistical test, and the methodology models the likelihood or power as a function of the informative variable. When the informative variable does not affect the likelihood or power, it reduces to the traditional q-value based FDR methodology. The functional FDR methodology induces and employs the functional q-value. At a pre-specified threshold for the functional q-value (which is in general equivalent to a nominal FDR level), a null hypothesis is rejected if its associated functional q-value is less than or equal to the threshold.

An overview of the functional FDR methodology can be found in the preprint ``The Functional False Discovery Rate with Applications to Genomics'' authored by Xiongzhi Chen, David G. Robinson and John D. Storey and available at [bioRxiv](https://doi.org/10.1101/241133).  We recommend reading the package vignette for users of the fFDR package who want a quick start.


Installation and documentation
----------------------------------

To install, open R and type:

```R 
install.packages("devtools")
library("devtools")
devtools::install_github("StoreyLab/fFDR")
```

The vignette can be viewed by typing:

```R
browseVignettes(package = "fFDR")
```

Package overview
--------

### Functions
* `fqvalue`: Estimate functional q-value based on p-value and an informative variable.
* `estimate_fpi0`: Estimate the functional proportion; a tuning parameter lambda will be used for this purpose. 
* `plot`: Plot the estimated functional proportion and/or construct a scatter plot of p-values and quantiles of the informative variable


### Quick start guide

First, the built-in function `simulate_t_tests` generates m=2000 hypotheses, their associated p-values (as p.value) that are induced by t-tests, and (as z) quantile transformed sample sizes of the Normal observations. Note that the power of a t-test is affected by the sample size, which in this setting is an informative variable.

```R
library(fFDR)
ttests_example <- simulate_t_tests(m = 2000)
p.value <- ttests_example$p.value
z <- ttests_example$n
```

Once p-values and observations from an informative are available, `fqvalue` implements the functional FDR methodology, and creates a fqvalue object fq.

```R
fq <- fqvalue(p.value, z)
```

The functional q-value associated with each hypothesis testing can be extract by:

```R
fqvalues <- fq$table$fq.value
```

With the individual functional q-values, at FDR level 0.05, null hypotheses that are rejected and the number of rejections can be found by:

```R
fqvalues <- fq$table$fq.value
which(fqvalues < .05)
sum(fqvalues < .05)
```

The fqvalue object fq can be visualized via 2 panels in a vertical layout: the top panel is a scatter plot of p-values against quantiles of the informative variable; the bottom shows the estimated functional proportion for different values of the tuning parameter lambda, for which the solid curve is the estimated corresponding to the optimal choice of lambda.

```R
p <- plot(fq)
grid::grid.draw(p)
```

For additional details, please see the package vignette.
