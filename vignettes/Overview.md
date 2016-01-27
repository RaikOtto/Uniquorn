---
title: "Overview Vignette for the Uniquorn Package"
author: "Raik Otto"
date: "January 14, 2016"
output: 
  rmarkdown::html_vignette: 
    toc: true

vignette: >
  %\VignetteIndexEntry{Overview Vignette for Uniquorn Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## Introduction

Yet to comne

## Data Origin


```r
list.files(system.file("extdata", package = "Uniquorn"))
```

```
## character(0)
```

## Adding Cancer cell line training sets

## Querying for contained Cancer cell lines and mutations

## Working with the full CCLE dataset


## Future directions

Yet to come

## Session Info

```r
   sessionInfo() 
```

```
## R version 3.2.1 (2015-06-18)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.2 (unknown)
## 
## locale:
## [1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] reshape2_1.4.1        CancerCellLines_0.4.9 shiny_0.12.2         
##  [4] scales_0.3.0          tidyr_0.3.1           readr_0.2.2          
##  [7] readxl_0.1.0          dplyr_0.4.3           Uniquorn_0.99.0      
## [10] RSQLite_1.0.0         DBI_0.3.1             devtools_1.9.1       
## [13] stringr_1.0.0         caret_6.0-64          ggplot2_2.0.0        
## [16] lattice_0.20-33       e1071_1.6-7           ROCR_1.0-7           
## [19] gplots_2.17.0        
## 
## loaded via a namespace (and not attached):
##  [1] gtools_3.5.0       splines_3.2.1      colorspace_1.2-6  
##  [4] htmltools_0.3      stats4_3.2.1       mgcv_1.8-10       
##  [7] nloptr_1.0.4       foreach_1.4.3      plyr_1.8.3        
## [10] MatrixModels_0.4-1 munsell_0.4.2      gtable_0.1.2      
## [13] caTools_1.17.1     codetools_0.2-14   memoise_0.2.1     
## [16] evaluate_0.8       knitr_1.11         WriteXLS_4.0.0    
## [19] SparseM_1.7        httpuv_1.3.3       quantreg_5.19     
## [22] pbkrtest_0.4-4     parallel_3.2.1     curl_0.9.4        
## [25] class_7.3-14       Rcpp_0.12.2        xtable_1.8-0      
## [28] KernSmooth_2.23-15 formatR_1.2.1      gdata_2.17.0      
## [31] mime_0.4           lme4_1.1-10        digest_0.6.8      
## [34] stringi_1.0-1      grid_3.2.1         tools_3.2.1       
## [37] bitops_1.0-6       magrittr_1.5       car_2.1-1         
## [40] MASS_7.3-45        Matrix_1.2-3       assertthat_0.1    
## [43] minqa_1.2.4        roxygen2_5.0.1     httr_1.0.0        
## [46] iterators_1.0.8    R6_2.1.1           nnet_7.3-11       
## [49] nlme_3.1-122
```
