KM plots based on TBXT mRNA levels
================
Helena
08/09/2020

## Before starting the analysis

Load the data and store it in an object ‘samples’. The data I am using
was downloaded from the UCSC Xena project. It combines RNASeq data from
TCGA (cancer tissues) and GTex (normal tissues) via a common
bioinformatics pipeline - so you can directly compare values from both
studies.

``` r
samples <- 
  read.delim("~/Desktop/crc_prognosis/tcga_target_gtex_colorectal.tsv", 
             header=TRUE)
```

Load required packages

``` r
library(survival)
library(survminer)
library(dplyr)
library(forcats)
```

## Do some preliminery analyses

Summarise numbers of samples of each type in your dataset. We are only
interested in ‘Primary tumor’(which comes from the TCGA dataset) and
‘Normal Tissue’ (GTex) for the following analyses.

(I have decided not to include the ‘Solid Tissue Normal’ samples in the
following analyses. These are from the TCGA and are normal samples
adjacent to some of the tumor samples as opposed to from non-cancer
patients and therefore not independent data points. I should do a
boxplot of these vs their primary tumor counterparts though.)

``` r
count(samples, X_sample_type)
```

    ##         X_sample_type   n
    ## 1          Metastatic   1
    ## 2       Normal Tissue 308
    ## 3       Primary Tumor 380
    ## 4     Recurrent Tumor   2
    ## 5 Solid Tissue Normal  51

This is a boxplot summarising TBXT mRNA levels in normal colon tissue
(GTex) vs primary colorectal adenoma.

``` r
normal <- filter(samples, X_sample_type == "Normal Tissue")
cancer <- filter(samples, X_sample_type == "Primary Tumor")
boxplot(normal$T, cancer$T,
        at = c(1,2),
        names = c("Normal", "Cancer"),
        xlab = "Tissue type",
        ylab = "T expression (log2(norm_count+1))")
```

![](kmplots_tbxt_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
