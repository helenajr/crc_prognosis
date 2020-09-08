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

## Including Plots

You can also embed plots, for example:

![](kmplots_tbxt_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
