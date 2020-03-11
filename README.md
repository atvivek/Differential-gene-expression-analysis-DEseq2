# Differential-gene-expression-analysis-DEseq2

###*Installation Requirements*

    Download the most recent versions of R and RStudio for your laptop:

    R
    RStudio

    Install the following packages using the instructions provided below.

    NOTE:  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable but know that it can take awhile.

(a) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them one by one:

install.packages("insert_first_package_name_in_quotations")
install.packages("insert__second_package_name_in_quotations")
& so on ...

Packages to install from CRAN (note that these package names are case sensitive!):

    BiocManager
    RColorBrewer
    pheatmap
    ggrepel
    devtools
    tidyverse

(b) Install the below packages from Bioconductor, using BiocManager::install() function 7 times for the 7 packages:

BiocManager::install("insert_first_package_name_in_quotations")
BiocManager::install("insert_second_package_name_in_quotations") 

Packages to install from Bioconductor (note that these package names are case sensitive!):

    DESeq2
    clusterProfiler
    DOSE
    org.Hs.eg.db
    pathview
    DEGreport
    EnsDb.Hsapiens.v86
    AnnotationHub
    ensembldb

    Finally, please check that all the packages were installed successfully by loading them one at a time using the library() function.

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tidyverse)
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(ensembldb)

    Once all packages have been loaded, run sessionInfo().

sessionInfo()
