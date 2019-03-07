# JCM

An R package to implement [Joint Clustering and Matching (JCM)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0100334)

### Installation

To install run
```r
#install.packages("devtools")
devtools::github_install("frycast/JCM")
```

You may need to install the geneplotter package from Bioconductor
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("geneplotter", version = "3.8")
```

