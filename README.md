# cloneCallR

Tools for clone calling using `demuxmix` posteriors and a naive Bayes classifier.

## Installation

```r
# install.packages("remotes")
remotes::install_github("samleenz/cloneCallR")
```

## Usage

```r
library(cloneCallR)

# run demuxmix + naive Bayes classification
# sce: SingleCellExperiment with barcode counts in an altExp
# clone_mat: binary barcode-by-clone matrix
sce <- runNBayesClassification(
  sce = sce,
  clone_mat = cloneCallR::bsvtk_clone_matrix,
  altExp_assay = "lego"
)

# classification columns are added to colData with the "cc_" prefix
```

Please cite demuxmix if you use this package

```r
citation("demuxmix")
```

