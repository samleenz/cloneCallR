#' Extract (feature barcode) counts from a SingleCellExperiment
#'
#' @param sce A SingleCellExperiment with feature counts of barcodes.
#' @param altexp The name of the feature type containing the barcode counts
#' @param altexp_name The name for the resulting altExp
#' @param aggregate Whether to aggregate counts across features per barcode tag.
#' @return Modified SingleCellExperiment with barcodes counts stored as an altExp.
#' @export
extract_lego_counts <- function(sce, altexp = "Antibody Capture", altexp_name = "lego", aggregate = TRUE){
  sce <- SingleCellExperiment::splitAltExps(
    sce,
    SummarizedExperiment::rowData(sce)$Type
  )

  # get the lego altExp
  if (!(altexp %in% SingleCellExperiment::altExpNames(sce))) {
    stop(paste("Alternative experiment", altexp, "not found in sce."))
  }

  bar_sce <- SingleCellExperiment::altExp(sce, altexp)


  colnames(bar_sce) <- bar_sce$Barcode

  if(aggregate){
    # sum across the features per tag
    bar_sce <- scater::aggregateAcrossFeatures(
      bar_sce,
      ids = sub(
        "_.*$",
        "",
        rownames(bar_sce)
      ),
      use.assay.type = "counts"
    )
  }
  
  
  SummarizedExperiment::colData(bar_sce) <- SummarizedExperiment::colData(sce)
  SingleCellExperiment::altExp(sce, "Antibody Capture") <- bar_sce
  
  i <- which(SingleCellExperiment::altExpNames(sce) == altexp)
  SingleCellExperiment::altExpNames(sce)[i] <- altexp_name

  sce
  }
