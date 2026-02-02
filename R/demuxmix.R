#' A wrapper to run `demuxmix` on the barcode counts in the altExp of a SCE object
#'
#' @param sce a sce object with an altExp containig barcode counts 
#' @param altExp_assay name of the altExp assay containing the barcode counts
#' @param barcodes which barcodes in the legoteknik counts to use
#' @param return return the dmm
#'
#' @returns 
#' @export
#'
#' @examples
runDemuxmix <- function(
    sce, 
    altExp_assay = "lego",
    barcodes = c("eBFP2", "Sapphire", "Venus", "Tomato", "Katushka")
){
  require(demuxmix)
  require(SingleCellExperiment)
  # sanity check args
  if(! altExp_assay %in% altExpNames(sce)){
    stop(paste0('altExp assay "', altExp_assay, '" not found in `sce`'))
  }
  
  if(is.null(sce$detected)){
    require(scuttle)
    df <- scuttle::perCellQCMetrics(sce)
    detected <- df$detected
  } else{
    detected <- sce$detected
  }
  
  # run demuxmix
  dmm <- demuxmix::demuxmix(
    counts(altExp(sce, altExp_assay))[barcodes, ], 
    rna = detected
  )
  
  
  return(dmm)
  
}