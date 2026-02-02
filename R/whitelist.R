#' Create a binary barcode-by-clone code matrix
#' 
#' The default LeGOteknik clone matrix is available as the dataset `bsvtk_clone_matrix`.
#'
#' @param barcodes character vector of barcode names (order defines row order).
#' @returns A numeric matrix with `length(barcodes)` rows and `2^length(barcodes) - 1`
#' columns. Rows are barcodes and columns are clone IDs.
#' @export

#'
#' @examples
#' createCodeMatrix(c("eBFP2", "Sapphire", "Venus", "Tomato", "Katushka"))
createCodeMatrix <- function(barcodes) {
  n <- length(barcodes)
  n_clones <- 2^n - 1
  clones <- seq_len(n_clones)

  clone_mat <- sapply(clones, function(i) {
    as.integer(intToBits(i))[seq_len(n)]
  })

  rownames(clone_mat) <- barcodes
  colnames(clone_mat) <- as.character(clones)

  clone_mat
}
