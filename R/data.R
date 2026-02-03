#' BSVTK LeGOteknik clone matrix
#'
#' Binary barcode-by-clone matrix for the canonical BSVTK LeGOteknik panel
#' (eBFP2, Sapphire, Venus, Tomato, Katushka).
#'
#' @format A numeric matrix with 5 rows (barcodes) and 31 columns (clones).
#' Row names are barcode names and column names are clone IDs.
#'
#' @source Generated via `createCodeMatrix()` in `data-raw/bsvtk_clone_matrix.R`.
"bsvtk_clone_matrix"

#' CloneCallR palette
#'
#' A list of named color values and ggplot2 scales used for clone plotting.
#'
#' @format A named list with elements `hex`, `colour`, and `fill`. `hex` is a
#' named character vector of hex colors keyed by clone ID; `colour` and `fill`
#' are ggplot2 scale objects built from `hex`.
#'
#' @source Generated in `data-raw/clonepalette.R`.
"clonepalette"
