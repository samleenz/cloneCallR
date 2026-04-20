#' Plot a heatmap of naive Bayes posterior probabilities
#'
#' Draw a `ComplexHeatmap` heatmap of clone posterior probabilities stored in
#' the `nbayes_posteriors` alternative experiment of a
#' `SingleCellExperiment`. Rows are annotated with clone assignments and any
#' additional columns requested from `colData(sce)`.
#'
#' @param sce A `SingleCellExperiment` containing an alternative experiment
#'   named `nbayes_posteriors`. The function also expects `colData(sce)` to
#'   contain the requested annotation columns and `cc_call` for row splitting.
#' @param annots Character vector of additional annotation column names from
#'   `colData(sce)` to display alongside the clone assignment.
#' @param clone_col Name of the clone assignment column in `colData(sce)`.
#'
#' @returns Invisibly returns the drawn `ComplexHeatmap` object.
#' @export
#'
#' @examples
#' # TODO: add example usage
plotPosteriorPHeatmap <- function(sce,
                                  annots = c("scDblFinder.score"),
                                  clone_col = "cc_clone") {
    
  mat <- t(counts(altExp(sce, "nbayes_posteriors")))
  ann_df <- colData(sce)
  annot_cols <- unique(c(clone_col, annots))
  missing_cols <- setdiff(c(annot_cols, "cc_call"), colnames(ann_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in colData(sce): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  row_annot_df <- as.data.frame(ann_df[, annot_cols, drop = FALSE])
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    df = row_annot_df,
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    show_legend = FALSE
  )
  
  ## ---- scico colour function (lajolla)
  col_fun <- circlize::colorRamp2(
    seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 256),
    # scico::scico(256, palette = "batlow")
    viridisLite::plasma(256)
  )
  
  ## ---- heatmap
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "pProb",
    left_annotation = row_ha,
    cluster_columns = FALSE,
    row_split = ann_df[["cc_call"]],
    cluster_row_slices = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    col = col_fun,
    heatmap_legend_param = list(
      legend_direction = "horizontal", 
      legend_width = unit(5, "cm")
    )
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "right",
    show_annotation_legend = FALSE
  )
  
  invisible(ht)
}
