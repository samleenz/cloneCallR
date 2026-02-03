plotPosteriorPHeatmap <- function(sce, annots = c("cc_clone", "scDblFinder.score")) {
    
  mat <- t(counts(altExp(sce, "nbayes_posteriors")))
  ann_df <- colData(sce)
  
  
  
  dbl_col_fun <- circlize::colorRamp2(
    c(0, 1),
    c("white", "firebrick3")
  )
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    clone     = ann_df$cc_clone,
    Doublet   = ann_df$scDblFinder.score,
    col = list(
      clone     = clonepalette$hex[unique(sce$cc_clone)],
      # codeword  = laf::plotting_scales$pal_all[unique(bayes_df$codeword)],
      # dmm_clone = laf::plotting_scales$pal_all[unique(sce$dmm_clone)],
      Doublet = dbl_col_fun
    ),
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
    row_split = ann_df$cc_call,
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