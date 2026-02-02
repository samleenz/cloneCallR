#' Naive Bayes clone assignment from demuxmix posteriors
#'
#' @param dmm demuxmix model object with a posteriorProb slot.
#' @param clone_mat binary barcode-by-clone matrix (barcodes in rows, clones in columns).
#' @param min_conf minimum posterior probability to call a clone "assigned".
#' @param eps small value used to clamp probabilities away from 0 and 1 for stability.
#'
#' @returns A data frame with columns `cell`, `clone`, `prob`, `margin`, and `call`.
#'
#' @examples
#' # TODO: add example usage
calculateNBayesClassification <- function(dmm, clone_mat, min_conf = 0.9, eps = 1e-12){
  # extract the posterior probabilities
  Ppos <- slot(dmm, "posteriorProb")
  
  stopifnot(is.matrix(Ppos))
  H <- nrow(Ppos) # barcodes
  N <- ncol(Ppos) # cells
  
  # shift probabilities away from 0/1 for numerical stability
  P <- pmin(pmax(Ppos, eps), 1 - eps)
  
  K <- ncol(clone_mat)
  
  # log P(clone k | data_i) up to a constant:
  # sum_h [ c_h * log(p_h) + (1-c_h)*log(1-p_h) ]
  logP  <- t(clone_mat) %*% log(P) + t(1 - clone_mat) %*% log(1 - P)  # K x N
  
  # Convert to normalized posteriors across the clone_mat (softmax per column)
  # subtract max for stability
  maxlp <- apply(logP, 2, max) ## max logP per cell
  w <- exp(logP - matrix(maxlp, nrow = K, ncol = N, byrow = TRUE))
  post <- w / matrix(colSums(w), nrow = K, ncol = N, byrow = TRUE)  # K x N
  
  best_k <- max.col(t(post), ties.method = "first")
  best_name <- colnames(clone_mat)[best_k]
  best_prob <- post[cbind(best_k, seq_len(N))]
  
  # optional: margin to 2nd best
  ord2 <- apply(post, 2, function(x) sort(x, decreasing = TRUE)[2])
  second_prob <- ord2
  margin <- best_prob - second_prob
  
  ret <- data.frame(
    cell = colnames(Ppos),
    clone = best_name,
    prob = best_prob,
    margin = margin,
    call = ifelse(best_prob >= min_conf, "assigned", "uncertain"),
    stringsAsFactors = FALSE
  )

  return(ret)
}

#' Run demuxmix and naive Bayes clone classification on a SCE
#'
#' @param sce a SingleCellExperiment with barcode counts in an altExp.
#' @param clone_mat binary barcode-by-clone matrix (barcodes in rows, clones in columns).
#' @param min_conf minimum posterior probability to call a clone "assigned".
#' @param eps small value used to clamp probabilities away from 0 and 1 for stability.
#' @param altExp_assay name of the altExp assay containing the barcode counts.
#' @param barcodes which barcodes in the legoteknik counts to use.
#'
#' @returns The input `sce` with classification columns appended to `colData`.
#'
#' @examples
#' # TODO: add example usage
runNBayesClassification <- function(
    sce,
    clone_mat,
    min_conf = 0.9,
    eps = 1e-12,
    altExp_assay = "lego",
    barcodes = c("eBFP2", "Sapphire", "Venus", "Tomato", "Katushka")
){
  dmm <- runDemuxmix(
    sce = sce,
    altExp_assay = altExp_assay,
    barcodes = barcodes
  )
  ret <- calculateNBayesClassification(
    dmm = dmm,
    clone_mat = clone_mat,
    min_conf = min_conf,
    eps = eps
  )

  ret$clone_filtered <- ifelse(
    ret$call == "assigned",
    ret$clone,
    "uncertain"
  )

  idx <- match(colnames(sce), ret$cell)
  if(any(is.na(idx))){
    stop("Not all SCE cell names found in classification output.")
  }
  ret <- ret[idx, , drop = FALSE]
  ret$cell <- NULL
  colnames(ret) <- paste0("cc_", colnames(ret))



  colData(sce) <- cbind(colData(sce), ret)
  return(sce)
}
