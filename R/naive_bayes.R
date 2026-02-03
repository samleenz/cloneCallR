#' calculate the clone possterior pobabilities using a naive Bayes approach
#'
#' Given a matrix of posterior probabilities for barcodes (rows) across cells (columns),
#' and a codebook matrix of barcodes (rows) by clones (columns), calculate the posterior
#' probabilities for each clone in each cell using a naive Bayes approach.
#' 
#' @param ppos posterior probability matrix, e.g., `dmm@posteriorProb`.
#' @param clone_mat binary barcode-by-clone matrix (barcodes in rows, clones in columns).
#' @param eps small value used to clamp probabilities away from 0 and 1 for stability.
#'
#' @returns a matrix of normalised posterior probabilities (clones x cells).
#' @export
#' 
#' @examples
#' # TODO: add example usage
calculateClonePosteriors <- function(ppos, clone_mat, eps = 1e-12){
  stopifnot(is.matrix(ppos))
  stopifnot(all(rownames(ppos) == rownames(clone_mat)))

  H <- nrow(ppos) # barcodes
  N <- ncol(ppos) # cells
  
  # shift probabilities away from 0/1 for numerical stability
  P <- pmin(pmax(ppos, eps), 1 - eps)
  
  K <- ncol(clone_mat)
  
  # log P(clone k | data_i) up to a constant:
  # sum_h [ c_h * log(p_h) + (1-c_h)*log(1-p_h) ]
  logP  <- t(clone_mat) %*% log(P) + t(1 - clone_mat) %*% log(1 - P)  # K x N
  
  # Convert to normalized posteriors across the clone_mat (softmax per column)
  # subtract max for stability
  maxlp <- apply(logP, 2, max) ## max logP per cell
  w <- exp(logP - matrix(maxlp, nrow = K, ncol = N, byrow = TRUE))
  post <- w / matrix(colSums(w), nrow = K, ncol = N, byrow = TRUE)  # K x N

  return(post)
}

#' Naive Bayes clone assignment from demuxmix posteriors
#'
#' @param post matrix of posterior probabilties (clones x cells), probably from `calculateClonePosteriors()`.
#' @param min_conf minimum posterior probability to call a clone "assigned".
#'
#' @returns A data frame with columns `cell`, `clone`, `prob`, `margin`, and `call`.
#' @export
#' 
#' @examples
#' # TODO: add example usage
calculateNBayesClassification <- function(post, min_conf = 0.9){
  N <- ncol(post)

  best_k <- max.col(t(post), ties.method = "first")
  best_clone <- rownames(post)[best_k]
  best_prob <- post[cbind(best_k, seq_len(N))]
  
  # margin between best and second best probabilities
  ord2 <- apply(post, 2, function(x) sort(x, decreasing = TRUE)[2])
  second_prob <- ord2
  margin <- best_prob - second_prob
  
  ret <- data.frame(
    cell = colnames(post),
    clone = best_clone,
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
#' @returns The input `sce` with classification columns appended to `colData` and the posterior probabilities in an altExp.
#' @export
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

  post <- calculateClonePosteriors(
    ppos = dmm@posteriorProb,
    clone_mat = clone_mat,
    eps = eps
  )

  ret <- calculateNBayesClassification(
    post = post,
    min_conf = min_conf
  )

  ret$clone_filtered <- ifelse(
    ret$call == "assigned",
    ret$clone,
    "0"
  )

  idx <- match(colnames(sce), ret$cell)
  if(any(is.na(idx))){
    stop("Not all SCE cell names found in classification output.")
  }
  ret <- ret[idx, , drop = FALSE]
  ret$cell <- NULL
  colnames(ret) <- paste0("cc_", colnames(ret))

  altExp(sce, "nbayes_posteriors") <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = post
    )
  )

  assay(altExp(sce, altExp_assay), "demuxmix_posterior", withDimnames = FALSE) <- matrix(NA, nrow = nrow(altExp(sce, altExp_assay)), ncol = ncol(sce))
  assay(altExp(sce, altExp_assay), "demuxmix_posterior")[barcodes, ] <- dmm@posteriorProb[barcodes, ]

  SummarizedExperiment::colData(sce) <- cbind(colData(sce), ret)
  return(sce)
}
