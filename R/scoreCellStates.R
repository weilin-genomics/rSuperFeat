#' the R implemention of a computational framework to evaluate the cellular states
#' @description score the cell states on the single-cell RNA-seq digital expression profiles
#' @param mat a features x barcodes expression matrix
#' @param state cell state(s) to evaluate. For now, Optional states: Exhaustion, EMT, CellCycle, Hypoxia, Angiogenesis, Differentiation, Inflammation, Quiescent,
#'  MacM1Polarization, MacM2Polarization, apCAFSignature, iCAFSignature, mCAFSignature, vCAFSignature, progenitorCAF in human, progenitorCAF in mouse.
#' @param species human or mouse for now
#' @param n.chunks split these cells into n.chunks batch
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @import Matrix
#' @return a scores matrix
#' @export
#' @references
#' Quantitative Learning of Cellular Features From Single-cell Transcriptomics Data Facilitates Effective Drug Repurposing
#' @examples
#' \dontrun{
#' library(Seurat)
#' myscores <- scoreStates(pbmc_small@assays$RNA@counts, state = "CellCycle")
#' pbmc_small = AddMetaData(pbmc_small, metadata = myscores)
#' }
scoreStates <- function(mat,
                        state = c("Exhaustion",'EMT','CellCycle','Hypoxia'),
                        species = "human",
                        n.chunks = min(10,ncol(mat))){
  features = switch (species, "human" = features.hs, "mouse" = features.mm)
  bias = switch (species, "human" = bias.hs, "mouse" = bias.mm)
  wt = switch (species, "human" = wt.hs, "mouse" = wt.mm)

  is = unique(state[state %in% colnames(wt)])
  ns = state[!state %in% colnames(wt)]
  if(length(ns) > 0) message("these states not available: ", paste0(ns, collapse = " "))
  if(length(is) == 0) stop(" No states to score.")
  state = is
  wt = wt[,state,drop=F]

  ig <- intersect(features, rownames(mat))
  if(length(ig)/length(features) < 2/3){
    warning('Less than 2/3 required features. Scoring might be biased.')
  }
  og <- setdiff(features, ig)
  # complete the data matrix
  imat <- mat[ig,,drop=F]
  imat[imat > 0] = 1
  omat <- matrix(0, nrow=length(og), ncol = ncol(mat), dimnames = list(og, colnames(mat)))
  m <- rbind(imat, omat)
  m <- m[features,,drop=F]

  chunk <- function(x, n) split(x, sort(rank(x) %% n))
  cks <- chunk(colnames(m), n.chunks)
  pb <- txtProgressBar(min = 0, max = length(cks), style = 3, char = '=')
  probs <- lapply(1:length(cks), FUN = function(i){
    mi = m[,cks[[i]],drop=F]
    s <- Matrix::t(wt[,state,drop=F]) %*% mi
    s <- s+bias[,state]
    setTxtProgressBar(pb, i)
    return(s)
  })
  close(pb)
  out = as.data.frame(Matrix::t(do.call(cbind, probs)))
  return(out)
}

