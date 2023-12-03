#' the R implemention of a computational framework to evaluate the cellular states
#' @description score the cell states based on your model
#' @param mat a features x barcodes expression matrix
#' @param w1_file weight file downloaded from the online training website
#' @param b1_file bias file downloaded from the online training website
#' @param add_state_name  cell state name
#' @param species support human or mouse for now
#' @param n.chunks split these cells into n.chunks batch
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#' @import Matrix
#' @return a scores matrix
#' @export
#' @references
#' Quantitative Learning of Cellular Features From Single-cell Transcriptomics Data Facilitates Effective Drug Repurposing
#' @examples
#' \dontrun{
#' library(Seurat)
#' myscores <- scoreStates_selftrain(pbmc_small@assays$RNA@counts,
#'     w1_file = system.file('extdata','w1.csv',package = 'rSuperFeat'),
#'     b1_file = system.file('extdata','b1.csv',package = 'rSuperFeat'),
#'     species = "human")
#' pbmc_small = AddMetaData(pbmc_small, metadata = myscores)
#' # or on your model
#' myscores <- scoreStates_selftrain(pbmc_small@assays$RNA@counts,
#'     w1_file = "w1.csv", b1_file = "b1.csv", species = "human")
#' pbmc_small = AddMetaData(pbmc_small, metadata = myscores)
#'
#' }

scoreStates_selftrain <- function(mat,
                                  w1_file,
                                  b1_file,
                                  add_state_name = "mystate",
                                  species = "human",
                                  n.chunks = min(10,ncol(mat))){
  w1 = read.csv(w1_file, header = FALSE, check.names = F, stringsAsFactors = F)
  b1 = read.csv(b1_file, header = FALSE, check.names = F, stringsAsFactors = F)

  features = switch (species, "human" = features.hs, "mouse" = features.mm)

  ig <- intersect(features, rownames(mat))
  if(length(ig)/length(features) < 2/3){
    warning('Less than 2/3 query features exist in database. Scoring might be biased.')
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
    s <- Matrix::t(Matrix::t(w1) %*% mi)
    s <- s + b1$V1
    setTxtProgressBar(pb, i)
    return(s)
  })
  close(pb)
  out = do.call(rbind, probs)
  colnames(out) = add_state_name
  return(as.data.frame(out))
}
