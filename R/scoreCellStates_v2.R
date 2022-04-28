#' score the cell states in cells in tumor micro-environment based on the single-cell RNA-seq digital expression profiles
#' @param mat a features x barcodes expression matrix
#' @param w1_file weight file downloaded from the online training website
#' @param b1_file bias file downloaded from the online training website
#' @param cell_state_name  cell state name
#' @param normalize whether to normalize the matrix to 0-1 range
#' @param n.chunks split these cells into n.chunks batch
#' @importFrom utils txtProgressBar setTxtProgressBar read.csv
#' @return a scores matrix
#' @export
#' @references
#' SuperFeat: A Framework of Quantitative Feature Learning and Assessment from Single-cell Transcriptomics Data
#' @examples
#' \dontrun{
#' myscores <- scoreStates_selftrain(pbmc_small@assays$RNA@data,
#'     w1_file = system.file('extdata','w1.csv',package = 'rSuperFeat'),
#'     b1_file = system.file('extdata','b1.csv',package = 'rSuperFeat'))
#' }

scoreStates_selftrain <- function(mat, w1_file = "extdata/w1.csv", b1_file = "extdata/b1.csv", cell_state_name = NULL, normalize = TRUE, n.chunks = min(100,ncol(mat))){
  w1 = read.csv(w1_file, header = FALSE)
  b1 = read.csv(b1_file, header = FALSE)

  ig <- intersect(features, rownames(mat))
  if(length(ig)/length(features) < 2/3){
    warning('Less than 2/3 query features exist in database. Scoring might be biased.')
  }
  og <- setdiff(features, ig)

  # complete the data matrix
  imat <- mat[ig,,drop=F]
  omat <- matrix(0, nrow=length(og), ncol = ncol(mat), dimnames = list(og, colnames(mat)))
  m <- rbind(imat, omat)
  m <- m[features,]

  chunk <- function(x, n) split(x, sort(rank(x) %% n))
  cks <- chunk(colnames(m), n.chunks)
  pb <- txtProgressBar(min = 0, max = 100, style = 3, char = '=')
  probs <- lapply(1:length(cks), FUN = function(i){
    mi = m[,cks[[i]],drop=F]
    if(normalize) mi <- apply(mi, 2, .normalize)

    s <- t(as.matrix(t(w1)) %*% mi)
    s <- s + b1$V1

    setTxtProgressBar(pb, i)
    return(s)
  })
  close(pb)
  out = do.call(rbind, probs)
  colnames(out) = cell_state_name
  return(as.data.frame(out))
}
