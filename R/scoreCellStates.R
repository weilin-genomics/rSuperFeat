#' score the cell states in cells in tumor micro-environment based on the single-cell RNA-seq digital expression profiles
#' @param mat a features x barcodes expression matrix
#' @param normalize whether to normalize the matrix to 0-1 range
#' @param state cell state(s) to evaluate
#' @param n.chunks split these cells into n.chunks batch
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return a scores matrix
#' @export
#' @references
#' SuperFeat: A Framework of Quantitative Feature Learning and Assessment from Single-cell Transcriptomics Data
#' @examples
#' \dontrun{
#' myscores <- scoreStates(pbmc_small@assays$RNA@data)
#' pbmc_small = AddMetaData(pbmc_small, metadata = myscores)
#' }
scoreStates <- function(mat, state = c("Exhaustion",'EMT','CellCycle','Hypoxia',"Cytotoxic"), normalize = TRUE, n.chunks = min(100,ncol(mat))){
  ig <- intersect(features, rownames(mat))
  if(length(ig)/length(features) < 2/3){
    warning('Less than 2/3 required features. Scoring might be biased.')
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

    s <- t(as.matrix(t(wt[,state,drop=F])) %*% mi)
    s <- s + bias[,state,drop=F]

    setTxtProgressBar(pb, i)
    return(s)
  })
  close(pb)
  out = do.call(rbind, probs)
  return(as.data.frame(out))
}

