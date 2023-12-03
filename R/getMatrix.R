#' prepare the matrix and cell labels for training
#' @description prepare the matrix and cell labels for training
#' @param obj SeuratObject
#' @param column which column of cell labels
#' @param state1 target cell type
#' @param state0 non-target cell state
#' @param prefix prefix of file to save data
#' @param species human or mouse
#' @importFrom utils write.csv
#' @return
#' save matrix and cell information as prefix_data.csv and prefix_info.csv
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' getMatrix(pbmc_small,column = 'groups', state1 = "g2", state0 = "g1")
#' }
getMatrix <- function(obj, column, state1, state0, prefix = "./train", species = "human"){
  meta = obj@meta.data[,column,drop=F]
  names(meta) = 'group'
  meta$Cell = rownames(meta)
  meta = meta[which(meta$group %in% c(state1,state0)),,drop = F]
  meta[which(meta$group == state1),"group"] = "state1"
  meta[which(meta$group == state0),"group"] = "state0"
  if(nrow(meta) == 0) stop("No available cells.")
  meta = meta[,c("Cell","group")]
  cnt = table(meta$group)
  message("the cell number of state1 is ", cnt[["state1"]])
  message("the cell number of state0 is ", cnt[["state0"]])

  mat = obj@assays$RNA@counts[,rownames(meta),drop=F]
  features = switch (species, "human" = features.hs, "mouse" = features.mm)

  ig <- intersect(features, rownames(mat))
  if(length(ig)/length(features) < 2/3){
    warning('Less than 2/3 required features. Scoring might be biased.')
  }
  og <- setdiff(features, ig)

  imat <- mat[ig,,drop=F]
  # binarize the matrix
  imat[imat > 0] = 1

  # complete the data matrix
  omat <- matrix(0, nrow=length(og), ncol = ncol(mat), dimnames = list(og, colnames(mat)))
  m <- rbind(imat, omat)
  m <- m[features,,drop=F]

  # save data
  dataN = paste0(prefix,'_data.csv')
  infoN = paste0(prefix,'_info.csv')

  write.csv(m, dataN, quote = F)
  write.csv(meta, infoN, row.names = F, quote = F)

  message("save data to ", dataN, "  ",infoN)
}

