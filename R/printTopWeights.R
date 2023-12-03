#' print top features for specific state
#' @description print top features for specific state
#' @param stateName cell state to print
#' @param w1_file if set, stateName will be ignored
#' @param myStateName state name of w1_file
#' @param species human or mouse for now
#' @param showN show n features
#' @import dplyr
#' @import ggplot2
#' @import cowplot
#' @return a list of positive and bottom features
#' @export
#' @examples
#' \dontrun{
#' printTopWeights(stateName = "EMT",showN=20)
#' printTopWeights(w1_file="w1.csv",myStateName = "newState")
#' }
printTopWeights <- function(stateName = "EMT", w1_file = NULL, myStateName = "newState",species="human", showN = 50, print=T){
  if(is.null(w1_file)){
    wt = switch (species, "human" = wt.hs, "mouse" = wt.mm)
    wti = as.data.frame(wt[,stateName,drop=F])
    name = stateName
  }else{
    wti = read.csv(w1_file, header = FALSE, check.names = F, stringsAsFactors = F)
    rownames(wti) = switch (species, "human" = features.hs, "mouse" = features.mm)
    name = myStateName
  }
  names(wti) = "weight"
  wti$stateName = name
  wti$gene = rownames(wti)
  wti$absW = abs(wti[,1])
  wti = wti[order(wti$absW, decreasing = T),, drop = F]
  wti$direct = ifelse(wti[,1] >= 0, "pos","neg")

  a = wti %>% group_by(direct) %>% top_n(showN,absW) %>% mutate(rank = order(absW, decreasing = T))
  b = split(a, f=a$direct)
  b$pos$gene = factor(b$pos$gene, levels = rev(b$pos$gene))
  b$neg$gene = factor(b$neg$gene, levels = b$neg$gene)
  g1=ggplot(b$pos, aes(direct, gene, fill = absW)) +
    geom_point(aes(size = absW), shape = 21) + theme_classic() +
    scale_fill_distiller(palette = 'Reds', direction = 1) +
    guides(fill = 'none', size='none') + ylab(NULL) + xlab(NULL) +
    theme(axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          plot.subtitle = element_text(size = 8), axis.text.y = element_text(face = 'italic')) +
    scale_x_discrete(expand = c(0,0)) + scale_size_continuous(range = c(0,4)) +
    labs(title = name,subtitle = paste0('', paste0("Positive Weights: ",round(range(b$pos$weight), 2)[1],"~",round(range(b$pos$weight), 2)[2])))
  g2=ggplot(b$neg, aes(direct, gene, fill = absW)) +
    geom_point(aes(size = absW), shape = 21) + theme_classic() +
    scale_fill_distiller(palette = 'Greens', direction = 1) +
    guides(fill = 'none', size='none') + ylab(NULL) + xlab(NULL) +
    scale_y_discrete(position = "right") +
    scale_x_discrete(expand = c(0,0)) + scale_size_continuous(range = c(0,4)) +
    theme(axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          plot.subtitle = element_text(size = 8), axis.text.y = element_text(face = 'italic')) +
    labs(title = "",subtitle = paste0('', paste0("Negative Weights: ",round(range(b$neg$weight), 2)[1],"~",round(range(b$neg$weight), 2)[2])))
  gg = plot_grid(g1,g2,ncol = 2)
  print(gg)
  if(print){
    return(list(pos = as.data.frame(b$pos), neg = as.data.frame(b$neg)))
  }else{
    return(NULL)
  }
}

