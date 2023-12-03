## rSuperFeat 
A computational framework to evaluate known cellular states that are oftentimes associated with the critical clinical decisions and drug search. Also, you can train your cell state model to score new datasets. 

## Installation
```{r}
if (! requireNamespace("devtools")) 
      install.packages("devtools")
devtools::install_github('weilin-genomics/rSuperFeat')
```

## Getting started
For now, several cell states are supported: Exhaustion, EMT, CellCycle, Hypoxia, Angiogenesis, Differentiation, Inflammation, Quiescent, MacM1Polarization, MacM2Polarization, apCAFSignature, iCAFSignature, mCAFSignature, vCAFSignature, progenitorCAF in human, progenitorCAF in mouse. 

**1.score cell states on your scRNA-seq data**

rSuperFeat takes binary count matrix as input and score cell states.
```{r}
library(rSuperFeat)
library(Seurat)

myscores <- scoreStates(pbmc_small@assays$RNA@data, state = c("CellCycle","EMT","Exhaustion","Hypoxia"))

# set n.chunks if large amount of cells
myscores <- scoreStates(pbmc_small@assays$RNA@data, state = c("CellCycle","EMT","Exhaustion","Hypoxia"), n.chunks = 20)

# add myscores to seurat object to visualize
pbmc_small = AddMetaData(pbmc_small, metadata = myscores)

# show top features for specific state
printTopWeights(stateName = "EMT", showN = 150, print = T)
```
drugs search using top weighted genes
```{r}
library(signatureSearch)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# LINCS L1000
eh <- ExperimentHub()
lincs_path <- eh[['EH3226']]
geneList = printTopWeights(stateName = "EMT",showN=150)
upset = mapIds(org.Hs.eg.db, as.character(geneList$pos$gene), 'ENTREZID', 'SYMBOL') 
dnset = mapIds(org.Hs.eg.db, as.character(geneList$neg$gene), 'ENTREZID', 'SYMBOL') 

qsig_lincs <- qSig(query = list(upset=upset, downset=dnset), 
                   gess_method="LINCS", refdb=lincs_path)
lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE, workers=1)
result_lincs<-result(lincs) 
result_lincs_filter = result_lincs[result_lincs$trend == 'down',]
result_lincs_filter = result_lincs_filter[order(result_lincs_filter$WTCS),]
drugs <- unique(result_lincs_filter$pert[1:10])
celltypes <- unique(result_lincs_filter$cell[1:10])
g = gess_res_vis(result_lincs, drugs = drugs, col = "WTCS") 
ggplot(g$data,aes(x=pert,y=WTCS,color=cell_type,shape=cell_class)) +
  geom_point() + xlab(NULL) + ylab("Weighted Connectivity Score") +
  theme_bw(base_size = 16) + 
  theme(plot.margin = unit(c(2,1,0,1),'cm'), legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 315, vjust = 1, hjust = 0)) + 
  ggsci::scale_color_ucscgb() + 
  guides(color = guide_legend(title = "cell line", keyheight = unit(0.2,"cm"), ncol = 2, override.aes = aes(size = 2.5)), 
         shape = guide_legend(title = "cell class", override.aes = aes(size = 2.5))) +
  scale_shape_manual(values = c("normal"=16,"tumor"=17))
ggsave('results/dot.lincs.pdf', width = 7, height = 4)

# CMap
eh <- ExperimentHub()
cmap_path <- eh[["EH3223"]]
qsig_cmap <- qSig(query = list(upset=upset, downset=dnset), 
                  gess_method="CMAP", refdb=cmap_path)
cmap <- gess_cmap(qSig=qsig_cmap, chunk_size=5000, workers=1)
result_cmap<-result(cmap) 
result_cmap_filter = result_cmap[result_cmap$trend == 'down',]
result_cmap_filter = result_cmap_filter[order(result_cmap_filter$raw_score),]
drugs <- unique(result_cmap_filter$pert[1:10])
g = gess_res_vis(result_cmap, drugs = drugs, col = "raw_score") 
ggplot(g$data,aes(x=pert,y=raw_score,color=cell_type,shape=cell_class)) +
  geom_point() + xlab(NULL) + ylab("raw score") +
  theme_bw(base_size = 16) + 
  theme(plot.margin = unit(c(2,1,0,1),'cm'), legend.background = element_rect(fill="transparent"),
        axis.text.x = element_text(angle = 315, vjust = 1, hjust = 0)) + 
  ggsci::scale_color_ucscgb() + 
  guides(color = guide_legend(title = "cell line", keyheight = unit(0.2,"cm"), ncol = 1, override.aes = aes(size = 2.5)), 
         shape = guide_legend(title = "cell class", override.aes = aes(size = 2.5))) +
  scale_shape_manual(values = c("normal"=16,"tumor"=17))
ggsave('results/dot.cmap.pdf', width = 6, height = 4)
```

**2.train a new model for your interested cell state.**

(1) prepare training data. This step will fetch cells of two states in group column in which g2 is target cell population and g1 is non-target cell population and data required to train the model will be saved to train_data.csv and train_info.csv. Then use train_new_model/SuperFeat_trainingCode.py to train your model in python console. The weights and bias will be saved to models folder when the traning process finished.
```{r}
getMatrix(pbmc_small, column = "groups", state1 = "g2", state0 = "g1", prefix = "./data/train")
```
(2) use your model to score new scRNA-seq dataset
```{r}
myscores <- scoreStates_selftrain(seuObj@assays$RNA@data, w1_file = "models/w1_binarized.csv", b1_file = "models/b1_binarized.csv")
seuObj = AddMetaData(seuObj, metadata = myscores)
```
(3) display some top weighted genes and further do drug search using above-mentioned code
```{r}
geneList = printTopWeights(w1_file="./inst/extdata/w1.csv",myStateName = "newState")
```

## References
[1] Xie Peng and Gao Mingxuan (2019). SuperCT: a supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles. https://doi.org/10.1093/nar/gkz116. Nucleic Acids Research.

[2] https://github.com/weilin-genomics/SuperCT

[3] Quantitative Learning of Cellular Features From Single-cell Transcriptomics Data Facilitates Effective Drug Repurposing

## Contact
If you have any suggestion, questions and bugs report, feel free to contact weilin.baylor@gmail.com.
