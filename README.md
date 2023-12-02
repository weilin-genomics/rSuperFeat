## rSuperFeat 
A computational framework to evaluate known cellular states that are oftentimes associated with the critical clinical decisions and drug search. Also, you can train your cell state model to score new datasets. 

## Installation
```{r}
if (! requireNamespace("devtools")) 
      install.packages("devtools")
devtools::install_github('weilin-genomics/rSuperFeat')
```

## Getting started
For now, four cell states are supported Exhaustion, EMT, CellCycle and Hypoxia. 

**1.score cell states on your scRNA-seq data**

rSuperFeat takes count matrix as input then converted to a binary matrix
```{r}
library(rSuperFeat)
library(Seurat)

myscores <- scoreStates(pbmc_small@assays$RNA@data)

# set n.chunks if large amount of cells
myscores <- scoreStates(pbmc_small@assays$RNA@data, n.chunks = 20)

# add myscores to seurat object to visualize
pbmc_small = AddMetaData(pbmc_small, metadata = myscores)

# print top features for specific state
printTopWeights(stateName = "EMT", posN = 150, negN = 150)
```

**.train new models for your interested cell states.**

(1) prepare training data
```{r}
# this function will fetch Exhaustion and nonExhaustion cells in group column. binarized_data.csv and binarized_info.csv will be saved to data folder. Then use SuperFeat_trainingCode.py to train your model in python console. the weights and bias, model will be saved to models folder when traning finished.
getMatrix(seurat, column = "group", state1 = "Exhaustion", state0 = "nonExhaustion", prefix = "data/binarized")
```
(2) use your model to score new scRNA-seq dataset
```{r}
myscores <- scoreStates_selftrain(seuObj@assays$RNA@data, w1_file = "models/w1_binarized.csv", b1_file = "models/b1_binarized.csv")
seuObj = AddMetaData(seuObj, metadata = myscores)
```

## References
[1] Xie Peng and Gao Mingxuan (2019). SuperCT: a supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles. https://doi.org/10.1093/nar/gkz116. Nucleic Acids Research.

[2] https://github.com/weilin-genomics/SuperCT

[3] Quantitative Learning of Cellular Features From Single-cell Transcriptomics Data Facilitates Effective Drug Repurposing

## Contact
If you have any suggestion, questions and bugs report, feel free to contact weilin.baylor@gmail.com.
