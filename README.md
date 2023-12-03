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

**2.train a new model for your interested cell state.**

(1) prepare training data
```{r}
# this step will fetch cells of two states in group column in which g2 is target cell population and g1 is non-target cell population and data required to train the model will be saved to train_data.csv and train_info.csv. Then use SuperFeat_trainingCode.py to train your model in python console. The weights and bias will be saved to models folder when the traning process finished.
getMatrix(pbmc_small, column = "groups", state1 = "g2", state0 = "g1", prefix = "./data/train")
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
