## rSuperCT 
the R implemention of a computational framework to evaluate the cellular states that are oftentimes associated with the critical clinical decisions and quickly but thoroughly evaluate the known cellular states/features that are critical to the study of the samples. 
For more information, see:

  SuperFeat: A Framework of Quantitative Feature Learning and Assessment from Single-cell Transcriptomics Data

## Installation
```{r}
if (! requireNamespace("devtools")) 
      install.packages("devtools")
devtools::install_github('weilin-genomics/rSuperFeat')
```

## Getting started
Follow below steps to score five cell states (Exhaution, EMT, Hypoxia, CellCycle, Cytotoxic) for your scRNA-seq data. For now, only homo sapiens are supported. Take *pbmc_small* in Seurat package as example.

1) To score cell states on your scRNA-seq data, library-size normalized data is recommended such as data slot in pbmc_small@assays$RNA
```{r}
library(rSuperCT)
library(Seurat)
# default to score all cell states for your input data
myscores <- scoreStates(pbmc_small@assays$RNA@data)

# set larger n.chunks to split cells into more batches 
myscores <- scoreStates(pbmc_small@assays$RNA@data, n.chunks = 100)

# then myscores can be added to seurat object to visualize and so on
pbmc_small = AddMetaData(pbmc_small, metadata = myscores)
```

2) To score target cell states trained online on your scRNA-seq data,
```{r}
path_to_model_weight_file <- "w1.csv" 
path_to_model_bias_file <- "b1.csv"
# both can be downloaded from online website after self-defined cell states training

myscores <- scoreStates_selftrain(seuObj@assays$RNA@data,
     w1_file = path_to_model_weight_file,
     b1_file = path_to_model_bias_file)
     
# or get through with example model file:
myscores <- scoreStates_selftrain(seuObj@assays$RNA@data,
     w1_file = system.file('extdata','w1.csv',package = 'rSuperFeat'),
     b1_file = system.file('extdata','b1.csv',package = 'rSuperFeat'),
     cell_state_name = "test")
```

## References
[1] Xie Peng and Gao Mingxuan (2019). SuperCT: a supervised-learning framework for enhanced characterization of single-cell transcriptomic profiles. https://doi.org/10.1093/nar/gkz116. Nucleic Acids Research.

[2] https://github.com/weilin-genomics/SuperCT

## Contact
If you have any suggestion, questions and bugs report, feel free to contact weilin.baylor@gmail.com.
