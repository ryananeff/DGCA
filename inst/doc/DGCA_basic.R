## ---- results="hide", warning=FALSE, message=FALSE-----------------------
library(DGCA, quietly = TRUE)
data(darmanis)
data(design_mat)

## ---- message = FALSE----------------------------------------------------
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"),
  adjust = "none", nPerm = 0, nPairs = 100)
head(ddcor_res)

## ---- fig.width = 7, fig.height = 7, message = FALSE, warning = FALSE----
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"),
  adjust = "perm", nPerm = 5, splitSet = "RTN4")
head(ddcor_res)

## ---- fig.width = 7.5, fig.height = 4, message = FALSE, warning = FALSE----
plotCors(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"), geneA = "RTN4", geneB = "COX6A1")

## ---- fig.width = 6, fig.height = 5, message = FALSE, warning = FALSE----
plotVals(inputMat = darmanis, design = design_mat,
  compare = c("oligodendrocyte", "neuron"), gene = "RTN4") 

