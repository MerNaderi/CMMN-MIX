# CMMN-MIX
Parsimonious mixtures of contaminated mean-mixture of normal distributions

Please copy the files to the "current working directory" of the R package.

Use your own pathname in PATH().

./Functions

contains 

    (1) Mix-CMMNE.r: EM-based estimating program for the CMMNE-MIX;

    (2) Mix-CMMNW.r: EM-based estimating program for the CMMNW-MIX;

    (3) Mix-CrSN.r: EM-based estimating program for the CrSN-MIX;

    (4) Additional functions.r: function required for random data generation and EM implementation;

  
./Code

   contains

    (1) Fig 1: to reproduce figure 1,

    (2) Fig 2-RMSE.r: to reproduce figure 2 for RMSE;

    (3) Fig 2-RMSE.r: to reproduce figure 2 for Bias;

    (4) Fig 3.r: to reproduce figure 3;

    (5) Fig 4.r: to reproduce figure 3;

    (6) wine test.r: to reproduce table 7;

    (7) wine.r: to reproduce table 8 and fitting models to the wine dataset;

    
Notes:

  1. Some R packages are required to get the results. These include "mvtnorm", "kmed", "trimcluster", 
  "dendextend", "mcclust", "matrixcalc", "gclus", "GGally", "ggpubr", "gridExtra", "ggplot2",
