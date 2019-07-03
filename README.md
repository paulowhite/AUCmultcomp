# Appendix to the manuscript "Closed tests for multiple comparisons of areas under the ROC curve"

This repository contains data and R code to supplement the manuscript entitled "Closed tests for multiple comparisons of areas under the ROC curve", by Blanche et al.

The R codes to reproduce the results of the manuscript for the example settings A, B, C and D are provided in the following scripts.

 - ExampleA.R
 - ExampleB.R
 - ExampleC.R
 - ExampleD.R

Each of the above script calls R functions defined in the following R scripts

 - PlotClosedTest.R
 - SingleStep.R

The first script contains the main function to produce "pedagogical" plots such as those of Figures 2, 4 and 5 of the manuscript. The second contains a function which is called by the main function.

The plot of the four ROC curves at time t=5 years is produced in the script ExampleA.R.

The data (paquid.csv) are similar to some others which are available from the timeROC and riskRegression packages, which are avaialble on the CRAN. See the manuscript for details.