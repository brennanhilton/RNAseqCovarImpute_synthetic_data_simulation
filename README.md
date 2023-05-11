RNAseqCovarImpute_synthetic_data_simulation
================
2023-05-11

Code for analyses in the manuscript “RNAseqCovarImpute: a multiple
imputation procedure that outperforms complete case and single
imputation differential expression analysis”

## Installation of RNAseqCovarImpute package

[Detailed instructions
here](https://brennanhilton.github.io/Impute_Covariate_Data_in_RNA_sequencing_Studies.html)

``` r
# Install from github
library(devtools)
install_github("brennanhilton/RNAseqCovarImpute")

# Install from Bioconductor (not available yet on bioconductor)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RNAseqCovarImpute")
```

## simulate data and counts.Rmd does the following:

A third set of simulations utilized completely synthetic covariate data
with modified RNA-seq counts and explored scenarios with varying
starting sample sizes (500, 200, 100, and 50). At each target sample
size, we generated five random normal variables (V1-V5) with means of 0
and standard deviations of 1. We defined one variable, V1, as the main
predictor of interest. The remaining variables, V2-V5, were defined as
confounders of the association between V1 and gene expression. We set
correlations between V1 and V2-V5 at Spearman’s rank correlation
coefficient of 0.1, but no correlations among V2-V5. For 2,000 randomly
sampled genes, we used the seqgendiff package to modify the real RNA seq
data to adjust gene expression associations to match coefficients
randomly drawn from a gamma distribution with an 82.5% null association
rate.

All code in this file cannot be run, as the participant RNA-seq counts
used to generate the synthetic counts cannot be readily shared. However,
the modified counts are provided so that anyone can run the analysis
below.

## simulate missingness and compare methods.Rmd does the following:

We first conducted differential expression analysis using the limma-voom
pipeline on the entire set of observations with their complete covariate
data (hereinafter “full data”). These models estimated the effect of V1
on gene expression while controlling for V2-V5. Genes significantly
associated with the predictor of interest at FDR \<0.05 in these full
data models were defined as true differentially expressed genes (DEGs).
Missingness was then simulated to emulate a common situation in
scientific research where an investigator has complete data for a
predictor of interest (V1), but may have missing data for other
important covariates (V2-V5). Therefore, missingness was only induced in
covariates and not the predictor of interest. We explored scenarios with
various levels of missing data ranging from 5-55% of participants having
at least one missing data point, and under two missingness mechanisms:
missing completely at random (MCAR) and missing at random (MAR). We
simulated ten datasets for each missingness mechanism at each level of
missingness before applying RNAseqCovarImpute, SI, and CC methods and
comparing the results with the full data model. Missingness was
simulated using the ampute function from the mice package. The
limma-voom pipeline was applied for CC and SI as described for the full
data model.

## plot error and bias.Rmd plots the following:

True DEGs from the full data model that were also identified as
significant by a given method were defined as true positives. We report
the percent of true DEGs identified as significant for each method out
of the total number of true DEGs from the full data model. Genes
erroneously identified as significant by a given method that were not
true DEGs from the full data analysis were defined as false positives.
We report the percent of false positives out of the total number of
significant results for each method as the FDR. Bias was calculated as
the difference between coefficients from the full data model and
coefficients produced by each method, expressed as a percent of the full
data model coefficient. We report the average bias per gene per
missingness scenario across all 10 simulations.
