# OPLSDA_TeA
These scripts accompany the paper **Machine learning for data classification, regression and toxicokinetic modelling of mycotoxin biomarkers of exposure in a human intervention trial**; they describe the OPLS-DA modelling and chemometric filtering employed to screen the untargeted data obtained through the analysis of urinary samples of subjects exposed to the mycotoxin tenuazonic acid (TeA) and control subjects.

## Provided data
The 'data' repository contains the input file (HUIRBJFRHRRE) used during the study to build and cross-validate the OPLS-DA model and determine which are the metabolic features that were most influenced by the exposure to TeA. The dataset consisted of 69 samples collected from 9 volunteers divided into a control group (5 subjects, i.e., those volunteers taking a placebo solution) and a TeA group (4 subjects, i.e., those volunteers taking TeA at the TTC level). Further details about the collected data are availbale in the related paper **Unraveling biomarkers of exposure for tenuazonic acid through urinary metabolomics**. 

## Getting started
Ensure the following R-packages (including their dependencies) are installed and loaded correctly in your environment. The computations were carried out using the R version 4.2.3 (2023-03-15 ucrt) -- "Shortstop Beagle" and the RStudio version RStudio 2023.12.1+402.

The required R packages are:

* structToolbox: extensive set of data (pre-)processing and analysis methods and tools for metabolomics and other omics;
* ropls: package specifically developed for OPLS-DA of metabolic data characterized by multi-collinearity among variables;
* ggplo2: a well-known library for the creation of straightforward graphs and plots;
* gridExtra: a package to be used in combination with ggplot2 to arrange multiple grid-based plots and to prepare tables.
* cowplot: provides various features that help with creating publication-quality figures;
* openxlsx: simplifies the creation of .xlsx files;
* caret: is a framework for building machine learning models.

Additional information can be found on the website **Bioconductor - ropls** (https://bioconductor.org/packages/release/bioc/html/ropls.html) that released the package ropls used for this data elaboration.
  
## Scripts
Under the Scripts directory, script "OPLS-DA.R" is shared. The script containes the R codes require to build the OPLS-DA model and performe the chemometric filtering of the features. 

## Figures
All the plots generated by the script are reported in the Figures directory.

## Usage
To run the code, simply source the scripts in RStudio (after specifying the working directory within the script) – or work following the procedural part of the paper step-by-step.