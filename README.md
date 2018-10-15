# A Visually-Interpretable, Dictionary-Based Approach to Radiogenomics

This repo contains code for the experiments in the 2018 paper *A Visually-Interpretable, Dictionary-Based Approach to Imaging-Genomic Modeling, with Low-Grade Glioma as a Case Study.* 
 
## Requirements

All experiments are written in the MATLAB programming language and require [VLFeat](http://www.vlfeat.org/overview/hog.html) for representation extraction and [DICTOL](https://github.com/tiepvupsu/DICTOL) for dictionary learning.

## Execution

To perform the experiments, run either the "mainHOG" header file (for use with the Histogram of Gradients representation), the "mainSIFT" header file (for use with the Scale Invariant Feature Transform representation), or the "mainPIXEL" header file (for use with the trivial representation of the raw pixels). Each header file requires two arguments - the *method* argument should be passed 1 for use with the DLSI learning algorithm, 2 for the COPAR learning algorithm, or 3 for the FDDL learning algorithm, and the *seq* argument should be passed 1 for the Flair modality, 2 for the T1 modality, 3 for the T1C modality, or 4 for the T2 modality. Each header file will return a 1x3 vector *mean_aucs* containing the mean AUCs for prediction of IDH1 Status, Tumor Grade, and Codeletion Status, respectively.
