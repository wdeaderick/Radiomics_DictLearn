# A Visually-Interpretable, Dictionary-Based Approach to Radiogenomics

This repo contains code for the experiments in the 2018 paper *A Visually-Interpretable, Dictionary-Based Approach to Imaging-Genomic Modeling, with Low-Grade Glioma as a Case Study.* 
 
## Requirements

All exeriments are written in the MATLAB programming and require [VLFeat](http://www.vlfeat.org/overview/hog.html) for representation extraction and [DICTOL](https://github.com/tiepvupsu/DICTOL) for dictionary learning.

## Execution

To perform the experiments, run either the "mainHOG" header file (for use with the Histogram of Gradients representation), the "mainSIFT" header file (for use with the Scale Invariant Feature Transform representation), or the "mainPIXEL" header file (for use with the trivial representation of the raw pixels). Each header file requires two arguments - 
