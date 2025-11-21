# Burn-blood-cell-type-timeline

This repository contains the R code and input data used to analyze
**cell type–specific transcriptomic changes** in peripheral blood after burn injury.  
The files support the following analysis modules: differential expression analysis,
MFuzz clustering, transcription factor analysis, and construction of an
early prognostic model for burn patients.

---

### **code/**
R scripts for the different analysis modules:

- **2559celltimeDEG.R**  
  Performs time-series differential expression analysis for multiple cell types
  using a relaxed threshold (|logFC| > 1, adj. p value < 0.05), yielding 2,559
  DEGs across time points.

- **bulk_model.R**  
  Modeling based on bulk (whole-blood) transcriptomic data.

- **bulk_time_DE.R**  
  Time-dependent differential expression analysis in bulk samples.

- **Lasso_model.R**  
  Uses LASSO regression to construct a group-pattern model and to build models
  linking cell fractions with gene expression.

- **Mfuzz_finalplot.R**  
  Performs MFuzz clustering and visualizes temporal gene expression trajectories.

- **TF_analysis.R**  
  Scripts for transcription factor (TF)–related regulatory analysis.

---

### **data/**
Processed datasets used in the analyses above:

- **CIBERSORTx_941_high_resolution_output/**  
  High-resolution CIBERSORTx deconvolution results (cell fractions and related
  outputs) based on 941 genes with |logFC| > 1.5.

- **GSE182616_high_resolution/**  
  High-resolution expression matrices derived from the GSE182616 cohort.

- **cell_time_DEGS_logFC.RData**  
  Log fold-change matrices for time-series DEGs across different cell types.

- **RFmodel.RData**  
  Trained random forest model object used for prediction and feature evaluation.

- **g1.zip, g2.zip, g3.zip**  
  Additional CIBERSORTx high-resolution outputs (cell-type deconvolution results
  and related files) based on 941 genes with |logFC| > 1.

---

This repository contains only the analysis code and related data objects
required to reproduce the main computational steps of the study.
