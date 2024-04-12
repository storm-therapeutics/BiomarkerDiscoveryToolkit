# BiomarkerDiscoveryToolkit

An R package for biomarker discovery.

## Background

A biomarker is defined as "a measurable indicator of some biological state or condition" ([Wikipedia](https://en.wikipedia.org/wiki/Biomarker)).
In medicine, important application areas for biomarkers include disease *diagnosis*, *prognosis* of disease progression, *prediction* of therapy response, and *pharmacodynamic* assessments.
In drug discovery and development, biomarkers that can predict drug efficacy in patients are immensely valuable, and may even be required for patient selection for clinical trials.
Earlier in the drug discovery process, identifying biomarkers for drug response *in vitro* (cell lines) or *in vivo* (animal models) can help to elucidate mechanisms of action or generate hypotheses for patient biomarkers.

## Nomenclature

In a bioinformatic context, a biomarker is a (molecular) *feature* that can be used to predict a biomedical property of interest (*response*) - either by itself or in combination with other biomarkers as part of a *biomarker signature*.
Biomarker discovery thus means evaluating available features for associations with the response; either one by one (*univariate* analysis) or in combinations (*multivariate* analysis).

Examples of *responses* include:
- Genetic dependency from CRISPR/RNAi experiments
- Drug sensitivity from *in vitro* or *in vivo* experiments
- Patient response in a clinical trial

Examples of *features* include:
- Cell line or tumour lineage/tissue of origin
- Omics data:
  - Gene expression
  - Mutations
  - Copy numbers
  - Protein expression

## Current functionality

### Drug dose responses (`curve_fit.R`)
- Dose-response curve fitting (using the [dr4pl](https://cran.r-project.org/package=dr4pl) package)
- AUC calculation
- Visualization of dose-response curves (incl. IC50 and AUC annotations)

### Univariate analyses (`univariate.R`)

#### 1. Correlation-based analysis (`correlations.R`)
- For continuous response and continuous features
- Convenient calculation of correlation coefficients for many features
- Generation of null distributions by permutation
- Significance estimates
- Visualization of distributions and top hits

#### 2. Categorical analysis
- For categorical response and continuous features
- Compare feature values grouped by response using statistical tests
- Visualization of top hits

#### 3. Analysis of mutation data
- For continuous response and binary features
- Compare response values grouped by features using statistical tests
- Visualization of top hits

### Enrichment analyses (`enrichment.R`)
- Gene set enrichment analyses (GSEA; using the [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) package and others)
- Gene Ontology, Reactome and MSigDB gene sets
- Visualization using dot plots (with indication of directionality)

## Planned functionality
- Multivariate prediction models with measures of feature importance (e.g. random forests)
