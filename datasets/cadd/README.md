# Combined Annotation-Dependent Depletion (CADD)

CADD is a method for integrating genomic variant (SNV and small indel) annotations into a single metric (C score) for each position in the human genome. 

## Features



## Imputation and data preprocessing for modeling


## Data sources:
  * [Publication](http://www.ncbi.nlm.nih.gov/pubmed/24487276)
  * [Supplemental file](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3992975/bin/NIHMS555958-supplement-1.pdf)
  * [Website](http://cadd.gs.washington.edu/) - main website with information, files and online submission for scoring (and annotating) up to 100,000 variants.
  * [Data source 1](http://krishna.gs.washington.edu/members/mkircher/download/CADD/) - files server for the CADD website.
  * [Data source 2](http://krishna.gs.washington.edu/members/mkircher/download/CADD/)

### Data source 2

A file server specifically linked to the first author. This source contains the updated (CADD v1.3) files, fully annotated testing and training datasets as well as the simulated and observed variant sites in vcf format. The testing and training files available are:
  * test.csv.gz - File contains the "imputed" data matrix with XX rows and XX columns. The documentation
  * train.csv.gz - File contains the "imputed" data matrix with XX rows and XX columns. 
  * training_data.tsv.gz - File contains the "non-imputed" data matrix with all the training sites and v1.3 featureset. There are 35043061 rows and XX columns.
