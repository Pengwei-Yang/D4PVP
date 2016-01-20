# PVP

Documentation, scripts and data for pathogenic variant prediction modeling.

## CADD

Data and documentation sources:
  * [A general framework for estimating the relative pathogenicity of human genetic variants.](http://www.ncbi.nlm.nih.gov/pubmed/24487276)
  * Publication [website](cadd.gs.washington.edu)
  * Data storage website [here](http://krishna.gs.washington.edu/download/CADD/) and [here](http://krishna.gs.washington.edu/members/mkircher/download/CADD/).

## DANN

DANN aimed to improve upon the CADD pathogenic variant prediction using a "non-linear" modeling approach using neural nets. The same dataset used by CADD was used by DANN for training, validation and testing. They also analyze real data from ClinVar and ESP, similar to the CADD paper.

Data and documentation sources:
  * [DANN: a deep learning approach for annotating the pathogenicity of genetic variants](http://www.ncbi.nlm.nih.gov/pubmed/25338716)
  * Publication [website](https://cbcl.ics.uci.edu/public_data/DANN/).
  * Github [website](https://github.com/uci-cbcl/DeepCADD)

### Publication website datasets

After reading the DANN paper and analyzing the publication website it was still unclear as to whether the data used for the reported analysis for available for use in either of the identified data sources. I will try to lay out all the convolutions and my conclusions to moving forward. First, the DANN publication cites the use of 10,000 pathogenic and benign variants each for their ROC curve. The publication website contains multiple datasets from the publication:
  1. testing.X.npz, testing.svmlight.gz, testing.y.npy
  2. training.X.npz, training.svmlight.gz, training.y.npy
  3. validation.X.npz, validation.svmlight.gz, validation.y.npy
  4. ClinVar_ESP.X.npz, ClinVar_ESP.svmlight.gz, ClinVar_ESP.y.npy

These four datasets are stored in numpy / scipy formats (.npz and .npy) and a custom svmlight gzip format. Hopefully the naming of these four datasets indicate their purpose. Unfortunately, the numpy / scipy formats are not easily readable without first loading them in python and then rewriting them into a data matrix or tab-delimited format. There is code in the other directories on the publication website to quickly load these datasets. I will review the ClinVar and ESP data, since that is of interest at the moment. In python (assuming the data are stored in the same directory):

```
import numpy
import scipy.sparse

# Load the feature set
X_ClinVar_ESP = numpy.load('ClinVar_ESP.X.npz')  
X_ClinVar_ESP = scipy.sparse.csr_matrix((X_ClinVar_ESP['data'], X_ClinVar_ESP['indices'], X_ClinVar_ESP['indptr']), shape=X_ClinVar_ESP['shape'])

# Load the binary response / target data.
y_ClinVar_ESP = numpy.load('ClinVar_ESP.y.npy')

# Expand sparse matrix into a dense matrix to identify shape
xDense = X_ClinVar_ESP.dense

# Show matrix dimensions - it should be (61406, 949)
xDense.shape 
```

After loading the numpy datasets in python, it is clear that the number of features (n=949) matches the number reported from the CADD paper (after expanding the categorical variables and imputing - see above). The number of rows in the matrix is a bit confusing since you would expect 20,000 rows as cited by the paper from their analysis (n_clinvar = 10,000 and n_esp = 10,000). Since the data had been encoded in numpy format there is not meta-information, like a header, for these datasets, and the data has been expanded into numerical values for input into the model already. To try to resolve what exactly was in this dataset, I sought to compare it to the data available on the github website.

### Github website datasets

The github page has a nice README on it to go over the models and how to run the code with the datasets. However, the finer details of the origins of the ClinVar and ESP datasets were a bit lost. First there is a Testing directory which contains ESP and ClinVar dataset files, but there is also a Jul2 directory that contains the same ClinVar and ESP dataset files and some imputation scripts. The files available are:
  1. clinvar.vcf.gz
  2. clinvar_20140303_pathogenic.anno_all.tsv.gz
  3. clinvar_CADD.tsv.gz
  4. clinvar_imputed.csv
  5. ESP6500.vcf.gz
  6. ESP6500SI.V2.MAF5.anno_all.tsv.gz
  7. ESP6500_CADD.tsv.gz
  8. ESP6500_imputed.csv
  9. impute2csv_mod.py
  10. impute_mod.py

According the README in the DeepCADD directory here is how the datasets were formatted:

> After constantly e-mailing the original CADD authors, they finally gave me some scripts to convert the output from their CADD online service into the proper format. All files are in the Testing/Jul2 folder. The svmlight file you need for testing is Jul2_testing.svmlight. Below are all the command lines I used to generate the file in case you want to follow:
>
> ```
> awk '{print $1"\t"$2"\t.\t"$3"\t"$5}' ESP6500SI.V2.MAF5.anno_all.tsv > ESP6500.vcf
> awk '{print $1"\t"$2"\t.\t"$3"\t"$5}' clinvar_20140303_pathogenic.anno_all.tsv > clinvar.vcf
> #submitted both vcf files to the CADD online service to get _CADD.tsv files
> cat clinvar_CADD.tsv | python impute_mod.py | python impute2csv_mod.py > clinvar_imputed.csv
> cat ESP6500_CADD.tsv | python impute_mod.py | python impute2csv_mod.py > ESP6500_imputed.csv
> python csv2svmlight.py  
> ```

This is great; however, it is still unclear as to where the two initial files originated (ESP6500SI.V2.MAF5.anno_all.tsv and clinvar_20140303_pathogenic.anno_all.tsv), and is further confused by the fact that these files contain the same columns as the "annotated CADD" files (clinvar_CADD.tsv and ESP6500_CADD.tsv). It seems likely that these two initial files were originally pulled from the CADD data files in some manner that is not fully documented. Nonethless, when you look at the contents of the .tsv files here is the summary:

| File name | Number of rows | Number of columns   | Number of unique genomic positions   | Number of header rows   |
|---|---|---|---|---|
| clinvar.vcf.gz | 29477 | 5  | 15202 | 0 |
| clinvar_20140303_pathogenic.anno_all.tsv.gz | 29479 | 90 | 15202 | 2 |
| clinvar_CADD.tsv.gz | 29317 | 90 | 15193 | 2 |
| clinvar_imputed.csv | 29316 | 950 | 15193 | 1 |
| ESP6500.vcf.gz | 34608 | 5 | 18852 | 0 |
| ESP6500SI.V2.MAF5.anno_all.tsv.gz | 34610 | 90 | 18852 | 2 |
| ESP6500_CADD.tsv.gz | 32093 | 90 | 18852 | 2 |
| ESP6500_imputed.csv | 32092 | 950 | 18852 | 1 |

If you recall the dataset that was available on the publication website in numpy / scipy format had a dimension of 61406 x 949, which appears to align exactly with the clinvar_imputed.csv and ESP6500_imputed.csv datasets together (29315 + 32091 = 61406 rows and 950 - 1 (target variable) = 949 feature columns). A couple of hanging questions that come out of looking at these files:
  1. Why does the number from the original files (clinvar_20140303_pathogenic.anno_all.tsv.gz and ESP6500SI.V2.MAF5.anno_all.tsv.gz) reduce after going through the CADD annotation process?
    * I don't know off hand, but I suspect the particular variant, meaning the reference-alternate allele combination at a particular position is not in the CADD dataset.
  2. Are the pathogenic variants in the ClinVar file curated by reviewer status at all, and are the ESP benign variants matched by allele frequency or consequence?
  3. Why are there so many duplicate genomic positions in the files?
    * This was not readily apparent after looking at the available files and the CADD online service that supposedly generated the annotations. First, the CADD online scoring service will only report a single line of annotation and score for any given position that is input. In the Information section of the [website](http://cadd.gs.washington.edu/info) it states:

> Given a set of variants, we use Ensembl Variant Effect Predictor (VEP) on the set of variants to obtain the gene model based predictions. We run VEP with the --per_gene option, which will return a "representative" transcript with the most severe effect for a certain gene. If a position overlaps multiple genes, it will return multiple annotations for this variant.

  Given this information, there should only be multiple lines for a single variant position if it lands in multiple genes, which is not the case for a majority of the duplicate entries in the files provided on the DANN github website. The only resolution that makes sense is that the raw data files provided by the CADD authors on their [download page](http://cadd.gs.washington.edu/download) were used to extract these variant annotations, since these files contain multiple entries for many genomic positions. This also begs the question whether the CADD and DANN authors used the multi-annotation per variant position files for their model testing or collapsed them down to a single entry per variant position. Based on communication with the DANN first author, Daniel Quanq, the files provided were the data used in the analysis from the paper and 10,000 were sampled from the ClinVar and 10,000 sampled from ESP datasets.

### Summary for DANN ClinVar-ESP dataset

The best source of the ClinVar-ESP dataset in its (semi) original form as provided by the DANN author(s) are in the files clinvar_CADD.tsv.gz and ESP6500_CADD.tsv.gz. The "imputed" feature files (clinvar_imputed.csv and ESP6500_imputed.csv) represent what was originally done to format the original CADD v1.0 feature set to input into a model. These contain 61406 total variants (29315 pathogenic and 32091 benign variants) and 949 expanded features.

## SSCM



## Datasets

| Data source | File short name | File description | Number of samples | Number of features | File location |
|---|---|---|---|---|---|
| ClinVar | clinvar_pathogenic_cadd_annov1 | The original ClinVar pathogenic variant dataset generated and analyzed by the CADD and DANN publications (note they sampled from this set of variants). There are multiple entries for many of the variant positions and there are no protein domain features in this version (CADD v1.0). The file is presumably extracted from the CADD preprocessed files provided on their website, but it was obtained from the DANN website. | 29315 | 90 | data/clinvar_CADD.tsv.gz |
| ClinVar | clinvar_pathogenic_cadd_annov1_imputed | The processed / imputed file generated from the clinvar_pathogenic_cadd_annov1 file. See CADD notes above and/or CADD publication supplementary information as to how the missing data was imputed and categorical variables were expanded in to N-binary variables. The scripts to perform the formatting are available in scripts/imput_mod.py and scripts/impute2csv_mod.py. This file is ready to input modeling. See clinvar_pathogenic_cadd_annov1 file description for the source of this data. | 29315 | 949 | data/clinvar_imputed.csv |
| ESP | esp_benign_cadd_annov1 | The original ESP benign variant dataset generated and analyzed by CADD and DANN publications. There are multiple entries for many of the variant positions and there are no protein domain features in this version (CADD v1.0). The file is presumably extracted from the CADD preprocessed files provided on their website, but it was obtained from the DANN website. | 32091 | 90 | data/ESP6500_CADD.tsv.gz |
| ESP | esp_benign_cadd_annov1_imputed | The processed / imputed file generated from the esp_benign_cadd_annov1 file. See CADD notes above and/or CADD publication supplementary information as to how the missing data was imputed and categorical variables were expanded in to N-binary variables. The scripts to perform the formatting are available in scripts/cadd_annofeatures_impute.py and scripts/cadd_annofeatures_impute2csv.py. This file is ready to input modeling. See esp_benign_cadd_annov1 file description for the source of this data. | 32091 | 949 | data/ESP6500_imputed.csv |

