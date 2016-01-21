# D4PVP

A couple of methods have recently been published aiming to provide metrics or predictions to determine the pathogenicity of genomic mutations. This repository contains documentation of the available datasets from these published methods and scripts to help obtain the data and code (if available) to recreate the previous analyses.

### Previously published methods and datasets

| Dataset | Description | Number of samples (path / benign) | Number of features | Number of expanded features |
|---|---|---|---|---|
| [ClinVar-ESP6500SI.V2 CADDv1 features](https://github.com/ryanabo/D4PVP/blob/master/datasets/clinvar_esp_caddv1/README.md) | The original ClinVar-ESP variant dataset generated and analyzed by the CADD and DANN publications. | 61406 (29315 / 32091) | 90 | 949 |
| [CADDv1.3 training](https://github.com/ryanabo/D4PVP/blob/master/datasets/cadd/README.md) | An updated dataset with CADDv1.3 features, which have added new features like protein domains. | 35M (17521530 / 17521530) | 102 | 1063 |
| [CADDv1.3 testing](https://github.com/ryanabo/D4PVP/blob/master/datasets/cadd/README.md) | An updated testing set with CADDv1.3 features | 350051 (174800 / 175251) | 102 | 1063 |  
| [DANN training, testing, validation CADDv1 features](https://github.com/ryanabo/D4PVP/blob/master/datasets/dann/README.md) | A reuse of the CADD training and testing datasets for a different model. The ClinVar-ESP data is also reused from the CADD authors. | | | |
| [SSCM training, testing, validation](https://github.com/ryanabo/D4PVP/blob/master/datasets/sscm/README.md) | 
