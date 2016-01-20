# PVP

CADD



DANN: a deep learning approach for annotating the pathogenicity of genetic variants

DANN aimed to improve upon the CADD pathogenic variant prediction using a "non-linear" modeling approach using neural nets. The same dataset used by CADD was used by DANN for training, validation and testing. They also analyze real data from ClinVar and ESP, similar to the CADD paper.

Data and documentation sources:
  * Publication [website](https://cbcl.ics.uci.edu/public_data/DANN/).
  * Github [website](https://github.com/uci-cbcl/DeepCADD)

After reading the DANN paper and analyzing the publication website it was still unclear as to whether the data used for the reported analysis for available for use in either of the identified data sources. I will try to lay out all the convolutions and my conclusions to moving forward. First, the DANN publication cites the use of 10,000 pathogenic and benign variants each for their ROC curve. The publication website contains multiple datasets from the publication:
  1. testing.X.npz, testing.svmlight.gz, testing.y.npy
  2. training.X.npz, training.svmlight.gz, training.y.npy
  3. validation.X.npz, validation.svmlight.gz, validation.y.npy
  4. ClinVar_ESP.X.npz, ClinVar_ESP.svmlight.gz, ClinVar_ESP.y.npy

These four datasets are stored in numpy / scipy formats (.npz and .npy) and a custom svmlight gzip format. Hopefully the naming of these four datasets indicate their purpose. Unfortunately, the numpy / scipy formats are not easily readable without first loading them in python and then rewriting them into a data matrix or tab-delimited format. There is code in the other directories on the publication website to quickly load these datasets. In python (assuming the data are stored in the same directory):

`import numpy
import scipy.sparse

X_ClinVar_ESP = numpy.load('ClinVar_ESP.X.npz')  
X_ClinVar_ESP = scipy.sparse.csr_matrix((X_ClinVar_ESP['data'], X_ClinVar_ESP['indices'], X_ClinVar_ESP['indptr']), shape=X_ClinVar_ESP['shape'])
y_ClinVar_ESP = numpy.load('ClinVar_ESP.y.npy')`

Data sets:
  * CADD
    * Testing and training - 
  * DANN
    * ClinVar-ESP testing dataset

