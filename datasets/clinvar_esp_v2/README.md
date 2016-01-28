Data generated using:

Jan. 2016 ClinVar dataset:
See clinvar directory for details on how this was pulled from ClinVar and formatted. The most likely pathogenic variants from ClinVar based on review status were extracted based on the following criteria:
1. Review status had multiple submitters, reviewed by expert panel or practice guideline.
2. No conflicting clinical significance statuses.
3. Pathogenic clinical significance status.

Set 3 = expert reviewed, practice guidelines and multiple submitters and no conflicts
`awk -F"\t" '($16==1 && $17==0) && ($9=="criteria provided, multiple submitters, no conflicts" || $9=="reviewed by expert panel" || $9=="practice guideline") { print $0 }' clinvar.combined.dedup.tsv > clinvar_20160125.set3.tsv`

The ClinVar pathogenic variant sites were then annotated using the CADD annotations:
1. scripts/get_cadd.batch.py - uses tabix to extract individual variant sites from the whole genome annotated CADD files on their webserver. Note that this is not fully complete as it does not account for the specific variant site ref and alt alleles and tabix will take a genomic position as a query and return all the results that match that position. There are often multiple entries in the CADD annotated files for a single variant site due to different alternate alleles or multiple transcript annotations or overlapping genes at the same location. This STILL needs to be resolved.
2. The CADD annotated file (clinvar_CADDv1.3_annotations.tsv.gz) was then formatted for input into a modeling framework (clinvar_CADDv1.3_annotations.imputed.csv)
  * The formatting scripts are scripts/cadd_impute_v1.3.py and scripts/cadd_impute_v1.3.tocsv.py
     * `gzip -c clinvar_CADDv1.3_annotations.tsv.gz | python scripts/cadd_impute_v1.3.py -y <STATUS=0,1> | python scripts/cadd_impute_v1.3.tocsv.py > clinvar_CADDv1.3_annotations.imputed.csv`
3. Alternatively, to avoid having multiple entries for a single variant position, the ClinVar positions were submitted to the CADD scoring [webservice](http://cadd.gs.washington.edu/score) with the option to include underlying annotation in ouput and v1.3 selected.
  * These were then processed similarly to the annotated file in (2) and output to file clinvar_CADDv1.3_singleannotations.imputed.csv

ESP6500 dataset
The ESP6500 variant sites that had been previously selected for the CADD and DANN publications were used to select a "matching" number of benign sites for the updated ClinVar variant sites selected in this dataset version. See the clinvar-esp-caddv1 directory for details on the origins of this file.
Steps to format this data:
1. Format the ESP6500_CADD.tsv file by removing the two header rows to import into R
2. R
r <- read.table("ESP6500_CADD.tsv", header=T, sep="\t")
chrom.pos <- paste(r$Chrom, r$Pos, sep=":")
r$cp <- chrom.pos
rsorted <- r[sort.list(r$cp),]

nondupIdx <- which(!duplicated(rsorted$cp))
r <- rsorted[nondupIdx,]

tg.cols <- grep("TG", colnames(r))
esp.cols <- grep("ESP", colnames(r))

min.af.cols <- apply(r[,c(tg.cols, esp.cols)], 1, min, na.rm=T)
pdf("ESP6500_CADD_Cscore.hist.pdf")
hist(r$Cscore, breaks=100)
hist(r$Cscore[which(min.af.cols >= 0.05)], breaks=100, col='red', add=T)
dev.off()

rr <- r[which(min.af.cols > 0.05), ]
print(nrow(rr)) # 10552 values left

# Sample 2853 from 10552
rsample <- rr[sample(seq(1,nrow(rr)),2853), ]

write.table(rsample, file="ESP6500_clinvar_set3_match.rsample.tsv", sep="\t", quote=F, col.names=F, row.names=F)

# Use ESP6500_clinvar_set3_match.rsample.tsv to get new CADD annotations

After these steps, the same process as described above for the ClinVar data was used to annotate and format the ESP data.

Finally, the "imputed" files ready for modeling were merged using cat:
1. Using tabix to query CADD annotated files (i.e., multiple rows per variant site) - `cat clinvar_CADDv1.3_annotations.imputed.csv <(tail -n +2 ESP6500_CADDv1.3_annotations.imputed.csv) > clinvar_ESP6500_merged.CADDv1.3_annotations.imputed.csv`
2. Using CADD scoring and annotation webservice - `cat clinvar_CADDv1.3_singleannotations.imputed.csv <(tail -n +2 ESP6500_CADDv1.3_singleannotations.imputed.csv) > clinvar_ESP6500_merged.CADDv1.3_singleannotations.imputed.csv`
