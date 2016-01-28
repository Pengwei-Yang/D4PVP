Data generated using:

Jan. 2016 ClinVar dataset:
See clinvar directory for details on how this was pulled from ClinVar and formatted. The most likely pathogenic variants from ClinVar based on review status were extracted based on the following criteria:
1. Review status has submission criteria, reviewed by expert panel or practice guideline.
2. No conflicting clinical significance statuses.
3. Pathogenic clinical significance status and "missense" for molecular consequence.

To generate the clinvar missense data set from the full clinvar flat file:
`awk -F"\t" '($16==1 && $17==0) && ($9 ~ /submitter/ || $9=="reviewed by expert panel" || $9=="practice guideline") && ($15 ~ /missense/) { print $0 }' clinvar.combined.dedup.tsv > clinvar_missense.set2.tsv`

The ClinVar pathogenic variant sites were then annotated using the CADD annotations:
1. The ClinVar positions were submitted to the CADD scoring [webservice](http://cadd.gs.washington.edu/score) with the option to include underlying annotation in ouput and v1.3 selected.
  * Get vcf file `awk -F"\t" '( $1"\t"$2"\t.\t"$3"\t"$5 ) clinvar_missense.set2.tsv > clinvar_missense.set2.vcf`
  * The CADD annotated file was formatted for input into a modeling framework (clinvar_missense_CADDv1.3_singleannotations.imputed.csv)
     * The formatting scripts are scripts/cadd_impute_v1.3.py and scripts/cadd_impute_v1.3.tocsv.py
     * `gzip -c clinvar_CADDv1.3_annotations.tsv.gz | python scripts/cadd_impute_v1.3.py -y <STATUS=0,1> | python scripts/cadd_impute_v1.3.tocsv.py > clinvar_CADDv1.3_annotations.imputed.csv`

ESP6500 dataset
1. Download CADD ESP6500 annotated file.
2. R

r <- read.table("ESP6500SI_inclAnno.tsv", header=T, skip=1, sep="\t", stringsAsFactors=F, comment.char="")
r.maf.5_10 <- r[which(r$ESP_AF >= 0.05 & r$ESP_AF <= 0.1),]
chrom.pos <- paste(r.maf.5_10$X.Chrom, r.maf.5_10$Pos, sep=":")
r.maf.5_10$cp <- chrom.pos
rsorted <- r.maf.5_10[sort.list(r.maf.5_10$cp),]
nondupIdx <- which(!duplicated(rsorted$cp))
nondups <- rsorted[nondupIdx,]
nonbenign.idx <- which(nondups$SIFTcat == "deleterious" | nondups$PolyPhenCat == "possibly_damaging" | nondups$PolyPhenCat == "probably_damaging")
nrow(nondups[-nonbenign.idx,])
esp.benign <- nondups[-nonbenign.idx,]

select same number to match clinvar set - 6491. tail -n +2 <file> | gshuf -n 6491 > ESP6500_maf5_10.tsv

These positions were extracted into a vcf formatted file and submitted to the CADD webservice and formatted the same as the clinvar file. 

Note that there is also a merged ClinVar and ESP6500 file, where the first column is the status variable (-1 for ESP6500 benign and 1 for ClinVar pathogenic_.
