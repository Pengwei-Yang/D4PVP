See scripts/get_clinvar.sh

Final formatted ClinVar table is clinvar.combined.dedup.tsv

Steps in obtaining and formatting ClinVar data from NCBI ftp site:
1. Download full data set in XML format and flat txt file (these contain different values).
2. Format and merge the XML and flat txt files to generate a single flat tsv file with as much information as desired from both files.
  * clinvar_format.py -d <output directory> -f <fasta file> -s <Rscript for joining files>
  * Note the caviats of using this script to process the ClinVar data (i.e. it will not keep any variant without a sequence location).
  * Generates a clinvar.combined.dedup.tsv file.
3. Summarize the number of variants based on different criteria - see scripts/get_clinvar.sh - a work in progress...
4. Subset the clinvar dataset based on review status, clinical significance, molecular consequence...
  * This can be done using awk on the file, such as:

Obtain the pathogenic variants that are most likely pathogenic based on review status:
`awk -F"\t" '($16==1 && $17==0) && ($9=="criteria provided, multiple submitters, no conflicts" || $9=="reviewed by expert panel" || $9=="practice guideline") { print $0 }' clinvar.combined.dedup.tsv > clinvar_pathogenic.set3.tsv`

Obtain the pathogenic variants that have some criteria for submission, no conflicting submissions, and may be missense consequence based on a transcript annotation:
`awk -F"\t" '($16==1 && $17==0) && ($9 ~ /submitter/ || $9=="reviewed by expert panel" || $9=="practice guideline") && ($15 ~ /missense/) { print $0 }' clinvar.combined.dedup.tsv > clinvar_missense.set2.tsv`
