#! /bin/bash

# 1. Get current clinvar positions selected for XYZ
# 2. Get CADD annotations on clinvar positions.
# 3. Format CADD output file for input into ML modle.

# 1. Get ESP6500 positions.
# 2. Get CADD annotations on ESP6500 positions.
# 3. Format ESP6500 CADD output file for input into ML model.

# ClinVar current positions for pathogenic, no conflict.

# Jan. 2016 ClinVar has 134821 variant entries in total - with 112859 unique variant positions
# Normalization results in ref == alt: 164, wrong ref: 33, invalid nucleotide: 9 

# Final processed file contains 112455 variant entries and 106295 unique variant positions.
# 40925 where "pathogenic" in clinical significance description field, 279 where "pathogenic" and "benign" in clinical significance submissions
# 1450 where review status is "criteria provided, multiple submitters, no conflicts" and pathogenic and no conflicts.
#

# output clinvar flat file.
python ~/ryanabo_gh_repos/scripts/variant_annotation/clinvar_formatting/clinvar_format.py -d $PWD -f $PWD/human_g1k_v37_decoy.fasta -s ~/ryanabo_gh_repos/scripts/variant_annotation/clinvar_formatting/joinData.R

cat clinvar.combined.dedup.tsv | cut -f1,2,8,9,15,16,17 > clinvar.combined.dedup.reduced.tsv

> table(r$pathogenic, r$conflicted)

#                    conflict
#                    0     1
# benign       0 71529     0
# pathogenic   1 40646   279

> missense <- grep("missense", r$all_molecular_consequences)
> length(missense)
# 60507

> rr <- r[which(r$conflicted==0),]
> table(rr$review_status, rr$pathogenic)

#                                                          pathogenic & no conflicts (clinical sig has path and benign)
#                                                           0     1
#  criteria provided, conflicting interpretations        1967   729
#  criteria provided, multiple submitters, no conflicts  5804  1450
#  criteria provided, single submitter                  39230 13097
#  no assertion criteria provided                        9797 23801
#  no assertion for the individual variant                 31   155
#  no assertion provided                                12747    11
#  practice guideline                                       0    23
#  reviewed by expert panel                              1953  1380
# 1. Pathogenic no conflicts = 39917
# 2. Pathogenic no conflicts, with criteria = 14547
# 3. Pathogenic no conflicts, with multiple submitters, practice guideline, expert panel = 2853

> noconflict.missense.idx <- grep("missense", rr$all_molecular_consequences)
> length(noconflict.missense.idx)
# 60252

> nconflict.missense.df <- rr[noconflict.missense.idx,]
> table(noconflict.missense.df$review_status, noconflict.missense.df$pathogenic)

#                                                           0     1
#  criteria provided, conflicting interpretations        1294   651
#  criteria provided, multiple submitters, no conflicts  3091   633
#  criteria provided, single submitter                  21902  5629
#  no assertion criteria provided                        5361 13590
#  no assertion for the individual variant                 27   133
#  no assertion provided                                 6993     9
#  practice guideline                                       0     8
#  reviewed by expert panel                               710   221

# 1. Pathogenic no conflicts = 633+5629+13590+133+9+8+221 = 20223
# 2. Pathogenic no conflicts, with criteria = 633+5629+8+221 = 6491
# 3. Pathogenic no conflicts, with multiple submitters, practice guideline, expert panel = 633+8+221 = 862

# Get clinvar pathogenic selection set 1
awk -F"\t" '$15==1 && $16==0 { print $0 }' clinvar.combined.dedup.tsv | grep -v "conflicting interpretations" > clinvar_20160122.set1.tsv

# Get clinvar pathogenic select set 3

awk -F"\t" '($15==1 && $16==0) && ($9=="criteria provided, multiple submitters, no conflicts" || $9=="reviewed by expert panel" || $9=="practice guideline") { print $0 }' clinvar.combined.dedup.tsv > clinvar_20160125.set3.tsv
