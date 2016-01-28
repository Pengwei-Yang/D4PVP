#!/usr/bin/env Rscript

options(stringsAsFactors=F)

# load what we've extracted from the XML so far
xml_extract = read.table('clinvar_table.dedup.normalized.tsv',sep='\t',comment.char='',quote='',header=T)

# load the tab-delimited summary
# > head(r)
#  X.AlleleID                      Type                                                                                       Name GeneID GeneSymbol
#1      15041                     indel NM_014855.2(AP5Z1):c.80_83delGGATinsTGCTGTAAACTGTAACTGTAAA (p.Arg27_Ala362delinsLeuLeuTer)   9907      AP5Z1
#2      15041                     indel NM_014855.2(AP5Z1):c.80_83delGGATinsTGCTGTAAACTGTAACTGTAAA (p.Arg27_Ala362delinsLeuLeuTer)   9907      AP5Z1
#3      15042                  deletion                            NM_014855.2(AP5Z1):c.1413_1426delGGACCTGCCCTGCT (p.Leu473Glyfs)   9907      AP5Z1
#4      15042                  deletion                            NM_014855.2(AP5Z1):c.1413_1426delGGACCTGCCCTGCT (p.Leu473Glyfs)   9907      AP5Z1
#5      15043 single nucleotide variant                                               NM_014630.2(ZNF592):c.3136G>A (p.Gly1046Arg)   9640     ZNF592
#6      15043 single nucleotide variant                                               NM_014630.2(ZNF592):c.3136G>A (p.Gly1046Arg)   9640     ZNF592
#    ClinicalSignificance RS...dbSNP. nsv..dbVar. RCVaccession TestedInGTR                                     PhenotypeIDs   Origin Assembly
#1             Pathogenic   397704705           - RCV000000012           N MedGen:C3150901,OMIM:613647,Orphanet:ORPHA306511 germline   GRCh37
#2             Pathogenic   397704705           - RCV000000012           N MedGen:C3150901,OMIM:613647,Orphanet:ORPHA306511 germline   GRCh38
#3             Pathogenic   397704709           - RCV000000013           N MedGen:C3150901,OMIM:613647,Orphanet:ORPHA306511 germline   GRCh37
#4             Pathogenic   397704709           - RCV000000013           N MedGen:C3150901,OMIM:613647,Orphanet:ORPHA306511 germline   GRCh38
#5 Uncertain significance   150829393           - RCV000000014           N  MedGen:C1847114,OMIM:606937,Orphanet:ORPHA83472 germline   GRCh37
#6 Uncertain significance   150829393           - RCV000000014           N  MedGen:C1847114,OMIM:606937,Orphanet:ORPHA83472 germline   GRCh38
#  Chromosome    Start     Stop Cytogenetic                   ReviewStatus                                            HGVS.c..
#1          7  4820844  4820847      7p22.1 no assertion criteria provided NM_014855.2:c.80_83delGGATinsTGCTGTAAACTGTAACTGTAAA
#2          7  4781213  4781216      7p22.1 no assertion criteria provided NM_014855.2:c.80_83delGGATinsTGCTGTAAACTGTAACTGTAAA
#3          7  4827366  4827379      7p22.1 no assertion criteria provided            NM_014855.2:c.1413_1426delGGACCTGCCCTGCT
#4          7  4787735  4787748      7p22.1 no assertion criteria provided            NM_014855.2:c.1413_1426delGGACCTGCCCTGCT
#5         15 85342440 85342440       15q25 no assertion criteria provided                               NM_014630.2:c.3136G>A
#6         15 84799209 84799209     15q25.3 no assertion criteria provided                               NM_014630.2:c.3136G>A
#                                   HGVS.p.. NumberSubmitters LastEvaluated Guidelines                         OtherIDs VariantID ReferenceAllele
#1 NP_055670.1:p.Arg27_Ala362delinsLeuLeuTer                1  Jun 29, 2010          - OMIM Allelic Variant:613653.0001         2            GGAT
#2 NP_055670.1:p.Arg27_Ala362delinsLeuLeuTer                1  Jun 29, 2010          - OMIM Allelic Variant:613653.0001         2            GGAT
#3                 NP_055670.1:p.Leu473Glyfs                1  Jun 29, 2010          - OMIM Allelic Variant:613653.0002         3  GGACCTGCCCTGCT
#4                 NP_055670.1:p.Leu473Glyfs                1  Jun 29, 2010          - OMIM Allelic Variant:613653.0002         3  GGACCTGCCCTGCT
#5                  NP_055445.2:p.Gly1046Arg                1  Jun 29, 2015          - OMIM Allelic Variant:613624.0001         4               G
#6                  NP_055445.2:p.Gly1046Arg                1  Jun 29, 2015          - OMIM Allelic Variant:613624.0001         4               G
#         AlternateAllele SubmitterCategories ChromosomeAccession
#1 TGCTGTAAACTGTAACTGTAAA                   1        NC_000007.13
#2 TGCTGTAAACTGTAACTGTAAA                   1        NC_000007.14
#3                      -                   1        NC_000007.13
#4                      -                   1        NC_000007.14
#5                      A                   1         NC_000015.9
#6                      A                   1        NC_000015.10
#
txt_download = read.table('variant_summary.txt.gz',sep='\t',comment.char='',quote='',header=T)

# subset the tab-delimited summary to desired rows and cols
colnames(txt_download) = gsub('\\.','_',tolower(colnames(txt_download)))
desired_columns = c('variantid','genesymbol','clinicalsignificance','reviewstatus','hgvs_c__','hgvs_p__')
txt_extract = subset(txt_download, assembly == 'GRCh37', select=desired_columns)
colnames(txt_extract) = c('measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p')

# join on measureset_id
combined = merge(xml_extract, txt_extract, by='measureset_id')

# re-order the columns
combined = combined[,c('chrom','pos','ref','alt','mut','measureset_id','symbol','clinical_significance','review_status','hgvs_c','hgvs_p','all_submitters','all_traits','all_pmids', 'all_molecular_consequences')]

# add some layers of interpretation on top of this
# note: we are trying to get the "overall" interpretation that is displayed in the upper right of the clinvar web pages but
# it is not in any of the available FTP downloads, so this is a stopgap
combined$pathogenic = as.integer(grepl('athogenic',combined$clinical_significance)) # 1 if at least one submission says path or likely path, 0 otherwise
combined$conflicted = as.integer(grepl('athogenic',combined$clinical_significance) & grepl('enign',combined$clinical_significance)) # 1 if at least one submission each of [likely] benign and [likely] pathogenic

write.table(combined,'clinvar.combined.tsv',sep='\t',row.names=F,col.names=T,quote=F)