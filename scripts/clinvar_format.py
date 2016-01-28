#! /usr/bin/env python

'''
'''

import gzip
import re
import sys
import os
import argparse
from collections import defaultdict
import xml.etree.ElementTree as ET
import pysam


def parse_clinvar_tree(handle, dest=sys.stdout, verbose=True):
    '''
    '''

    mentions_pubmed_regex = '(?:PubMed|PMID)(.*)' # group(1) will be all the text after the word PubMed or PMID
    extract_pubmed_id_regex = '[^0-9]+([0-9]+)[^0-9](.*)' # group(1) will be the first PubMed ID, group(2) will be all remaining text
    dest.write(('\t'.join( ['chrom', 'pos', 'ref', 'alt', 'mut', 'measureset_id', 'all_submitters', 'all_traits', 'all_pmids', 'all_molecular_consequences'] ) + '\n').encode('utf-8'))
    counter = 0
    skipped_counter = defaultdict(int)
    for event, elem in ET.iterparse(handle):
        if event == 'end' and elem.tag == 'ClinVarSet':
            # find the GRCh37 VCF representation
            sequence_locations = elem.findall('.//SequenceLocation')
            grch37 = None
            for sequence_location in sequence_locations:
                if sequence_location.attrib.get('Assembly') == 'GRCh37':
                    if all(entry is not None for entry in [sequence_location.attrib.get(key) for key in ['referenceAllele','alternateAllele','start','Chr']]):
                        grch37 = sequence_location
            if grch37 is None:
                skipped_counter['missing SequenceLocation'] += 1
                continue # don't bother with variants that don't have a VCF location
            else:
                chrom = grch37.attrib['Chr']
                pos = grch37.attrib['start']
                ref = grch37.attrib['referenceAllele']
                alt = grch37.attrib['alternateAllele']
            measureset = elem.findall('.//MeasureSet')
            if measureset is None:
                skipped_counter['missing MeasureSet'] += 1
                continue # skip variants without a MeasureSet ID
            measureset_id = measureset[0].attrib['ID']
            mutant_allele = 'ALT' # default is that each entry refers to the alternate allele

            # Old code
            # attributes = elem.findall('.//Attribute')
            # molConsequences = []
            # for attribute in attributes:
            #     attribute_type = attribute.attrib.get('Type')
            #     if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
            #         if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
            #             mutant_allele = 'REF' # that is their funny way of saying this assertion refers to the reference allele

            attributeSets = elem.findall('.//AttributeSet')
            molConsequences = []
            for attributeSet in attributeSets:
                attributes = attributeSet.findall('.//Attribute')
                for attribute in attributes:
                    attribute_type = attribute.attrib.get('Type')
                    if attribute_type is not None and "HGVS" in attribute_type and "protein" not in attribute_type: # if this is an HGVS cDNA, _not_ protein, annotation:
                        if attribute.text is not None and "=" in attribute.text: # and if there is an equals sign in the text, then
                            mutant_allele = 'REF' # that is their funny way of saying this assertion refers to the reference allele
                    if attribute_type == "MolecularConsequence":
                        xrefs = attributeSet.findall('.//XRef')
                        for xref in xrefs:
                            if xref.attrib.get("ID").find("SO:") == -1:
                                molConsequences.append(attribute.text + "(" + xref.attrib.get("ID") + ")")
            # find all the Citation nodes, and get the PMIDs out of them
            citations = elem.findall('.//Citation')
            pmids = []
            for citation in citations:
                pmids += [id_node.text for id_node in citation.findall('.//ID') if id_node.attrib.get('Source')=='PubMed']
            # now find the Comment nodes, regex your way through the comments and extract anything that appears to be a PMID
            comments = elem.findall('.//Comment')
            comment_pmids = []
            for comment in comments:
                mentions_pubmed = re.search(mentions_pubmed_regex,comment.text)
                if mentions_pubmed is not None and mentions_pubmed.group(1) is not None:
                    remaining_text = mentions_pubmed.group(1)
                    while True:
                        pubmed_id_extraction = re.search(extract_pubmed_id_regex,remaining_text)
                        if pubmed_id_extraction is None:
                            break
                        elif pubmed_id_extraction.group(1) is not None:
                            comment_pmids += [pubmed_id_extraction.group(1)]
                            if pubmed_id_extraction.group(2) is not None:
                                remaining_text = pubmed_id_extraction.group(2)
            all_pmids = list(set(pmids + comment_pmids))
            # now find any/all submitters
            submitters = []
            submitter_nodes = elem.findall('.//ClinVarSubmissionID')
            for submitter_node in submitter_nodes:
                if submitter_node.attrib is not None and submitter_node.attrib.has_key('submitter'):
                    submitters.append(submitter_node.attrib['submitter'])
            # now find the disease(s) this variant is associated with
            traitsets = elem.findall('.//TraitSet')
            all_traits = []
            for traitset in traitsets:
                trait_type = ''
                trait_values = []
                if traitset.attrib is not None:
                    trait_type = str(traitset.attrib.get('Type'))
                    disease_name_nodes = traitset.findall('.//Name/ElementValue')
                    for disease_name_node in disease_name_nodes:
                        if disease_name_node.attrib is not None:
                            if disease_name_node.attrib.get('Type') == 'Preferred':
                                trait_values.append(disease_name_node.text)
                all_traits += trait_values

            # now we're done traversing that one clinvar set. print out a cartesian product of accessions and pmids
            # outLst = [chrom, pos, ref, alt, mutant_allele, measureset_id, ';'.join(submitters), ';'.join(all_traits), ','.join(all_pmids), ','.join(molConsequences)]
            # print 'Output list of values', outLst
            dest.write(('\t'.join([chrom, pos, ref, alt, mutant_allele, measureset_id, ';'.join(submitters), ';'.join(all_traits), ','.join(all_pmids), ','.join(molConsequences)]) + '\n').encode('utf-8'))
            counter += 1
            if counter % 100 == 0:
                dest.flush()
            if verbose:
                sys.stderr.write("{0} entries completed, {1}, {2} total \r".format(
                    counter, ', '.join('%s skipped due to %s' % (v, k) for k, v in skipped_counter.items()),
                    counter + sum(skipped_counter.values())))
                sys.stderr.flush()
            elem.clear()

def get_handle(path):
    if path[-3:] == '.gz':
        handle = gzip.open(path)
    else:
        handle = open(path)
    return (handle)

def get_clinvar_data_files(outDir):
    xmlFn = 'ClinVarFullRelease_00-latest.xml.gz'
    summaryFn = 'variant_summary.txt.gz'
    clinVarDataDict = {'xmlFn': os.path.join(outDir, xmlFn),
                       'summaryFn': os.path.join(outDir, summaryFn)}
    if not os.path.isfile(clinVarDataDict['xmlFn']):
        cmd1 = 'wget -P %s ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/%s' % (outDir, xmlFn)
        os.system(cmd1)
    if not os.path.isfile(clinVarDataDict['summaryFn']):
        cmd2 = 'wget -P %s ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/%s' % (outDir, summaryFn)
        os.system(cmd2)
    return clinVarDataDict

def sort_dedup_clinvar_records(inFn, outFn):
    '''Run through a sorted clinvar_table.tsv file from the parse_clinvar_xml script, and make it unique on CHROM POS REF ALT
    '''

    outF = open(outFn, 'w')
    tmpFn = os.path.join(os.path.dirname(inFn), 'tmpSort.tsv')

    # Sort header row
    sorts = 'head -n 1 %s > %s; ' % (inFn, tmpFn)
    # Numerically sort chroms 1-22
    sorts += 'cat %s | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> %s; ' % (inFn, tmpFn)
    # Lexicographically sort non-numerical chroms at end
    sorts += 'cat %s | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> %s' % (inFn, tmpFn)
    os.system(sorts)
    # cat clinvar_table_raw.tsv | head -1 > clinvar_table_sorted.tsv # header row
    # cat clinvar_table_raw.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv 
    # cat clinvar_table_raw.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_table_sorted.tsv 
    tmpF = open(tmpFn, 'r')
    # Dedup - python dedup_clinvar.py < clinvar_combined_sorted.tsv > clinvar_combined_sorted_dedup.tsv
    header = tmpF.readline()
    outF.write(header)
    columnNames = header.strip('\n').split('\t')
    firstDataLine = tmpF.readline()
    lastData = dict(zip(columnNames, firstDataLine.strip('\n').split('\t')))
    lastUniqueId = '_'.join([lastData['chrom'], str(lastData['pos']), lastData['ref'], lastData['alt']])
    counter = 0
    for line in tmpF.readlines():
        data = dict(zip(columnNames, line.strip('\n').split('\t')))
        uniqueId = '_'.join([data['chrom'], str(data['pos']), data['ref'], data['alt']])
        if uniqueId == lastUniqueId:
            data = dedup_records(data, lastData)
        else:
            # note that using a comprehension instead of just data.values() preserves column order
            outF.write('\t'.join([lastData[colname] for colname in columnNames])+'\n')
        lastData = data
        lastUniqueId = uniqueId
        counter += 1
    outF.write('\t'.join([lastData[colName] for colName in columnNames])+'\n')
    outF.close()

def dedup_records(data1, data2):
    '''
    De-duplicate two ClinVar records
    '''
    # if one is REF and one is ALT, discard the REF
    if data1['mut'] == 'REF' and data2['mut'] == 'ALT':
        return data2
    elif data1['mut'] == 'ALT' and data2['mut'] == 'REF':
        return data1
    else:
        combined_data = data1 # this sets defaults, now we fix it:
    # discard one MeasureSet ID if there are two
    if data1['measureset_id'] != data2['measureset_id']:
        combined_data['measureset_id'] = str(min(map(int,[data1['measureset_id'],data2['measureset_id']])))
    combined_data['all_pmids'] = ','.join(set(data1['all_pmids'].split(',') + data2['all_pmids'].split(',')))
    combined_data['all_submitters'] = ';'.join(set(data1['all_submitters'].split(';') + data2['all_submitters'].split(';')))
    combined_data['all_traits'] = ';'.join(set(data1['all_traits'].split(';') + data2['all_traits'].split(';')))    
    return combined_data

'''
This script is a python implementation of the algorithm for variant
normalization described by Tan et al 2015:
Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants.
Bioinformatics. 2015 Jul 1;31(13):2202-4. doi: 10.1093/bioinformatics/btv112.
Epub 2015 Feb 19. PubMed PMID: 25701572.
The authors have made a C++ implementation available in vt as vt normalize
And the source code is viewable here: https://github.com/atks/vt
For our purposes, we wanted a Python implementation so that we could
build end-to-end scripts in Python.
If you use this, please cite Tan et al 2015.
A note about when this is useful. In VCFs generated with GATK (or probably
other tools) from short read sequencing data, variants are already left-aligned
but may be non-minimal to the extent that indels overlap with other variants.
For those cases, minimal_representation.py is sufficient to convert variants
to minimal representation. However, genomic coordinates converted from HGVS
(we have encountered this when parsing the ClinVar XML dump) may be not only
non-minimal but also right-aligned rather than left-aligned, and may contain
hyphens. For those situations, use this script (or just run vt normalize).
Usage: normalize.py -R $b37ref < bad_file.txt > good_file.txt
'''

class RefEqualsAltError(Exception):
    '''
    An Error class for rare cases where REF == ALT (seen in ClinVar XML)
    '''
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InvalidNucleotideSequenceError(Exception):
    '''
    An Error class for REF or ALT values that are not valid nucleotide sequences
    '''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class WrongRefError(Exception):
    '''
    An Error class for variants where the REF does not match the reference genome
    '''
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
   

def normalize(pysam_fasta, chrom, pos, ref, alt):
    '''
    Accepts a pysam FastaFile object pointing to the reference genome, and
    chrom, pos, ref, alt genomic coordinates, and normalizes them.
    '''
    pos = int(pos) # make sure position is an integer
    ref = ref.upper()
    alt = alt.upper()

    # Remove variants that contain invalid nucleotides
    if any(letter not in ['A', 'C', 'G', 'T', 'N', '-'] for letter in ref + alt):
        raise InvalidNucleotideSequenceError('Invalid nucleotide sequence: %s %s %s %s'%(chrom, pos, ref, alt))
    # use blanks instead of hyphens
    if ref == '-':
        ref = ''
    if alt == '-':
        alt = ''
    # check whether the REF is correct
    true_ref = pysam_fasta.fetch(chrom, pos - 1, pos - 1 + len(ref))
    if ref != true_ref:
        raise WrongRefError('Incorrect REF value: %s %s %s %s (actual REF should be %s)'%(chrom, pos, ref, alt, true_ref))
    # Prevent infinte loops in cases where REF == ALT.
    # We have encountered this error in genomic coordinates from the ClinVar XML file
    if ref == alt:
        raise RefEqualsAltError('The REF and ALT allele are the same: %s %s %s %s'%(chrom, pos, ref, alt))
    # Time-saving shortcut for SNPs that are already minimally represented
    if len(ref) == 1 and len(alt) == 1 and ref in ['A','C','G','T'] and alt in ['A','C','G','T']:
        return chrom, pos, ref, alt
    # This first loop left-aligns and removes excess nucleotides on the right.
    # This is Algorithm 1 lines 1-6 from Tan et al 2015
    keep_working = True
    while keep_working:
        keep_working = False
        if len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            keep_working = True
        if len(ref) == 0 or len(alt) == 0:
            preceding_base = pysam_fasta.fetch(chrom, pos-2, pos-1)
            ref = preceding_base + ref
            alt = preceding_base + alt
            pos = pos - 1
            keep_working = True
    # This second loop removes excess nucleotides on the left. This is Algorithm 1 lines 7-8.
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos = pos + 1
    return chrom, pos, ref, alt

def normalize_tab_delimited_file(inFn, outFn, refFn, verbose=True):
    '''
    This function takes a tab-delimited file with a header line containing columns
    named chrom, pos, ref, and alt, plus any other columns. It normalizes the
    chrom, pos, ref, and alt, and writes all columns out to another file.
    '''
    inFile = open(inFn, 'r')
    outFile = open(outFn, 'w')
    fastaFile = pysam.FastaFile(refFn)
    header = inFile.readline() # get header of input file
    columns = header.strip('\n').split('\t')  # parse col names 
    outFile.write('\t'.join(columns) + '\n') # write header line plus the CpG col to be generated
    counter = 0
    ref_equals_alt = 0
    wrong_ref = 0
    invalid_nucleotide = 0
    for line in inFile.readlines():
        data = dict(zip(columns, line.strip('\n').split('\t')))
        # fill the data with blanks for any missing data
        for column in columns:
            if column not in data.keys():
                data[column] = ''
        pos = int(data['pos'])
        try:
            data['chrom'], pos, data['ref'], data['alt'] = normalize(fastaFile, data['chrom'], pos, data['ref'], data['alt'])
        except RefEqualsAltError as e:
            sys.stderr.write('\n' + str(e) + '\n')
            ref_equals_alt += 1
            continue
        except WrongRefError as e:
            sys.stderr.write('\n' + str(e) + '\n')
            wrong_ref += 1
            continue
        except InvalidNucleotideSequenceError as e:
            sys.stderr.write('\n' + str(e) + '\n')
            invalid_nucleotide += 1
            continue
        data['pos'] = str(pos)
        outFile.write('\t'.join([data[column] for column in columns]) + '\n')
        counter += 1
        if verbose:
            sys.stderr.write("\r%s records processed\n"%(counter))
    outFile.write('\n\n')
    outFile.close()
    if verbose:
        sys.stderr.write("Final counts of variants discarded:\nREF == ALT: %s\nWrong REF: %s\nInvalid nucleotide: %s\n"%(ref_equals_alt, wrong_ref, invalid_nucleotide))

def test_normalize(pysam_fasta):
    '''
    Battery of test cases for normalize
    '''
    sys.stdout.write(str(normalize(pysam_fasta, '7', 117199646, 'CTT', '-'))+'\n') # HGVS translation of CFTR p.F508del, should be ('7', 117199644, 'ATCT', 'A')
    sys.stdout.write(str(normalize(pysam_fasta, '13', 32914438, 'T', '-'))+'\n') # HGVS translation of a BRCA2 Ashkenazi founder variant, should be ('13', 32914437, 'GT', 'G')


PARSER = argparse.ArgumentParser(description='Program to format the ClinVar data files into a flat file.', usage='%(prog)s [options]', add_help=True)
PARSER.add_argument('-d', '--directory', dest='analysisDir', default=None, help='Full path to the analysis directory')
PARSER.add_argument('-f', '--reference_fasta', dest='referenceFastaFn', default=None, help='Full path to the reference fasta file (GRCh37)')
PARSER.add_argument('-s', '--join_script', dest='joinScriptFn', default=None, help='Full path to the R script to join tables.')

pArgs = PARSER.parse_args()

# Set output directory to current directory if not set in options.
outDir = os.getcwd() if pArgs.analysisDir is None else pArgs.analysisDir

# Get clinvar datasets in XML and flat file format.
clinVarDataDict = get_clinvar_data_files(outDir)

# Extract the GRCh37 coordinates, mutant alleles, measureset ID and pubmed IDs from it.
clinVarDataDict['rawTable'] = os.path.join(outDir, 'clinvar_table.raw.tsv')
parse_clinvar_tree(get_handle(clinVarDataDict['xmlFn']), dest=open(clinVarDataDict['rawTable'], 'w'))

# 4. De-duplicate records
clinVarDataDict['dedupTable'] = os.path.join(outDir, 'clinvar_table.dedup.tsv')
sort_dedup_clinvar_records(inFn=clinVarDataDict['rawTable'], outFn=clinVarDataDict['dedupTable'])

# 5. Normalize variants by converting to minimal representation and left-align
clinVarDataDict['dedupNormalizedTable'] = os.path.join(outDir, 'clinvar_table.dedup.normalized.tsv')
normalize_tab_delimited_file(inFn=clinVarDataDict['dedupTable'], outFn=clinVarDataDict['dedupNormalizedTable'], refFn=args.referenceFastaFn)

# 6. Join information from the tab-delimited summary to the normalized genomic coordinates.
clinVarDataDict['combined'] = os.path.join(outDir, 'clinvar.combined.tsv')
joinDataRscript = pArgs.joinScriptFn
joinDataCmd = 'Rscript joinData.R'
os.system(joinDataCmd)

# 7. Sort again by genomic coordinates
# cat clinvar_combined.tsv | head -1 > clinvar_combined_sorted.tsv # header row
# cat clinvar_combined.tsv | tail -n +2 | egrep -v "^[XYM]" | sort -k1,1n -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # numerically sort chroms 1-22
# cat clinvar_combined.tsv | tail -n +2 | egrep "^[XYM]" | sort -k1,1 -k2,2n -k3,3 -k4,4 >> clinvar_combined_sorted.tsv # lexicographically sort non-numerical chroms at end

# 8. De-duplicate again because the tab-delimited file contains duplicates.
clinVarDataDict['combinedDedup'] = os.path.join(outDir, 'clinvar.combined.dedup.tsv')
sort_dedup_clinvar_records(inFn=clinVarDataDict['combined'], outFn=clinVarDataDict['combinedDedup'])

# 9. Create a text file.
# cp clinvar_combined_sorted_dedup.tsv clinvar.tsv
# gzip -c clinvar.tsv > clinvar.tsv.gz  # create compressed version

# 10. Create a vcf formatted file.

# Clean up files.
