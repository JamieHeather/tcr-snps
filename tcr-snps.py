#!/usr/bin/env python

"""
tcr-snps.py: Extract reference and variant sequence for SNPs and other polymorphisms in the human T-cell receptor loci
https://github.com/JamieHeather/tcr-snps
"""

import collections as coll
import os
import sys
import snpfunctions as fnx
import Levenshtein as lev

__version__ = '3.2.2'

def populate_empty_imgt_fields(gene, snp_id):
    """
    The IMGT data does not come with the annotations that the ExAC data does, and thus will have empty fields
    In order to get it to plot the same these need to be populated with empty strings
    """
    for field in ['position', 'type', 'freq', 'rsid']:
        snps[gene][snp_id][field] = ""
    return


if __name__ == '__main__':

    # Set up necessary data structures
    ids = coll.defaultdict(str)
    chromosomes = coll.defaultdict(int)
    snps = coll.defaultdict(fnx.double_nest)
    path_to_data = 'Raw_Files/'

    # Define pad distance, i.e. how many nt to take on either side of the variant position
    pad = 12

    # Read in gene level data
    with open("TCR_Gene_Details.csv") as inf:
        for line in inf:
            if "Symbol" not in line:
                bits = line.rstrip().split(",")
                gene = bits[0]
                ens_id = bits[1]
                chromosome = bits[2]
                ids[ens_id] = gene
                chromosomes[gene] = chromosome

    # Get exac data
    exac_files = [path_to_data+x for x in os.listdir(path_to_data)
                  if x.startswith("exac_") is True and x.endswith(".csv") is True][::-1]

    # Get IMGT data
    imgt_files = [path_to_data+x for x in os.listdir(path_to_data)
                  if x.startswith("imgt_") is True and x.endswith(".fasta") is True]

    # Read in the SNP data (from ExAC files harvested from gnomAD), give every row a number,
    for fl in exac_files:
        print "Reading in file", str(exac_files.index(fl)+1), "out of", str(len(exac_files)) + ":", fl
        ens_id = fl.split("_")[2]
        with open(fl) as inf:
            cnt = 1
            for line in inf:
                if "Chrom" in line:
                    key = line
                else:
                    bits = line.replace("\"", "").rstrip().split(",")
                    snps[ids[ens_id]][fnx.tidy(cnt)]['position'] = int(bits[1])
                    snps[ids[ens_id]][fnx.tidy(cnt)]['chr'] = int(bits[0])
                    snps[ids[ens_id]][fnx.tidy(cnt)]['rsid'] = bits[2]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['from'] = bits[3]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['to'] = bits[4]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'] = bits[8]

                    # Get SNP, if we're dealing with a single nt change (that's a difference base!)
                    if len(snps[ids[ens_id]][fnx.tidy(cnt)]['from']) == 1 \
                            and len(snps[ids[ens_id]][fnx.tidy(cnt)]['to']) == 1 \
                            and 'dup' not in snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'] \
                            and 'del' not in snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'] \
                            and snps[ids[ens_id]][fnx.tidy(cnt)]['from'] != snps[ids[ens_id]][fnx.tidy(cnt)]['to']:
                        snps[ids[ens_id]][fnx.tidy(cnt)]['ref'], snps[ids[ens_id]][fnx.tidy(cnt)]['alt'] = \
                            fnx.get_snp(snps[ids[ens_id]][fnx.tidy(cnt)]['position'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['chr'], pad,
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['from'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['to'])

                    # Otherwise look for deletions...
                    elif 'del' in snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons']:
                        snps[ids[ens_id]][fnx.tidy(cnt)]['ref'], snps[ids[ens_id]][fnx.tidy(cnt)]['alt'] = \
                            fnx.get_del(snps[ids[ens_id]][fnx.tidy(cnt)]['position'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['chr'], pad,
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['from'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'])

                    # ... or duplications
                    elif 'dup' in bits[8]:
                        snps[ids[ens_id]][fnx.tidy(cnt)]['ref'], snps[ids[ens_id]][fnx.tidy(cnt)]['alt'] = \
                            fnx.get_dup(snps[ids[ens_id]][fnx.tidy(cnt)]['position'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['chr'], pad,
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['from'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['to'],
                                        snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'])

                    # Otherwise I don't think we can hazard an accurate guess of what it should be, so ignore
                    else:
                        snps[ids[ens_id]][fnx.tidy(cnt)]['ref'], snps[ids[ens_id]][fnx.tidy(cnt)]['alt'] = '', ''

                    snps[ids[ens_id]][fnx.tidy(cnt)]['type'] = bits[9]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['prot_consequence'] = bits[6]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['freq'] = float(bits[14])
                    cnt += 1

    # Loop through IMGT file, taking the prototypic allele (*01) as reference infer and record variant positions

    generic_error = "Double check input IMGT fasta file:", \
        "must contain ordered alleles (including *01), only for genes with polymorphisms"
    for readid, seq, qual in fnx.readfq(open(imgt_files[0], 'rU')):
        fullgene = readid.split('|')
        gene = fullgene[1].split('*')[0]
        allele = fullgene[1].split('*')[1]
        #
        if allele == '01':
            prototype = seq
            prototype_gene = gene
            prototype_fullgene = fullgene
        #
        else:
            # Basic check
            if prototype_gene != gene:
                print "Error: prototype gene is different to comparison gene."
                print generic_error
                sys.exit()
            #
            # Check whether sequences are same length
            if len(prototype) != len(seq):
                # Establish which end of which sequence is truncated, and truncate the other accordingly
                prototype_check, seq_check = fnx.trim_sequences(prototype, seq)
            else:
                prototype_check, seq_check = prototype, seq
            #
            # Use Levenshtein package to infer differences
            diffs = lev.editops(prototype_check, seq_check)
            if len(diffs) == 0:
                print "Error: no difference detected between", prototype_gene, "and", gene + "."
                print generic_error
                sys.exit()
            elif len(diffs) == 1 and 'replace' in diffs[0]:
                # If the polymorphism is less than one pad away from 5', just go from start of what we have
                # (Polymorphisms less than a pad from the 3' end will have similarly truncated sequence contexts)
                left_pad, right_pad = pad, pad
                position = diffs[0][1]
                if position < pad:
                    left_pad = position
                if len(prototype_check) - position < pad:
                    right_pad = len(prototype_check) - position
                snps[gene]['*'+allele+'-1']['ref'] = fnx.capitalise_position(
                    prototype_check[position-left_pad:position+right_pad+1], left_pad)
                snps[gene]['*'+allele+'-1']['alt'] = fnx.capitalise_position(
                    seq_check[position-left_pad:position+right_pad+1], left_pad)
                populate_empty_imgt_fields(gene, '*'+allele+'-1')
            elif len(diffs) > 1 and set([x[0] for x in diffs]) == set(['replace']):
                # If more than one difference relative to reference, output each separately unless their pad overlaps
                if diffs[0][1] < pad:
                    start = 0
                else:
                    start = diffs[0][1] - pad
                to_capitalise = [diffs[0][1]]
                for i in range(len(diffs)):
                    snp_id = '*'+allele+'-'+str(i+1)
                    # If have hit end, round it off
                    if i+1 == len(diffs):
                        end = diffs[i][1] + 1 + pad
                        snps[gene][snp_id]['ref'] = fnx.multiple_capitalise(prototype_check, to_capitalise)[start:end]
                        snps[gene][snp_id]['alt'] = fnx.multiple_capitalise(seq_check, to_capitalise)[start:end]
                        populate_empty_imgt_fields(gene, snp_id)
                    # Otherwise check whether the next polymorphism would fall in the pad of the current
                    elif diffs[i+1][1] - diffs[i][1] < pad:
                        to_capitalise.append(diffs[i+1][1])
                    else:
                        end = diffs[i][1] + 1 + pad
                        snps[gene][snp_id]['ref'] = fnx.multiple_capitalise(prototype_check, to_capitalise)[start:end]
                        snps[gene][snp_id]['alt'] = fnx.multiple_capitalise(seq_check, to_capitalise)[start:end]
                        populate_empty_imgt_fields(gene, snp_id)
                        start = diffs[i+1][1] - pad
                        to_capitalise = [diffs[i+1][1]]

    # TODO - add in option to get none-replace variants
    # TODO - check whether SNPs covered by ExAC? If so try to retain both infos. Then output all to same files

    # Write data out to those gene specific files
    print "Writing data out..."
    with open('TRAV.snps', 'w') as av_file, open('TRBV.snps', 'w') as bv_file, \
            open('TRAJ.snps', 'w') as aj_file, open('TRBJ.snps', 'w') as bj_file:

        # Add comments detailing version number
        for out_file in [av_file, bv_file, aj_file, bj_file]:
            summary_str = "# SNPs produced with tcr-snps.py, version " + str(__version__) + \
                          "\n# See https://github.com/JamieHeather/tcr-snps\n"
            out_file.write(summary_str)

        for g in snps:
            snp_ids = snps[g].keys()
            snp_ids.sort()
            for s in snp_ids:
                if snps[g][s]['ref'] not in [[], ''] and snps[g][s]['alt'] not in [[], '']:
                    outline = ','.join([g, str(s), str(snps[g][s]['position']),
                                        snps[g][s]['ref'], snps[g][s]['alt'],
                                        snps[g][s]['type'], str(snps[g][s]['freq']),
                                        snps[g][s]['rsid'], ]) + "\n"
                    if 'TRAV' in g or 'TRDV' in g:
                        av_file.write(outline)
                    elif 'TRBV' in g:
                        bv_file.write(outline)
                    elif 'TRAJ' in g:
                        aj_file.write(outline)
                    elif 'TRBJ' in g:
                        bj_file.write(outline)
                    else:
                        print "Error: missing chain info?"

        # TODO : write SNPs to own directory?


sys.exit()
# TODO use the following kind of approach to go through and check whether the IMGT-data is covered by the ExAC data - if it is just change the name of the ExAC data so all info is retained
for gg in snps.keys():
    for x in snps[gg]:
        for y in snps[gg]:
            if x != y:
                # if snps[gg][x]['ref'] != "" and snps[gg][y]['ref'] == "":
                    if snps[gg][x]['ref'] == snps[gg][y]['ref'] and snps[gg][x]['alt'] == snps[gg][y]['alt']:
                        print gg, x,y
