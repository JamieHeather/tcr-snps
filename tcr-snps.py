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

__version__ = '3.3.3'


def populate_empty_imgt_fields(gene, snp_id):
    """
    The IMGT data does not come with the annotations that the gnomAD data does, and thus will have empty fields
    In order to get it to plot the same these need to be populated with empty strings
    """
    for field in ['position', 'type', 'freq', 'rsid']:
        out_snps[gene][snp_id][field] = ""
    return


if __name__ == '__main__':

    # Set up necessary data structures
    ids = coll.defaultdict(str)
    chromosomes = coll.defaultdict(int)
    snps = coll.defaultdict(fnx.double_nest)
    path_to_data = 'Raw_Files/'
    out_dir = 'Output_SNP_Files/'

    # Define pad distance, i.e. how many nt to take on either side of the variant position
    pad = 12

    # Read in gene level data
    with open(path_to_data + "TCR_Gene_Details.csv") as inf:
        for line in inf:
            if "Symbol" not in line:
                bits = line.rstrip().split(",")
                gene = bits[0]
                ens_id = bits[1]
                chromosome = bits[2]
                ids[ens_id] = gene
                chromosomes[gene] = chromosome

    # Get gnomad data
    gnomad_files = [path_to_data+x for x in os.listdir(path_to_data)
                  if x.startswith("gnomad_") is True and x.endswith(".csv") is True][::-1]

    # Get IMGT data
    imgt_files = [path_to_data+x for x in os.listdir(path_to_data)
                  if x.startswith("imgt_") is True and x.endswith(".fasta") is True]

    # Read in the SNP data (from gnomAD files harvested from gnomAD), give every row a number,
    for fl in gnomad_files:
        print "Reading in file", str(gnomad_files.index(fl)+1), "out of", str(len(gnomad_files)) + ":", fl
        ens_id = fl.split("_")[2]
        with open(fl) as inf:
            cnt = 1
            for line in inf:
                bits = line.replace("\"", "").rstrip().split(",")
                if "Chrom" in line:
                    key = bits
                else:
                    snps[ids[ens_id]][fnx.tidy(cnt)]['position'] = int(bits[1])
                    snps[ids[ens_id]][fnx.tidy(cnt)]['chr'] = int(bits[0])
                    snps[ids[ens_id]][fnx.tidy(cnt)]['rsid'] = bits[2]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['from'] = bits[3]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['to'] = bits[4]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['t_cons'] = bits[10]

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

                    snps[ids[ens_id]][fnx.tidy(cnt)]['type'] = bits[11]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['prot_consequence'] = bits[9]
                    snps[ids[ens_id]][fnx.tidy(cnt)]['freq'] = float(bits[16])
                    cnt += 1

    # Loop through IMGT file, taking the prototypic allele (*01) as reference infer and record variant positions
    imgt_snps = coll.defaultdict(fnx.double_nest)
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
                imgt_snps[gene]['*'+allele+'-1']['ref'] = fnx.capitalise_position(
                    prototype_check[position-left_pad:position+right_pad+1], left_pad)
                imgt_snps[gene]['*'+allele+'-1']['alt'] = fnx.capitalise_position(
                    seq_check[position-left_pad:position+right_pad+1], left_pad)
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
                        imgt_snps[gene][snp_id]['ref'] = fnx.multiple_capitalise(prototype_check,
                                                                                 to_capitalise)[start:end]
                        imgt_snps[gene][snp_id]['alt'] = fnx.multiple_capitalise(seq_check, to_capitalise)[start:end]
                    # Otherwise check whether the next polymorphism would fall in the pad of the current
                    elif diffs[i+1][1] - diffs[i][1] < pad:
                        to_capitalise.append(diffs[i+1][1])
                    else:
                        end = diffs[i][1] + 1 + pad
                        imgt_snps[gene][snp_id]['ref'] = fnx.multiple_capitalise(prototype_check,
                                                                                 to_capitalise)[start:end]
                        imgt_snps[gene][snp_id]['alt'] = fnx.multiple_capitalise(seq_check, to_capitalise)[start:end]
                        start = diffs[i+1][1] - pad
                        to_capitalise = [diffs[i+1][1]]

    # See whether the IMGT polymorphism is already present in the gnomAD data
    # If it is, combine the two (keeping the name from IMGT and the rest of the data from gnomAD) in third dict
    out_snps = coll.defaultdict(fnx.double_nest)

    for gene in snps.keys():
        if gene in imgt_snps.keys():
            # If the gnomAD SNP is also present in the IMGT data, find which it is and combine them
            imgt_combos = [imgt_snps[gene][x]['ref']+'|'+imgt_snps[gene][x]['alt'] for x in imgt_snps[gene]]
            for snp_id in snps[gene]:
                if snps[gene][snp_id]['ref'] not in ['', []] and snps[gene][snp_id]['alt'] not in ['', []]:
                    this_combo = snps[gene][snp_id]['ref'] + '|' + snps[gene][snp_id]['alt']
                    if this_combo in imgt_combos:
                        for imgt_id in imgt_snps[gene]:
                            imgt_combo = imgt_snps[gene][imgt_id]['ref'] + '|' + imgt_snps[gene][imgt_id]['alt']
                            if imgt_combo == this_combo:
                                out_snps[gene][imgt_id] = snps[gene][snp_id]
                                imgt_snps[gene][imgt_id]['ignore'] = True
                                break
                    else:
                        out_snps[gene][snp_id] = snps[gene][snp_id]

            # Having added all the IMGT sequences that are also in gnomAD, add the remaining IMGT SNPs
            for imgt_id in imgt_snps[gene]:
                if imgt_snps[gene][imgt_id]['ignore'] is not True:
                    out_snps[gene][imgt_id] = imgt_snps[gene][imgt_id]
                    # Need to add empty values for IMGT-only SNPs for writing section to work
                    populate_empty_imgt_fields(gene, imgt_id)

        else:
            # If this particular gene shows no variants in IMGT, just take the dict at gene level
            out_snps[gene] = snps[gene]

    # TODO (potentially) - add in option to get non-'replace' variants? Only 1 or 2, maybe not worth its

    # Write data out to those gene specific files
    print "Writing data out..."
    with open(out_dir + 'TRAV.snps', 'w') as av_file, open(out_dir + 'TRBV.snps', 'w') as bv_file, \
            open(out_dir + 'TRAJ.snps', 'w') as aj_file, open(out_dir + 'TRBJ.snps', 'w') as bj_file:

        # Add comments detailing version number
        for out_file in [av_file, bv_file, aj_file, bj_file]:
            summary_str = "# SNPs produced with tcr-snps.py, version " + str(__version__) + \
                          "\n# See https://github.com/JamieHeather/tcr-snps\n"
            out_file.write(summary_str)

        # Write headers to comments
        out_file.write("#Position,Reference,Alternative,Type,Frequency,RSID\n")

        # Write out only those SNPs for which we could find appropriate info
        for g in out_snps:
            snp_ids = out_snps[g].keys()
            snp_ids.sort()
            for s in snp_ids:
                if out_snps[g][s]['ref'] not in [[], ''] and out_snps[g][s]['alt'] not in [[], '']:
                    if out_snps[g][s]['ignore'] is not True:
                        outline = ','.join([g, str(s), str(out_snps[g][s]['position']),
                                            out_snps[g][s]['ref'], out_snps[g][s]['alt'],
                                            out_snps[g][s]['type'], str(out_snps[g][s]['freq']),
                                            out_snps[g][s]['rsid'], ]) + "\n"
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
