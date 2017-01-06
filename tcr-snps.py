#!/usr/bin/env python

"""
tcr-snps.py: Extract reference and variant sequence for SNPs and other polymorphisms in the human T-cell receptor loci
https://github.com/JamieHeather/tcr-snps
"""

import collections as coll
import os
import urllib2
import string

__version__ = '3.0.1'


def rev_comp(seq):
    """
    Reverse complement a DNA strand (using only a default python module)
    """
    return seq.translate(string.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


def nest():
    """
    Use to generate a nested default dict
    """
    return coll.defaultdict(list)


def double_nest():
    """
    Use to generate a doubly nested default dict
    """
    return coll.defaultdict(nest)


def tidy(number):
    """
    Poorly named, but exists to give counts with leading zeroes while saving space on busy lines
    """
    return str(number).zfill(4)


def get_hg19_seq(chrm, seq_from, seq_to):
    """
    Given an hg19 range, return the nucleotides via the UCSC DAS server.
    Used to extract the reference nt sequence, from which the alternative variant can be determined.
    """
    base_url = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr"
    page = urllib2.urlopen(base_url + str(chrm) + ":" + str(seq_from) + "," + str(seq_to))
    contents = []
    for line in page:
        if "<" not in line:
            contents.append(line.rstrip())
    full_seq = "".join(contents).upper()
    return full_seq


def get_snp(site, chrm, padding, alt_from, alt):
    """
    Takes a genomic position, a chromosome number, a value by which to pad the position, the original sequence
     and the variant nucleotide.
    Returns two values: the reference and alternative sequences
    """
    full_seq = get_hg19_seq(int(chrm), seq_from=int(site)-padding, seq_to=int(site)+padding)
    ref_seq = full_seq[:padding].lower() + full_seq[padding:-padding].upper() + full_seq[-padding:].lower()
    if full_seq[padding:-padding].upper() != alt_from.upper():
        print "WEPAKRelkae/wjrlkewa"
        return '', ''
    alt_seq = full_seq[:padding].lower() + alt.upper() + full_seq[-padding:].lower()
    return ref_seq, alt_seq


def get_del(site, chrm, padding, alt_from, cons):
    """
    Takes a genomic position, chromosome number, padding value, the alternative 'from' sequence and
     the transcript consequence in order to deduce the deleted sequence.
    Then returns the reference and alternative sequences (with the appropriate nucleotides deleted).
    """
    full_seq = get_hg19_seq(int(chrm), seq_from=int(site) - padding, seq_to=int(site) + padding)
    left = full_seq[:padding].lower()
    middle = full_seq[padding:padding+len(alt_from)].upper()
    right = full_seq[padding+len(alt_from):].lower()
    midright = middle + right
    ref = left + middle + right
    # Extract the deleted string from the 'consequence'
    to_del = cons[cons.index('del')+3:]
    if to_del.upper() in midright.upper():
        del_site = midright.upper().index(to_del.upper())
    # Try to reverse complement the deletion, if it's long enough and can't find the stated sequence
    elif len(to_del) >= 3 and rev_comp(to_del).upper() in midright.upper():
        del_site = midright.upper().index(rev_comp(to_del).upper())
    else:
        # Failed, deleted residues not found in sequence
        return "", ""
    deleted = left + midright[:del_site] + midright[del_site+len(to_del):]
    return ref, deleted


def get_dup(site, chrm, padding, alt_from, alt_to, cons):
    """
    Takes a genomic position, chromosome number, padding value, the alternative 'from' and 'to' sequences and
     the transcript consequence in order to deduce/verify the duplicated sequence.
    Then returns the reference and alternative sequences (with the appropriate nucleotides duplicated).
    """
    full_seq = get_hg19_seq(int(chrm), seq_from=int(site) - padding, seq_to=int(site) + padding)
    left = full_seq[:padding].lower()
    middle = full_seq[padding:padding+len(alt_from)].upper()
    right = full_seq[padding+len(alt_from):].lower()
    midright = middle + right
    # Check the 'to' variant lies at the right position
    if midright.upper().index(alt_to.upper()) != 0:
        return "", ""
    ref = left + middle + right
    # Extract the deleted string from the 'consequence'
    to_dup = cons[cons.index('dup')+3:]
    # Check that the stated duplicated sequence matches what the consequence said, and is in the right place
    if alt_from + to_dup != alt_to or right.upper().index(to_dup.upper()) != 0:
        return "", ""
    duplicated = left + middle + to_dup + right
    return ref, duplicated


if __name__ == '__main__':

    # Set up necessary data structures
    ids = coll.defaultdict(str)
    chromosomes = coll.defaultdict(int)
    snps = coll.defaultdict(double_nest)
    path_to_data = 'Raw_Files/'

    # Define pad distance, i.e. how many nt to take on either side of the variant position
    pad = 12

    # Read in gene level data
    with open("TCR_Gene_Details.csv") as inf:
        for line in inf:
            # print line
            if "Symbol" not in line:
                bits = line.rstrip().split(",")
                gene = bits[0]
                ens_id = bits[1]
                chromosome = bits[2]
                ids[ens_id] = gene
                chromosomes[gene] = chromosome

    # Get exac data
    exac_files = [path_to_data+x for x in os.listdir(path_to_data) \
                  if x.startswith("exac_") == True and x.endswith(".csv") == True][::-1]

    # Read in the SNP data (from ExAC files harvested from gnomAD), give every row a number,
    for fl in exac_files:
        ens_id = fl.split("_")[1]
        with open(fl) as inf:
            cnt = 1
            for line in inf:
                if "Chrom" in line:
                    key = line
                else:
                    bits = line.replace("\"", "").rstrip().split(",")
                    snps[ids[ens_id]][tidy(cnt)]['position'] = int(bits[1])
                    snps[ids[ens_id]][tidy(cnt)]['chr'] = int(bits[0])
                    snps[ids[ens_id]][tidy(cnt)]['rsid'] = bits[2]
                    snps[ids[ens_id]][tidy(cnt)]['from'] = bits[3]
                    snps[ids[ens_id]][tidy(cnt)]['to'] = bits[4]
                    snps[ids[ens_id]][tidy(cnt)]['t_cons'] = bits[8]

                    # Get SNP, if we're dealing with a single nt change (that's a difference base!)
                    if len(snps[ids[ens_id]][tidy(cnt)]['from']) == 1 and len(snps[ids[ens_id]][tidy(cnt)]['to']) == 1 \
                            and 'dup' not in snps[ids[ens_id]][tidy(cnt)]['t_cons'] \
                            and 'del' not in snps[ids[ens_id]][tidy(cnt)]['t_cons'] \
                            and snps[ids[ens_id]][tidy(cnt)]['from'] != snps[ids[ens_id]][tidy(cnt)]['to']:
                        snps[ids[ens_id]][tidy(cnt)]['ref'], snps[ids[ens_id]][tidy(cnt)]['alt'] = \
                            get_snp(snps[ids[ens_id]][tidy(cnt)]['position'], snps[ids[ens_id]][tidy(cnt)]['chr'], pad,
                                    snps[ids[ens_id]][tidy(cnt)]['from'], snps[ids[ens_id]][tidy(cnt)]['to'])

                    # Otherwise look for deletions...
                    elif 'del' in snps[ids[ens_id]][tidy(cnt)]['t_cons']:
                        snps[ids[ens_id]][tidy(cnt)]['ref'], snps[ids[ens_id]][tidy(cnt)]['alt'] = \
                            get_del(snps[ids[ens_id]][tidy(cnt)]['position'], snps[ids[ens_id]][tidy(cnt)]['chr'], pad,
                                snps[ids[ens_id]][tidy(cnt)]['from'], snps[ids[ens_id]][tidy(cnt)]['t_cons'])

                    # ... or duplications
                    elif 'dup' in bits[8]:
                        snps[ids[ens_id]][tidy(cnt)]['ref'], snps[ids[ens_id]][tidy(cnt)]['alt'] = \
                             get_dup(snps[ids[ens_id]][tidy(cnt)]['position'], snps[ids[ens_id]][tidy(cnt)]['chr'], pad,
                                 snps[ids[ens_id]][tidy(cnt)]['from'], snps[ids[ens_id]][tidy(cnt)]['to'],
                                 snps[ids[ens_id]][tidy(cnt)]['t_cons'])

                    # Otherwise I don't think we can hazard an accurate guess of what it should be, so ignore
                    else:
                        snps[ids[ens_id]][tidy(cnt)]['ref'], snps[ids[ens_id]][tidy(cnt)]['alt'] = '', ''

                    snps[ids[ens_id]][tidy(cnt)]['type'] = bits[9]
                    snps[ids[ens_id]][tidy(cnt)]['prot_consequence'] = bits[6]
                    snps[ids[ens_id]][tidy(cnt)]['freq'] = float(bits[14])
                    cnt += 1

    # Write the data out to gene-specific files
    with open('TRAV.snps', 'w') as av_file, open('TRBV.snps', 'w') as bv_file, \
            open('TRAj.snps', 'w') as aj_file, open('TRBJ.snps', 'w') as bj_file:
        for g in snps:
            print g
            # for s in range(len(snps[g])):
            # if snps[g][s+1]['snp_from'] and snps[g][s+1]['position'] != []:
            #   both = get_snp(snps[g][s+1]['position'], snps[g][s+1]['chr'], 10, snps[g][s+1]['snp_to'])
            #   outline = ','.join([g, str(s+1), str(snps[g][s+1]['position']), \
            #     both[0], both[1], \
            #     snps[g][s+1]['type'], str(snps[g][s+1]['freq']), snps[g][s+1]['rsid'], ]) + "\n"
            #   if 'TRAV' in g or 'TRDV' in g:
            #     a_file.write(outline)
            #   elif 'TRBV' in g:
            #     b_file.write(outline)
            #   else:
            #     print "Error: missing chain info?"

    # TODO finish converting the output code from v2 to v3
