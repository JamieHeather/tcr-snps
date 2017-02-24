#!/usr/bin/env python

"""
snp-functions.py: The user-defined functions to accompany 
https://github.com/JamieHeather/tcr-snps
"""

import collections as coll
import string
import urllib2
import Levenshtein as lev

__version__ = '1.1.2'


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
        return '', ''
    alt_seq = full_seq[:padding].lower() + alt.upper() + full_seq[-padding:].lower()
    return ref_seq, alt_seq


def get_del(site, chrm, padding, alt_from, cons):
    """
    Takes a genomic position, chromosome number, padding value, the alternative 'from' sequence and
     the transcript consequence in order to deduce the deleted sequence.
    Then returns the reference and alternative sequences (with the appropriate nucleotides deleted).
    """
    # Extract the deleted string from the 'consequence'
    to_del = cons[cons.index('del')+3:]
    full_seq = get_hg19_seq(int(chrm), seq_from=int(site) - padding, seq_to=int(site) + padding + len(to_del))
    left = full_seq[:padding].lower()
    middle = full_seq[padding:padding+len(alt_from)].upper()
    right = full_seq[padding+len(alt_from):].lower()
    midright = middle + right
    ref = left + middle + right
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
    # Extract the duplicated string from the 'consequence'
    to_dup = cons[cons.index('dup')+3:]
    full_seq = get_hg19_seq(int(chrm), seq_from=int(site) - padding, seq_to=int(site) + padding + len(to_dup))
    left = full_seq[:padding].lower()
    middle = full_seq[padding:padding+len(alt_from)].upper()
    right = full_seq[padding+len(alt_from):].lower()
    midright = middle + right
    # Check the 'to' variant lies at the right position
    if midright.upper().index(alt_to.upper()) != 0:
        return "", ""
    ref = left + middle + right
    # Check duplicated string is in string; if not, try reverse complement
    if to_dup.upper() not in midright.upper():
        to_dup = rev_comp(to_dup)
        if to_dup.upper() not in midright.upper():
            return '', ''
    # Check that the stated duplicated sequence matches what the consequence said, and is in the right place
    if alt_from + to_dup != alt_to or right.upper().index(to_dup.upper()) != 0:
        return "", ""
    duplicated = left + middle + to_dup + right
    return ref[:-len(to_dup)], duplicated[:-len(to_dup)]


def readfq(fp):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def trim_sequences(ref, alt):
    """
    Given two sequences (one reference and one alternative allele) of different length, return trimmed to align.
    Uses the Levenshtein distance to determine which end to truncate
    """
    if len(ref) > len(alt):
        len_diff = len(ref) - len(alt)
        return pick_truncation(ref, alt, len_diff), alt
    #
    else:
        len_diff = len(alt) - len(ref)
        return ref, pick_truncation(alt, ref, len_diff)


def pick_truncation(long_seq, short_seq, length):
    """
    Given two sequencese (from trim_sequences), truncate the longer at both the 5' or 3' ends
    and return the sequence that matches shorter best
    """
    truncate_3 = long_seq[:-length]
    truncate_5 = long_seq[length:]
    distance_3 = lev.distance(truncate_3, short_seq)
    distance_5 = lev.distance(truncate_5, short_seq)
    if distance_3 < distance_5:
        return truncate_3
    elif distance_5 < distance_3:
        return truncate_5
    else:
        # Cannot determine which end to chop from, i.e. each truncation is equally similar to short sequence (unlikely!)
        return ""


def capitalise_position(seq, position):
    """
    Given an string, return said string with the character at a certain position capitalised.
    Note that it is blind to the pre-existing case (so that it can be applied multiple times to same string)
    """
    return seq[:position] + seq[position].upper() + seq[position+1:]


def multiple_capitalise(seq, list_positions):
    """
    Given a sequence and a list of positions, return the string having capitalised all the letters at those positions
    """
    outstr = seq
    for pos in list_positions:
        outstr = capitalise_position(outstr, pos)
    return outstr
