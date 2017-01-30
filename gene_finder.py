# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    nucArray = "ACGT"
    reverseIndex = nucArray.find(nucleotide)
    return nucArray[len(nucArray)-reverseIndex-1]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverseComplement = ""

    for i in range(len(dna)):
        reverseComplement += get_complement(dna[i])
    return reverseComplement[::-1]

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stopCodons = ['TAG', 'TAA', 'TGA']
    for i in stopCodons:
        for j in range(len(dna)//3):
            jIndex = 3*j
            if dna[jIndex] == i[0] and dna[jIndex+1] == i[1] and dna[jIndex+2] == i[2]:
                return dna[:jIndex]

    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    orf = []
    i = 0
    while i < (len(dna)//3):
            index = 3*i
            currentCodon = [dna[index], dna[index+1], dna[index+2]]
            if currentCodon[0] == 'A' and currentCodon[1] =='T' and currentCodon[2] == 'G':
                foundOrf = rest_of_ORF(dna[index:])
                orf.append(foundOrf)
                i += len(foundOrf)//3
            i += 1
    return orf


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    orfs = []
    for i in range(3):
        res = find_all_ORFs_oneframe(dna[i:])
        for i in res:
            orfs.append(i)
    return orfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allOrfs = []
    revComp = get_reverse_complement(dna)
    orfsRev = find_all_ORFs(revComp)
    orfsFor = find_all_ORFs(dna)

    for i in orfsFor:
        allOrfs.append(i)
    for i in orfsRev:
        allOrfs.append(i)

    return allOrfs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence

        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
"""
print(get_complement("T"))
print(get_reverse_complement("AACCGGTT"))
print(rest_of_ORF("ATGAGATAGG"))
print(find_all_ORFs_oneframe("ATGCATGAATGTTAGATAGATGTGCCC"))
print(find_all_ORFs("ATGCATGAATGTAG"))
print(find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA"))
"""
