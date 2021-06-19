"""
This module contains all functions necessary for generating your own
BLOSUM-inspired subsitution matrix for copy number profile alignments
"""
import pandas as pd
import math
import json

from .format import Profile
from .utils import get_chrom_order
from .CGHcall import convert_CGHcall


def get_sequence_dataframe(data, bin_size):
    """Converts dictionary of segmented profile dataframes into a dataframe of
    samples (columns) bins (indices) and values (G, N, L) for all profiles.

    Args:
        data (dictionary): Dictionary of keys (sample names).
        and values (Clonality segment datataframes).

    Returns:
        pd.DataFrame with columns (sample_names), rows (bins)
        and values (G, N, L).
    """
    # Create dictionary of patient + sample (keys) and
    # states per bin (values)
    all_seqs = {}
    # For all patients, find unique samples
    for exp in data:
        samples = data[exp]['ID'].unique()
        # Subset all unique samples and convert to sequences
        for sample in samples:
            df = data[exp][data[exp]['ID'] == sample]
            S = Profile(sample, bin_size, df)
            arms = S.get_arms()
            # Format CNA information correctly
            temp = {}
            for arm in arms:
                for i in arm.bins:
                    temp[i.chrom+'_'+str(i.start)+'_'+str(i.end)] = i.state
            all_seqs[exp+'_'+sample] = temp
    return pd.DataFrame(all_seqs)


def blosum(dictionary, bin_size, states):
    """Creates a BLOSUM-like subsittution matrix based on
    copy number alterations instead of nucleotides.
    Based upon: https://tinyurl.com/x5tcksxf

    Args:
        data (dictionary): Dictionary with sample names.
        (keys) and segmented copy number profile dataframes (values).
        bin_size (int): Size of bins
    """
    # Convert to dataframe representation
    profiles = get_sequence_dataframe(dictionary, bin_size).dropna()
    # Re-convert to dictionary format
    pattern_dict = profiles.to_dict(orient='index')

    # For every position, find and count all possible
    # combinations of two samples (add one to prevent zero-division errors in
    # case of unobserved combinations)
    pair_counts = {}
    for i in states:
        for j in states:
            pair_counts[i+j] = 1

    for pos in pattern_dict:
        for sample1 in pattern_dict[pos]:
            for sample2 in pattern_dict[pos]:
                if sample1 != sample2:
                    pair_counts[pattern_dict[pos][sample1] +
                                pattern_dict[pos][sample2]] += 1

    # Calculate observed probability of occurrence pairs
    d = len(pattern_dict[pos])
    w = len(pattern_dict)
    n = ((w*d) * (d - 1)) / 2
    pair_probabilities = {k: v/n for k, v in pair_counts.items()}

    # Calculate probability of occurrence for every alteration
    cnv_probabilities = {}
    for i in states:
        cnv_probabilities[i] = 0

    for i in cnv_probabilities:
        numer = 0
        for j in pair_counts:
            if j.count(i) == 2:
                numer += pair_counts[j]
            elif j.count(i) == 1:
                numer += pair_counts[j] / 2

        # Divide by two because we're counting both (e.g. NG/GN, NL/LN)
        cnv_probabilities[i] = (numer / n) / 2

    # Calculate expected probability of ocurrence pairs
    exp_pair_probabilities = {}
    for cnv1 in cnv_probabilities:
        for cnv2 in cnv_probabilities:
            if cnv1 == cnv2:
                exp_pair_probabilities[cnv1+cnv2] = \
                    cnv_probabilities[cnv1] * cnv_probabilities[cnv2]
            else:
                exp_pair_probabilities[cnv1+cnv2] = \
                    2 * cnv_probabilities[cnv1] * cnv_probabilities[cnv2]

    # Create final matrix
    blosum_matrix = {}
    for i in pair_probabilities:
        blosum_matrix[i] = 2 * math.log((pair_probabilities[i] /
                                         exp_pair_probabilities[i]), 2)

    # Format final matrix as nested dict
    nested_dict = {}
    for i in states:
        nested_dict[i] = {}
        for j in states:
            nested_dict[i][j] = blosum_matrix[i+j]
    return nested_dict


def create_blosum_matrices(paths, bin_size, per_arm=False, outfile=None,
                           CGHcall=False, split=False):
    """
    Args:
        paths (dictionary): Dictionary with sample names.
        (keys) and paths (values).
        bin_size (int): Size of bins.
        per_arm (bool): Calculate BLOSUM per arm. If false,
        same matrix is returned for every arm based on
        genome-wide results.

    Optional:
        outfile (str): Full path for output file (including .json)
    """
    # Load all files
    segments = {}
    for sample in paths:
        segments[sample] = pd.read_csv(paths[sample])
        if convert_CGHcall:
            segments[sample] = convert_CGHcall(segments[sample], bin_size)
    order = get_chrom_order(split)
    matrices = {}
    matrix = blosum(segments, bin_size, states=['A', 'G', 'N', 'L', 'D'])
    # Insert same matrix for every arm
    for i in order:
        matrices[i] = matrix
    # Save matrices
    if outfile is not None:
        with open(outfile, 'w') as json_file:
            json.dump(matrices, json_file)
    return matrices
