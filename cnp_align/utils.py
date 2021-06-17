"""
This module contains a variety of utility functions that do not belong
in one particular module.
"""
import pandas as pd
from statistics import mean, median


def get_chrom_order():
    """Returns ordered list of chromosomes in format:
            ['chr1p', 'chr1q', ... 'chr22q']
       Excluding sex chromosomes!
    """
    order = []
    for i in range(1, 23):
        if len(str(i)) == 1:
            order.append('chr0'+str(i)+'p')
            order.append('chr0'+str(i)+'q')
        else:
            order.append('chr'+str(i)+'p')
            order.append('chr'+str(i)+'q')
    return order


def get_chrom_proportions(chrom):
    """Returns the proportion of every chromosome arms contribution to genome
    size."""
    total = 0
    for i in chrom:
        for j in chrom[i]['seq1']:
            if j != '-':
                total += 1
    prop = []
    for i in chrom:
        subtotal = 0
        for j in chrom[i]['seq1']:
            if j != '-':
                subtotal += 1
        prop.append(subtotal/total)
    return prop


def find_missing_chrom(dictionary, verbose=True):
    """Looks for missing chroms in a dictionary and prints the missing values."""
    order = get_chrom_order()
    missing = [i for i in order if i not in dictionary.keys()]
    if verbose:
        print('Detected', len(missing), 'missing arms')
        print('Missing arms: ' + ' '.join(missing))
    return missing


def get_chrom_sizes(df):
    """Returns dictionary of chromosome information used to
    re-create accurate bins."""
    chrom = {}
    arms = df['chrom'].unique()
    for i in arms:
        chrom[i] = {}
        chrom[i]['start'] = df[df['chrom'] == i]['loc.start'].min()
        chrom[i]['end'] = df[df['chrom'] == i]['loc.end'].max()
    return chrom


def get_state_ranges(states):
    """From a list of characters, finds all regions with similar values.
    Returns regions with respective values as:
        {(start, end):value, (start, end):value}
    """
    ranges = {}
    values = []
    for i in range(len(states)):
        if i == 0:
            # Start position
            values.append(i)
        elif i == len(states)-1:
            # End position, stop range
            values.append(i)
            ranges[(min(values), i)] = states[i-1]
        elif states[i] == states[i-1]:
            # Extend range
            values.append(i)
            continue
        elif states[i] != states[i-1]:
            # Stop existing range and begin new range
            if len(values) == 0:
                ranges[i-1, i] = states[i-1]
            else:
                ranges[(min(values), i)] = states[i-1]
            values = []
    return ranges


def autoresolve(dataframe, bin_size):
    """Fills in missing areas in segment dataframe by
    intra and/or extrapoloation. Fills in left value if only one
    bin size is missing.

    Example (numbers are missing bins):
        NNNNNN-1-2-3-NNNNN. -> 1:'N', 2:'N', 3:'N'
        NNNNNN-1-2-3-GGGGG. -> 1:'N', 2:'N', 3:'G'
        NNNNNN-1-2-3-4-GGGGG. -> 1:'N', 2:'N', 3:'G', 4:'G'
        NNNNNN-1-GGGGG. -> 1:'N'

    Args:
        df (pd.DataFrame): Pandas Dataframe of Clonality output. Should contain
        the columns 'chrom', 'loc.start', 'loc.end' and 'state'.
        bin_size (int): Size of segment bins
    """
    df = dataframe.copy()
    dfs = []
    # Loop over all chromosomes
    for chrom in get_chrom_order():
        subset = df[df['chrom'] == chrom]
        cnvs = subset.to_dict('records')
        new_entries = []

        # Loop over all segments
        for i in enumerate(cnvs):
            if i[0] == 0:
                continue

            # Evaluate whether there is a gap between the segments
            elif cnvs[i[0]]['loc.start'] != cnvs[i[0]-1]['loc.end'] + bin_size:

                # Calculate size difference
                size = cnvs[i[0]]['loc.start'] - cnvs[i[0]-1]['loc.end']
                n_bins = (size - 2 * bin_size) / bin_size

                if n_bins == 0:
                    # The missing field cannot be inserted because
                    # there is a gap of exactly 2 bins. Because
                    # the start position is equal to the last position + bin
                    # and the end the next start - bin this is not possible.
                    # Extend previous one by one bin
                    cnvs[i[0]-1]['loc.end'] += bin_size
                    cnvs[i[0]-1]['fill'] = 'modified'
                    continue

                if n_bins == 1:
                    # Fill in as last observed
                    entry = {'ID': cnvs[i[0]]['ID'],
                             'chrom': cnvs[i[0]]['chrom'],
                             'loc.start': cnvs[i[0]-1]['loc.end'] + bin_size,
                             'loc.end': cnvs[i[0]]['loc.start'] - bin_size,
                             'num.mark': None,
                             'seg.mean': None,
                             'state': cnvs[i[0]-1]['state'],
                             'fill': 'previous'}
                    new_entries.append(entry)
                    continue

                # Uneven number of bins (fill in as left)
                if n_bins % 2 != 0:
                    first = (n_bins + 1) / 2
                else:
                    first = n_bins / 2

                # First half
                entry = {'ID': cnvs[i[0]]['ID'],
                         'chrom': cnvs[i[0]]['chrom'],
                         'loc.start': cnvs[i[0]-1]['loc.end'] + bin_size,
                         'loc.end': cnvs[i[0]-1]['loc.end'] + first * bin_size,
                         'num.mark': None,
                         'seg.mean': None,
                         'state': cnvs[i[0]-1]['state'],
                         'fill': 'first_half'}
                new_entries.append(entry)
                # Second half
                entry = {'ID': cnvs[i[0]]['ID'],
                         'chrom': cnvs[i[0]]['chrom'],
                         'loc.start': cnvs[i[0]-1]['loc.end'] +
                         first * bin_size + bin_size,
                         'loc.end': cnvs[i[0]]['loc.start'] - bin_size,
                         'num.mark': None,
                         'seg.mean': None,
                         'state': cnvs[i[0]]['state'],
                         'fill': 'second_half'}
                new_entries.append(entry)
        # Store modified segments
        if len(new_entries) > 0:
            cnvs.extend(new_entries)
        dfs.append(pd.DataFrame.from_records(cnvs))
    # Re-generate segment dataframe
    df = pd.concat(dfs)
    df = df.astype({"loc.start": int, "loc.end": int})
    return df


def format_alignment_results(alignment, null_scores=None):
    """Returns a dataframe with rows (arms) and
    columns (alignment features).
    """
    all_scores = {'arm': [], 'score': [], 'adjusted_score': []}
    if null_scores is not None:
        all_scores['match_proba'] = []
        all_scores['mismatch_proba'] = []
    # Retrieve scores for every arm
    for arm in alignment:
        all_scores['arm'].append(arm)
        all_scores['score'].append(alignment[arm]['score'])
        all_scores['adjusted_score'].append(alignment[arm]['adjusted_score'])
        if null_scores is not None:
            all_scores['match_proba'].append(alignment[arm]['match_proba'])
            all_scores['mismatch_proba'].append(alignment[arm]['mismatch_proba'])
    return pd.DataFrame(all_scores)


def summarize_results(dfs, exp_id, null_scores=None):
    """For a dictionary of dataframes created by format_alignment_results,
    summarize the statistics using the mean, median and if possible probability
    cut-offs."""
    results = []
    for i in dfs:
        regular_scores = dfs[i]['score']
        adjusted_scores = dfs[i]['adjusted_score']
        result = {'id': exp_id,
                  'pair': i,
                  'total_score': sum(regular_scores),
                  'mean_score': mean(regular_scores),
                  'median_score': median(regular_scores),
                  'total_adjusted_score': sum(adjusted_scores),
                  'mean_adjusted_score': mean(adjusted_scores),
                  'median_adjusted_score': median(adjusted_scores)}
        if null_scores is not None:
            match_probas = dfs[i]['match_proba']
            mismatch_probas = dfs[i]['mismatch_proba']
            result['match_proba_5'] = len([i for i in match_probas if i <= 0.05])
            result['match_proba_10'] = len([i for i in match_probas if i <= 0.1])
            result['match_proba_20'] = len([i for i in match_probas if i <= 0.2])
            result['mismatch_proba_60'] = len([i for i in mismatch_probas if i <= 0.6])
            result['mismatch_proba_50'] = len([i for i in mismatch_probas if i <= 0.5])
            result['mismatch_proba_40'] = len([i for i in mismatch_probas if i <= 0.4])
        results.append(result)
    return pd.DataFrame(results)


def format_subst_matrix(subst_matrix):
    output = {}
    for i in subst_matrix:
        for j in subst_matrix[i]:
            if (i, j) not in output:
                output[(i, j)] = subst_matrix[i][j]
    return output
