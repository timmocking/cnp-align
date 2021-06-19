"""
This module contains all functions necessary for parsing
CGHcall segment output to a format similar to Clonality's
"""
from .utils import get_chrom_order
import pandas as pd


def convert_to_char(col):
    """Converts numerical output of CGHcall
    to characters.
    """
    vals = []
    for i in col:
        if i == 0.0:
            vals.append('N')
        elif i == -1.0:
            vals.append('L')
        elif i == -2.0:
            vals.append('D')
        elif i == 1.0:
            vals.append('G')
        elif i == 2.0:
            vals.append('A')
        else:
            vals.append(i)
    return vals


def add_chr_str(col):
    """Adds chr or chr0 before integer chromosome
    names.
    """
    new_values = []
    for i in col:
        if len(str(i)) == 1:
            new_values.append('chr0'+str(i))
        else:
            new_values.append('chr'+str(i))
    return new_values


def autoresolve_CGHcall(dataframe, bin_size):
    """Modified version of blosum.autoresolve.
    Will be refractored in future version!
    """
    df = dataframe.copy()
    dfs = []
    # Loop over all chromosomes
    for chrom in get_chrom_order(split=False):
        subset = df[df['chrom'] == chrom]
        cnvs = subset.to_dict('records')
        new_entries = []

        # Loop over all segments
        for i in enumerate(cnvs):
            if i[0] == 0:
                continue

            # Evaluate whether there is a gap between the segments
            elif cnvs[i[0]]['loc.start'] + 1 != cnvs[i[0]-1]['loc.end']:

                # Calculate size difference
                size = cnvs[i[0]]['loc.start'] - 1 - cnvs[i[0]-1]['loc.end']
                n_bins = size / bin_size
                
                if n_bins == 1:
                    # Fill in as last observed
                    entry = {'ID': cnvs[i[0]]['ID'],
                             'chrom': cnvs[i[0]]['chrom'],
                             'loc.start': cnvs[i[0]-1]['loc.end'] + 1,
                             'loc.end': cnvs[i[0]-1]['loc.end'] + bin_size,
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
                         'loc.start': cnvs[i[0]-1]['loc.end'] + 1,
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
                         first * bin_size + 1,
                         'loc.end': cnvs[i[0]]['loc.start'] - 1,
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


def convert_CGHcall(df, bin_size, autoresolve=True):
    """
    Converts CGHcall format to segment dataframe
    """
    # Name index
    df = df.rename(columns={'Unnamed: 0':'bin'})

    # Convert CGHcall format (i.e. -2 is double loss) to characters
    for col in df.columns:
        df[col] = convert_to_char(df[col])
        
    # Pre-process for collapse
    bins = df['bin']
    df['chrom'] = [i.split(':')[0] for i in bins]
    df['loc.start'] = [int(i.split(':')[1].split('-')[0]) for i in bins]
    df['loc.end'] = [int(i.split(':')[1].split('-')[1]) for i in bins]

    # Add 'chr' or 'chr0' before chrom
    df['chrom'] = add_chr_str(df['chrom'])
    
    # Collapse row-based format into segment based format
    samples = list(df.columns[1:-3])
    segments = []
    for chrom in df['chrom'].unique():
        chrom_subset = df[df['chrom']==chrom]
        
        tracker_dict = {}
        for i in samples:
            tracker_dict[i] = {}
        
        cnvs = chrom_subset.to_dict('records')
        for sample in samples:
            for i in range(len(cnvs)):
                # Initialize         
                if i == 0:
                    temp_dict = {'ID':sample,
                                'chrom':chrom,
                                'loc.start': cnvs[i]['loc.start'],
                                'loc.end': cnvs[i]['loc.end'],
                                'state':cnvs[i][sample]}
                    tracker_dict[sample] = temp_dict
                
                elif tracker_dict[sample]['loc.end'] + 1 != cnvs[i]['loc.start']:
                    # Missing bin identified (start new one)
                    segments.append(tracker_dict[sample])

                    temp_dict = {'ID':sample,
                    'chrom':chrom,
                    'loc.start': cnvs[i]['loc.start'],
                    'loc.end': cnvs[i]['loc.end'],
                    'state':cnvs[i][sample]}
                    tracker_dict[sample] = temp_dict
                
                elif i == len(cnvs)- 1:
                    # Reached end 
                    segments.append(tracker_dict[sample])

                # Check if current value equal to last one
                elif tracker_dict[sample]['state'] == cnvs[i][sample]:
                    # Increase loc.end
                    tracker_dict[sample]['loc.end'] = cnvs[i]['loc.end']
                elif cnvs[i][sample] != tracker_dict[sample]['state']:
                    # Create new segment input
                    segments.append(tracker_dict[sample])

                    # Create new segment
                    temp_dict = {'ID':sample,
                    'chrom':chrom,
                    'loc.start': cnvs[i]['loc.start'],
                    'loc.end': cnvs[i]['loc.end'],
                    'state':cnvs[i][sample]}
                    tracker_dict[sample] = temp_dict
    segments = pd.DataFrame(segments)

    if autoresolve:
        segments = autoresolve_CGHcall(segments, bin_size)

    return segments