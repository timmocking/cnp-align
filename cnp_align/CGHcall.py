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


def convert_CGHcall_probs(df, bin_size, autoresolve=True):
    """
    Converts CGHcall format to segment dataframe
    """
    # Name index
    df = df.rename(columns={'Unnamed: 0':'bin'})
        
    # Pre-process for collapse
    bins = df['bin']
    df['chrom'] = [i.split(':')[0] for i in bins]
    df['loc.start'] = [int(i.split(':')[1].split('-')[0]) for i in bins]
    df['loc.end'] = [int(i.split(':')[1].split('-')[1]) for i in bins]

    # Add 'chr' or 'chr0' before chrom
    df['chrom'] = add_chr_str(df['chrom'])
    
    # Add a new column 'sample' designating the sample
    samples = []
    for i in df.columns:
        if i.startswith('PROBAMP'):
            sample = i.split('_')[1:]
            if i[8:] not in samples:
                samples.append(i[8:])
    dfs = []
    for sample in samples:
        cols = ['bin', 'chrom', 'loc.start', 'loc.end']
        for column in df.columns:
            if sample in column:
                cols.append(column)
        # Subset
        subset = df[cols]
        subset['sample'] = sample
        # Rename columns
        proba_cols = ['CALL', 'PROBAMP', 'PROBGAIN', 'PROBNORM', 'PROBLOSS', 'PROBDEL']
        for i in proba_cols:
            subset = subset.rename(columns={i+'_'+sample: i})
        dfs.append(subset)
    df = pd.concat(dfs, axis=0)

    # df = df[df['chrom']=='chr01']
    # df = df[df['loc.end']<= 1700001]    
    # df = df[df['sample']=='FR14066719_S74_R1_001']

    samples = df['sample'].unique()
    
    segments = []
    for chrom in df['chrom'].unique():
        chrom_subset = df[df['chrom']==chrom]
        
        tracker_dict = {}
        for i in samples:
            tracker_dict[i] = {}

        for sample in samples:
            cnvs = chrom_subset[chrom_subset['sample']==sample].to_dict('records')

            for i in range(len(cnvs)):
                # print(sample)
                # print('pos:')
                # print(cnvs[i])
                # print('tracker:')
                # print(segments)
                # print('-------------------------------------------------')

                # Initialize         
                if i == 0:
                    
                    ######
                    temp_dict = {'ID':sample,
                                'chrom':chrom,
                                'loc.start': cnvs[i]['loc.start'],
                                'loc.end': cnvs[i]['loc.end']}
                    for proba_col in proba_cols:
                        temp_dict[proba_col] = cnvs[i][proba_col]
                    tracker_dict[sample] = temp_dict
                    
                    ######
            
                elif tracker_dict[sample]['loc.end'] + 1 != cnvs[i]['loc.start']:
                    # Missing bin identified (start new one)
                    segments.append(tracker_dict[sample])
                    ######
                    temp_dict = {'ID':sample,
                                'chrom':chrom,
                                'loc.start': cnvs[i]['loc.start'],
                                'loc.end': cnvs[i]['loc.end']}
                    for proba_col in proba_cols:
                        temp_dict[proba_col] = cnvs[i][proba_col]
                    tracker_dict[sample] = temp_dict
                    ######

                elif i == len(cnvs) - 1:

                    # Reached end 
                    segments.append(tracker_dict[sample])
                    
                # Check if probabilities equal to last bin
                elif check_proba_similarity(proba_cols, tracker_dict[sample], cnvs[i]):
                    # Increase loc.end
                    tracker_dict[sample]['loc.end'] = cnvs[i]['loc.end']
                
                else:
                    # new segment
                    segments.append(tracker_dict[sample])

                    ######
                    temp_dict = {'ID':sample,
                                'chrom':chrom,
                                'loc.start': cnvs[i]['loc.start'],
                                'loc.end': cnvs[i]['loc.end']}
                    for proba_col in proba_cols:
                        temp_dict[proba_col] = cnvs[i][proba_col]
                    tracker_dict[sample] = temp_dict
                    ######
    segments = pd.DataFrame(segments)

    if autoresolve:
        segments = autoresolve_CGHcall_proba(segments, bin_size, proba_cols)
    return segments


def check_proba_similarity(proba_cols, dictionary1, dictionary2):
    dict1 = {key: value for key, value in dictionary1.items() if key in proba_cols}
    dict2 = {key: value for key, value in dictionary1.items() if key in proba_cols}
    return list(dict1.values()) == list(dict2.values())


def autoresolve_CGHcall_proba(dataframe, bin_size, proba_cols):
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
                             'fill': 'previous'}
                    for proba_col in proba_cols:
                        entry[proba_col] = cnvs[i[0]-1][proba_col]
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
                         'fill': 'first_half'}
                for proba_col in proba_cols:
                        entry[proba_col] = cnvs[i[0]-1][proba_col]
                
                new_entries.append(entry)
                # Second half
                entry = {'ID': cnvs[i[0]]['ID'],
                         'chrom': cnvs[i[0]]['chrom'],
                         'loc.start': cnvs[i[0]-1]['loc.end'] +
                         first * bin_size + 1,
                         'loc.end': cnvs[i[0]]['loc.start'] - 1,
                         'num.mark': None,
                         'seg.mean': None,
                         'fill': 'second_half'}
                for proba_col in proba_cols:
                        entry[proba_col] = cnvs[i[0]][proba_col]
                new_entries.append(entry)
                 
        # Store modified segments
        if len(new_entries) > 0:
            cnvs.extend(new_entries)
        dfs.append(pd.DataFrame.from_records(cnvs))
    # Re-generate segment dataframe
    df = pd.concat(dfs)
    df = df.astype({"loc.start": int, "loc.end": int})
    return df



