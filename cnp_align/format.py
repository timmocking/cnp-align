"""
This file contains all classes necessary for storing and processing segmented
copy number profiles
"""
from .utils import get_chrom_sizes


class Bin():
    """Class for storing all bin-level segment information.

    Atrributes:
        chrom (str): Chromosome (arm).
        start (int): Start position of bin.
        end (int): End position of bin.
        state (str): Copy number state (either 'G', 'N' or 'L').
    """

    def __init__(self, chrom, start, end, state):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.state = state


class Arm():
    """Class for storing all chromosome arm level segment information.

    Attributes:
        chrom (str): Chromosome (arm).
        chrom_sizes (dict): Nested dictionary with chrom info.
        bin_size (int): Size of segment bins
        data (Pd.DataFrame): Pandas Dataframe of Clonality subsetted
        for chromosme arm and patient.
        bin_size

    Methods:
        get_sequence: Return string sequence representation of aberrations for
        all chromosome bins (e.g. NNNNGNN).
    """

    def __init__(self, chrom, bin_size, chrom_sizes, data):
        self.chrom = chrom
        self.data = data
        self.bins = []
        # Re-create bins of correct size
        rng = range(chrom_sizes[chrom]['start'], chrom_sizes[chrom]['end'],
                    bin_size)
        cnvs = data.to_dict('records')
        # Find states for every bin
        unassigned = []
        for i in rng:
            assigned = False
            for j in cnvs:
                if i >= j['loc.start'] and i <= j['loc.end']:
                    self.bins.append(Bin(self.chrom, i, i+bin_size,
                                         j['state'][0]))
                    assigned = True
            # Store unassigned bins
            if not assigned:
                unassigned.append(self.chrom+'_'+str(i))
        # Print statistics on unassigned bins
        if len(unassigned) != 0:
            print('Warning: detected', len(unassigned), 'missing bins on '
                  + chrom)

    def get_sequence(self):
        """Returns string sequence of alterations in bin.
        """
        seq = ''
        for i in self.bins:
            seq += i.state
        return seq


class Profile():
    """Class for storing all profile level segment information. Should be
    initialized with a segment dataframe for that given profile.

    Args:
        id (str): Sample id.
        bin_size (int): Size of segment bins.
        df (pd.DataFrame): Pandas Dataframe of Clonality output. Should contain
        only one patient and the columns 'chrom', 'loc.start', 'loc.end' and
        'state'.
    """

    def __init__(self, sample_id, bin_size, df):
        self.id = sample_id
        self.arms = []
        # Get chromosome sizes. Used to re-create the bins from the segments
        chroms = get_chrom_sizes(df)
        # Create Arm class for every arm
        for chrom in chroms:
            subset = df[df['chrom'] == chrom][['chrom', 'loc.start', 'loc.end',
                                               'state']]
            self.arms.append(Arm(chrom, bin_size, chroms, subset))

    def get_dict(self):
        """Returns dictionary of chromsome arms (keys) and
        CNAs as one long string (values)."""
        arms = {}
        for i in self.arms:
            arms[i.chrom] = i.get_sequence()
        return arms

    def get_arms(self):
        return self.arms
