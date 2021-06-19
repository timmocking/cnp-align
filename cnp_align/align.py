"""
This module contains the Alignment class which is used for performing and
interacting with CNP-ALIGN's "Alignment" class objects.
"""
from Bio import pairwise2

from .plot import plot_alignment
from .utils import get_chrom_order, find_missing_chrom, format_subst_matrix


class Alignment():
    """Class for performing and visualizing copy-number alignment Biopython's
    pairwise2 global sequence alignment algorithm.

    Args:
        profile1 (Sample): First Sample class object.
        profile2 (Sample): Second Sample class object.
    """

    def __init__(self, profile1, profile2):
        self.profile1 = profile1
        self.profile2 = profile2
        self.alignment = None

    def align(self, sub_matrix, gap_open=10000, gap_extend=10000,
              null_scores=None):
        """Performs copy-number alignment using Biopython's pairwise2 global
        sequence alignment and returns various alignment metrics.

        Note: By default opening and extension gap penalties are set at 10000
        in order to enforce a gapless alignment. Perform gapped alignments
        at own risk!

        Args:
            sub_matrix (dict): Nested dictionary representing
            substitution values for every possible sequence pair.

        Optional:
            gap_open (float): Opening gap penalty. Default=10000
            gap_extend (float): Extending gap penalty. Default=10000
            null_scores (dict): Dictionary containing list of alignment
            scores for a given population.

        Returns:
            dict: Nested dict of keys (arms) and alignment
            output (values). Output is formatted as follows:
                seq1 (str): aligned sequence 1 (incl. gaps).
                seq2 (str): aligned sequence 2 (incl. gaps).
                score (float): Alignment score.
                adjusted_score (float): Alignment score divided by seq. length.
                seq1_gaps (int): Number of introduced gaps in seq1.
                seq2_gaps (int): Number of introduced gaps in seq2.

            Optional output:
                match_proba (float): Proportion of null_scores
                scoring equal or better than adjusted_score
                mismatch_proba (float): Proportion of null_scores
                scoring worse than adjusted_score
        """
        # Get sequences (string format) from Sample class
        profile1_dict = self.profile1.get_dict()
        profile2_dict = self.profile2.get_dict()
        # Check if both profiles have same amount of chromosome arms
        assert len(profile1_dict) == len(profile2_dict)
        # Perform alignment per chromosome arm
        # Return alignments as nested dictionary with keys (arms) and alignment
        # output (values)
        alignments = {}
        for arm in profile1_dict:
            # Convert substitution matrix to Biopython format
            formatted_matrix = format_subst_matrix(sub_matrix[arm])
            # Perform alignment
            result = pairwise2.align.globalds(profile1_dict[arm],
                                              profile2_dict[arm],
                                              formatted_matrix,
                                              gap_open,
                                              gap_extend)
            # Format output of alignment
            alignments[arm] = {'seq1': result[0].seqA,
                               'seq2': result[0].seqB,
                               'score': result[0].score,
                               'adjusted_score': result[0].score /
                               len(result[0].seqA),
                               'seq1_gaps': result[0].seqA.count('-'),
                               'seq2_gaps': result[0].seqB.count('-')}
            # Add match/mismatch probability based on panel of null-scores
            if null_scores is not None:
                bigger = len([i for i in null_scores[arm]
                              if i >= alignments[arm]['adjusted_score']])
                smaller = len([i for i in null_scores[arm]
                              if i < alignments[arm]['adjusted_score']])
                alignments[arm]['match_proba'] = bigger / len(null_scores[arm])
                alignments[arm]['mismatch_proba'] = smaller / len(null_scores[arm])
        # Store and return alingments
        self.alignment = alignments
        return self.alignment

    def plot(self, chromosomes=None, null_scores=None, match_thresh=0.1,
             mismatch_thresh=0.5, save=False, figname=None, split=True):
        """Return and/or save visualisation of profiles with
        alignment features."""
        if self.alignment is None:
            raise Exception('Please perform alignment before visualization.')

        if chromosomes is None:
            # Plot all chromosomes
            # Arrange chromosome order and remove missing arms
            chromosome_order = get_chrom_order(split)
            missing = find_missing_chrom(self.alignment, split, verbose=False)
            for i in missing:
                chromosome_order.remove(i)
        else:
            chromosome_order = chromosomes

        if split:
            pos = chromosome_order.index('chr11p')
        else:
            pos = chromosome_order.index('chr12')

        # Plot first half
        plot_alignment(self.alignment, chromosome_order[:pos], null_scores,
                       match_thresh, mismatch_thresh, save, figname+'_chr1-chr10')

        # Plot second half
        plot_alignment(self.alignment, chromosome_order[pos:], null_scores,
                       match_thresh, mismatch_thresh, save, figname+'_chr10-chr22')
