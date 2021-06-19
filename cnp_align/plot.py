"""
This module contains functions for the visualisation of copy number profiles
and their CNP-ALIGN alignment features.
"""
import warnings
import matplotlib.pyplot as plt
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

from .config import Config
from .utils import get_chrom_proportions, get_state_ranges


def plot_aberration(ax, chrom):
    """Plots copy numbers as lines at levels +1 (G), 0 (N) and -1 (L)."""
    # Remove all gaps from alignment
    chrom = chrom.replace('-', '')
    y = []
    for i in range(len(chrom)):
        if chrom[i] == 'A':
            y.append(2)
        elif chrom[i] == 'G':
            y.append(1)
        elif chrom[i] == 'L':
            y.append(-1)
        elif chrom[i] == 'Z':
            y.append(-2)
        else:
            y.append(0)
    # Collapse continuous lines (plotting speed optimalization)
    lines = get_state_ranges(y)
    for i in lines:
        ax.plot(list(i), [lines[i], lines[i]], color=Config().line_col,
                lw=Config().line_width)


def find_concordance(chrom1, chrom2):
    """Returns concordance (list of colors) of two CNPs.
    Refer to config.py for setting background colors."""
    # Remove gaps from alignments
    chrom1= chrom1.replace('-', '')
    chrom2= chrom2.replace('-', '')
    conc = []
    for i in range(len(chrom1)):
        try:
            if chrom1[i] == chrom2[i]:
                if chrom1[i] == 'N':
                    conc.append(Config().normal_col)
                else:
                    conc.append(Config().match_col)
            else:
                conc.append(Config().mismatch_col)
        except KeyError:
            # Most likely uneven amount of bins (autoresolve disabled)
            conc.append(Config().error_col)
    return conc


def plot_profile_concordance(ax, conc):
    """Plots background colour based on whether two profiles share CNAs"""
    # Collapse continuous background colors (plotting speed optimalization)
    states = get_state_ranges(conc)
    for i in states:
        ax.axvspan(i[0], i[1], facecolor=states[i], alpha=0.5)


def format_box(ax):
    ax.set_ylim(0, 2)
    ax.set_xlim(0, 2)
    plt.yticks([])
    plt.xticks([])
    ax.tick_params(axis='y', length=0)


def plot_bin_box(ax, conc, box_type, plot_row_name):
    """Plots box with number of match/mismatching bins."""
    format_box(ax)
    if box_type == 'matches':
        n = conc.count(Config().match_col)
    else:
        n = conc.count(Config().mismatch_col)
    plt.text(0.5, 0.5, str(n), horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,
             fontsize=Config().font_size)
    if plot_row_name:
        plt.yticks([1], ['Number of '+box_type],
                   fontsize=Config().font_size)


def plot_score_box(ax, chrom, plot_row_name):
    """Plots box with adjusted alignment score."""
    format_box(ax)
    if plot_row_name:
        plt.yticks([1], ['Adjusted alignment score'],
                   fontsize=Config().font_size)
    # Plot score
    score = str(round(chrom['adjusted_score'], 1))
    if score in ('0.0', '-0.0'):
        score = '0'
    plt.text(0.5, 0.5, score, horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,
             fontsize=Config().font_size)


def plot_proba_box(ax, chrom, significant, box_type, plot_row_name):
    """Plots box with match/mismatching probability."""
    format_box(ax)
    if box_type == 'match':
        if significant:
            ax.axvspan(0, 2, facecolor=Config().match_col, alpha=0.5)
        if plot_row_name:
            plt.yticks([1], ['Match probability'], fontsize=Config().font_size)
    if box_type == 'mismatch':
        if significant:
            ax.axvspan(0, 2, facecolor=Config().mismatch_col, alpha=0.5)
        if plot_row_name:
            plt.yticks([1], ['Mismatch probability'], fontsize=Config().font_size)
    probability = str(round(chrom[box_type+'_proba'], 3))
    plt.text(0.5, 0.5, probability, horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes,
             fontsize=Config().font_size)


def plot_alignment(alignment, order, null_scores, match_thresh=0.1,
                   mismatch_thresh=0.1, save=False, figname=None):
    """Visualizes similarity of two copy number profiles along with a number
    of alignment-based features.

    Args:
        alignment (dict): Dictionary of alignment features.
        order (list): List (ordered) of chromosome arms to be visualized.
        null_scores (dict): Dictionary containing list of alignment
        scores for a given population.
        match_thresh (float): Threshold for designating match probability as
        significant.
        mismatch_thresh (float): Threshold for designating match probability as
        significant.
        save (bool): Save figure or not.
        figname (str): Filename for saving figure (if save=True).
    """
    fig = plt.figure(figsize=(45, 6))

    # Create Gridspec frame for visualising subplots
    n_rows = 5
    width_ratios = get_chrom_proportions(alignment)
    height_ratios = [1, 1, 0.25, 0.25, 0.25]
    if null_scores is not None:
        # Create additional boxes for visualizing null_scores and
        # match/mismatch probabilities
        n_rows = 7
        height_ratios = [1, 1, 0.25, 0.25, 0.25, 0.25, 0.25]
    gs = fig.add_gridspec(n_rows, len(alignment), width_ratios=width_ratios,
                          height_ratios=height_ratios)

    # Loop over all selected chromosomes in defined order
    for i in enumerate(order):
        #  Create subplots from top to bottom
        for j in enumerate(['seq1', 'seq2']):

            # Plot aberrations
            # Aberrations are visualized as lines with 3 levels:
            # +1 (Gain), 0 (Normal), -1 (Loss)
            ax = fig.add_subplot(gs[j[0], i[0]])
            plot_aberration(ax, alignment[i[1]][j[1]])
            # Plot name of chromosome
            if j[0] == 0:
                plt.title(i[1], rotation=45, fontsize=Config().font_size)
            # Plot sequence name
            if i[0] == 0:
                plt.yticks([0], [j[1]], fontsize=Config().font_size)
            else:
                plt.yticks([])
            # Remove superfluous ticks and labels
            plt.xticks([])
            ax.tick_params(axis='y', length=0)
            plt.grid(False)
            ax.set_ylim(-2.5, 2.5)

            # Set length of x-axis equal to sequence
            # Removes possible gaps from the alignment in order to align seqs
            ax.set_xlim(0, len(alignment[i[1]][j[1]].replace('-', '')))

            # Identify & plot background colour based on CNA concordance
            conc = find_concordance(alignment[i[1]]['seq1'],
                                    alignment[i[1]]['seq2'])
            plot_profile_concordance(ax, conc)

        # Plot number of matching/mismatching bins
        ax = fig.add_subplot(gs[2, i[0]])
        plot_bin_box(ax, conc, 'matches', plot_row_name=i[0] == 0)
        ax = fig.add_subplot(gs[3, i[0]])
        plot_bin_box(ax, conc, 'mismatches', plot_row_name=i[0] == 0)

        # Plot alignment score
        ax = fig.add_subplot(gs[4, i[0]])
        plot_score_box(ax, alignment[i[1]], plot_row_name=i[0] == 0)

        if null_scores is not None:
            # Plot match/mismatching probabilities
            ax = fig.add_subplot(gs[5, i[0]])
            plot_proba_box(ax, alignment[i[1]],
                           alignment[i[1]]['match_proba'] <= match_thresh,
                           'match', plot_row_name=i[0] == 0)
            ax = fig.add_subplot(gs[6, i[0]])
            plot_proba_box(ax, alignment[i[1]],
                           alignment[i[1]]['mismatch_proba'] <= mismatch_thresh,
                           'mismatch', plot_row_name=i[0] == 0)

    # Remove superfluous spacing between subplots
    plt.subplots_adjust(wspace=0, hspace=0)
    if save:
        plt.savefig(figname+'.png', bbox_inches='tight', facecolor='w',
                    dpi=Config().dpi)
