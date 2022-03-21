"""
This module contains enables running CNP-ALIGN from the command-line.
"""
import os
import argparse
import pandas as pd
import json

from cnp_align.utils import autoresolve, format_alignment_results, summarize_results
from cnp_align.align import Alignment
from cnp_align.format import Profile
from cnp_align.CGHcall import convert_CGHcall_probs


def parse_args():
    """Parse all command-line arguments."""
    parser = argparse.ArgumentParser(prog='python3 core.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='')
    parser.add_argument('-file', dest='file',
                        help='Input file (.csv)')
    parser.add_argument('-format', dest='format', default='clonality',
                        help='Format of input file')
    parser.add_argument('-id', dest='exp_id', default='PATIENT',
                        help='Identifier of experiment/patient/sample')
    parser.add_argument('-out', dest='out', default='results',
                        help='Path of output folder')
    parser.add_argument('-autoresolve', dest='autoresolve', default=True,
                        help='Resolve missing bins')
    parser.add_argument('-split', dest='split', default='True',
                        help='Format of chroms (whole or split)')
    parser.add_argument('-sub_mat', dest='sub_mat',
                        default='data/general.blosum.json',
                        help='Path of subst. matrix (.json)')
    parser.add_argument('-null_scores', dest='null_scores',
                        default=None,
                        help='Path of null scores (.json)')
    parser.add_argument('-bin_size', dest='bin_size', default=100000,
                        help='Bin size (bp)')
    parser.add_argument('-gap_open', dest='gap_open', default=-100000,
                        help='Opening gap penalty')
    parser.add_argument('-gap_extend', dest='gap_extend', default=-100000,
                        help='Extension gap penalty')
    parser.add_argument('-match_thresh', dest='match_thresh', default=0.10,
                        help='Probability threshold for matching arms')
    parser.add_argument('-mismatch_thresh', dest='mismatch_thresh', default=0.50,
                        help='Probability threshold for mismatching arms')
    parser.add_argument('-verbose', dest='verbose', default=True,
                        help='Print progress (rec. for large nr. profiles).')
    args = parser.parse_args()
    return args


def main(args=False):
    if not args:
        args = parse_args()

    # Load data
    data = pd.read_csv(args.file)

    # Convert CGHcall data segment format
    if args.format == 'CGHcall':
        data = convert_CGHcall_probs(data, args.bin_size, args.autoresolve)

    with open(args.sub_mat) as f:
        matrix = json.load(f)
    null_scores = args.null_scores
    if args.null_scores is not None:
        with open(args.null_scores) as f:
            null_scores = json.load(f)

    # Prepare output folder
    if not os.path.exists(f'{args.out}'):
        os.makedirs(f'{args.out}')

    # Create alignments/visualizations for all unique subcombinations
    # of patients in the segment files
    samples = list(data['ID'].unique())
    observed = []
    dfs = {}
    for s1 in samples:
        for s2 in samples:
            if s1 != s2 and {s1, s2} not in observed:
                if args.verbose:
                    print("Processing pair " + s1 + '/' + s2)
                sample1 = data[data['ID'] == s1]
                sample2 = data[data['ID'] == s2]
                # Resolve missing bins
                if args.autoresolve and args.format != 'CGHcall':
                    sample1 = autoresolve(sample1, args.bin_size)
                    sample2 = autoresolve(sample2, args.bin_size)
                    if args.verbose:
                        print("Resolved all missing areas")
                # Convert to Profile class objects
                profile1 = Profile(s1, args.bin_size, sample1, format=args.format)
                profile2 = Profile(s2, args.bin_size, sample2, format=args.format)
                # Perform alignments
                if args.verbose:
                    print("Aligning pair " + s1 + '/' + s2 + '...')
                A = Alignment(profile1, profile2)
                results = A.align(matrix, format=args.format, 
                                gap_open=args.gap_open, gap_extend=args.gap_extend,
                                null_scores=null_scores, )
                # Dump alignment to json
                outfile = f'{args.out}/{s1}_{s2}_alignment.json'
                with open(outfile, 'w') as json_file:
                    json.dump(results, json_file)
                # Save alignment features split by chromosome arm
                df = format_alignment_results(results, null_scores)
                df.to_csv(f'{args.out}/{s1}_{s2}_split_alignment_stats.csv')
                # Visualisation
                if args.verbose:
                    print("Visualizing pair " + s1 + '/' + s2 + '...')
                A.plot(null_scores=null_scores, save=True,
                       split=eval(args.split),
                       figname=f'{args.out}/{s1}_{s2}',
                       match_thresh=float(args.match_thresh),
                       mismatch_thresh=float(args.mismatch_thresh))
                dfs[f'{s1}_{s2}'] = df
                observed.append({s1, s2})

    # Save alignment stats for all possible subcombinations
    # Summarize statistics using mean and median + probability cut-offs
    df = summarize_results(dfs, args.exp_id, null_scores)
    df.to_csv(f'{args.out}/{args.exp_id}_alignment_stats.csv')


if __name__ == '__main__':
    main()
