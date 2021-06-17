# CNP-ALIGN

CNP-ALIGN: A Python package for assessing the similarity of copy number profiles based on global sequence alignment. 

To assess the similarities between copy number profiles in order to distinguish, we developed a Python package named CNP-ALIGN. CNP-ALIGN uses global sequence alignment in association with (BLOSUM-inspired) substitution matrices in order to identify and quantify the similarities of two copy number profiles in their character representation. Similarity metrics can be used to assess clonal relationships between two or more tumors. 



## Installation

CNP-ALIGN is fully build around Python 3. 

```bash
# Clone the repository
git clone https://github.com/timmocking/cnp-align.git

# Create a new conda environment
conda create -n cnp-align pip
conda activate cnp-align

# Install packages
pip install -r requirements.txt
```



## Quick start

CNP-ALIGN requires information about segmented copy number profiles in a comma separated (.csv) format. It is recommended to obtain these through the one-step circular binary segmentation algorithm of the R package [Clonality](https://www.bioconductor.org/packages/release/bioc/html/Clonality.html)

The input dataframes should at least have the following columns and follow this format:

| **ID**   | **chrom** | **loc.start** | **loc.end** | **state** |
| -------- | --------- | ------------- | ----------- | --------- |
| Sample.1 | chr1p     | 1             | 5000001     | Normal    |
| Sample.1 | chr1p     | 5100001       | 10000001    | Loss      |
| Sample.2 | chr1p     | 1             | 3000001     | Normal    |
| Sample.2 | chr1p     | 3100001       | 7000001     | Gain      |
| Sample.2 | chr1p     | 7100001       | 10000001    | Loss      |
| ...      | ...       | ...           | ...         | ...       |

CNP-ALIGN uses a substitution matrix based on clinical data from 80 patients with one or multiple tumors in combination with a list of null-scores from 100 known non-clonal relationships in order to calculate different features. These are enabled by default.

To run CNP-ALIGN:

```bash
python cnp-align.py -id <EXPERIMENT_NAME> -file <SEGMENT_FILE>
```

For every unique pair of IDs in the segment file, 3 different files are generated:

* A .json file containing the "raw" output for every chromosome arm.
* A .csv file containing all features for every chromosome arm.
* A .png visualisation of the alignment, which includes the fatures for every chromosome arm.

Additionally, a summary file is generated for all the unique pairs. 



## Command line

You can directly type `python cnp-align.py --help` and will find all the options below:

```bash
usage: python3 core.py [-h] [-file FILE] [-id EXP_ID] [-out OUT] [-autoresolve AUTORESOLVE] [-sub_mat SUB_MAT]
                       [-null_scores NULL_SCORES] [-bin_size BIN_SIZE] [-gap_open GAP_OPEN] [-gap_extend GAP_EXTEND]
                       [-match_thresh MATCH_THRESH] [-mismatch_thresh MISMATCH_THRESH] [-verbose VERBOSE]

optional arguments:
  -h, --help            show this help message and exit
  -file FILE            Clonality or DNAcopy segment file (.csv)
  -id EXP_ID            Identifier of experiment/patient/sample
  -out OUT              Path of output folder
  -autoresolve AUTORESOLVE
                        Resolve missing bins
  -sub_mat SUB_MAT      Path of subst. matrix (.json)
  -null_scores NULL_SCORES
                        Path of null scores (.json)
  -bin_size BIN_SIZE    Bin size (bp)
  -gap_open GAP_OPEN    Opening gap penalty
  -gap_extend GAP_EXTEND
                        Extension gap penalty
  -match_thresh MATCH_THRESH
                        Probability threshold for matching arms
  -mismatch_thresh MISMATCH_THRESH
                        Probability threshold for mismatching arms
  -verbose VERBOSE      Print progress (rec. for large nr. profiles).
```
