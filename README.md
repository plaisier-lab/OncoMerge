# OncoMerge
Software to integrate somatic protein affecting mutations (PAMs) gene fusions, and copy number alterations (CNAs) for downstream computational analyses.

# Installation
Installation requires installing Python dependencies. Most of the dependencies are met by installing Anaconda (https://www.anaconda.com/). 

## Dependencies
In Anaconda:
 - pandas
 - numpy
From PyPI only:
 - statsmodels (2.1.0)
 - tqdm

### Command to install dependencies:
```
pip install pandas numpy statsmodels==0.10.0 tqdm
```
## Test our installed packages:
```
python
import pandas
import numpy
import statsmodels
import tqdm
```

# Input
- PAM:
  - Somatic mutation matrix
  - MutSigCV2 result file
- Fusion:
  - PRADA output file
- CNA:
  - GISTIC2 output files:
    - Amp file
    - Del file
    - all_data_by_genes file

# Usage
From the command line -h or --help will provide the usage. Parameters can be given either as a config file formatted in JSON, or as named comand line parameters:
```
usage: oncoMerge.py [-h] [-cf CONFIG_FILE] [-gp GISTIC_PATH] [-df DEL_FILE] [-af AMP_FILE] [-gdf GENE_DATA_FILE]
                    [-aaf ALTERNATE_ANNOTATION_FILE] [-ln LABEL_NAME] [-tf THRESH_FILE] [-gt GISTIC_THRESHOLD] [-pam PAM_FILE]
                    [-mscv MUTSIG2CV_FILE] [-fus FUSIONS_FILE] [-op OUTPUT_PATH] [-mmf MIN_MUT_FREQ] [-pq PERM_QV] [-sp]
                    [-lp LOAD_PERMUTATION] [-mlg MIN_LOCI_GENES] [-mpf MIN_PAM_FREQ] [-tcga TCGA] [-bl BLOCKLIST]

OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.

optional arguments:
  -h, --help            show this help message and exit
  -cf CONFIG_FILE, --config_file CONFIG_FILE
                        Path to JSON encoded configuration file, overrides command line parameters
  -gp GISTIC_PATH, --gistic_path GISTIC_PATH
                        Path to GISTIC output folder
  -df DEL_FILE, --del_file DEL_FILE
                        Path to GISTIC deletion file (default = del_genes.conf_99.txt)
  -af AMP_FILE, --amp_file AMP_FILE
                        Path to the GISTIC amplification file (default = amp_genes.conf_99.txt)
  -gdf GENE_DATA_FILE, --gene_data_file GENE_DATA_FILE
                        Path to the GISTIC gene data file (default = all_data_by_genes.txt)
  -aaf ALTERNATE_ANNOTATION_FILE, --alternate_annotation_file ALTERNATE_ANNOTATION_FILE
                        Supply alternate annotation file to convert gene symbols to Entrez IDs (default does not import
                        alternate annotation file and instead uses conversion embedded in GISTIC output files).
  -ln LABEL_NAME, --label_name LABEL_NAME
                        Label for Entrez ID column in GISTIC gene data file (default = 'Gene ID')
  -tf THRESH_FILE, --thresh_file THRESH_FILE
                        Path to the GISTIC all_thresholded file (default = all_thresholded.by_genes.txt)
  -gt GISTIC_THRESHOLD, --gistic_threshold GISTIC_THRESHOLD
                        Cutoff value for amplifications and deletions applied to the GISTIC all_thresholded file (default = 2)
  -pam PAM_FILE, --pam_file PAM_FILE
                        Path to the protein affecting mutation (PAM) file (CSV matrix where columns are patients and genes are
                        rows) [0 = not mutated, and 1 = mutated]
  -mscv MUTSIG2CV_FILE, --mutsig2cv_file MUTSIG2CV_FILE
                        Path to a MutSig2CV output file
  -fus FUSIONS_FILE, --fusions_file FUSIONS_FILE
                        Path to the gene fusions file (CSV matrix where columns are patients and genes are rows) [0 = not
                        fused, and 1 = fused]
  -op OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path you would like to output OncoMerged files (default = current directory)
  -mmf MIN_MUT_FREQ, --min_mut_freq MIN_MUT_FREQ
                        Minimum frequency of mutation (range = 0-1; default = 0.05)
  -pq PERM_QV, --perm_qv PERM_QV
                        Permuted p-value FDR BH corrected cutoff (default = 0.1)
  -sp, --save_permutation
                        Run and save out permutation analysis to be used for comparability in another OncoMerge run (default
                        off)
  -lp LOAD_PERMUTATION, --load_permutation LOAD_PERMUTATION
                        Do not run permutation anlaysis and load permutation anlaysis from previous run (default off)
  -mlg MIN_LOCI_GENES, --min_loci_genes MIN_LOCI_GENES
                        Minimum number of genes in loci to apply maximum final frequency filter (default = 10)
  -mpf MIN_PAM_FREQ, --min_pam_freq MIN_PAM_FREQ
                        Minimum PAM frequency (default = 0.01)
  -tcga TCGA, --tcga TCGA
                        Clip gistic TCGA names.
  -bl BLOCKLIST, --blocklist BLOCKLIST
                        List of patients (one per line) to exclude for frequency calculations.
 ```

# Output
There are three output files from OncoMerge:
 1. background.json - set of genes tested for hypergeometric/Fisher's exact tests
 2. oncoMerge_ActLofPermPV.csv - permuted p-values
 3. oncomerge_CNA_loci.csv - a file which holds the loci included in the run
 4. oncoMerge_mergedMuts.csv - integrated somatic mutation matrix with patients as columns and genes are rows
 5. oncoMerge_summaryMatrix.csv - collection of all information integrated to generate the integrated somatic mutation matrices
 
 ### Optional output
 To facilitate comparisons the permutation data can be saved out using the --save_permutaton parameters, and loaded back up as needed with --load_permutation options. Which will generate:
 - oncomerge_ampPerm.npy
 - oncomerge_delPerm.npy
