[![DOI](https://zenodo.org/badge/197803088.svg)](https://zenodo.org/badge/latestdoi/197803088)

# OncoMerge
Software to integrate somatic protein affecting mutations (PAMs) gene fusions, and copy number alterations (CNAs) for downstream computational analyses.

## Installation
Installation requires installing Python dependencies. Most of the dependencies are met by installing Anaconda (https://www.anaconda.com/). 

### Dependencies
In Anaconda:
 - pandas
 - numpy

From PyPI only:
 - statsmodels (2.1.0)
 - tqdm

#### Command to install dependencies:
```
pip install pandas numpy statsmodels==0.10.0 tqdm
```
#### Test out installed packages:
```
python
import pandas
import numpy
import statsmodels
import tqdm
```

#### If any dependency issues arise please consult dependency doucmentation to determine how to address these issues.

# OncoMerge input
## Summary
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

### Inputs to OncoMerge
An essential first step in OncoMerge is loading up and binarizing the somatic mutation data. The somatic mutation data comprised of four primary matrices: 1) PAMs, 2) fusions, 3) CNA amplifications (CNAamps), and 4) CNA deletions (CNAdels) (Figure 1). In addition, two derivative matrices Act and LoF are created by merging the PAM with the CNAamps or CNAdels matrices, respectively (Figure 1). All files are formatted as comma-separated values (CSV) files with genes as rows and patients as columns unless otherwise noted.
 - PAM matrix - The matrix values are [0 or 1]:  zero indicates the gene is not mutated in a patient tumor, and one indicates the gene is mutated in a patient tumor.
 - Fusion matrix - The matrix values are [0 or 1]:  zero indicates no gene fusion in a patient tumor, and one indicates the gene fused to another genomic locus in a patient tumor.
 - CNAamp and CNAdel matrices – The all_thresholded_by_genes.csv GISTIC output file is used to populate the CNAamp and CNAdel matrices. The all_thresholeded_by_genes matrix values range from -2 and have no positive bound, and the values indicate the copy number relative to the background. A cutoff of greater than or equal to 2 was used to identify deep amplifications and less than or equal to -2 for deep deletions. Only deep amplifications or deletions were included in these studies due to heterogeneity of cell types and tumor biopsy purity. Oncomerge allows this threshold to be modified through a command line parameter ('-gt' or '--gistic-threshold').
    - CNAamp matrix – The matrix values are [0 or 1]:  zero indicates a gene is not amplified in a patient tumor, and one indicates the gene is amplified in a patient tumor.
    - CNAdel matrix – The matrix values are [0 or 1]:  zero indicates a gene is not deleted in a patient tumor, and one indicates a gene is deleted in a patient tumor.
 - Act matrix – The Act matrix is the bitwise OR combination of the PAM, Fusion, and CNAamp matrices. The Act matrix has genes as rows and patients as columns. The matrix values are [0 or 1]: zero indicates the gene is not mutated or amplified in a patient tumor, and one indicates the gene is either mutated, fused, amplified, or some combination in a patient tumor.
 - LoF matrix – The LoF matrix is the bitwise OR combination of the PAM, Fusion, and CNAdel matrices. The LoF matrix has genes as rows and patients as columns. The matrix values are [0 or 1]:  zero indicates the gene is not mutated or deleted in a patient tumor, and one indicates the gene is either mutated, fused, deleted, or some combination in a patient tumor.\
 
#### Example input data
The [TCGA OncoMerge input data](https://doi.org/10.6084/m9.figshare.21760964.v1) can be downloaded from Figshare.

### Seeding OncoMerge with putative somatic mutations
OncoMerge focuses on likely causal somatic mutations by considering only somatic mutations that were statistically shown to be mutated more often than expected by chance alone. These statistically significant mutations were used as seeds for OncoMerge integration. Somatic PAMs used as seeds were identified with [MutSigCV2](https://pubmed.ncbi.nlm.nih.gov/25770567/) q-values less than or equal to 0.1 and a mutation frequency greater than 5%. Gene fusions used as seeds were identified as significant in [PRADA](https://pubmed.ncbi.nlm.nih.gov/29099951/) and a mutation frequency greater than 5%. CNAamps or CNAdels used as seeds were identified as significantly amplified or deleted from the amplified genes (amp_genes) or deleted genes (del_genes) GISTIC output files with residual q-values less than or equal to 0.05. CNAs from sex chromosomes (X and Y) were excluded. Genes from sex chromosomes can enter OncoMerge as seeds from PAMs or fusions. These seed genes become the starting point of the OncoMerge integration. Subsequent steps determine if Act or LoF merged mutation profiles or their component PAM, Fusion, CNAamp, or CNAdel mutation roles are the most appropriate integration model for a gene.

### Merging somatic mutations in OncoMerge
For each seed gene the mutation role is assigned based on the following criteria:
 - If Act frequency (PAM+Fusion+CNAamp) > PAM+Fusion frequency and the Act frequency ≥ 5% then the mutation role is set to Act.
 - Else LoF frequency (PAM+Fusion+CNAdel) > PAM+Fusion frequency and the LoF frequency ≥ 5% then the mutation role is set to LoF.
 - Else if the gene mutation role is not set to Act or LoF:
    - If the gene is a PAM seed gene (MutSig2CV q-value ≤ 0.1 and frequency ≥ 5%) and has a frequency greater than Fusion, CNAamp, and CNAdel, then the mutation role is set to PAM.
    - Else if the gene is a Fusion seed gene (TumorFusion.org frequency ≥ 5%) and has a frequency greater than PAM, CNAamp, and CNAdel, then the mutation role is set to Fusion.
    - Else if the gene CNAamp frequency ≥ 5% and has a frequency greater than PAM, Fusion, and CNAdel, then the mutation role is set to CNAamp.
    - Else if the gene CNAdel frequency ≥ 5% and has a frequency greater than PAM, Fusion, and CNAamp, then the mutation role is set to CNAdel.

### Permuted q-value (PQ) filter
For putative Act and LoF mutations, a permuted q-value is computed by randomizing the order of rows in the PAM, Fusion, and CNA mutation matrices' and then calculating the randomized frequency distribution for Acts and LoFs. The observed frequency for an Act or Lof mutation is then be compared to the randomized frequency distribution to compute the permuted p-value. Permuted p-values are corrected into q-values using the multiple-test Benjamini-Hochberg FDR-based correction method. Only Acts or LoFs that had a permuted q-value ≤ 0.1 were retained. Any Act or LoF with a permuted q-value > 0.1 was set to the mutation role of either PAM, Fusion, CNAamp, or CNAdel based on which mutation role had the highest frequency. The permuted q-value cutoff can be set through a command line parameter ('-pq', --perm_qv').

### Minimum final frequency (MFF) filter
A low-pass genomic filter was applied to each CNA locus if the CNA locus had ≥ 10 underlying genes. The number of genes underlying a CNA locus can be set through a command line parameter ('-mlg', --min_loci_genes'). The filter keeps only the gene(s) with the maximum mutation frequency, and all genes with the maximum mutation frequency are kept for ties.

### Microsatellite hypermutation censoring (MHC) filter
The TCGA tumors used in this study have been characterized for both [MSI](https://pubmed.ncbi.nlm.nih.gov/29850653/) and [hypermutation](https://pubmed.ncbi.nlm.nih.gov/29625053/). The tumors with MSI or hypermutation are loaded as a blocklist of patient IDs through a command line parameter ('-bl' or '--blocklist'). All tumors in the blocklist are excluded from consideration by the PQ and MFF filters while determining the genes to include in the final somatic mutation matrix. The mutation status for blocklist tumors are included in the final integrated mutation matrix.

## Usage
From the command line -h or --help will provide the usage. Parameters can be given either as a config file formatted in JSON as entries with command line argument names as keys in a dictionary, or as named comand line parameters:
```
usage: oncoMerge.py [-h] [-cf CONFIG_FILE] [-gp GISTIC_PATH] [-df DEL_FILE] [-af AMP_FILE] [-gdf GENE_DATA_FILE]
                    [-aaf ALTERNATE_ANNOTATION_FILE] [-ln LABEL_NAME] [-tf THRESH_FILE] [-gt GISTIC_THRESHOLD] [-pam PAM_FILE]
                    [-mscv MUTSIG2CV_FILE] [-fus FUSIONS_FILE] [-op OUTPUT_PATH] [-mmf MIN_MUT_FREQ] [-pq PERM_QV] [-sp]
                    [-lp LOAD_PERMUTATION] [-mlg MIN_LOCI_GENES] [-mpf MIN_PAM_FREQ] [-tcga TCGA] [-bl BLOCKLIST]

OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.

optional arguments:
  -h, --help            show this help message and exit
  -cf CONFIG_FILE, --config_file CONFIG_FILE
                        Path to JSON encoded configuration file, overrides command line parameters using argument names as entries in a JSON formatted dictionary
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
 
## Output
OncoMerge provides four output files that provide valuable information about the integration process and the final integrated mutation matrix that can be used in downstream studies. Here is a brief description of each file and its contents:
 - oncoMerge_mergedMuts.csv – The integrated mutation matrix is comprised of genes (rows) by patient tumors (columns) of mutation status after integration by OncoMerge. The matrix values are [0 or 1]:  zero indicates that the gene is not mutated in a patient tumor, and one indicates that the gene was mutated in a patient tumor.
 - oncoMerge_CNA_loci.csv – A list of the genes mapping to each CNAamp or CNAdel locus included in the OncoMerge integrated mutation matrix.
 - oncoMerge_ActLofPermPV.csv – List of all significant Act and LoF genes, their OncoMerge mutation role, frequency, empirical p-value, and empirical q-value. This output is before the application of the low-pass frequency filter.
 - oncoMerge_summaryMatrix.csv – Matrix of genes (rows) by all information gathered by OncoMerge.

### Optional output
To aid in comparisons between runs, we provide the save permutation option ('-sp' or '--save_permutation') to output permutation results so that the same permuted distribution can be used with different parameters in separate runs. We also provide the load permutation option ('-lp' or '--load_permutation') to load up the permuted distribution from a previous run. The permuted distributions are saved in the following files if requested:
 - oncomerge_ampPerm.npy, oncomerge_delPerm.npy – Snapshot of the non-deterministic permutation results from combining PAM, Fusion, and CNAamp or PAM, Fusion, and CNAdel frequencies, respectively.
