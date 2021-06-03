# OncoMerge
Software to integrate somatic protein affecting mutations (PAMs) and copy number alterations (CNAs) for downstream computational analyses.

# Usage
From the command line help. The Oncomerge program requires files from a GISTIC run, somatic mutation matrix, and MutSig2CV output file. The program has two different cutoff parameters:
 - Minimum mutation frequency - the minimum somatic mutation frequency that will be allowed into the output file
 - Permuted p-value cutoff - alpha value cutoff for the FDR BH corrected q-value from the permuted p-values
 
 Finally, there are two parameters which can be used to label a run:
 - Tumor type
 - Run name
```
usage: oncoMerge.py [-h] [-tt TUMOR_TYPE] [-rn RUN_NAME] [-mmf MIN_MUT_FREQ]
                    [-dp DEL_PATH] [-ap AMP_PATH] [-gdp GENE_DATA_PATH]
                    [-ln LABEL_NAME] [-tp THRESH_PATH] [-smp SOM_MUT_PATH]
                    [-mscv MUTSIG2_CV] [-pp PERM_PV] [-op OUTPUT_PATH]

Used to choose a cancer type

optional arguments:
  -h, --help            show this help message and exit
  -tt TUMOR_TYPE, --tumor_type TUMOR_TYPE
                        descriptive three/four letter code (e.g. TCGA cancer
                        codes)
  -rn RUN_NAME, --run_name RUN_NAME
                        Descriptive name the run
  -mmf MIN_MUT_FREQ, --min_mut_freq MIN_MUT_FREQ
                        Minimum frequency of mutation (range = 0-1; default =
                        0.05)
  -dp DEL_PATH, --del_path DEL_PATH
                        Path to GISTIC deletion file (e.g.
                        del_genes.conf_99.txt)
  -ap AMP_PATH, --amp_path AMP_PATH
                        Path to the GISTIC amplification file
                        (amp_genes.conf_99.txt)
  -gdp GENE_DATA_PATH, --gene_data_path GENE_DATA_PATH
                        Path to the GISTIC gene data file
                        (all_data_by_genes.txt)
  -ln LABEL_NAME, --label_name LABEL_NAME
                        Label for Entrez ID column in GISTIC gene data file
                        (default = 'Gene ID')
  -tp THRESH_PATH, --thresh_path THRESH_PATH
                        Path to the GISTIC all_thresholded file
                        (all_thresholded.by_genes.txt)
  -smp SOM_MUT_PATH, --som_mut_path SOM_MUT_PATH
                        Path to the somatic mutations file (CSV matrix where
                        columns are patients and genes are rows) [0 = not
                        mutated, and 1 = mutated]
  -mscv MUTSIG2_CV, --mutsig2_cv MUTSIG2_CV
                        Path to a MutSig2CV2 output file
  -pp PERM_PV, --perm_pv PERM_PV
                        Permuted p-value FDR BH corrected cutoff (default =
                        0.1)
  -op OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path you would like to output OncoMerged files
                        (default = current directory)
```

# Output
There are three output files from OncoMerge:
 1. *_CNA_loci.csv - a file which holds the loci included in the run
 2. *_LoF_GoF_sig.csv - a file which holds the permutation analysis from the LoFs and GoFs
 3. *_finalMutFile_deep_filtered_mf_*.csv - the merged PAM and CNA matrix file (patients are the columns, and the rows are the merged somatic mutations)
