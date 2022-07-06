##########################################################
## OncoMerge:  oncoMerge.py                             ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @github: https://github.com/plaisier-lab/OncoMerge   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################


#####################
## Import packages ##
#####################

import json
import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests # pip install statsmodels==0.10.0
import itertools
import sys
import os
from scipy.stats import hypergeom
from tqdm import tqdm
#import mygene
#mg = mygene.MyGeneInfo() # pull in function to map genes

##################################
## Read in command line options ##
##################################

parser = argparse.ArgumentParser(description='OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.')
parser.add_argument('-cf', '--config_file', help='Path to JSON encoded configuration file, overrides command line parameters', type = str)
parser.add_argument('-gp', '--gistic_path', help='Path to GISTIC output folder', type = str)
parser.add_argument('-df', '--del_file', help='Path to GISTIC deletion file (default = del_genes.conf_99.txt)', type = str, default = 'del_genes.conf_99.txt')
parser.add_argument('-af', '--amp_file', help='Path to the GISTIC amplification file (default = amp_genes.conf_99.txt)', type = str, default = 'amp_genes.conf_99.txt')
parser.add_argument('-gdf', '--gene_data_file', help='Path to the GISTIC gene data file (default = all_data_by_genes.txt)', type = str, default = 'all_data_by_genes.txt')
parser.add_argument('-aaf', '--alternate_annotation_file', help='Supply alternate annotation file to convert gene symbols to Entrez IDs (default does not import alternate annotation file and instead uses conversion embedded in GISTIC output files).', type = str)
parser.add_argument('-ln', '--label_name', help='Label for Entrez ID column in GISTIC gene data file (default = \'Gene ID\')', type = str, default = 'Gene ID')
parser.add_argument('-tf', '--thresh_file', help='Path to the GISTIC all_thresholded file (default = all_thresholded.by_genes.txt)', type = str, default = 'all_thresholded.by_genes.txt')
parser.add_argument('-gt', '--gistic_threshold', help='Cutoff value for amplifications and deletions applied to the GISTIC all_thresholded file (default = 2)', type = int, default = 2)
parser.add_argument('-pam', '--pam_file', help='Path to the protein affecting mutation (PAM) file (CSV matrix where columns are patients and genes are rows) [0 = not mutated, and 1 = mutated]', type = str)
parser.add_argument('-mscv', '--mutsig2cv_file', help='Path to a MutSig2CV output file', type = str)
parser.add_argument('-fus', '--fusions_file', help='Path to the gene fusions file (CSV matrix where columns are patients and genes are rows) [0 = not fused, and 1 = fused]', type = str, default='none')
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default = '.')
parser.add_argument('-mmf', '--min_mut_freq', help='Minimum frequency of mutation (range = 0-1; default = 0.05)', type = float, default = 0.05)
parser.add_argument('-pq', '--perm_qv', help='Permuted p-value FDR BH corrected cutoff (default = 0.1)', type = float, default = 0.1)
parser.add_argument('-sp', '--save_permutation', help='Run and save out permutation analysis to be used for comparability in another OncoMerge run (default off)', action='store_true')
parser.add_argument('-lp', '--load_permutation', help='Do not run permutation anlaysis and load permutation anlaysis from previous run (default off)', type = str, default = None)
parser.add_argument('-mlg', '--min_loci_genes', help='Minimum number of genes in loci to apply maximum final frequency filter (default = 10)', type = int, default = 10)
parser.add_argument('-mpf', '--min_pam_freq', help='Minimum PAM frequency (default = 0.01)', type = float, default = 0.01)
parser.add_argument('-tcga', '--tcga', help='Clip gistic TCGA names.', type = bool, default = False)
parser.add_argument('-bl', '--blocklist', help='List of patients (one per line) to exclude for frequency calculations.', type = str, default='none')
args = parser.parse_args()


#######################
## Define parameters ##
#######################

params = args.__dict__
if args.config_file:
    with open(args.config_file, "r") as cfg:
        tmp = json.loads(cfg.read())
        for i in tmp:
            params[i] = tmp[i]
if (not params['gistic_path']) or (not params['pam_file']) or (not params['mutsig2cv_file']) or (params['save_permutation'] and not params['load_permutation']==None):
    parser.print_help()
    sys.exit(1)


##################
## Load up data ##
##################

print('Loading data...')

# Create conversion series for gene symbol to Entrez ID
if not params['alternate_annotation_file']:
    n1 = pd.read_csv(params['gistic_path']+'/'+params['gene_data_file'],index_col=0,sep='\t',usecols=[0,1])
    n1.index = [i.split('|')[0] for i in n1.index]
    n1 = n1.drop_duplicates()
    n1 = n1[params['label_name']]
else:
    n1 = pd.read_csv(params['alternate_annotation_file'],index_col=0)[params['label_name']].apply(int)

# load up significantly mutated genes
mutSig2CV = pd.read_csv(params['mutsig2cv_file'],index_col=1)
mutSig2CV = mutSig2CV.loc[mutSig2CV.index.map(lambda x: x in n1.index)]
mutSig2CV.index = mutSig2CV.index.map(lambda x: n1.loc[x])
mutSig2CV = mutSig2CV.loc[~mutSig2CV.index.duplicated(keep='first')]
sigPAMs = list(mutSig2CV.index[mutSig2CV['q']<=0.1]) # Set at 0.1 based on Lawerence et al., Nature 2013 (PMID = 23770567)

# Get list of significantly CNA amplified genes
ampLoci = {}
ampLoci_qv = {}
amp1 = pd.read_csv(params['gistic_path']+'/'+params['amp_file'],index_col=0,sep='\t')
for col1 in amp1.columns:
    if float(amp1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
        ampLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
        ampLoci_qv[col1] = float(amp1[col1]['residual q value'])

# Get list of significantly CNA deleted genes
delLoci = {}
delLoci_qv = {}
del1 = pd.read_csv(params['gistic_path']+'/'+params['del_file'],index_col=0,sep='\t')
for col1 in del1.columns:
    if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
        delLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
        delLoci_qv[col1] = float(del1[col1]['residual q value'])

# Set up background gene numbers for gold-standard enrichment analysis
backgrounds = {'Activating':[], 'LossOfFunction':[], 'Aggregate':[]}
actTmp = [gene for locus in ampLoci for gene in ampLoci[locus]]
backgrounds['Activating'] += actTmp
backgrounds['Aggregate'] += actTmp
delTmp = [gene for locus in delLoci for gene in delLoci[locus]]
backgrounds['LossOfFunction'] += delTmp
backgrounds['Aggregate'] += delTmp
backgrounds['Aggregate'] += sigPAMs

# Load up somatically mutated genes
somMuts = pd.read_csv(params['pam_file'],index_col=0,header=0)
somMuts[somMuts != 0] = 1
if not somMuts.index.dtype=='int64':
    somMuts = somMuts.loc[somMuts.index.map(lambda x: x in n1.index)]
    somMuts.index = somMuts.index.map(lambda x: n1.loc[x])

somMuts = somMuts.loc[~somMuts.index.duplicated(keep='first')]

# Load up gene fusions
if not params['fusions_file']=='none':
    fusions = pd.read_csv(params['fusions_file'],index_col=0,header=0)
    fusions[fusions != 0] = 1
    if not fusions.index.dtype=='int64':
        fusions = fusions.loc[somMuts.index.map(lambda x: x in n1.index)]
        fusions.index = fusions.index.map(lambda x: n1.loc[x])
    fusions = fusions.loc[~fusions.index.duplicated(keep='first')]

# Read in gistic2 all_data_by_genes file
with open(params['gistic_path']+'/'+params['thresh_file'],'r') as inFile:
    tmp = inFile.readline().strip().split('\t')
    numCols1 = len(tmp)

d1 = pd.read_csv(params['gistic_path']+'/'+params['thresh_file'],index_col=0,sep='\t').drop(tmp[1], axis = 1)
if params['tcga']:
    d1.columns = [i[:12] for i in d1.columns]

d1 = d1.loc[d1.index.map(lambda x: x.split('|')[0] in n1.index)]
d1.index = d1.index.map(lambda x: n1.loc[x.split('|')[0]])
d1.index.name = 'Locus ID'

# Removing sex chromosomes (issues in CNA analysis) from d1
lociThresh = d1['Cytoband']
include = []
for i in lociThresh:
    if not(i[0]=='X' or i[0]=='Y'):
        include.append(True)
    else:
        include.append(False)

d1 = d1.loc[lociThresh[include].index].drop('Cytoband', axis = 1)

# Make sure somMuts and gistic have same samples
if params['fusions_file']=='none':
    pats = list(set(d1.columns).intersection(somMuts.columns))
else:
    pats = list((set(d1.columns).intersection(somMuts.columns)).intersection(fusions.columns))
    fusions = fusions[pats]

somMuts = somMuts[pats]
d1 = d1[pats]

# Get rid of duplicated rows
d1 = d1[~d1.index.duplicated(keep='first')]

## Fill out summary matrix
cols = ['Symbol', 'PAM_freq', 'MutSig2CV_qvalue']
if not params['fusions_file'] == 'none':
    cols += ['Fusion_freq']
cols += ['CNA_freq', 'CNA_locus', 'CNA_type', 'GISTIC_residual_q_value', 'Act_freq', 'LoF_freq', 'OM_type_selected', 'OM_empirical_p_value', 'OM_empirical_q_value', 'Genes_in_locus',  'Final_mutation_type', 'Final_freq', 'Delta_over_PAM']
summaryMatrix = pd.DataFrame(index= list(set([gene for locus in ampLoci.values() for gene in locus] + [gene for locus in delLoci.values() for gene in locus] + list(somMuts.index))), columns = cols)

# Add gene symbols
toMatch = [i for i in summaryMatrix.index if i in n1.values]
summaryMatrix.loc[toMatch,'Symbol'] = [n1.index[n1==i][0] for i in toMatch]

# Add MutSig2CV q-values
toMatch = [i for i in summaryMatrix.index if i in mutSig2CV.index]
summaryMatrix.loc[toMatch,'MutSig2CV_qvalue'] = mutSig2CV.loc[toMatch,'q']

# Load up black list of patients if given one
blacklist = []
if not params['blacklist']=='none':
    with open(params['blacklist'], 'r') as inFile:
        while 1:
            line = inFile.readline()
            if not line:
                break
            blacklist.append(line.strip().split(',')[0])

#Print out some useful information
print('\tSize of somatic mutation matrix: '+str(somMuts.shape))
if not params['fusions_file']=='none':
    print('\tSize of fusions matrix: '+str(fusions.shape))
print('\tSize of CNA matrix: '+str(d1.shape))
print('Finished loading data.')

# Make the output directory if it doesn't exists already
if not os.path.exists(params['output_path']):
    os.mkdir(params['output_path'])


########################
## Begin OnocoMerging ##
########################

# Cutoff somatic mutations based on the minimum mutation frequency (mf)
freq1 = somMuts[somMuts.columns.difference(blacklist)].sum(axis=1)/len(list(somMuts.columns.difference(blacklist)))
somMutPoint = freq1[freq1>=params['min_mut_freq']].index
tmp1 = [i for i in summaryMatrix.index if i in freq1.index]
summaryMatrix.loc[tmp1, 'PAM_freq'] = freq1.loc[tmp1]

# Cutoff fusions based on the minimum mutation frequency (mf)
if params['fusions_file']:
    freq2 = fusions[fusions.columns.difference(blacklist)].sum(axis=1)/len(list(fusions.columns.difference(blacklist)))
    fusionsMut = freq2[freq2>=params['min_mut_freq']].index
    tmp1 = [i for i in summaryMatrix.index if i in freq2.index]
    summaryMatrix.loc[tmp1, 'Fusion_freq'] = freq2.loc[tmp1]

# Precompute positive and negative dichotomized matrices
print('Precomputing dichotomized matrices...')
posdicot = (lambda x: 1 if x>=params['gistic_threshold'] else 0)
posD1 = d1.applymap(posdicot)
posFreq = posD1[posD1.columns.difference(blacklist)].mean(axis=1)
ampGenes = {j:i for i in ampLoci for j in ampLoci[i]}
negdicot = (lambda x: 1 if x<=(-params['gistic_threshold']) else 0)
negD1 = d1.applymap(negdicot)
negFreq = negD1[negD1.columns.difference(blacklist)].mean(axis=1)
delGenes = {j:i for i in delLoci for j in delLoci[i]}
print('Finished precomputing dichotomized matrices.')

# Add back genes that are high frequency amplification or deletion (20%)
for gene1 in list([i for i in posD1.index if posFreq[i]>=0.2]):
    if (not gene1 in ampGenes) and isinstance(lociThresh[gene1], str) and (lociThresh[gene1] in ampLoci_qv) and (gene1 in posD1.index):
        ampGenes[gene1] = lociThresh[gene1]
        ampLoci[lociThresh[gene1]].append(gene1)
        #print('Pos', gene1, posFreq[gene1], lociThresh[gene1])

for gene1 in list([i for i in negD1.index if negFreq[i]>=0.2]):
    if (not gene1 in delGenes) and isinstance(lociThresh[gene1], str) and (lociThresh[gene1] in delLoci_qv) and (gene1 in negD1.index):
        delGenes[gene1] = lociThresh[gene1]
        delLoci[lociThresh[gene1]].append(gene1)
        #print('Neg', gene1, negFreq[gene1], lociThresh[gene1])

print('Generating CNA summaryMatrix data...')
somMuts1 = list(set(list(ampGenes.keys())+list(delGenes.keys())).intersection(list(d1.index)))
delGenesSet = set(delGenes.keys())
ampGenesSet = set(ampGenes.keys())
bothGenesSet = (ampGenesSet.intersection(delGenesSet)).intersection(somMuts1)
onlyDel = delGenesSet.difference(ampGenesSet).intersection(somMuts1)
onlyAmp = ampGenesSet.difference(delGenesSet).intersection(somMuts1)
with tqdm(total=len(bothGenesSet)+len(onlyDel)+len(onlyAmp)) as pbar:
    for both1 in bothGenesSet:
        # Choose amplification if is higher frequency, otherwise choose deletion
        if posFreq[both1]>negFreq[both1]:
            summaryMatrix.loc[both1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Amp', ampGenes[both1], ampLoci_qv[ampGenes[both1]], posFreq[both1]]
        # Otherwise choose deletion (deletion becomse default for equal frequency; unlikely to happen)
        else:
            summaryMatrix.loc[both1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Del', delGenes[both1], delLoci_qv[delGenes[both1]], negFreq[both1]]
        pbar.update(1)
    # Straight up deletion
    for del1 in onlyDel:
        summaryMatrix.loc[del1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Del', delGenes[del1], delLoci_qv[delGenes[del1]], negFreq[del1]]
        pbar.update(1)
    # Straight up amplification
    for amp1 in onlyAmp:
        summaryMatrix.loc[amp1,['CNA_type','CNA_locus','GISTIC_residual_q_value','CNA_freq']] = ['Amp', ampGenes[amp1], ampLoci_qv[ampGenes[amp1]], posFreq[amp1]]
        pbar.update(1)
print('Done.')

# Bundle together loci for deletions and amplifications that are synonymous
lociCNAgenes = {}
lociCNA = pd.DataFrame(columns=d1.columns)
print('Bundling amplification loci...')
for loci1 in ampLoci:
    # Get matrix of CNAs for genes in loci
    dt = posD1.loc[set(posD1.index).intersection(ampLoci[loci1])]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAamp'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]

print('Bundling deletion loci...')
for loci1 in delLoci:
    # Get matrix of CNAs for genes in loci
    dt = negD1.loc[set(posD1.index).intersection(delLoci[loci1])]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAdel'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]

# Make combined matrix
# LoF = deletions + somatic point mutations (+ fusions if have data)
# Act = amplifications + somatic point mutations (+ fusions if have data)
print('Starting somatic mutations...')
pamLofAct = {}
freq = {}
potMuts = somMutPoint
if not params['fusions_file'] == 'none':
    potMuts = potMuts.append(fusionsMut)
for s1 in potMuts:
    if s1>0:
        if not str(s1) in pamLofAct:
            pamLofAct[str(s1)] = {}
        if s1 in somMuts.index:
            tmpSom = somMuts.loc[s1]
            tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
            # If potential PAM, store PAM
            if not (str(s1)+'_PAM' in pamLofAct[str(s1)] or sum(tmpSom)==0):
                pamLofAct[str(s1)][str(s1)+'_PAM'] = tmpSom
        if (not params['fusions_file'] == 'none') and (s1 in fusions.index):
            tmpFusion = fusions.loc[s1]
            if not (str(s1)+'_Fusion' in pamLofAct[str(s1)] or sum(tmpFusion)==0):
                pamLofAct[str(s1)][str(s1)+'_Fusion'] = tmpFusion
            if s1 in somMuts.index:
                tmpSom = tmpSom.add(tmpFusion).clip(0,1)
            else:
                tmpSom = tmpFusion.clip(0,1)
        if (s1 in negD1.index and s1 in posD1.index):
            tmpNeg = negD1.loc[s1]
            tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index].clip(0,1)
            tmpPos = posD1.loc[s1]
            tmpAct = tmpSom.add(tmpPos)[tmpPos.index].clip(0,1)
            if not s1 in freq:
                if (not params['fusions_file'] == 'none') and (s1 in fusions.index):
                    freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                else:
                    freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
        else:
            if not s1 in freq:
                if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                    freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}
                else:
                    freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':0,'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}

print('Starting amplifications...')
for loci1 in ampLoci:
    for s1 in set(ampLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Act
            if s1 in somMuts.index and s1 in posD1.index:
                tmpSom = somMuts.loc[s1]
                tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
                if not params['fusions_file'] == 'none' and (s1 in fusions.index):
                    tmpFusion = fusions.loc[s1]
                    tmpSom = tmpSom.add(tmpFusion).clip(0,1)
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index].clip(0,1)
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index].clip(0,1)
                if not s1 in freq:
                    if params['fusions_file'] == 'none':
                        freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                    else:
                       freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                # Store Act
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_Act' in pamLofAct[str(s1)] or tmpAct.equals(tmpSom) or tmpAct.equals(tmpPos)):
                    pamLofAct[str(s1)][str(s1)+'_Act'] = tmpAct
                    pamLofAct[str(s1)][str(s1)+'_CNAamp'] = tmpPos
    for s1 in set(ampLoci[loci1]).difference(somMuts.index):
        if s1>0:
            tmpNeg = negD1.loc[s1]
            tmpPos = posD1.loc[s1]
            if not s1 in freq:
                if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                    freq[str(s1)] = {'PAM':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                else:
                    freq[str(s1)] = {'PAM':np.nan,'Fusion':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
            # Store Amp
            if not str(s1) in pamLofAct:
                pamLofAct[str(s1)] = {}
            if not str(s1)+'_CNAamp' in pamLofAct[str(s1)]:
                pamLofAct[str(s1)][str(s1)+'_CNAamp'] = tmpPos

print('Starting deletions...')
for loci1 in delLoci:
    for s1 in set(delLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Lof
            if s1 in somMuts.index and s1 in negD1.index:
                tmpSom = somMuts.loc[s1]
                tmpSomMean = tmpSom[tmpSom.index.difference(blacklist)].mean()
                if not params['fusions_file'] == 'none' and (s1 in fusions.index):
                    tmpFusion = fusions.loc[s1]
                    tmpSom = tmpSom.add(tmpFusion).clip(0,1)
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index].clip(0,1)
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index].clip(0,1)
                if not s1 in freq:
                    if params['fusions_file'] == 'none':
                        freq[str(s1)] = {'PAM':tmpSomMean,'Fusion':tmpFusion[tmpFusion.index.difference(blacklist)].mean(),'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                    else:
                       freq[str(s1)] = {'PAM':tmpSomMean,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':tmpLoF[tmpLoF.index.difference(blacklist)].mean(),'Act':tmpAct[tmpAct.index.difference(blacklist)].mean()}
                # Store LoF
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_LoF' in pamLofAct[str(s1)] or tmpLoF.equals(tmpSom) or tmpLoF.equals(tmpNeg)):
                    pamLofAct[str(s1)][str(s1)+'_LoF'] = tmpLoF
                    pamLofAct[str(s1)][str(s1)+'_CNAdel'] = tmpNeg
    for s1 in set(delLoci[loci1]).difference(somMuts.index):
        if s1>0:
            tmpNeg = negD1.loc[s1]
            tmpPos = posD1.loc[s1]
            if not s1 in freq:
                if (params['fusions_file'] == 'none') and (s1 in fusions.index):
                    freq[str(s1)] = {'PAM':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
                else:
                    freq[str(s1)] = {'PAM':np.nan,'Fusion':np.nan,'CNAdel':tmpNeg[tmpNeg.index.difference(blacklist)].mean(),'CNAamp':tmpPos[tmpPos.index.difference(blacklist)].mean(),'LoF':np.nan,'Act':np.nan}
            # Store Del
            if not str(s1) in pamLofAct:
                pamLofAct[str(s1)] = {}
            if not str(s1)+'_CNAdel' in pamLofAct[str(s1)]:
                pamLofAct[str(s1)][str(s1)+'_CNAdel'] = tmpNeg

# Decide which mutations will be tested in permutation analysis
print('Screening for frequency...')
keepPAM = []
keepFusion = []
keepDel = []
keepAmp = []
keepers = {}
calcSig = []
for s1 in pamLofAct:
    if s1 in freq:
        freqPAM = freq[s1]['PAM']
        if 'Fusion' in freq[s1]:
            freqFusion = freq[s1]['Fusion']
        freqPos = freq[s1]['CNAamp']
        freqNeg = freq[s1]['CNAdel']
        freqAct = freq[s1]['Act']
        summaryMatrix.loc[int(s1), 'Act_freq'] = freq[s1]['Act']
        freqLoF = freq[s1]['LoF']
        summaryMatrix.loc[int(s1), 'LoF_freq'] = freq[s1]['LoF']
        if freqLoF>=params['min_mut_freq'] or freqAct>=params['min_mut_freq'] or freqPAM>=params['min_mut_freq'] or freqPos>=params['min_mut_freq'] or freqNeg>=params['min_mut_freq'] or ('Fusion' in freq[s1] and freqFusion>=params['min_mut_freq']):
            name1 = 'Unknown'
            if sum(n1.isin([int(s1)]))==1:
                name1 = n1.index[n1==int(s1)][0]
            if 'Fusion' in freq[s1]:
                print('\t'+''.join([str(i) for i in [name1+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqFusion: ', round(freqFusion,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
            else:
                print('\t'+''.join([str(i) for i in [name1+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
        # Add PAM
        if freqPAM>0 and freqPAM>=params['min_mut_freq'] and int(s1) in somMutPoint and int(s1) in sigPAMs:
            keepers[str(s1)+'_PAM'] = pamLofAct[str(s1)][str(s1)+'_PAM']
            keepPAM.append(str(s1)+'_PAM')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'PAM'
        # Add Fusion
        if 'Fusion' in freq[s1]:
            if freqFusion>0 and freqFusion>=params['min_mut_freq'] and int(s1) in fusionsMut:
                keepers[str(s1)+'_Fusion'] = pamLofAct[str(s1)][str(s1)+'_Fusion']
                keepFusion.append(str(s1)+'_Fusion')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'Fusion'
        # Add Act
        if str(s1)+'_Act' in pamLofAct[str(s1)] and freqAct>freqPAM and freqAct>=params['min_mut_freq'] and freqAct>freqLoF:
            if freqPAM>=params['min_pam_freq']:
                keepers[str(s1)+'_Act'] = pamLofAct[str(s1)][str(s1)+'_Act']
                calcSig.append(str(s1)+'_Act')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'Act'
        # Add LoF
        if str(s1)+'_LoF' in pamLofAct[str(s1)] and freqLoF>freqPAM and freqLoF>=params['min_mut_freq']and freqLoF>freqAct:
            if freqPAM>=params['min_pam_freq']:
                keepers[str(s1)+'_LoF'] = pamLofAct[str(s1)][str(s1)+'_LoF']
                calcSig.append(str(s1)+'_LoF')
                summaryMatrix.loc[int(s1),'OM_type_selected'] = 'LoF'
        # Add CNAamp
        if (str(s1)+'_CNAamp' in pamLofAct[str(s1)]) and (freqPos>=params['min_mut_freq']) and (not ((summaryMatrix.loc[int(s1),'OM_type_selected']=='Act') or (summaryMatrix.loc[int(s1),'OM_type_selected']=='LoF'))):
            keepers[str(s1)+'_CNAamp'] = pamLofAct[str(s1)][str(s1)+'_CNAamp']
            keepAmp.append(str(s1)+'_CNAamp')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'CNAamp'
        # Add CNAdel
        if ((str(s1)+'_CNAdel' in pamLofAct[str(s1)]) and (freqNeg>=params['min_mut_freq'])) and (not ((summaryMatrix.loc[int(s1),'OM_type_selected']=='Act') or (summaryMatrix.loc[int(s1),'OM_type_selected']=='LoF'))):
            keepers[str(s1)+'_CNAdel'] = pamLofAct[str(s1)][str(s1)+'_CNAdel']
            keepDel.append(str(s1)+'_CNAdel')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'CNAdel'

##################################
## Conduct permutation analysis ##
##################################
print('Permutation anlaysis...')
numPermutes = 1000

# Permute to get frequency
def singlePermute(somMutsMF, somFusionMF, somCNAsMF):
    # Get common genes
    if type(somFusionMF)==pd.core.series.Series:
        subset1 = list((set(somMutsMF.index).intersection(somCNAsMF.index)).intersection(somFusionMF.index))
    else:
        subset1 = list(set(somMutsMF.index).intersection(somCNAsMF.index))
    # Subset genes
    pam = somMutsMF.loc[subset1,:]
    cna = somCNAsMF.loc[subset1,:]
    if type(somFusionMF)==pd.core.frame.DataFrame:
        fus = somFusionMF.loc[subset1,:]
        return list((pam.sample(frac=1)+cna.sample(frac=1)+fus.sample(frac=1)).clip(0,1).mean(axis=1))
    else:
        return list((pam.sample(frac=1)+cna.sample(frac=1)).clip(0,1).mean(axis=1))

permMF_neg = []
permMF_pos = []
if params['load_permutation']==None:
    ## Compute permutations if not loading from previous run
    # Deletions
    print('\tPermuting deletions...')
    somMutsMF = somMuts[somMuts.columns.difference(blacklist)]
    somFusionMF = 'none'
    if not params['fusions_file'] == 'none':
        somFusionMF = fusions[fusions.columns.difference(blacklist)]
    somCNAsMF = negD1[negD1.columns.difference(blacklist)]
    with tqdm(total=numPermutes) as pbar:
        for i in range(numPermutes):
            permMF_neg += singlePermute(somMutsMF, somFusionMF, somCNAsMF)
            pbar.update(1)

    # Amplifications
    print('\tPermuting amplifications...')
    somMutsMF = somMuts[somMuts.columns.difference(blacklist)]
    somFusionMF = 'none'
    if not params['fusions_file'] == 'none':
        somFusionMF = fusions[fusions.columns.difference(blacklist)]
    somCNAsMF = posD1[posD1.columns.difference(blacklist)]
    with tqdm(total=numPermutes) as pbar:
        for i in range(numPermutes):
            permMF_pos += singlePermute(somMutsMF, somFusionMF, somCNAsMF)
            pbar.update(1)

    # Change precision so that is the same as what will be written out,
    # - Fixes bug where precision change leads to different behavior from freshly run (then saved) and loaded permutation data
    permMF_neg = [float(str(i)) for i in permMF_neg]
    permMF_pos = [float(str(i)) for i in permMF_pos]

    # If requested to save out permutations
    if params['save_permutation']==True:
        print('\tSaving permutations...')
        np.save(params['output_path']+'/oncomerge_delPerm', permMF_neg)
        np.save(params['output_path']+'/oncomerge_ampPerm', permMF_pos)
elif not params['load_permutation']==None:
    ## Load up permutations from previous run
    print('\tLoading previous permutations...')
    permMF_neg = np.load(params['load_permutation']+'/oncomerge_delPerm.npy')
    permMF_pos = np.load(params['load_permutation']+'/oncomerge_ampPerm.npy')

# Write Permutation Analysis file Lof_Act_sig
lofActSig = pd.DataFrame(columns = ['Symbol', 'Type','Freq','Emp.p_value'], index = calcSig)
for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_LoF'))][0]
        lofActSig['Type'].loc[sig1] = 'LoF'
        lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_LoF')]['LoF']
    elif sig1.find('Act')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1==int(sig1.rstrip('_Act'))][0]
        lofActSig['Type'].loc[sig1] = 'Act'
        lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_Act')]['Act']

# Precalculate the permuted p-values for each frequency
permdict1 = {}
permdict1['LoF'] = {}
permdict1['Act'] = {}
oMfreqs = [freq[sig1.split('_')[0]][sig1.split('_')[1]] for sig1 in calcSig]
for f in set(oMfreqs):
    permdict1['LoF'][f] = float(len([i for i in permMF_neg if i >= f]))/len(permMF_neg)
    permdict1['Act'][f] = float(len([i for i in permMF_pos if i >= f]))/len(permMF_pos)

# Add permuted p-values to summary matrix and permutation summary
for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['LoF'][freq[sig1.rstrip('_LoF')]['LoF']]
        summaryMatrix.loc[int(sig1.rstrip('_LoF')),'OM_empirical_p_value'] = permdict1['LoF'][freq[sig1.rstrip('_LoF')]['LoF']]
    elif sig1.find('Act')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['Act'][freq[sig1.rstrip('_Act')]['Act']]
        summaryMatrix.loc[int(sig1.rstrip('_Act')),'OM_empirical_p_value'] = permdict1['Act'][freq[sig1.rstrip('_Act')]['Act']]

# Filter LoF and Act based on permuted p-values
if len(lofActSig)>0:
    lofActSig['q_value'] = multipletests(lofActSig['Emp.p_value'], 0.05, method='fdr_bh')[1]
    summaryMatrix.loc[[int(i.split('_')[0]) for i in lofActSig.index],'OM_empirical_q_value'] = list(lofActSig['q_value'])
    lofActSig.sort_values('q_value').to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    # Screen out LoF and Act that don't meet significance cutoffs
    keepLofAct0 = list(lofActSig.index[lofActSig['q_value']<=params['perm_qv']])
else:
    # No LoF or Act to filter
    lofActSig.to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    keepLofAct0 = []

# Function to map mutation to locus
def findLoci(mutation, ampLoci, delLoci):
    gene, mutType = mutation.split('_')
    if (mutType == 'Act') or (mutType == 'CNAamp'):
        loci = ampLoci
    if (mutType == 'LoF') or (mutType == 'CNAdel'):
        loci = delLoci
    return [locus1 for locus1 in loci.keys() if int(gene) in loci[locus1]]

# Tabulate the number of genes per locus
lofActLoci = {}
for mut in keepLofAct0:
    mutLoci = findLoci(mut, ampLoci, delLoci)
    for locus1 in mutLoci:
        if locus1 not in lofActLoci.keys():
            lofActLoci[locus1] = []
        lofActLoci[locus1].append(mut)

# Tabulate the number of genes per locus
combinedLoci = lofActLoci
for mut in keepDel+keepAmp:
    mutLoci = findLoci(mut, ampLoci, delLoci)
    for locus1 in mutLoci:
        if not locus1 in lofActLoci:
            if locus1 not in combinedLoci.keys():
                combinedLoci[locus1] = []
            combinedLoci[locus1].append(mut)

## Decide whether to apply the maximum final frequency filter
keepLofAct1 = []
for locus in combinedLoci.keys():
    # Don't filter further with maximum final frequency filter
    if len(combinedLoci[locus]) < params['min_loci_genes']:
        keepLofAct1 += combinedLoci[locus]
        for mut in combinedLoci[locus]:
            gene, mutType = mut.split('_')
            summaryMatrix.loc[int(gene), 'Final_mutation_type'] = mutType
            if mutType=='Act' or mutType=='LoF':
                summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), mutType+'_freq']
                summaryMatrix.loc[int(gene), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene), mutType+'_freq'] - summaryMatrix.loc[int(gene), 'PAM_freq']
            else:
                summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), 'CNA_freq']
            summaryMatrix.loc[int(gene), 'Genes_in_locus'] = len(combinedLoci[locus])
    # Filter with maximum final frequency filter
    else:
        maxFF = max([summaryMatrix.loc[int(mut.split('_')[0]), mut.split('_')[1]+'_freq'] for mut in combinedLoci[locus]])
        gl1 = len([mut for mut in combinedLoci[locus] if summaryMatrix.loc[int(mut.split('_')[0]), mut.split('_')[1]+'_freq']==maxFF])
        for mut in combinedLoci[locus]:
            gene, mutType = mut.split('_')
            if summaryMatrix.loc[int(gene), mutType+'_freq']==maxFF:
                keepLofAct1.append(mut)
                summaryMatrix.loc[int(gene), 'Genes_in_locus'] = gl1
                summaryMatrix.loc[int(gene), 'Final_mutation_type'] = mutType
                if mutType=='Act' or mutType=='LoF':
                    summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), mutType+'_freq']
                    summaryMatrix.loc[int(gene), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene), mutType+'_freq'] - summaryMatrix.loc[int(gene), 'PAM_freq']
                else:
                    summaryMatrix.loc[int(gene), 'Final_freq'] = summaryMatrix.loc[int(gene), 'CNA_freq']

keepLofAct = list(set(keepLofAct1))

# Screen out PAMs that are LoF/Act
newKeepPAM = []
for pam1 in keepPAM:
    found = 0
    tmp1 = pam1.split('_')[0]
    for lofAct in keepLofAct:
        if tmp1==lofAct.split('_')[0]:
            found = 1
    if found==0:
        newKeepPAM.append(pam1)
        summaryMatrix.loc[int(tmp1), 'Final_mutation_type'] = 'PAM'
        summaryMatrix.loc[int(tmp1), 'Final_freq'] = summaryMatrix.loc[int(tmp1), 'PAM_freq']
        summaryMatrix.loc[int(tmp1), 'Delta_over_PAM'] = float(0)

# Screen out PAMs that are LoF/Act
for fus1 in keepFusion:
    found = 0
    tmp1 = fus1.split('_')[0]
    for lofAct in keepLofAct:
        if tmp1==lofAct.split('_')[0]:
            found = 1
    if found==0:
        newKeepPAM.append(fus1)
        summaryMatrix.loc[int(tmp1), 'Final_mutation_type'] = 'Fusion'
        summaryMatrix.loc[int(tmp1), 'Final_freq'] = summaryMatrix.loc[int(tmp1), 'Fusion_freq']
        summaryMatrix.loc[int(tmp1), 'Delta_over_PAM'] = summaryMatrix.loc[int(tmp1), 'Fusion_freq']-summaryMatrix.loc[int(tmp1), 'PAM_freq']

## Screen out loci that have a representative gene
# Mutations that are at or above minimum mutation frequency cutoff
highFreqLoci = lociCNA.loc[lociCNA.mean(axis=1)>=params['min_mut_freq']]

# Figure out what loci are explained by current Act or LoF genes
explainedLoc = []
for locus1 in highFreqLoci.index:
    genesInLocus = [i for i in lociCNAgenes[locus1] if (str(i)+'_Act' in keepLofAct or str(i)+'_LoF' in keepLofAct or str(i)+'_CNAamp' in keepLofAct or str(i)+'_CNAdel' in keepLofAct)]
    if len(genesInLocus)>0:
        explainedLoc.append(locus1.split('_')[0])

# Screen out all other loci in that region
keepLoc = []
for locus1 in highFreqLoci.index:
    if not locus1.split('_')[0] in explainedLoc:
        keepLoc.append(locus1)

keepLoc_dict = {}
for locus1 in set(['_'.join([i.split('_')[0],i.split('_')[-1]]) for i in keepLoc]):
    locus2 = locus1.split('_')[0]
    if locus2 in ampLoci.keys():
        keepLoc_dict[locus1] = posD1.loc[ampLoci[locus2]]
    if locus2 in delLoci.keys():
        keepLoc_dict[locus1] = negD1.loc[delLoci[locus2]]

def mode2(pat1col):
    tmpser = pat1col.mode()
    if len(tmpser)>1:
        tmp10 = 1
    else:
        tmp10 = tmpser.iloc[0]
    return tmp10

# Use most common value to determine the mutational profile of a loci
keepLoc_df = pd.DataFrame(columns = d1.columns)
tmpMode = pd.Series(index = d1.columns)
for locus1 in keepLoc_dict.keys():
    for pat1 in keepLoc_dict[locus1].columns:
         tmpMode.loc[pat1] = mode2(keepLoc_dict[locus1][pat1])
    if tmpMode.mean()>=0.05:
        keepLoc_df.loc[locus1] = tmpMode

keepLoc_df = keepLoc_df.applymap(int)

####################################
## Compile OncoMerge output files ##
####################################
finalMutFile = pd.concat([pd.DataFrame(keepers).transpose().loc[newKeepPAM].sort_index(), pd.DataFrame(keepers).transpose().loc[keepLofAct].sort_index(), keepLoc_df.sort_index()], sort=True)

# Rename all loci with only one gene
ind_list = finalMutFile.index.tolist()
for locus1 in keepLoc_df.index:
    splitUp = locus1.split('_')
    if splitUp[-1]=='CNAamp' and len(ampLoci[splitUp[0]])==1:
        idx = ind_list.index(locus1)
        ind_list[idx] = str(ampLoci[splitUp[0]][0]) + '_CNAamp'
        summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Final_mutation_type'] = 'CNAamp'
        summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Final_freq'] = summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'CNA_freq']
        summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'Delta_over_PAM'] = summaryMatrix.loc[int(ampLoci[splitUp[0]][0]), 'CNA_freq']

    if splitUp[-1]=='CNAdel' and len(delLoci[splitUp[0]])==1:
        idx = ind_list.index(locus1)
        ind_list[idx] = str(delLoci[splitUp[0]][0]) + '_CNAdel'
        summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Final_mutation_type'] = 'CNAdel'
        summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Final_freq'] = summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'CNA_freq']
        summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'Delta_over_PAM'] = summaryMatrix.loc[int(delLoci[splitUp[0]][0]), 'CNA_freq']

finalMutFile.index = ind_list

# Write out final mutation matrix to csv file
finalMutFile.to_csv(params['output_path']+'/oncoMerge_mergedMuts.csv')

# Write out summary matrix to csv file
summaryMatrix.to_csv(params['output_path']+'/oncoMerge_summaryMatrix.csv')

## Write out loci information
# Prepare for writing out
writeLoci = ['Locus_name,Genes']
for locus1 in keepLoc_df.index:
    splitUp = locus1.split('_')
    if splitUp[-1]=='CNAamp':
        writeLoci.append(locus1+','+' '.join([str(i) for i in ampLoci[splitUp[0]]]))
    if splitUp[-1]=='CNAdel':
        writeLoci.append(locus1+','+' '.join([str(i) for i in delLoci[splitUp[0]]]))

# Write out loci file
with open(params['output_path']+'/oncoMerge_CNA_loci.csv','w') as outFile:
    outFile.write('\n'.join(writeLoci))

# Write out background information for hypergeometric analysis
backgrounds = {i:[int(j) for j in backgrounds[i]] for i in backgrounds}
with open(params['output_path']+'/background.json', 'w') as outFile:
    json.dump(backgrounds, outFile)

print('Done.')

