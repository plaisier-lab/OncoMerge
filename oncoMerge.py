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
## @Author:  Chris Plaisier, Rob Schultz                ##
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
from statsmodels.stats.multitest import multipletests
import itertools
import sys
import os
from scipy.stats import hypergeom


##################################
## Read in command line options ##
##################################

parser = argparse.ArgumentParser(description='OncoMerge merges patient Protein Affecting Mutations (PAMs) and Copy Number Alterations (CNAs) into a unified mutation matrix.')
parser.add_argument('-cf', '--config_path', help='Path to JSON encoded configuration file, overrides command line parameters', type = str)
parser.add_argument('-gp', '--gistic_path', help='Path to GISTIC output folder', type = str)
parser.add_argument('-dp', '--del_path', help='Path to GISTIC deletion file (default = del_genes.conf_99.txt)', type = str, default = 'del_genes.conf_99.txt')
parser.add_argument('-ap', '--amp_path', help='Path to the GISTIC amplification file (default = amp_genes.conf_99.txt)', type = str, default = 'amp_genes.conf_99.txt')
parser.add_argument('-gdp', '--gene_data_path', help='Path to the GISTIC gene data file (default = all_data_by_genes.txt)', type = str, default = 'all_data_by_genes.txt')
parser.add_argument('-cfp', '--conversion_file_path', help='Supply hand annotated file to convert gene symbols to Entrez IDs (default does not import hand annotated conversion file and instead uses conversion embedded in GISTIC output files).', type = str)
parser.add_argument('-ln', '--label_name', help='Label for Entrez ID column in GISTIC gene data file (default = \'Gene ID\')', type = str, default = 'Gene ID')
parser.add_argument('-tp', '--thresh_path', help='Path to the GISTIC all_thresholded file (default = all_thresholded.by_genes.txt)', type = str, default = 'all_thresholded.by_genes.txt')
parser.add_argument('-smp', '--som_mut_path', help='Path to the somatic mutations file (CSV matrix where columns are patients and genes are rows) [0 = not mutated, and 1 = mutated]', type = str)
parser.add_argument('-mscv', '--mutsig2cv_path', help='Path to a MutSig2CV output file', type = str)
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default = '.')
parser.add_argument('-mmf', '--min_mut_freq', help='Minimum frequency of mutation (range = 0-1; default = 0.05)', type = float, default = 0.05)
parser.add_argument('-pq', '--perm_qv', help='Permuted p-value FDR BH corrected cutoff (default = 0.1)', type = float, default = 0.1)
parser.add_argument('-sp', '--save_permutation', help='Run and save out permutation analysis to be used for comparability in another OncoMerge run (default off)', action='store_true')
parser.add_argument('-lp', '--load_permutation', help='Do not run permutation anlaysis and load permutation anlaysis from previous run (default off)', type = str, default = None)
parser.add_argument('-mcr', '--min_coincidence_rate', help='Minimum coincidence rate for the shallow coincidence filter (default = 0.001)', type = float, default = 0.001)
parser.add_argument('-lg', '--min_loci_genes', help='Minimum number of genes in loci to apply shallow coincidence filter (default = 2)', type = int, default = 2)
args = parser.parse_args()


#######################
## Define parameters ##
#######################

params = args.__dict__
if args.config_path:
    with open(args.config_path, "r") as cfg:
        tmp = json.loads(cfg.read())
        for i in tmp:
            params[i] = tmp[i]
if (not params['gistic_path']) or (not params['som_mut_path']) or (not params['mutsig2cv_path']) or (params['save_permutation'] and not params['load_permutation']==None):
    parser.print_help()
    sys.exit(1)


###############
## Functions ##
###############

def enrichment(overlap1, bakcground1, set1, set2):
    x = len(overlap1)
    M = len(background1)
    n = len(set1)
    N = len(set2)
    return 1-hypergeom.cdf(x, M, n, N)


##################
## Load up data ##
##################

print('Loading data...')

# Create conversion series for gene symbol to Entrez ID
if not params['conversion_file_path']:
    n1 = pd.read_csv(params['gistic_path']+'/'+params['gene_data_path'],index_col=0,sep='\t',usecols=[0,1])
    n1.index = [i.split('|')[0] for i in n1.index]
    n1 = n1.drop_duplicates()
    n1 = n1[params['label_name']]
else:
    n1 = pd.read_csv(params['conversion_file_path'],index_col=0)[params['label_name']].apply(int)


# load up significantly mutated genes
mutSig2CV = pd.read_csv(params['mutsig2cv_path'],index_col=1)
mutSig2CV = mutSig2CV.loc[mutSig2CV.index.map(lambda x: x in n1.index)]
mutSig2CV.index = mutSig2CV.index.map(lambda x: n1.loc[x])
mutSig2CV = mutSig2CV.loc[~mutSig2CV.index.duplicated(keep='first')]
sigPAMs = list(mutSig2CV.index[mutSig2CV['q']<=0.05])

# Get list of significantly CNA amplified genes
ampLoci = {}
ampLoci_qv = {}
amp1 = pd.read_csv(params['gistic_path']+'/'+params['amp_path'],index_col=0,sep='\t')
for col1 in amp1.columns:
    if float(amp1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
        ampLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
        ampLoci_qv[col1] = float(amp1[col1]['residual q value'])

# Get list of significantly CNA deleted genes
delLoci = {}
delLoci_qv = {}
del1 = pd.read_csv(params['gistic_path']+'/'+params['del_path'],index_col=0,sep='\t')
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

# Load up enrichment sets
"""enrichmentSets = {}
if 'enrichment_sets' in params:
    for set1 in params['enrichment_sets']:
        enSet1 = pd.read_csv(params['enrichment_sets'][set1], header=0)
        enrichmentSets[set1] = {}
        enrichmentSets[set1]['Aggregate'] = {'Activating': [], 'LossOfFunction': []}
        for cancer in set(enSet1['Cancer']):
            enrichmentSets[set1][cancer] = {}
            for type1 in ['Activating', 'LossOfFunction']:
                enrichmentSets[set1][cancer][type1] = list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])
                enrichmentSets[set1]['Aggregate'][type1] += list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])
"""

# Load up somatically mutated genes
somMuts = pd.read_csv(params['som_mut_path'],index_col=0,header=0)
if not somMuts.index.dtype=='int64':
    somMuts = somMuts.loc[somMuts.index.map(lambda x: x in n1.index)]
    somMuts.index = somMuts.index.map(lambda x: n1.loc[x])
somMuts = somMuts.loc[~somMuts.index.duplicated(keep='first')]

# Read in gistic all_data_by_genes file
with open(params['gistic_path']+'/'+params['thresh_path'],'r') as inFile:
    tmp = inFile.readline().strip().split('\t')
    numCols1 = len(tmp)
d1 = pd.read_csv(params['gistic_path']+'/'+params['thresh_path'],index_col=0,sep='\t').drop(tmp[1], axis = 1)
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
somMuts = somMuts[list(set(d1.columns).intersection(somMuts.columns))]
d1 = d1[list(set(d1.columns).intersection(somMuts.columns))]

# Get rid of duplicated rows
d1 = d1[~d1.index.duplicated(keep='first')]

## Fill out summary matrix
summaryMatrix = pd.DataFrame(index= list(set([gene for locus in ampLoci.values() for gene in locus] + [gene for locus in delLoci.values() for gene in locus] + list(somMuts.index))), columns = ['Symbol', 'PAM_freq', 'MutSig2CV_qvalue', 'CNA_freq', 'CNA_locus', 'CNA_type', 'GISTIC_residual_q_value', 'Act_freq', 'Act_shallow_coincidence_rate', 'LoF_freq', 'LoF_shallow_coincidence_rate', 'OM_type_selected', 'OM_empirical_p_value', 'OM_empirical_q_value', 'Genes_in_locus',  'Final_mutation_type', 'Final_freq', 'Delta_over_PAM'])

# Add gene symbols
toMatch = [i for i in summaryMatrix.index if i in n1.values]
summaryMatrix.loc[toMatch,'Symbol'] = [n1.index[n1==i][0] for i in toMatch]

# Add MutSig2CV q-values
toMatch = [i for i in summaryMatrix.index if i in mutSig2CV.index]
summaryMatrix.loc[toMatch,'MutSig2CV_qvalue'] = mutSig2CV.loc[toMatch,'q']

# Add CNA amp locus and GISTIC q-value
for locus1 in ampLoci:
    for gene1 in ampLoci[locus1]:
        summaryMatrix.loc[gene1,'CNA_type'] = 'Amp'
        summaryMatrix.loc[gene1,'CNA_locus'] = locus1
        summaryMatrix.loc[gene1,'GISTIC_residual_q_value'] = ampLoci_qv[locus1]

# Add CNA del locus and GISTIC q-value
for locus1 in delLoci:
    for gene1 in delLoci[locus1]:
        summaryMatrix.loc[gene1,'CNA_type'] = 'Del'
        summaryMatrix.loc[gene1,'CNA_locus'] = locus1
        summaryMatrix.loc[gene1,'GISTIC_residual_q_value'] = delLoci_qv[locus1]

# Load enrichment 


# Print out some useful information
print('\tSize of CNA matrix: '+str(d1.shape))
print('\tSize of somatic mutation matrix: '+str(somMuts.shape))
print('Finished loading data.')

# Make the output directory if it doesn't exists already
if not os.path.exists(params['output_path']):
    os.mkdir(params['output_path'])


########################
## Begin onocoMerging ##
########################

# Cutoff somatic mutations based on the minimum mutation frequency (mf)
freq1 = somMuts.sum(axis=1)/len(list(somMuts.columns))
somMutPoint = freq1[freq1>=params['min_mut_freq']].index

# Precompute positive and negative dichotomized matrices
posdicot = (lambda x: 1 if x>=2 else 0)
posD1 = d1.applymap(posdicot)
summaryMatrix.loc[list(summaryMatrix[summaryMatrix['CNA_type']=='Amp'].index),'CNA_freq'] = posD1.loc[list(summaryMatrix[summaryMatrix['CNA_type']=='Amp'].index)].mean(axis=1)
negdicot = (lambda x: 1 if x<=(-2) else 0)
negD1 = d1.applymap(negdicot)
summaryMatrix.loc[list(summaryMatrix[summaryMatrix['CNA_type']=='Del'].index),'CNA_freq'] = negD1.loc[list(summaryMatrix[summaryMatrix['CNA_type']=='Del'].index)].mean(axis=1)
print('Finished precomputing dichotomized matrices.')

# Merge loci for deletions and amplifications
lociCNAgenes = {}
lociCNA = pd.DataFrame(columns=d1.columns)
print('Combining amplification loci...')
for loci1 in ampLoci:
    # Get matrix of CNAs for genes in loci
    dt = posD1.loc[ampLoci[loci1]]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAamp'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]

print('Combining deletion loci...')
for loci1 in delLoci:
    # Get matrix of CNAs for genes in loci
    dt = negD1.loc[delLoci[loci1]]
    # Get unique rows
    dedup = dt.drop_duplicates(keep='first')
    # Get genes which match and add to output dictionaries
    for i in range(len(dedup.index)):
        cnaName = loci1+'_'+str(i)+'_CNAdel'
        lociCNA.loc[cnaName] = dedup.iloc[i]
        lociCNAgenes[cnaName] = [j for j in dt.index if dedup.iloc[i].equals(dt.loc[j])]

# Calculate shallow coincidence
shallowCoincidence = {'Act':{},'LoF':{}}
totalPats = len(d1.columns)

# Activating shallow coincidence
for loci1 in ampLoci:
    for ampGene in ampLoci[loci1]:
        if ampGene in somMuts.index and ampGene in d1.index:
            coincPats = set(d1.columns[list(d1.loc[ampGene].ge(1))]).intersection(somMuts.columns[list(somMuts.loc[ampGene]==1)])
            shallowCoincidence['Act'][ampGene] = float(len(coincPats))/float(totalPats)
            summaryMatrix.loc[ampGene, 'Act_shallow_coincidence_rate'] = float(len(coincPats))/float(totalPats)

# Loss of Function shallow coincidence
for locus1 in delLoci:
    for delGene in delLoci[locus1]:
        if delGene in somMuts.index and delGene in d1.index:
            coincPats = set(d1.columns[list(d1.loc[delGene].le(-1))]).intersection(somMuts.columns[list(somMuts.loc[delGene]==1)])
            shallowCoincidence['LoF'][delGene] = float(len(coincPats))/float(totalPats)
            summaryMatrix.loc[delGene, 'LoF_shallow_coincidence_rate'] = float(len(coincPats))/float(totalPats)

# Make combined matrix
# LoF = deletions + somatic point mutations
# Act = amplifications + somatic point mutations
print('Starting somatic mutations...')
pamLofAct = {}
freq = {}
for s1 in somMutPoint:
    if s1>0:
        if not str(s1) in pamLofAct:
            pamLofAct[str(s1)] = {}
        tmpSom = somMuts.loc[s1]
        # If potential PAM, store PAM
        if not (str(s1)+'_PAM' in pamLofAct[str(s1)] or sum(tmpSom)==0):
            pamLofAct[str(s1)][str(s1)+'_PAM'] = tmpSom
        if (s1 in negD1.index and s1 in posD1.index):
            tmpNeg = negD1.loc[s1]
            tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
            tmpLoF[tmpLoF==2] = 1
            tmpPos = posD1.loc[s1]
            tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
            tmpAct[tmpAct==2] = 1
            if not s1 in freq:
                freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
        else:
            if not s1 in freq:
                freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':0,'CNAamp':0,'LoF':0,'Act':0}

print('Starting amplifications...')
for loci1 in ampLoci:
    for s1 in set(ampLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Act
            if s1 in somMuts.index and s1 in posD1.index:
                tmpSom = somMuts.loc[s1]
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                tmpLoF[tmpLoF==2] = 1
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                tmpAct[tmpAct==2] = 1
                if not s1 in freq:
                    freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
                # Store Act
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_Act' in pamLofAct[str(s1)] or tmpAct.equals(tmpSom) or tmpAct.equals(tmpPos)):
                    pamLofAct[str(s1)][str(s1)+'_Act'] = tmpAct

print('Starting deletions...')
for loci1 in delLoci:
    for s1 in set(delLoci[loci1]).intersection(somMuts.index):
        if s1>0:
            # If potential Lof
            if s1 in somMuts.index and s1 in negD1.index:
                tmpSom = somMuts.loc[s1]
                tmpNeg = negD1.loc[s1]
                tmpLoF = tmpSom.add(tmpNeg)[tmpNeg.index]
                tmpLoF[tmpLoF==2] = 1
                tmpPos = posD1.loc[s1]
                tmpAct = tmpSom.add(tmpPos)[tmpPos.index]
                tmpAct[tmpAct==2] = 1
                if not s1 in freq:
                    freq[str(s1)] = {'PAM':tmpSom.mean(),'CNAdel':tmpNeg.mean(),'CNAamp':tmpPos.mean(),'LoF':tmpLoF.mean(),'Act':tmpAct.mean()}
                # Store LoF
                if not str(s1) in pamLofAct:
                    pamLofAct[str(s1)] = {}
                if not (str(s1)+'_LoF' in pamLofAct[str(s1)] or tmpLoF.equals(tmpSom) or tmpLoF.equals(tmpNeg)):
                    pamLofAct[str(s1)][str(s1)+'_LoF'] = tmpLoF

# Decide which mutations will be tested in permutation analysis
print('Screening for frequency...')
keepPAM = []
keepers = {}
calcSig = []
for s1 in pamLofAct:
    if s1 in freq:
        freqPAM = freq[s1]['PAM']
        summaryMatrix.loc[int(s1), 'PAM_freq'] = freq[s1]['PAM']
        freqPos = freq[s1]['CNAamp']
        freqNeg = freq[s1]['CNAdel']
        if summaryMatrix.loc[int(s1),'CNA_type']=='Del':
            summaryMatrix.loc[int(s1), 'CNA_freq'] = freq[s1]['CNAdel']
        elif summaryMatrix.loc[int(s1),'CNA_type']=='Amp':
            summaryMatrix.loc[int(s1), 'PAM_freq'] = freq[s1]['CNAamp']
        freqAct = freq[s1]['Act']
        summaryMatrix.loc[int(s1), 'Act_freq'] = freq[s1]['Act']
        freqLoF = freq[s1]['LoF']
        summaryMatrix.loc[int(s1), 'LoF_freq'] = freq[s1]['LoF']
        if freqLoF>=0.05 or freqAct>=0.05 or freqPAM>=0.05:
            print('\t'+''.join([str(i) for i in [n1.index[n1==int(s1)][0]+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
        if freqPAM>0 and freqPAM>=params['min_mut_freq'] and int(s1) in somMutPoint and int(s1) in sigPAMs:
            keepers[str(s1)+'_PAM'] = pamLofAct[str(s1)][str(s1)+'_PAM']
            keepPAM.append(str(s1)+'_PAM')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'PAM'
        if str(s1)+'_Act' in pamLofAct[str(s1)] and freqAct>freqPAM and freqAct>=params['min_mut_freq']:
            keepers[str(s1)+'_Act'] = pamLofAct[str(s1)][str(s1)+'_Act']
            calcSig.append(str(s1)+'_Act')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'Act'
        if str(s1)+'_LoF' in pamLofAct[str(s1)] and freqLoF>freqPAM and freqLoF>=params['min_mut_freq']:
            keepers[str(s1)+'_LoF'] = pamLofAct[str(s1)][str(s1)+'_LoF']
            calcSig.append(str(s1)+'_LoF')
            summaryMatrix.loc[int(s1),'OM_type_selected'] = 'LoF'


##################################
## Conduct permutation analysis ##
##################################
print('Permutation anlaysis...')
numPermutes = 1000

# Permute to get frequency
def singlePermute(somMutsMF, somCNAsMF):
    tmp1 = pd.Series(np.random.permutation(somMutsMF), index=somMutsMF.index)
    tmp2 = pd.Series(np.random.permutation(somCNAsMF), index=somCNAsMF.index)
    subset1 = set(somMutsMF.index).intersection(somCNAsMF.index)
    return list(tmp1.loc[subset1]+tmp2.loc[subset1])

# If not loading up from previous run
permMF_neg = []
permMF_pos = []
if params['load_permutation']==None:
    # Deletions
    print('\tPermuting deletions...')
    somMutsMF = somMuts.transpose().mean()
    somCNAsMF = negD1.transpose().mean()
    for i in range(numPermutes):
        permMF_neg += singlePermute(somMutsMF, somCNAsMF)

    # Amplifications
    print('\tPermuting amplifications...')
    somMutsMF = somMuts.transpose().mean()
    somCNAsMF = posD1.transpose().mean()
    for i in range(numPermutes):
        permMF_pos += singlePermute(somMutsMF, somCNAsMF)
    
    # Change precision so that is the same as what will be written out,
    # - Fixes bug where precision change leads to different behavior from freshly run (then saved) and loaded permutation data
    permMF_neg = [float(str(i)) for i in permMF_neg]
    permMF_pos = [float(str(i)) for i in permMF_pos]

    # If requested to save out permutations
    if params['save_permutation']==True:
         print('\tSaving permutations...')
         with open(params['output_path']+'/oncomerge_delPerm.txt','w') as delPermFile, open(params['output_path']+'/oncomerge_ampPerm.txt','w') as ampPermFile:
            delPermFile.write('\n'.join([str(i) for i in permMF_neg]))
            ampPermFile.write('\n'.join([str(i) for i in permMF_pos]))
# Load up permutations from previous run
elif not params['load_permutation']==None:
    print('\tLoading previous permutations...')
    with open(params['load_permutation']+'/oncomerge_delPerm.txt','r') as delPermFile, open(params['load_permutation']+'/oncomerge_ampPerm.txt','r') as ampPermFile:
        permMF_neg = [float(i.strip()) for i in delPermFile.readlines()]
        permMF_pos = [float(i.strip()) for i in ampPermFile.readlines()]

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

permdict1 = {}
permdict1['LoF'] = {}
permdict1['Act'] = {}
oMfreqs = [freq[sig1[:-4]][sig1.split('_')[1]] for sig1 in calcSig]
for f in set(oMfreqs):
    permdict1['LoF'][f] = float(len([i for i in permMF_neg if i >= f]))/len(permMF_neg)
    permdict1['Act'][f] = float(len([i for i in permMF_pos if i >= f]))/len(permMF_pos)

for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['LoF'][freq[sig1.rstrip('_LoF')]['LoF']]
        summaryMatrix.loc[int(sig1.rstrip('_LoF')),'OM_empirical_p_value'] = permdict1['LoF'][freq[sig1.rstrip('_LoF')]['LoF']]
    elif sig1.find('Act')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['Act'][freq[sig1.rstrip('_Act')]['Act']]
        summaryMatrix.loc[int(sig1.rstrip('_Act')),'OM_empirical_p_value'] = permdict1['LoF'][freq[sig1.rstrip('_Act')]['Act']]

if len(lofActSig)>0:
    lofActSig['q_value'] = multipletests(lofActSig['Emp.p_value'], 0.05, method='fdr_bh')[1]
    summaryMatrix.loc[[int(i.split('_')[0]) for i in lofActSig.index],'OM_empirical_q_value'] = list(lofActSig['q_value'])
    lofActSig.sort_values('q_value').to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    # Screen out LoF and Act that don't meet significance cutoffs
    keepLofAct0 = list(lofActSig.index[lofActSig['q_value']<=params['perm_qv']])
else:
    lofActSig.to_csv(params['output_path']+'/oncoMerge_ActLofPermPV.csv')
    keepLofAct0 = []

# Filter passenger mutations from keepLofAct using Shallow Coincidence
def FindLoci(inp, ampLoci, delLoci):
    gene1, mutType = inp.split('_')
    if mutType == 'Act':
        LociSet = ampLoci
    if mutType == 'LoF':
        LociSet = delLoci
    ret = []
    for locus1 in LociSet.keys():
        if int(gene1) in LociSet[locus1]:
            ret.append(locus1)
    return ret

LofActLoci_dict = {}
for mut in keepLofAct0:
    mutLoci1 = FindLoci(mut, ampLoci, delLoci)
    for locus1 in mutLoci1:
        if locus1 not in LofActLoci_dict.keys():
            LofActLoci_dict[locus1] = []
        LofActLoci_dict[locus1].append(mut)

keepLofAct1 = []
for locus1 in LofActLoci_dict.keys():
    if len(LofActLoci_dict[locus1]) < params['min_loci_genes']:
        keepLofAct1 += LofActLoci_dict[locus1]
        for mut in LofActLoci_dict[locus1]:
            gene1, mutType1 = mut.split('_')
            summaryMatrix.loc[int(gene1), 'Final_mutation_type'] = mutType1
            summaryMatrix.loc[int(gene1), 'Final_freq'] = summaryMatrix.loc[int(gene1), mutType1+'_freq']
            summaryMatrix.loc[int(gene1), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene1), mutType1+'_freq'] - summaryMatrix.loc[int(gene1), 'PAM_freq']
            summaryMatrix.loc[int(gene1), 'Genes_in_locus'] = len(LofActLoci_dict[locus1])
    else:
        for mut in LofActLoci_dict[locus1]:
            gene1, mutType1 = mut.split('_')
            summaryMatrix.loc[int(gene1), 'Genes_in_locus'] = len(LofActLoci_dict[locus1])
            if shallowCoincidence[mutType1][int(gene1)] < params['min_coincidence_rate']:
                print(mut + ' filtered with shallow coincidence = ' + str(shallowCoincidence[mutType1][int(gene1)]))
            else:
                keepLofAct1.append(mut)
                summaryMatrix.loc[int(gene1), 'Final_mutation_type'] = mutType1
                summaryMatrix.loc[int(gene1), 'Final_freq'] = summaryMatrix.loc[int(gene1), mutType1+'_freq']
                summaryMatrix.loc[int(gene1), 'Delta_over_PAM'] = summaryMatrix.loc[int(gene1), mutType1+'_freq'] - summaryMatrix.loc[int(gene1), 'PAM_freq']
                summaryMatrix.loc[int(gene1), 'Genes_in_locus'] = len(LofActLoci_dict[locus1])


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

## Screen out loci that have a representative gene
# Mutations that are at or above minimum mutation frequency cutoff
highFreqLoci = lociCNA[lociCNA.mean(axis=1)>=params['min_mut_freq']]

# Figure out what loci are explained by current Act or LoF genes
explainedLoc = []
for locus1 in highFreqLoci.index:
    genesInLocus = [i for i in lociCNAgenes[locus1] if (str(i)+'_Act' in keepLofAct or str(i)+'_LoF' in keepLofAct)]
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
finalMutFile = pd.concat([pd.DataFrame(keepers).transpose().loc[newKeepPAM].sort_index(), pd.DataFrame(keepers).transpose().loc[keepLofAct].sort_index(), keepLoc_df.sort_index()])

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

## Write out information for hypergeometric analysis
# Can be used to load in gene information from final mutation file:
#    [[int(j.split('_')[0]),j.split('_')[1]] for j in finalMutFile.index[[not i for i in list(finalMutFile.index.str.contains('p|q'))]]]

# Write out background information for hypergeometric analysis
backgrounds = {i:[int(j) for j in backgrounds[i]] for i in backgrounds}
with open(params['output_path']+'/background.json', 'w') as outFile:
    json.dump(backgrounds, outFile)

print('Done.')

