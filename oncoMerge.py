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

import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import itertools
import json
import argparse
import sys
import dill
# Read in command line options
parser = argparse.ArgumentParser(description='Used to choose a cancer type')
parser.add_argument('-tt', '--tumor_type', help='descriptive three/four letter code (e.g. TCGA cancer codes)', type = str)
parser.add_argument('-rn', '--run_name', help='Descriptive name the run', type = str)
parser.add_argument('-mmf', '--min_mut_freq', help='Minimum frequency of mutation (range = 0-1; default = 0.05)', type = float, default='0.05')
parser.add_argument('-dp', '--del_path', help='Path to GISTIC deletion file (e.g. del_genes.conf_99.txt)', type = str)
parser.add_argument('-ap', '--amp_path', help='Path to the GISTIC amplification file (amp_genes.conf_99.txt)', type = str)
parser.add_argument('-gdp', '--gene_data_path', help='Path to the GISTIC gene data file (all_data_by_genes.txt)', type = str)
parser.add_argument('-ln', '--label_name', help='Label for Entrez ID column in GISTIC gene data file (default = \'Gene ID\')', type = str, default='Gene ID')
parser.add_argument('-tp', '--thresh_path', help='Path to the GISTIC all_thresholded file (all_thresholded.by_genes.txt)', type = str)
parser.add_argument('-smp', '--som_mut_path', help='Path to the somatic mutations file (CSV matrix where columns are patients and genes are rows) [0 = not mutated, and 1 = mutated]', type = str)
parser.add_argument('-mscv', '--mutsig2_cv', help='Path to a MutSig2CV output file', type = str)
parser.add_argument('-pp', '--perm_pv', help='Permuted p-value FDR BH corrected cutoff (default = 0.1)', type = float, default='0.1')
parser.add_argument('-op', '--output_path', help='Path you would like to output OncoMerged files (default = current directory)', type = str, default='.')
args = parser.parse_args()

'''args.tumor_type = 'CHOL'
args.run_name = 'Test1'
args.min_mut_freq = 0.05
args.label_name = 'Locus ID'
args.perm_pv = 0.1
if args.tumor_type == 'SKCM':
    args.del_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt'
    args.amp_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt'
    args.gene_data_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt'
    args.thresh_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt'
else:
    args.del_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt'
    args.amp_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt'
    args.gene_data_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt'
    args.thresh_path = 'GISTIC_99/gdac.broadinstitute.org_'+args.tumor_type+'-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt'

args.som_mut_path = 'mutations_entrez_mc3/'+args.tumor_type+'_somMutMC3.csv'
args.mutsig2_cv = 'sig2cv/'+args.tumor_type+'_sig2cv.csv'
args.output_path = 'output/oncoMerge_4_3'
args.output_path = 'output/oncoMerge_Changes'
'''
if (not args.tumor_type) or (not args.run_name) or (not args.del_path) or (not args.amp_path) or (not args.gene_data_path) or (not args.thresh_path) or (not args.som_mut_path) or (not args.mutsig2_cv):
    parser.print_help()
    sys.exit(1)

#For saving out each level of filtering:
in_genes = {}

# Read in gistic all_data_by_genes file
n1 = pd.read_csv(args.gene_data_path,index_col=0,sep='\t',usecols=[0,1])
n1.index = [i.split('|')[0] for i in n1.index]
n1 = n1.drop_duplicates()

# Should pull these from mutSig2CV file
mutSig2CV = pd.read_csv(args.mutsig2_cv,index_col=1)
sigPAMs = [n1.loc[i][0] for i in list(mutSig2CV.index[mutSig2CV['q']<=0.05]) if i in n1.index]

# Get list of significantly CNA amplified genes
ampGenes = []
ampLoci = {}
amp1 = pd.read_csv(args.amp_path,index_col=0,sep='\t')
for col1 in amp1.columns:
    if float(amp1[col1]['residual q value'])<=0.05:
        ampGenes += [i.lstrip('[').rstrip(']') for i in list(amp1[col1].dropna()[3:])]
        ampLoci[col1] = list(set([n1[args.label_name].loc[i] for i in [i.lstrip('[').rstrip(']') for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1[args.label_name].loc[i]>0]))

# Get list of significantly CNA deleted genes
delGenes = []
delLoci = {}
del1 = pd.read_csv(args.del_path,index_col=0,sep='\t')
for col1 in del1.columns:
    if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
        delGenes += [i.lstrip('[').rstrip(']') for i in list(del1[col1].dropna()[3:])]
        delLoci[col1] = list(set([n1[args.label_name].loc[i] for i in [i.lstrip('[').rstrip(']') for i in list(del1[col1].dropna()[3:])] if i in n1.index and n1[args.label_name].loc[i]>0]))

# Convert gene ids for amp and del genes
ampGenes = [n1[args.label_name].loc[i] for i in ampGenes if i in n1.index]
delGenes = [n1[args.label_name].loc[i] for i in delGenes if i in n1.index]

# Load up somatically mutated genes
somMuts = pd.read_csv(args.som_mut_path,index_col=0,header=0)

# Read in gistic all_data_by_genes file
with open(args.thresh_path,'r') as inFile:
    numCols1 = len(inFile.readline().strip().split('\t'))

d1 = pd.read_csv(args.thresh_path,index_col=0,sep='\t',usecols=[i for i in range(numCols1) if not i in [0,2]])
d1.columns = [i[:12] for i in d1.columns]

# Removing sex chromosomes (issues in CNA analysis) from d1
lociThresh = pd.read_csv(args.thresh_path,index_col=0,sep='\t',usecols=[1,2])
include = []
for i in lociThresh['Cytoband']:
    if not(i[0]=='X' or i[0]=='Y'):
        include.append(True)
    else:
        include.append(False)

d1 = d1.loc[lociThresh[include].index]

# Removing sex chromosomes from ampLoci
delMes = [i for i in ampLoci if i[0]=='X' or i[0]=='Y']
for delMe in delMes:
    del ampLoci[delMe]

# Removing sex chromosomes from delLoci
delMes = [i for i in delLoci if i[0]=='X' or i[0]=='Y']
for delMe in delMes:
    del delLoci[delMe]

# Make sure somMuts and gistic have same samples
somMuts = somMuts[list(set(d1.columns).intersection(somMuts.columns))]
d1 = d1[list(set(d1.columns).intersection(somMuts.columns))]

# Cutoff somatic mutations based on the minimum mutation frequency (mf)
freq1 = somMuts.sum(axis=1)/len(list(somMuts.columns))
somMutPoint = freq1[freq1>=args.min_mut_freq].index
print(somMuts.shape)

# Get rid of duplicated rows
d1 = d1[~d1.index.duplicated(keep='first')]
print(d1.shape)
print('Finished loading data.')

# Precompute positive and negative dichotomized matrices
posdicot = (lambda x: 1 if x>=2 else 0)
posD1 = d1.applymap(posdicot)
negdicot = (lambda x: 1 if x<=(-2) else 0)
negD1 = d1.applymap(negdicot)
print('Finished precomputing dichotomized matrices.')

# Merge loci for deletions and amplifications
lociCNAgenes = {}
lociCNA = pd.DataFrame(columns=d1.columns)
cna1 = []
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

print('Combining amplification loci..')
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

print(lociCNA.shape)
in_genes['d1'] = d1.index.values
in_genes['somMuts'] = somMuts.index.values
in_genes['posD1'] = posD1.index.values
in_genes['negD1'] = negD1.index.values
# Make combined matrix
# LoF = deletions + somatic point mutations
# Act = amplifications + somatic point mutations
print('Starting somatic mutations...')
pamLofAct = {}
freq = {}
in_genes['somMutPoint'] = [s1 for s1 in somMutPoint]
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

print('Starting deletions...')
in_genes['delLoci'] = [i for loci1 in delLoci for i in delLoci[loci1]]
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

print('Starting amplifications...')
in_genes['ampLoci'] = [i for loci1 in ampLoci for i in ampLoci[loci1]]
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

print('Screening for frequency...')
keepPAM = []
keepers = {}
calcSig = []
in_genes['pamLofAct'] = pamLofAct
for s1 in pamLofAct:
    if s1 in freq:
        freqPAM = freq[s1]['PAM']
        freqNeg = freq[s1]['CNAdel']
        freqLoF = freq[s1]['LoF']
        freqPos = freq[s1]['CNAamp']
        freqAct = freq[s1]['Act']
        if freqLoF>=0.05 or freqAct>=0.05 or freqPAM>=0.05:
            print(''.join([str(i) for i in [n1.index[n1[args.label_name]==int(s1)][0]+' ('+str(s1),') - FreqPAM: ', round(freqPAM,3), ' | FreqNeg: ', round(freqNeg,3), ' | FreqLoF: ', round(freqLoF,3), ' | FreqPos: ', round(freqPos,3),' | FreqAct: ', round(freqAct,3)]]))
        if freqPAM>0 and freqPAM>=args.min_mut_freq and int(s1) in somMutPoint and int(s1) in sigPAMs:
            keepers[str(s1)+'_PAM'] = pamLofAct[str(s1)][str(s1)+'_PAM']
            keepPAM.append(str(s1)+'_PAM')
        if str(s1)+'_LoF' in pamLofAct[str(s1)]  and freqLoF>freqPAM and freqLoF>=args.min_mut_freq and len(delGenes)>0:
            keepers[str(s1)+'_LoF'] = pamLofAct[str(s1)][str(s1)+'_LoF']
            calcSig.append(str(s1)+'_LoF')
        if str(s1)+'_Act' in pamLofAct[str(s1)] and freqAct>freqPAM and freqAct>=args.min_mut_freq and len(ampGenes)>0:
            keepers[str(s1)+'_Act'] = pamLofAct[str(s1)][str(s1)+'_Act']
            calcSig.append(str(s1)+'_Act')

## Write out LoF and Act significance file
print('Permutation anlaysis...')
numPermutes = 1000
# Permute to get frequency
def singlePermute(somMutsMF, somCNAsMF):
    tmp1 = pd.Series(np.random.permutation(somMutsMF), index=somMutsMF.index)
    tmp2 = pd.Series(np.random.permutation(somCNAsMF), index=somCNAsMF.index)
    subset1 = set(somMutsMF.index).intersection(somCNAsMF.index)
    return list(tmp1.loc[subset1]+tmp2.loc[subset1])

# Deletions
print('\tPermuting deletions...')
permMF_neg = []
somMutsMF = somMuts.transpose().mean()
somCNAsMF = negD1.transpose().mean()
for i in range(numPermutes):
    permMF_neg += singlePermute(somMutsMF, somCNAsMF)

# Amplifications
print('\tPermuting amplifications...')
permMF_pos = []
somMutsMF = somMuts.transpose().mean()
somCNAsMF = posD1.transpose().mean()
for i in range(numPermutes):
    permMF_pos += singlePermute(somMutsMF, somCNAsMF)

# Write out file
in_genes['calcSig'] = [i.split('_')[0] for i in calcSig]
lofActSig = pd.DataFrame(columns = ['Symbol', 'Type','Freq','Emp.p_value'], index = calcSig)
for sig1 in calcSig:
    if sig1.find('LoF')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1[args.label_name]==int(sig1.rstrip('_LoF'))][0]
        lofActSig['Type'].loc[sig1] = 'LoF'
        lofActSig['Freq'].loc[sig1] = freq[sig1.rstrip('_LoF')]['LoF']
    elif sig1.find('Act')>0:
        lofActSig['Symbol'].loc[sig1] = n1.index[n1[args.label_name]==int(sig1.rstrip('_Act'))][0]
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
    elif sig1.find('Act')>0:
        lofActSig.loc[sig1, 'Emp.p_value'] = permdict1['Act'][freq[sig1.rstrip('_Act')]['Act']]

if len(lofActSig)>0:
    lofActSig['q_value'] = multipletests(lofActSig['Emp.p_value'], 0.05, method='fdr_bh')[1]
    lofActSig.sort_values('q_value').to_csv(args.output_path+'/Lof_Act_sig/'+args.run_name+'_'+args.tumor_type+'_'+str(args.min_mut_freq)+'_LoF_Act_sig.csv')
    ## Screen out LoF and Act that don't meet significance cutoffs
    keepLofAct = list(lofActSig.index[lofActSig['q_value']<=args.perm_pv])
else:
    lofActSig.to_csv(args.output_path+'/Lof_Act_sig/'+args.run_name+'_'+args.tumor_type+'_'+str(args.min_mut_freq)+'_LoF_Act_sig.csv')
    keepLofAct = []

## Screen out PAMs that are LoF/Act
newKeepPAM = []
for pam1 in keepPAM:
    found = 0
    tmp1 = pam1.split('_')[0]
    for lofAct in keepLofAct:
        if tmp1==lofAct.split('_')[0]:
            found = 1
    if found==0:
        newKeepPAM.append(pam1)

## Screen out loci that have a representative gene
# Mutations that are at or above minimum mutation frequency cutoff
highFreqLoci = lociCNA[lociCNA.mean(axis=1)>=args.min_mut_freq]

# Figure out what loci are likely explained by current Act or LoF genes
explainedLoc = []
for locus1 in highFreqLoci.index:
    genesInLocus = [i for i in lociCNAgenes[locus1] if (str(i)+'_Act' in keepers.keys() or str(i)+'_LoF' in keepers.keys())]
    if len(genesInLocus)>0:
        explainedLoc.append(locus1.split('_')[0])

# Screen out all other loci in that region
keepLoc = []
for locus1 in highFreqLoci.index:
    if not locus1.split('_')[0] in explainedLoc:
        keepLoc.append(locus1)


in_genes['keepers'] = [i.split('_')[0] for i in keepers.keys()]
in_genes['keepLofAct'] = [i.split('_')[0] for i in keepLofAct]
in_genes['keepLoc'] = [i.split('_')[0] for i in keepLoc]
in_genes['newKeepPAM'] = newKeepPAM

## Write out OncoMerge output file
finalMutFile = pd.concat([pd.DataFrame(keepers).transpose().loc[newKeepPAM].sort_index(),pd.DataFrame(keepers).transpose().loc[keepLofAct].sort_index(),lociCNA.loc[keepLoc].sort_index()])
finalMutFile.to_csv(args.output_path+'/oncoMerged/'+args.run_name+'_'+args.tumor_type+'_oncoMerged_mmf_'+str(args.min_mut_freq)+'.csv')

# Pickle in_genes for QC
#in_genes['final'] = [i.split('_')[0] for i in finalMutFile.index]
picme = open(args.output_path+'/pickles/'+args.tumor_type+'_in_genes_No_Changes.pkl', 'wb')
dill.dump_session(picme)
picme.close()

## Write out loci
# Prepare for writing out
writeLoci = ['Locus_name,Genes']
for locus1 in lociCNAgenes:
    writeLoci.append(locus1+','+' '.join([str(i) for i in lociCNAgenes[locus1]]))

# Write out file
with open(args.output_path+'/CNA_loci/'+args.run_name+'_'+args.tumor_type+'_'+str(args.min_mut_freq)+'_CNA_loci.csv','w') as outFile:
    outFile.write('\n'.join(writeLoci))

print('Done.')

