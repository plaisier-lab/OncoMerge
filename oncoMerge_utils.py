##########################################################
## OncoMerge:  oncoMerge_utils.py                       ##
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

# oncoMerge_utils.py full of useful functions for oncoMerge and the accompanying scripts
import pandas as pd

def DefinePaths(Run_Name,mutation_frequency):
#Please use this section to define the paths to each of the imput files
#Each must be a string.
    mutation_frequency = str(mutation_frequency)
    PathTo = {}
    if Run_Name == 'SKCM':
        p = 'TM'
    else:
        p = 'TP'
    PathTo['symbol_conversion'] = 'OncoMerge_input_g2e_converter.csv'
    PathTo['Dels'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/del_genes.conf_99.txt'
    #Amplifications File:
    PathTo['Amps'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/amp_genes.conf_99.txt'
    # Gene Data File:
    PathTo['gene_data'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt'
    # Title of Column containing entrez ids in the gistic Gene Data file (all_data_by_genes.txt) should match title of entrez id column in symbol_conversion_path:
    PathTo['label_name'] = 'Locus ID'
    #all_Thresholded File:
    PathTo['thresh'] = 'GISTIC_99/gdac.broadinstitute.org_'+Run_Name+'-'+p+'.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt'
    # Somatic Mutations File:
    PathTo['somatic_mutations'] = 'mutations_mc3/'+Run_Name+'_somMutMC3.csv'
    # MutSig2CV File:
    PathTo['mutsig2CV'] = 'sig2cv/'+Run_Name+'_sig2cv.csv'
    # Where to store output:
    PathTo['output'] = 'output/oncoMerged_'+Run_Name
    PathTo['log'] = 'output/log'
    PathTo['CNALoci'] = PathTo['output']+'/CNA_loci_'+mutation_frequency+'.csv'
    PathTo['oncoMerged'] = PathTo['output']+'/oncoMerged_'+mutation_frequency+'.csv'
    PathTo['LofActSig'] = PathTo['output']+'/LoF_Act_sig_'+mutation_frequency+'.csv'
    return PathTo


def Load_oncoMerge_Input(Run_Name,mutation_frequency):
    PathTo = DefinePaths(Run_Name,mutation_frequency)
    label_name = PathTo['label_name']
    # Read in gistic all_data_by_genes file
    if PathTo['symbol_conversion'] == 'No File':
        n1 = pd.read_csv(PathTo['gene_data'],index_col=0,sep='\t',usecols=[0,1])
        n1.index = [i.split('|')[0] for i in n1.index]
        n1 = n1.drop_duplicates()
    else:
        n1 = pd.read_csv(PathTo['symbol_conversion'],index_col=0)[label_name].apply(int)
    # Should pull these from mutSig2CV file
    mutSig2CV = pd.read_csv(PathTo['mutsig2CV'],index_col=1)
    mutSig2CV = mutSig2CV.loc[mutSig2CV.index.map(lambda x: x in n1.index)]
    mutSig2CV.index = mutSig2CV.index.map(lambda x: n1.loc[x])
    mutSig2CV = mutSig2CV.loc[~mutSig2CV.index.duplicated(keep='first')]
    sigPAMs = list(mutSig2CV.index[mutSig2CV['q']<=0.05])
    # Get list of significantly CNA amplified genes
    ampGenes = []
    ampLoci = {}
    amp1 = pd.read_csv(PathTo['Amps'],index_col=0,sep='\t')
    for col1 in amp1.columns:
        if float(amp1[col1]['residual q value'])<=0.05:
            ampGenes += [i.lstrip('[').rstrip(']') for i in list(amp1[col1].dropna()[3:])]
            ampLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(amp1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
    # Get list of significantly CNA deleted genes
    delGenes = []
    delLoci = {}
    del1 = pd.read_csv(PathTo['Dels'],index_col=0,sep='\t')
    for col1 in del1.columns:
        if float(del1[col1]['residual q value'])<=0.05 and not (col1[0]=='X' or col1=='Y'):
            delGenes += [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])]
            delLoci[col1] = list(set([n1.loc[i] for i in [i.lstrip('[').rstrip(']').split('|')[0] for i in list(del1[col1].dropna()[3:])] if i in n1.index and n1.loc[i]>0]))
    # Convert gene ids for amp and del genes
    ampGenes = [n1.loc[i] for i in ampGenes if i in n1.index]
    delGenes = [n1.loc[i] for i in delGenes if i in n1.index]
    # Load up somatically mutated genes
    somMuts = pd.read_csv(PathTo['somatic_mutations'],index_col=0,header=0)
    somMuts = somMuts.loc[somMuts.index.map(lambda x: x in n1.index)]
    somMuts.index = somMuts.index.map(lambda x: n1.loc[x])
    somMuts = somMuts.loc[~somMuts.index.duplicated(keep='first')]
    # Read in gistic all_data_by_genes file
    with open(PathTo['thresh'],'r') as inFile:
        numCols1 = len(inFile.readline().strip().split('\t'))
    d1 = pd.read_csv(PathTo['thresh'],index_col=0,sep='\t').drop(label_name, axis = 1)#,usecols=[i for i in range(numCols1) if not i in [1,2]])
    d1.columns = [i[:12] for i in d1.columns]
    d1 = d1.loc[d1.index.map(lambda x: x.split('|')[0] in n1.index)]
    d1.index = d1.index.map(lambda x: n1.loc[x.split('|')[0]])
    d1.index.name = 'Locus ID'
    # Removing sex chromosomes (issues in CNA analysis) from d1
    lociThresh = d1['Cytoband']
    include = []
    for i in lociThresh:
        if not(i[0]=='X' or i[0]=='Y'):
            include.append(     True)
        else:
            include.append(False)
    d1 = d1.loc[lociThresh[include].index].drop('Cytoband', axis = 1)
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
    # Get rid of duplicated rows
    d1 = d1[~d1.index.duplicated(keep='first')]
    return n1, mutSig2CV, sigPAMs, amp1, ampGenes, ampLoci, del1, delGenes, delLoci, somMuts, d1, lociThresh, PathTo



def Load_oncoMerge_Results(Run_Name, mutation_frequency):
    PathTo = DefinePaths(Run_Name, mutation_frequency)
    loci = pd.read_csv(PathTo['CNALoci'],index_col = 0)
    omdf = pd.read_csv(PathTo['oncoMerged'],index_col = 0, header = 0)
    actlof_Sig = pd.read_csv(PathTo['LofActSig'], index_col = 0)
    return loci, omdf, actlof_Sig, PathTo
