##########################################################
## aggregateTCGA.py                                     ##
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

import pandas as pd
import json
from scipy.stats import hypergeom


###############
## Functions ##
###############

def enrichment(overlap1, background1, set1, set2):
    x = len(overlap1)
    M = len(background1)
    n = len(set1)
    N = len(set2)
    return [x, n, N, M, 1-hypergeom.cdf(x, M, n, N)]


######################
## Common variables ##
######################

tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
#tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM'] # , 'TCGT'
analyses = ['no_filter','pq','mcr','pq_mcr']
toLoad = { "TCGA_Consensus": "gold_standards/TCGA_Consensus/TCGA_consensus_gold_standard.csv", 
           "CGC": "gold_standards/CGC/CGC_gold_standard.csv",
           "OncogeneDB": "gold_standards/OncogeneDB/OncogeneDB_gold_standard.csv",
           "TSGeneDB": "gold_standards/TSGeneDB/TSGeneDB_gold_standard.csv",
           "Vogelstein": "gold_standards/Vogelstein/Vogelstein_gold_standard.csv" }

###############
## Load data ##
###############

# Load up enrichment sets
enrichmentSets = {}
for set1 in toLoad:
    enSet1 = pd.read_csv(toLoad[set1], header=0)
    enrichmentSets[set1] = {}
    enrichmentSets[set1]['Aggregate'] = {'Activating': [], 'LossOfFunction': [], 'All': []}
    for cancer in set(enSet1['Cancer']):
        enrichmentSets[set1][cancer] = {}
        for type1 in ['Activating', 'LossOfFunction', 'All']:
            enrichmentSets[set1][cancer][type1] = list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])
            enrichmentSets[set1]['Aggregate'][type1] += list(enSet1.loc[(enSet1['Cancer']==cancer) & (enSet1['Type']==type1),'Entrez'])

# Load up background information
background = {}
for tumor in tumors:
    with open('TCGA/output/'+tumor+'/no_filter/background.json', 'r') as bkgd:
            background[tumor] = json.loads(bkgd.read())

# Load up gene sets for each analysis
merged = {}
geneSets = {}
for tumor in tumors:
    merged[tumor] = {}
    geneSets[tumor] = {}
    for analysis in analyses:
        merged[tumor][analysis] = pd.read_csv('TCGA/output/'+tumor+'/'+analysis+'/oncoMerge_mergedMuts.csv', header = 0, index_col = 0)
        geneSets[tumor][analysis] = [[int(j.split('_')[0]),j.split('_')[1]] for j in merged[tumor][analysis].index[[not i for i in list(merged[tumor][analysis].index.str.contains('p|q'))]]]


####################
## Hypergeometric ##
####################

# Compute hypergeometric for cancer-specific
csp_hypgeoResults = []
for analysis in analyses:
    # TCGA_Consensus (Activating)
    gs_cancerSpecific_Act = [[i,j[0]] for i in geneSets for j in geneSets[i][analysis] if j[1]=='Act' or j[1]=='CNAamp']
    eSet_cancerSpecific_Act = [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating'] if not i in ['PANCAN','Aggregate']]
    """eSet_cancerSpecific_Act = []
    for tmp1 in [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['Activating'] if not i in ['PANCAN','Aggregate']]+[[i,j] for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['Activating'] if not i in ['PANCAN','Aggregate']]:
        if not tmp1 in eSet_cancerSpecific_Act:
            eSet_cancerSpecific_Act.append(tmp1)
    """
    ovlp_cancerSpecific_Act = [i for i in gs_cancerSpecific_Act if i in eSet_cancerSpecific_Act]
    bkgd_cancerSpecific_Act = [[i,j] for i in background for j in background[i]['Activating']]

    # TCGA_Consensus Loss of function
    gs_cancerSpecific_LoF = [[i,j[0]] for i in geneSets for j in geneSets[i][analysis] if j[1]=='LoF' or j[1]=='CNAdel']
    eSet_cancerSpecific_LoF = [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]
    """eSet_cancerSpecific_LoF = []
    for tmp1 in [[i,j] for i in enrichmentSets['TCGA_Consensus'] for j in enrichmentSets['TCGA_Consensus'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]+[[i,j] for i in enrichmentSets['CGC'] for j in enrichmentSets['CGC'][i]['LossOfFunction'] if not i in ['PANCAN','Aggregate']]:
        if not tmp1 in eSet_cancerSpecific_LoF:
            eSet_cancerSpecific_LoF.append(tmp1)
    """
    ovlp_cancerSpecific_LoF = [i for i in gs_cancerSpecific_LoF if i in eSet_cancerSpecific_LoF]
    bkgd_cancerSpecific_LoF = [[i,j] for i in background for j in background[i]['LossOfFunction']]

    # Oncogene (http://ongene.bioinfo-minzhao.org/ongene_human.txt)
    gs_cancerSpecific_og = [j[0] for i in geneSets for j in geneSets[i][analysis] if j[1]=='Act' or j[1]=='CNAamp']
    eSet_cancerSpecific_og = [i for i in enrichmentSets['OncogeneDB']['Any']['Activating']]
    ovlp_cancerSpecific_og = [i for i in gs_cancerSpecific_og if i in eSet_cancerSpecific_og]
    bkgd_cancerSpecific_og = [j for i in background for j in background[i]['Activating']]
    
    # Tumor suppressors (https://bioinfo.uth.edu/TSGene/Human_TSGs.txt)
    gs_cancerSpecific_ts = [j[0] for i in geneSets for j in geneSets[i][analysis] if j[1]=='LoF' or j[1]=='CNAdel']
    eSet_cancerSpecific_ts = [i for i in enrichmentSets['TSGeneDB']['Any']['LossOfFunction']]
    ovlp_cancerSpecific_ts = [i for i in gs_cancerSpecific_ts if i in eSet_cancerSpecific_ts]
    bkgd_cancerSpecific_ts = [j for i in background for j in background[i]['LossOfFunction']]
        
    # Aggregate
    #gs_cancerSpecific_Agg = gs_cancerSpecific_Act+gs_cancerSpecific_LoF
    #eSet_cancerSpecific_Agg = eSet_cancerSpecific_Act+eSet_cancerSpecific_LoF
    #ovlp_cancerSpecific_Agg = ovlp_cancerSpecific_Act+ovlp_cancerSpecific_LoF
    #bkgd_cancerSpecific_Agg = bkgd_cancerSpecific_Act+bkgd_cancerSpecific_LoF
 
    # All = CGC + TCGA_Consensus + Vogelstein
    gs_cancerSpecific_All = list(set([j[0] for i in geneSets for j in geneSets[i][analysis]]))
    eSet_cancerSpecific_All = list(set(enrichmentSets['CGC']['Any']['All']+enrichmentSets['TCGA_Consensus']['Any']['All']+enrichmentSets['Vogelstein']['Any']['All']))
    ovlp_cancerSpecific_All = [i for i in gs_cancerSpecific_All if i in eSet_cancerSpecific_All]
    bkgd_cancerSpecific_All = list(set([j for i in background for j in background[i]['Aggregate']]))

    # Compute hypergeo
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (Act)','Activating']+enrichment(ovlp_cancerSpecific_Act, bkgd_cancerSpecific_Act, gs_cancerSpecific_Act, eSet_cancerSpecific_Act)+[' '.join([';'.join([str(j) for j in i]) for i in ovlp_cancerSpecific_Act])])
    tmp1 = [i for i in gs_cancerSpecific_Act if i in eSet_cancerSpecific_LoF]
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (Act)','LossOfFunction']+enrichment(tmp1, bkgd_cancerSpecific_Act, gs_cancerSpecific_LoF, eSet_cancerSpecific_Act)+[' '.join([';'.join([str(j) for j in i]) for i in tmp1])])
    csp_hypgeoResults.append([analysis,'OncogeneDB','Activating']+enrichment(ovlp_cancerSpecific_og, bkgd_cancerSpecific_og, gs_cancerSpecific_og, eSet_cancerSpecific_og)+[' '.join([str(i) for i in ovlp_cancerSpecific_og])])
    csp_hypgeoResults.append([analysis,'OncogeneDB','LossOfFunction']+enrichment([i for i in gs_cancerSpecific_ts if i in eSet_cancerSpecific_og], bkgd_cancerSpecific_og, gs_cancerSpecific_ts, eSet_cancerSpecific_og)+[' '.join([str(i) for i in ovlp_cancerSpecific_ts])])
    tmp1 = [i for i in gs_cancerSpecific_Act if i in eSet_cancerSpecific_LoF] 
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (LoF)','Activating']+enrichment(tmp1, bkgd_cancerSpecific_LoF, gs_cancerSpecific_Act, eSet_cancerSpecific_LoF)+[' '.join([';'.join([str(j) for j in i]) for i in [i for i in gs_cancerSpecific_Act if i in tmp1]])])
    csp_hypgeoResults.append([analysis,'TCGA_Consensus (LoF)','LossOfFunction']+enrichment(ovlp_cancerSpecific_LoF, bkgd_cancerSpecific_LoF, gs_cancerSpecific_LoF, eSet_cancerSpecific_LoF)+[' '.join([';'.join([str(j) for j in i]) for i in ovlp_cancerSpecific_LoF])])
    #tmp1 = [i for i in gs_cancerSpecific_og if i in eSet_cancerSpecific_ts]
    #csp_hypgeoResults.append([analysis,'TSGeneDB','Activating']+enrichment(tmp1, bkgd_cancerSpecific_ts, gs_cancerSpecific_og, eSet_cancerSpecific_ts)+[' '.join([str(i) for i in tmp1])])
    #csp_hypgeoResults.append([analysis,'TSGeneDB','LossOfFunction']+enrichment(ovlp_cancerSpecific_ts, bkgd_cancerSpecific_ts, gs_cancerSpecific_ts, eSet_cancerSpecific_ts)+[' '.join([str(i) for i in ovlp_cancerSpecific_ts])])
    #csp_hypgeoResults.append([analysis,'Aggregate','All']+enrichment(ovlp_cancerSpecific_Agg, bkgd_cancerSpecific_Agg, gs_cancerSpecific_Agg, eSet_cancerSpecific_Agg)+[' '.join([';'.join([str(j) for j in i]) for i in ovlp_cancerSpecific_Agg])])
    csp_hypgeoResults.append([analysis,'All','All']+enrichment(ovlp_cancerSpecific_All, bkgd_cancerSpecific_All, gs_cancerSpecific_All, eSet_cancerSpecific_All)+[' '.join([str(i) for i in ovlp_cancerSpecific_All])])
    



with open('csp_hypgeoResults.csv','w') as outFile:
    outFile.write('Analysis,Enrichment_Source,OncoMerge_Direction,Overlap,OncoMerge_Discovered,Enrichment_Set,Background,P-Value,Combos\n')
    outFile.write('\n'.join([','.join([str(j) for j in i]) for i in csp_hypgeoResults]))

