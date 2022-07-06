##########################################################
## OncoMerge:  runAgainstTCGA.py                        ##
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


from subprocess import *
import os

# Run for each tumor type in TCGA
for tumor in ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']:
    print(tumor)
    # Make the output directory if it doesn't exists already
    if not os.path.exists('TCGA/output_5_4_2022/'+tumor):
        os.mkdir('TCGA/output_5_4_2022/'+tumor)

    # Run first run with no filters
    print('Running no_filter...')
    cmd1 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_5_4_2022/'+tumor+'/no_filter',
            '-pq 1',
            '-mlg 10000',
            '-mpf 0',
            '-sp',
            '-tcga True',
            '-bl TCGA/blocklist/blocklist_29850653_29625053.csv']
    proc1 = Popen(' '.join(cmd1), shell=True, stdin=PIPE)
    out = proc1.communicate()

    # Run with only permuted q-value
    print('\nRunning pq...')
    cmd2 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_5_4_2022/'+tumor+'/pq',
            '-pq 0.1',
            '-mlg 10000',
            '-mpf 0',
            '-lp TCGA/output_5_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blocklist/blocklist_29850653_29625053.csv']
    proc2 = Popen(' '.join(cmd2), shell=True, stdin=PIPE)
    out = proc2.communicate()

    # Run with maximum frequency filter
    print('\nRunning mff...')
    cmd3 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_5_4_2022/'+tumor+'/mff',
            '-pq 1',
            '-mlg 10',
            '-mpf 0',
            '-lp TCGA/output_5_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blocklist/blocklist_29850653_29625053.csv']
    proc3 = Popen(' '.join(cmd3), shell=True, stdin=PIPE)
    out = proc3.communicate()

    # Run with both filters
    print('\nRunning pq_mff...')
    cmd4 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_5_4_2022/'+tumor+'/pq_mff',
            '-pq 0.1',
            '-mlg 10',
            '-mpf 0',
            '-lp TCGA/output_5_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blocklist/blocklist_29850653_29625053.csv']
    proc4 = Popen(' '.join(cmd4), shell=True, stdin=PIPE)
    out = proc4.communicate()

