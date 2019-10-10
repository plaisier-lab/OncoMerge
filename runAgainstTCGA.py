from subprocess import *
import os

# Run for each tumor type in TCGA
for tumor in ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']:
    # Make the output directory if it doesn't exists already
    if not os.path.exists('TCGA/output/'+tumor):
        os.mkdir('TCGA/output/'+tumor)

    # Run first run with no filters
    print('Running no_filter...')
    cmd1 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-cfp TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-smp TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/no_filter',
            '-pq 1',
            '-mcr 0',
            '-sp']
    proc1 = Popen(' '.join(cmd1), shell=True, stdin=PIPE)
    out = proc1.communicate()

    # Run first run with only permuted q-value
    print('\nRunning pq...')
    cmd2 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-cfp TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-smp TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/pq',
            '-pq 0.1',
            '-mcr 0',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc2 = Popen(' '.join(cmd2), shell=True, stdin=PIPE)
    out = proc2.communicate()

    # Run first run with no filters
    print('\nRunning mcr...')
    cmd3 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-cfp TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-smp TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/mcr',
            '-pq 1',
            '-mcr 0.001',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc3 = Popen(' '.join(cmd3), shell=True, stdin=PIPE)
    out = proc3.communicate()

    # Run first run with no filters
    print('\nRunning pq_mcr...')
    cmd4 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-cfp TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-smp TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/pq_mcr',
            '-pq 0.1',
            '-mcr 0.001',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc4 = Popen(' '.join(cmd4), shell=True, stdin=PIPE)
    out = proc4.communicate()

