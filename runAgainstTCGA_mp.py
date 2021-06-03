from subprocess import *
import os
import multiprocessing as mp

#tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
tumors = ['SKCM']

# Run for each tumor type in TCGA
def runOncoMerge(tumor):
    print(tumor)
    # Make the output directory if it doesn't exists already
    if not os.path.exists('TCGA/output/'+tumor):
        os.mkdir('TCGA/output/'+tumor)
    
    # Run first run with no filters
    print('Running no_filter...')
    cmd1 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/no_filter',
            '-pq 1',
            '-mlg 10000',
            '-mpf 0',
            '-sp']
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
            '-op TCGA/output/'+tumor+'/pq',
            '-pq 0.1',
            '-mlg 10000',
            '-mpf 0',
            '-lp TCGA/output/'+tumor+'/no_filter']
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
            '-op TCGA/output/'+tumor+'/mff',
            '-pq 1',
            '-mlg 10',
            '-mpf 0',
            '-lp TCGA/output/'+tumor+'/no_filter']
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
            '-op TCGA/output/'+tumor+'/pq_mff',
            '-pq 0.1',
            '-mlg 10',
            '-mpf 0',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc4 = Popen(' '.join(cmd4), shell=True, stdin=PIPE)
    out = proc4.communicate()

    '''
    # Run with minimum PAM frequency filter
    print('\nRunning mpf...')
    cmd5 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/mpf',
            '-pq 1',
            '-mlg 10000',
            '-mpf 0.01',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc5 = Popen(' '.join(cmd5), shell=True, stdin=PIPE)
    out = proc5.communicate()

    # Run with both filters
    print('\nRunning pq_mpf...')
    cmd6 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/pq_mpf',
            '-pq 0.1',
            '-mlg 10000',
            '-mpf 0.01',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc6 = Popen(' '.join(cmd6), shell=True, stdin=PIPE)
    out = proc6.communicate()
    
    # Run with both filters
    print('\nRunning pq_mff_mpf...')
    cmd7 = ['python oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-op TCGA/output/'+tumor+'/pq_mff_mpf',
            '-pq 0.1',
            '-mlg 10',
            '-mpf 0.01',
            '-lp TCGA/output/'+tumor+'/no_filter']
    proc7 = Popen(' '.join(cmd7), shell=True, stdin=PIPE)
    out = proc7.communicate()
    '''

# Parallelize
pool = mp.Pool(3)
pool.map(runOncoMerge, tumors)
pool.close()


