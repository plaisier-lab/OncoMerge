from subprocess import *
import os
import multiprocessing as mp

tumors = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
#tumors = ['ACC']

# Run for each tumor type in TCGA
def runOncoMerge(tumor):
    print(tumor)
    # Make the output directory if it doesn't exists already
    if not os.path.exists('TCGA/output/'+tumor):
        os.mkdir('TCGA/output/'+tumor)

    # Run first run with no filters
    print('Running no_filter...')
    cmd1 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-pq 1',
            '-mlg 10000',
            '-sp',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
    print(' '.join(cmd1))
    proc1 = Popen(' '.join(cmd1), shell=True, stdout=PIPE, stderr=PIPE)
    out = proc1.communicate()

    # Run with only permuted q-value
    print('\nRunning pq...')
    cmd2 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq',
            '-pq 0.1',
            '-mlg 10000',
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
    print(' '.join(cmd2))
    proc2 = Popen(' '.join(cmd2), shell=True, stdout=PIPE, stderr=PIPE)
    out = proc2.communicate()

    # Run with maximum frequency filter
    print('\nRunning mff...')
    cmd3 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/mff',
            '-pq 1',
            '-mlg 10',
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
    print(' '.join(cmd3))
    proc3 = Popen(' '.join(cmd3), shell=True, stdout=PIPE, stderr=PIPE)
    out = proc3.communicate()

    # Run with both filters
    print('\nRunning pq_mff...')
    cmd4 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq_mff',
            '-pq 0.1',
            '-mlg 10',
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
    print(' '.join(cmd4))
    proc4 = Popen(' '.join(cmd4), shell=True, stdout=PIPE, stderr=PIPE)
    out = proc4.communicate()

    # Both filters + GISTIC threshold shallow
    print('\nRunning pq_mff (GISTIC shallow)...')
    cmd5 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq_mff_gt1',
            '-pq 0.1',
            '-mlg 10',
            '-gt 1',
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
    proc5 = Popen(' '.join(cmd5), shell=True, stdout=PIPE, stderr=PIPE)
    out = proc5.communicate()

    # Both filters PQ options
    for pq in [0.01, 0.05, 0.2]: # 0.1 is default
        print('\nRunning pq_mff (pq theshold test) '+str(pq)+' ...')
        cmd6 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq_mff_pq'+str(pq).replace('\.','_'),
            '-pq '+str(pq),
            '-mlg 10',
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
        proc6 = Popen(' '.join(cmd6), shell=True, stdout=PIPE, stderr=PIPE)
        out = proc6.communicate()

    # Both filters MMF options
    for mmf in [0.01, 0.1, 0.2]: # 0.05 is default
        print('\nRunning pq_mff (pq theshold test) '+str(pq)+' ...')
        cmd7 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq_mff_mmf'+str(mmf).replace('\.','_'),
            '-pq 0.1',
            '-mlg 10',
            '-mmf '+str(mmf),
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
        proc7 = Popen(' '.join(cmd7), shell=True, stdout=PIPE, stderr=PIPE)
        out = proc7.communicate()

    # Both filters MLG options
    for mlg in [1, 5, 15, 20]: # 10 is default
        print('\nRunning pq_mff (mlg test) '+str(mlg)+' ...')
        cmd8 = ['python3 oncoMerge.py',
            '-gp TCGA/GISTIC/'+tumor,
            '-aaf TCGA/OncoMerge_input_g2e_converter.csv',
            '-ln "Locus ID"',
            '-pam TCGA/PAM/'+tumor+'_somMutMC3.csv',
            '-mscv TCGA/MutSig2cv/'+tumor+'_sig2cv.csv',
            '-fus TCGA/FUSIONS/'+tumor+'_fusions.csv',
            '-op TCGA/output_11_4_2022/'+tumor+'/pq_mff_mlg'+str(mlg).replace('\.','_'),
            '-pq 0.1',
            '-mlg '+str(mlg),
            '-lp TCGA/output_11_4_2022/'+tumor+'/no_filter',
            '-tcga True',
            '-bl TCGA/blacklist/blacklist_29850653_29625053.csv']
        proc8 = Popen(' '.join(cmd8), shell=True, stdout=PIPE, stderr=PIPE)
        out = proc8.communicate()

# Parallelize
pool = mp.Pool(15)
pool.map(runOncoMerge, tumors)
pool.close()


