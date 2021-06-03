from subprocess import *

# Run first run with no filters
print('Running no_filter...')
proc1 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -sp -pq 1 -mcr 0 -op Bueno/no_filter', shell=True, stdin=PIPE)
out = proc1.communicate()

# Run first run with only permuted q-value
print('\nRunning pq...')
proc2 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -pq 0.1 -mcr 0 -op Bueno/pq', shell=True, stdin=PIPE)
out = proc2.communicate()

# Run first run with no filters
print('\nRunning mcr...')
proc3 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -pq 1 -mcr 0.001 -op Bueno/mcr', shell=True, stdin=PIPE)
out = proc3.communicate()

# Run first run with no filters
print('\nRunning pq_mcr...')
proc4 = Popen('python oncoMerge.py -cf Bueno/BUENO_config_deep.json -lp Bueno/no_filter -pq 0.1 -mcr 0.001 -op Bueno/pq_mcr', shell=True, stdin=PIPE)
out = proc4.communicate()

