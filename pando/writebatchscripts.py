import pandas as pd
from os import system
import glob


script_dir = '/Users/anbr3575/LRTAnalysis/pando/scripts/'
output_pattern = '/Users/anbr3575/LRTAnalysis/output/%s.csv'
exe_pattern = 'python /Users/anbr3575/LRTAnalysis/pando/pipelinepando.py %s %s\n' #exe = execute

analysis = pd.read_pickle('/Users/anbr3575/LRTAnalysis/analysis/analysis.p')
subanalysis = analysis.query('Graph_order==6')
query = 'order6'
degseqdp = '/Users/anbr3575/LRTAnalysis/degreesequences/'
# degseqV = glob.glob(degseqdp+'*.txt')
degseqV = [degseqdp+fn for fn in list(subanalysis.index)]
for i, fp in enumerate(degseqV):
    input_file = fp
    output = open(script_dir + str(i) + '.sh', 'w')
    output.write('#!/bin/bash\n')
    fn = fp.split('/')[-1]
    output_file = output_pattern %(fn+query)
    output.write(exe_pattern %(input_file, output_file))
    output.close()
"""
Write launch scripts
"""

def writelaunch(minscripts, maxscripts, filename):
    """
     Write launch script for some portion of the data
    """
    output = open(script_dir + 'launch%s.sh' %(filename), 'w')
    launch = \
    """#!/bin/bash
    #PBS -N lrtanalysis%s
    #PBS -joe
    #PBS -t %d-%d%s
    #PBS -q long8gb
    #PBS -l pmem=8gb
    #PBS -l nodes=1:ppn=1

    /Users/anbr3575/LRTAnalysis/pando/scripts/${PBS_ARRAYID}.sh\n""" % (minscripts,maxscripts, '%50')
    output.write(launch)
    output.close()

num_scripts = len(glob.glob(script_dir+'*.sh')) - len(glob.glob(script_dir+'launch%s*' %(query)))
batch_size = 1000
loop_batches = num_scripts/batch_size
if num_scripts == batch_size*loop_batches:
    tot_batches = loop_batches
else:
    tot_batches = loop_batches+1

for i in range(loop_batches):
    minscripts = i*batch_size
    maxscripts = (i+1)*batch_size-1
    filename = "%s_%d_of_%d" %(query, i+1, tot_batches)
    writelaunch(minscripts,maxscripts, filename)

if loop_batches < tot_batches:
    minscripts = batch_size*loop_batches
    maxscripts = num_scripts-1
    filename = "%s_%d_of_%d" %(query, tot_batches, tot_batches)
    writelaunch(minscripts,maxscripts, filename)

# Make everything executable
system('chmod u+x %s/*' % script_dir)
