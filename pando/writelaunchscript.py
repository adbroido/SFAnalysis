import pandas as pd
from os import system
import glob


def writelaunch(minscripts, maxscripts, filename):
    """
     Write launch script for some portion of the data
    """
    output = open(script_dir + 'launch%s.sh' %(filename), 'w')
    launch = \
    """#!/bin/bash
    #PBS -N pandotest
    #PBS -joe
    #PBS -t %d-%d%s
    #PBS -q long8gb
    #PBS -l pmem=8gb
    #PBS -l nodes=1:ppn=1

    /Users/anbr3575/LRTAnalysis/pando/scripts/${PBS_ARRAYID}.sh\n""" % (minscripts,maxscripts, '%20')
    output.write(launch)
    output.close()

script_dir = '/Users/anbr3575/LRTAnalysis/pando/scripts/'
script_dir = '/Users/annabroido/Dropbox/Research/LRTAnalysis/LRTAnalysis/pando/scripts/'
query = '23'
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
    filename = "%s_%d_of_%d" %(query, i, tot_batches)
    writelaunch(minscripts,maxscripts, filename)

if loop_batches < tot_batches:
    minscripts = batch_size*loop_batches
    maxscripts = num_scripts-1
    filename = "%s_%d_of_%d" %(query, tot_batches, tot_batches)
    writelaunch(minscripts,maxscripts, filename)

# Make everything executable
system('chmod u+x %s/*' % script_dir)
