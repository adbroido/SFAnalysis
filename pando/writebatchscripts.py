import pandas as pd
from os import system
import glob


script_dir = '/Users/anbr3575/LRTAnalysis/pando/scripts/'
output_pattern = '/Users/anbr3575/LRTAnalysis/output/%s.csv'
exe_pattern = 'python /Users/anbr3575/LRTAnalysis/pipelinepando.py %s %s\n' #exe = execute

analysis = pd.read_pickle('analysis.p')
subanalysis = analysis.query('ppl>0.1 & ntail>50 & Graph_order<4')
queryinfo = ''
degseqdp = '/Users/annabroido/Dropbox/Research/LRTAnalysis/degreesequences/'
# degseqV = glob.glob(degseqdp+'*.txt')
degseqV = [degseqdp+fn for fn in list(subanalysis.index)]
for i, fp in enumerate(degseqV):
    input_file = fp
    output = open(script_dir + str(i) + '.sh', 'w')
    output.write('#!/bin/bash\n')
    fn = fp.split('/')[-1]
    output_file = output_pattern %(fn+queryinfo)
    output.write(exe_pattern %(input_file, output_file))
    output.close()
