import numpy as np
from os import system

patients = ['CM031115',
            'CS033015',
            'DC021015',
            'EC012815',
            'EP010815',
            'JC032415',
            'JO041315',
            'LA041615',
            'MD012815',
            'MH021115',
            'MM040615',
            'PD021315',
            'WS020415']

script_dir = '/Users/sawa6416/Projects/Hearts/rewiring/scripts/'
output_pattern = '/Users/sawa6416/Projects/Hearts/rewiring/output/%s.csv'
input_pattern = '/Users/sawa6416/Projects/Hearts/HumanAF/%s/%s/KSG_MI.mat'
exe_pattern = '/Users/sawa6416/Projects/Hearts/rewiring.py -i %s -o %s -s %f -m %s -n %d\n'
sd_steps = [0, 0.5, 1.0]
methods = ['product', 'eig', 'high', 'random']
num_iters = [1, 1, 1, 100]

i = 0
for atrium in ['RA1', 'LA1']:
    for tag in patients:
        for s_ind, sd in enumerate(sd_steps):
            for m, method in enumerate(methods):
                output = open(script_dir + str(i) + '.sh', 'w')
                output.write('#!/bin/bash\n')
                input_file = input_pattern % (tag, atrium)
                for n in xrange(num_iters[m]):
                    output_file = output_pattern % ('_'.join([tag, atrium, str(s_ind), method, str(n)]))
                    output.write(exe_pattern % (input_file, output_file, sd, method, n))
                output.close()
                i += 1

output = open(script_dir + 'launch.sh', 'w')
