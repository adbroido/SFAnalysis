#!/bin/bash
#PBS -N pandotest
#PBS -joe
#PBS -t 0-1000%20
#PBS -q long8gb
#PBS -l pmem=8gb
#PBS -l nodes=1:ppn=1

/Users/anbr3575/LRTAnalysis/pando/scripts/${PBS_ARRAYID}.sh
