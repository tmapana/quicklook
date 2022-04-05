#!/usr/bin/env python

# Script to start Python bytecode version of G2 processor.
# Code: j.m.horrell 2001-02-13

# import python modules required for this file
import os, sys

# path to the python bytecode processor files (ends in slash)
G2path = '/usr/local/g2/'

# construct command to run the processor 
cmd = 'python '+G2path+'g2main.pyc'
for i in range(1,len(sys.argv)):
    cmd = cmd+' '+sys.argv[i]

# shell out and run processor
os.system(cmd)
