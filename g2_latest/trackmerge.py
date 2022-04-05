#!/usr/bin/env python

"""
Code to read in data from two text files, subract/add the values,
and write out to another file for plotting with gnuplot.
Part of analysis of along track motion for SASAR 2000-02-11.
"""
import sys, string

print '\n--------------------'
print 'Prog: trackmerge.py - Ver 0.1 (jmh)\n'
if len(sys.argv) < 4:
    print 'USAGE: trackmerge.py [track1] [track2] [OutFile]'
    sys.exit()

#Misc
Trk1FileName = sys.argv[1]
Trk2FileName = sys.argv[2]
OutFileName = sys.argv[3]

#Open files
Trk1File = open(Trk1FileName,'r')
Trk2File = open(Trk2FileName,'r')
OutFile = open(OutFileName,'w')

#Read data from the two track files
while 1:  # Repeat for all lines in the track files
    field1 = string.split(Trk1File.readline(),' ')
    if not field1 or len(field1) != 2:
        break
    field2 = string.split(Trk2File.readline(),' ')
    if not field2 or len(field2) != 2:
        break

    pri = string.atoi(field1[0])
    diff = string.atof(field1[1]) - string.atof(field2[1])

    OutFile.write(`pri`+' '+`diff`+'\n')
    print 'PRI: '+`pri`+' / diff: '+`diff`

Trk1File.close()
Trk2File.close()
OutFile.close()
	