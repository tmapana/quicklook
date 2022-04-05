#!/usr/bin/python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2unpk_dgps.py
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
File to parse DGPS file from OmniStar system (used in SASAR) 
Extracts only the part of the DGPS data which spans the data take
as defined by the unpacked IMU file. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-08-02
COMMENTS:
Initial version, modified from dgps_g2unpk.py

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-13
COMMENTS:
Cleaned up a bit and made more efficient. Replaced tabs with spaces.

=========================================================
"""

import sys, string

print '\n--------------------'
print 'Prog: g2unpk_dgps.py - Ver 0.2 (jmh)\n'
if len(sys.argv) < 5:
    print 'USAGE: g2unpk_dgps-0.2.py [DGPS_File] [IMU_unpkFile] [OutFile] '\
          '[LBR_PRI_Incr] [PRF_Hz]'
    sys.exit()

# Misc
DGPSFileName = sys.argv[1]
IMUFileName = sys.argv[2]
OutFileG2Name = sys.argv[3]
PRI_Incr = string.atoi(sys.argv[4])			#LBR PRI increment
PRF_Hz = string.atof(sys.argv[5])
PRI_secs = 1.0 / PRF_Hz #secs

# Open files     
DGPSfile = open(DGPSFileName, 'r')
print 'Opened input DGPS file: ' + DGPSFileName
IMUfile = open(IMUFileName, 'r')
print 'Opened input unpacked LBR file: '+ IMUFileName
OutFileG2 = open(OutFileG2Name,'w')
print 'Opened output G2-unpacked DGPS file: ' + OutFileG2Name  

# Read start and end LBR PRIID and UTC time from IMU records file
IMUfile.readline()   # skip header line
field = string.split(IMUfile.readline())
startLBR_PRI = string.atoi(field[0])
startLBR_UTC = string.atof(field[1])
print '\nUnpacked LBR input data:'
print 'Start: LBR PRI= '+`startLBR_PRI`+' / UTC= '+`startLBR_UTC`+' secs'
while 1: # find last record
    field = string.split(IMUfile.readline())
    if not field: break
    lastfield = field
endLBR_PRI = string.atoi(lastfield[0])
endLBR_UTC = string.atof(lastfield[1])	    
print 'End  : LBR PRI= '+`endLBR_PRI`+' / UTC= '+`endLBR_UTC`+' secs\n'
IMUfile.close()

### PROCESS RELEVANT LINES OF DGPS FILE ### 
lines = 0L
addDay = 0.0
startFound = 0      # of DGPS processing block
endFound = 0	    # of DGPS processing block
while 1:  # Repeat for all lines in the DGPS file
    field = string.split(DGPSfile.readline(),',')
    if not field or len(field) != 15:
        break
	
    # find UTC time	
    UTCHr = string.atoi(field[1][0:2])
    UTCMin = string.atoi(field[1][2:4])
    UTCSec = string.atof(field[1][4:])
    timeSecs = UTCHr*3600+UTCMin*60+UTCSec

    # Find start params
    if not startFound and (timeSecs >= startLBR_UTC): 
        startLBR_PRI_DGPS = PRI_Incr * int(startLBR_PRI/PRI_Incr +  \
                            (timeSecs-startLBR_UTC)*PRF_Hz)			
        startTimeDGPS = timeSecs
        startFound = 1   # only enters here once 

    # Process only if start found 
    if startFound:
        if timeSecs < startTimeDGPS: # if data straddles 00h00
            addDay = 86400.0         # num secs in day	
        LBR_PRI_DGPS = PRI_Incr * int(startLBR_PRI_DGPS/PRI_Incr + \
                       (timeSecs+addDay-startTimeDGPS)*PRF_Hz)					

        # if in range, perform calcs and write to file
        if LBR_PRI_DGPS < endLBR_PRI:

            #Find Latitude (deg, min)
            if field[3] == 'S': 
                sign = -1
            else:
                sign  = 1  
            dotIndex = string.find(field[2],'.')
            LatDeg = sign*string.atoi(field[2][0:dotIndex-2])
            LatMin = sign*string.atof(field[2][dotIndex-2:])
            LatInDeg = LatDeg+LatMin/60.0

            #Find Longitude (deg, min)
            if field[5] == 'W': 
                sign = -1
            else:
                sign  = 1  
            dotIndex = string.find(field[4],'.')
            LongDeg = sign*string.atoi(field[4][0:dotIndex-2])
            LongMin = sign*string.atof(field[4][dotIndex-2:])
            LongInDeg = LongDeg+LongMin/60.0

            #Find Altitude (AMSL)
            Alt = string.atof(field[9])

            # write to output file
            OutFileG2.write(`LBR_PRI_DGPS`+' '+`timeSecs`+' '+`LatInDeg`+' '+ \
                            `LongInDeg`+' '+`Alt`+'\n')		 
            lines = lines + 1L
        else:
            break		# if we are past end of data take
	
print 'DGPS records processed: ' + `lines`
		
# tidy up (returns zero on success)
DGPSfile.close()
OutFileG2.close()

sys.exit(0)