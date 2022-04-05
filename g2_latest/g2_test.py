#!/usr/bin/env python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group (1994-1999)
FILE NAME: g2.py
CODE CONTROLLER: Jasper Horrell (jasper@eng.uct.ac.za)
DESCRIPTION:
The G2 SAR processor top-level controlling module. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-07-30
COMMENTS:
Initial python version.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-04
COMMENTS:
Calc FFT sizes automagically. Minor usage changes.

VERSION/AUTHOR/DATE : 0.3 / Jasper Horrell / 1999-08-14
COMMENTS:
Mocomp IMU/DGPS data merge functionality added.

VERSION/AUTHOR/DATE : 0.4 / Jasper Horrell / 1999-08-19
COMMENTS:
Add mocomp LBR PRI offset.

VERSION/AUTHOR/DATE : 0.5 / Jasper Horrell / 1999-08-27
COMMENTS:
Moved over to cleaned up motion comp calc (mocomp ver 0.5).

VERSION/AUTHOR/DATE : 0.6 / Jasper Horrell / 1999-09-13
COMMENTS:
Allow specification of interference suppresion parameters
(no actual changes in this file). Add in start and end
rng bin in call to mocomp for geocoding.

VERSION/AUTHOR/DATE : 0.7 / Jasper Horrell / 1999-10-01
COMMENTS:
Allow notch interference suppression. Allow output endian and
orientation specification.

VERSION/AUTHOR/DATE : test / Jasper Horrell / 2000-01-18
COMMENTS:
A test version which calls g2dgps_smooth.py rather than
g2mocfilt.py.
    

=========================================================
"""

import string, sys, os, g2tmpl, g2tools, Gnuplot, time

print """
=========================
G2 Airborne SAR Processor (Ver. 0.6)
Code: J.M. Horrell (C) University of Cape Town (1994-1999)
"""

if (len(sys.argv) != 2):
    print 'USAGE: g2.py [proc cfg file]'
    print '(To create template radar and proc cfg files, type: "g2.py -tmpl")'
    sys.exit()
elif len(sys.argv) == 2 and sys.argv[1] == '-tmpl':
    g2tmpl.createRadarCfgTmpl('Radar.CFG')
    print 'Radar template config file "Radar.CFG" written!'
    g2tmpl.createProcCfgTmpl('Proc.CFG')
    print 'Processor template config file "Proc.CFG" written!'
    sys.exit()

def ErrorExit(cause):
    print '\nERROR - '+cause+'!!'
    print 'Check relevant log file for more information on error.'
    print '\nAbnormal exit from G2 Processor!!\n'
    sys.exit(1)

# Set up the processor configuration dictionary 
G2 = g2tools.G2Config()    # create a new proc config object
if G2.addCfgFile(sys.argv[1]) != 0: # add contents of proc config file to object 
    ErrorExit('in adding contents of '+sys.argv[1])
if G2.addCfgFile(G2.cfg['$RadarCfgFile']) != 0: # add contents of radar config file
    ErrorExit('in adding contents of '+G2.cfg['$RadarCfgFile'])
if G2.inferParams() !=0:          # infer extra parameters
    ErrorExit('in inferring processing parameters')
    
# Set up some file names
IDRoot =  G2.cfg['$RunID']
UnpkIMU = IDRoot+'_unpkIMU'
UnpkDGPS = IDRoot+'_unpkDGPS'
MocMerge = IDRoot+'_IMU+DGPS'
LBRFile = G2.cfg['$LBRFile']
DCFile = IDRoot+'_sniffDC'
MocLog = IDRoot+'_mlog'
MocFile = IDRoot+'_moc'
G2.cfg['$MocFile'] = MocFile
MocTxtFile = MocFile+'txt'
RncFile = IDRoot+'.rnc'
G2.cfg['$RncFile'] = RncFile
CorFile = IDRoot+'.cor'
G2.cfg['$CorFile'] = CorFile
RngCmd = IDRoot+'_rngcmd'
AzCmd = IDRoot+'_azcmd'
RngLog = IDRoot+'_rnglog'
G2.cfg['$RngLog'] = RngLog
AzLog = IDRoot+'_azlog'
G2.cfg['$AzLog'] = AzLog
ByteFile = IDRoot+'.byte'   # used in float to byte conversion
MagFile = IDRoot+'.iqmag'
ImLog = IDRoot+'.imlog'
TiffFile = IDRoot+'.tif'
if G2.cfg['$DetectMethod'] == 'cmplx':
    AzcFile = IDRoot+'.slc'
elif G2.cfg['$DetectMethod'] == 'mag':    
    AzcFile = IDRoot+'.mag'
elif G2.cfg['$DetectMethod'] == 'pow':
    AzcFile = IDRoot+'.pow'
elif G2.cfg['$DetectMethod'] == 'powdB':
    AzcFile = IDRoot+'.powdB'
else:
    print 'ERROR - unknown detect method '+G2.cfg['$DetectMethod']
    ErrorExit('in setting up file names')   
G2.cfg['$AzcFile'] = AzcFile


#---------------------------------
# Unpack all IMU records to ASCII
if G2.cfg['$EnableUnpackIMU'] == 'y':
    cmd = 'imu_unpack '+LBRFile+' OutF='+UnpkIMU
    if os.system(cmd) != 0:
        ErrorExit('in IMU unpack program') 

#----------------------------
# Unpack DGPS Motion Records corresponding to the LBR PRI ID range
if G2.cfg['$EnableUnpackDGPS'] == 'y':
    cmd = 'g2unpk_dgps.py '+G2.cfg['$DGPSFile']+' '+UnpkIMU+' '+\
          UnpkDGPS+' '+G2.cfg['$PRIIncrLBR']+' '+G2.cfg['$RadarPRF']
    if os.system(cmd) != 0:
        ErrorExit('in DGPS unpack program')

    
#----------------------------
# Merge IMU and DGPS data
if G2.cfg['$EnableMergeMocData'] == 'y':
    cmd = 'g2dgps_smooth.py '+UnpkDGPS+' '+UnpkIMU+' '+MocMerge+' '+\
          G2.cfg['$MocompLBRPRIOffset']
    if os.system(cmd) != 0:
        ErrorExit('in Mocomp DGPS and IMU merge program')

#-----------------------------------------------------------
# Run mocomp range shift calc prog. This also provides ave
# ground speed, needed for az compression.
if G2.cfg['$EnableMocompCalc'] == 'y':
    cmd = 'mocomp MotF='+MocMerge+' LogF='+MocLog+\
          ' OutF='+MocFile+' OutTxtF='+MocTxtFile+\
          ' Delay0='+G2.cfg['$RadarDelayToStartSample']+\
          ' RefBin='+G2.cfg['$MocompRefRngBin']+\
          ' A2DFreq='+G2.cfg['$RadarA2DFreq']+\
          ' PRF='+G2.cfg['$RadarPRF']+\
          ' DataStartLP='+G2.cfg['$DataStartLP']+\
          ' ProcStartLP='+G2.cfg['$ProcStartLP']+\
          ' ProcEndLP='+G2.cfg['$ProcEndLP']+\
          ' LPIncr='+G2.cfg['$PRIIncrLBR']+\
          ' TerAlt='+G2.cfg['$TerrainAlt']+\
          ' NearRngBin='+G2.cfg['$StartRngBinToProc']+\
          ' FarRngBin='+G2.cfg['$EndRngBinToProc']
    if os.system(cmd) != 0:
        ErrorExit('in motion compensation calculation program')

if G2.cfg['$EnablePlotMotionError'] == 'y':
    g = Gnuplot.Gnuplot()
    g.xlabel('G2 PRI ID')
    g.ylabel('Motion Error (m)')
    data = Gnuplot.File(MocTxtFile,with='lines')
    g.plot(data)
#    g.hardcopy(filename=MocTxtFile+'.ps',color=1)
    time.sleep(5.0)
    del g,data
#    print '(Motion error plot written to postscipt file)'
    
#-----------------------------
# Find DC offsets and I/Q ratio
if G2.cfg['$EnableSniffDC'] == 'y':
    SDCRowsToProc = int(0.2*int(G2.cfg['$InputPRIsToProc']))
    SDCStartCol = int(0.25*int(G2.cfg['$RngBinsToProc']))
    SDCColsToProc = int(0.5*int(G2.cfg['$RngBinsToProc']))
    cmd = 'sniffdc '+G2.cfg['$RawDataFile']+\
          ' '+G2.cfg['$RadarAzSamples']+' '+G2.cfg['$RadarRngBins']+\
          ' StartRow='+G2.cfg['$StartG2PRI']+' RowsToProc='+`SDCRowsToProc`+\
          ' StartCol='+`SDCStartCol`+' ColsToProc='+`SDCColsToProc`+\
          ' LogF='+DCFile
    if os.system(cmd) != 0:
        ErrorExit('in sniff DC offset program')

#--------------------------------
# Read DC offsets, I/Q ratio and ground speed from the log files generated 
if G2.cfg['$EnableParseFromLogs'] == 'y':
    if G2.parseParam(DCFile,'Average I value','$DCOffsetI') != 0:
        ErrorExit('in parsing DC offset I from log')
    if G2.parseParam(DCFile,'Average Q value','$DCOffsetQ') != 0:
        ErrorExit('in parsing DC offset Q from log')
    if G2.parseParam(DCFile,'Average I/Q ratio','$IQRatio') != 0:
        ErrorExit('in parsing IQ ratio from log')
    if G2.parseParam(MocLog,'Proc average ground speed',\
              '$ProcAveGrndSpeed') != 0:
        ErrorExit('in parsing average ground speed from log')      

#-----------------------
# Range Processing Stage
if G2.cfg['$EnableRngProc'] == 'y':
    G2.createRngProcCmdFile(RngCmd)
    if os.system('rngcom '+RngCmd) != 0:
        ErrorExit('in range compression program')
  
#------------
# Corner Turn
if G2.cfg['$EnableCornerTurn'] == 'y':
    cmd = 'corner '+RncFile+' '+CorFile+' '+G2.cfg['$RngBinsToProc']+' 0 '+\
          G2.cfg['$PreSummedPRIsToProcLessOne']+' 8 mem=50'
    if os.system(cmd) != 0:
        ErrorExit('in corner turn program')     

#-------------------------
# Azimuth Processing Stage
if G2.cfg['$EnableAzProc'] == 'y':
    G2.createAzProcCmdFile(AzCmd)
    if os.system('azcom '+AzCmd) != 0:
        ErrorExit('in azimuth compression program')
    
#------------------
# TIFF Output Stage
if G2.cfg['$EnableFloat2Tiff'] == 'y':
    if G2.cfg['$DetectMethod'] == 'cmplx':  # first run iq2mag
        cmd = 'iq2mag '+AzcFile+' '+MagFile+' '+G2.cfg['$RngBinsToProc']+' '+\
              G2.cfg['$AzProcOutPRIs']+' pow float' 
        if os.system(cmd) != 0:
            ErrorExit('in IQ to mag conversion program')
        flt2byteInFile = MagFile
    else:                            # data not complex, so skip iq2mag       
        flt2byteInFile = AzcFile

    cmd = 'flt2byte '+flt2byteInFile+' '+ByteFile+' '+G2.cfg['$RngBinsToProc']+\
          ' '+G2.cfg['$AzProcOutPRIs']+' '+G2.cfg['$Float2ByteScale']+' 0'
    if os.system(cmd) != 0:      
        ErrorExit('in float to byte conversion program')

    cmd = 'b2tif '+ByteFile+' '+TiffFile+' '+G2.cfg['$RngBinsToProc']+' '+\
           G2.cfg['$AzProcOutPRIs']+' orient=3'
    if os.system(cmd) != 0:
        ErrorExit('in byte to tiff conversion program')

#---------------
# Swap endian order for big endian output. Note that the bytes per value
# for the swapend prog is half the bytes per value used in the corner turn
# if the data is complex. The endian swap is performed in-place using the
# AzcFile.
if G2.cfg['$EnableEndianSwap'] == 'y' and G2.cfg['$OutputEndian'] == 'big':
    if G2.cfg['$DetectMethod'] == 'cmplx':
        ValsPerLine = 2*int(G2.cfg['$AzProcOutPRIs'])
        SwapBytesPerVal = int(int(G2.cfg['$OutBytesPerVal'])/2)
    else:
        ValsPerLine = int(G2.cfg['$AzProcOutPRIs'])
        SwapBytesPerVal = int(G2.cfg['$OutBytesPerVal'])
    cmd = 'swapend '+AzcFile+' '+AzcFile+' '+G2.cfg['$RngBinsToProc']+\
          ' '+`ValsPerLine`+' '+`SwapBytesPerVal`
    if os.system(cmd) != 0:
        ErrorExit('in endian swap of azimuth compression output')

#---------------
# Transpose for rngline output.
if G2.cfg['$EnableOrient'] == 'y' and G2.cfg['$OutputOrient'] == 'rngline':
    RngBinsToProcLessOne = int(G2.cfg['$RngBinsToProc'])-1
    cmd = 'corner '+AzcFile+' '+CorFile+' '+G2.cfg['$AzProcOutPRIs']+' 0 '+\
          `RngBinsToProcLessOne`+' '+G2.cfg['$OutBytesPerVal']+' mem=50'
    if os.system(cmd) != 0:
        ErrorExit('in corner turn after azimuth compression')
    cmd = 'mv '+CorFile+' '+AzcFile
    if os.system(cmd) != 0:
        ErrorExit('in copying over .azc file after orient corner turn')

#---------------
if G2.cfg['$EnableImageLog'] == 'y':
    print '\n------------------'
    print 'Image Log Creation'
    if G2.writeImageLog(ImLog,MocLog,RngLog,AzLog) != 0:
        ErrorExit('in writing the image log file')
    
#---------------
# Clean Up Stage
if G2.cfg['$EnableCleanUp'] == 'y':
    print '\n----------------------'
    print 'Cleaning of Temp Files'
    cmd = 'rm -f '+UnpkIMU+' '+UnpkDGPS+' '+MocMerge+' '+MocLog+' '+MocFile+\
          ' '+MocTxtFile+' '+DCFile+' '+RngCmd+' '+AzCmd+' '+RngLog+' '+\
          AzLog+' '+RncFile+' '+CorFile+' '+ByteFile+' '+MagFile 
    if os.system(cmd) != 0:
        ErrorExit('in cleaning up temporary files')

#---------------
# Done!!
print '\nG2 Processor done!!'

    