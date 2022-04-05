#!/usr/bin/env python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group (1994-1999)
FILE NAME: g2main.py (formerly g2.py)
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

VERSION/AUTHOR/DATE : 0.8 / Jasper Horrell / 1999-10-01
COMMENTS:
Change the swap end command so that it is not performed in place.
This due to a subtle bug discovered in that the file pointer does
not seem to point to the correct place during an in-place operation
with large'ish files.
Ver. 0.81 (2000-02-10) - make ps file of motion error.
Ver. 0.82 (2000-02-16) - work with updated Flt2Byte prog (back compat).
                       Allow EnableFloat2Tiff => only.
2000-03-07 - change xtitle of mocomp error plot.

VERSION/AUTHOR/DATE : 0.9 / Jasper Horrell / 2000-03-24
COMMENTS:
Changes to allow for incorporation of Richard Lord's stepped
frequency processing routines. This also implies changes to
both the radar config and proc config files (the prog should
be backward compatible with the old versions). New strategy
to include a "wrk" dictionary with the G2Config class which
holds params which change with various modules (e.g. delay0,
num of range samples, etc). Calculate or read from log file
the new params after relevant modules have been run.

Certain modules (not all - see template command file) may be
switched "off". This is different to the case of a module being
set to 'n' where the processor assumes that the module's
calculations have been performed and changes certain parameters
accordingly (input file name to next module, PRF, etc.).

Added -? option for general help.

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 2000-08-04
COMMENTS:
Documented processor. New help system.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 2001-02-13
COMMENTS:
Usage changes to allow this code to be run as Python bytecode.
Renamed this file to g2main.py to assist with this.

TO DO LIST (for future versions):
---------------------------------    
If MocMerge file exists, skip processing up to this step. Decide on
suitable name for all polarizations. Also allow forced regeneration.

Add in slant to ground range projection. Allow output image to be
larger in range pixels so as not to change geocoding info.
                       
=========================================================
"""

import string, sys, os, g2tmpl, g2tools, time

print (
'=========================\n',
'G2 Airborne SAR Processor (Ver. 1.0)\n',
'Code: J.M. Horrell (C) University of Cape Town (1994-2001)\n'
)

if (len(sys.argv) != 2) or (sys.argv[1] == '--help'):
    g2tmpl.printHelp()
    sys.exit()
elif len(sys.argv) == 2 and sys.argv[1] == '--templates':
    g2tmpl.createRadarCfgTmpl('Radar.CFG')
    print ('Radar template config file "Radar.CFG" written to current dir')
    g2tmpl.createProcCfgTmpl('Proc.CFG')
    print ('Processor template config file "Proc.CFG" written to current dir')
    sys.exit()
elif len(sys.argv) == 2 and sys.argv[1] == '--user-manual':
    g2tmpl.createHTMLUserManual('G2UserManual.html')
    print ('G2 Processor User Manual "G2UserManual.html" written to current dir')
    sys.exit()

def ErrorExit(cause):
    print ('\nERROR - '+cause+'!!')
    print ('Check relevant log file for more information on error.')
    print ('\nAbnormal exit from G2 Processor!!\n')
    sys.exit(1)

# path to G2 Python bytecode files
tmp = string.split(sys.argv[0],'g2main.p')
G2Path = tmp[0]
if tmp[1] == 'yc': PyExt = 'pyc'
elif tmp[1] == 'y': PyExt = 'py'
else: ErrorExit('unknown extension for python files: '+tmp[1])

# Set up the processor configuration dictionary by reading in parameters
# from the radar- and proc config files. Also, ensure backward
# compatibility with older command file versions. Also, act if an 'only' enable
# has been set and infer extra parameters.
G2 = g2tools.G2Config()    # create a new proc config object
if G2.addCfgFile(sys.argv[1]) != 0: # add contents of proc config file to object 
    ErrorExit('in adding contents of '+sys.argv[1])
if G2.addCfgFile(G2.cfg['$RadarCfgFile']) != 0: # add contents of radar config file
    ErrorExit('in adding contents of '+G2.cfg['$RadarCfgFile'])
if G2.versionControl() !=0:       # compatibilty with older cmd file versions.
    ErrorExit('in executing versionControl function')
if G2.checkForOnly() !=0:         # reset the enables if an only
    ErrorExit('in executing checkForOnly function')
if G2.inferParams() !=0:          # infer extra parameters
    ErrorExit('in inferring processing parameters')
if G2.initChecks() !=0:           # perform some checks
    ErrorExit('in initial parameter checks')   
# Set up some file names
IDRoot =  G2.cfg['$RunID']
UnpkIMU = IDRoot+'_unpkIMU'
UnpkDGPS = IDRoot+'_unpkDGPS'
MocMerge = IDRoot+'_IMU+DGPS'
LBRFile = G2.wrk['$LBRFile']
DCFile = IDRoot+'_sniffDC'
MocLog = IDRoot+'_mlog'
MocFile = IDRoot+'_moc'
G2.cfg['$MocFile'] = MocFile
MocTxtFile = MocFile+'txt'
MocPsFile = MocFile+'.ps'
RngCmd = IDRoot+'_rngcmd'
RngLog = IDRoot+'_rnglog'
G2.cfg['$RngLog'] = RngLog
RncFile = IDRoot+'.rnc'
G2.cfg['$RncFile'] = RncFile
StepFreqCmd = IDRoot+'_stepfcmd'
StepFreqUser = IDRoot+'_stepfuser'
StepFreqLog = IDRoot+'_stepflog'
G2.cfg['$StepFreqLog'] = StepFreqLog
StepFreqFile = IDRoot+'.stepf'
G2.cfg['$StepFreqFile'] = StepFreqFile
CorFile = IDRoot+'.cor'
G2.cfg['$CorFile'] = CorFile
AzCmd = IDRoot+'_azcmd'
AzLog = IDRoot+'_azlog'
G2.cfg['$AzLog'] = AzLog
if G2.cfg['$DetectMethod'] == 'cmplx':
    AzcFile = IDRoot+'.slc'
elif G2.cfg['$DetectMethod'] == 'mag':    
    AzcFile = IDRoot+'.mag'
elif G2.cfg['$DetectMethod'] == 'pow':
    AzcFile = IDRoot+'.pow'
elif G2.cfg['$DetectMethod'] == 'powdB':
    AzcFile = IDRoot+'.powdB'
else:
    print ('ERROR - unknown detect method '+G2.cfg['$DetectMethod'])
    ErrorExit('in setting up file names')   
G2.cfg['$AzcFile'] = AzcFile
ByteFile = IDRoot+'.byte'   # used in float to byte conversion
MagFile = IDRoot+'.iqmag'
ImLog = IDRoot+'.imlog'
TiffFile = IDRoot+'.tif'
SwpEndFile = IDRoot+'_swpend'


#---------------------------------
# Unpack all IMU records to ASCII
if G2.cfg['$EnableUnpackIMU'] == 'y':
    cmd = 'imu_unpack '+LBRFile+' OutF='+UnpkIMU
    if os.system(cmd) != 0:  # run the imu_unpack program
        ErrorExit('in IMU unpack program') 

#----------------------------
# Unpack DGPS Motion Records corresponding to the LBR PRI ID range
if G2.cfg['$EnableUnpackDGPS'] == 'y':
    cmd = "python "+G2Path+'g2unpk_dgps.'+PyExt+' '+G2.wrk['$DGPSFile']+' '+UnpkIMU+' '+\
          UnpkDGPS+' '+G2.cfg['$PRIIncrLBR']+' '+G2.cfg['$RadarPRF']
    if os.system(cmd) != 0:  # run the g2unpk_dgps.py script
        ErrorExit('in DGPS unpack program')

    
#----------------------------
# Merge IMU and DGPS data
if G2.cfg['$EnableMergeMocData'] == 'y':
    cmd = "python "+G2Path+'g2mocfilt.'+PyExt+' '+UnpkDGPS+' '+UnpkIMU+' '+MocMerge+' '+\
          G2.cfg['$MocompLBRPRIOffset']
    if os.system(cmd) != 0:  # run the g2mocfilt.py script
        ErrorExit('in Mocomp DGPS and IMU merge program')

#-----------------------------------------------------------
# Run mocomp range shift calc prog. This also writes ave ground
# speed to motion log file, needed for az compression.

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
    if os.system(cmd) != 0:  # run the mocomp program
        ErrorExit('in motion compensation calculation program')

if G2.cfg['$EnableMocompCalc'] != 'off':  # parse the average ground speed
    G2.cfg['$ProcAveGrndSpeed'] = ''
    if G2.parseParam(MocLog,'Proc average ground speed','$ProcAveGrndSpeed','cfg') != 0:
        ErrorExit('in parsing average ground speed from log')
else:
    if G2.cfg['$AveGroundSpeed'] != 'null':  # use radar config file value, if exists
        G2.cfg['$ProcAveGrndSpeed'] = G2.cfg['$AveGroundSpeed']
    else:
        ErrorExit('AveGroundSpeed must be specified if MocompCalc is off')
        
'''if G2.cfg['$EnablePlotMotionError'] == 'y': # make a pretty plot of the motion error
    g = Gnuplot.Gnuplot()
    g.xlabel('Track Time (secs)')
    g.ylabel('Motion Error (m)')
    data = Gnuplot.File(MocTxtFile,title='Motion error')#,with='lines')
    g.plot(data)
    g.hardcopy(filename=MocPsFile,color=1,fontsize=20)
    time.sleep(5.0)
    del g,data
    print ('Motion error plot written to postscipt file')'''
    
#-----------------------------
# Find DC offsets and I/Q ratio
if G2.cfg['$EnableSniffDC'] == 'y':
    SDCRowsToProc = int(0.2*int(G2.cfg['$InputPRIsToProc']))
    SDCStartCol = int(0.25*int(G2.cfg['$RngBinsToProc']))
    SDCColsToProc = int(0.5*int(G2.cfg['$RngBinsToProc']))
    cmd = 'sniffdc '+G2.wrk['$RawDataFile']+\
          ' '+G2.cfg['$RadarAzSamples']+' '+G2.cfg['$RadarRngBins']+\
          ' StartRow='+G2.cfg['$StartG2PRI']+' RowsToProc='+'SDCRowsToProc'+\
          ' StartCol='+'SDCStartCol'+' ColsToProc='+'SDCColsToProc'+\
          ' LogF='+DCFile
    if os.system(cmd) != 0:  # run the sniffdc program
        ErrorExit('in sniff DC offset program')

if G2.cfg['$EnableSniffDC'] != 'off':  # parse the sniffdc ouput
    if G2.parseParam(DCFile,'Average I value','$DCOffsetI','cfg') != 0:
        ErrorExit('in parsing DC offset I from log')
    if G2.parseParam(DCFile,'Average Q value','$DCOffsetQ','cfg') != 0:
        ErrorExit('in parsing DC offset Q from log')
    if G2.parseParam(DCFile,'Average I/Q ratio','$IQRatio','cfg') != 0:
        ErrorExit('in parsing IQ ratio from log')
else:                             # use the radar config file values, if exist
    if G2.cfg['$DCOffsetI'] == 'null' or G2.cfg['$DCOffsetQ'] == 'null' or \
       G2.cfg['$IQRatio'] == 'null':
        ErrorExit('DC offsets and IQ ratio must be specified if SniffDC off')
        
#-----------------------
# Range Processing Stage

G2.wrk['$Delay0'] = G2.cfg['$RadarDelayToStartSample']
G2.wrk['$PRF'] = G2.cfg['$RadarPRF']
G2.wrk['$AzSamples'] = int(int(G2.cfg['$InputPRIsToProc'])/\
                                  int(G2.cfg['$ProcPresumRatio']))
G2.wrk['$RngSamples'] = G2.cfg['$RngBinsToProc']
G2.wrk['$A2DFreq'] = G2.cfg['$RadarA2DFreq']

if G2.cfg['$RawDataType'] == 'byte':
    G2.wrk['$InputDataType'] = '0'
elif G2.cfg['$RawDataType'] == 'float':
    G2.wrk['$InputDataType'] = '3'

if G2.cfg['$EnableStepFreqProc'] == 'off':
    G2.wrk['$StepFreqUserFile'] = 'null'
elif G2.cfg['$StepFreqMode'] == 'normal':
    G2.wrk['$StepFreqUserFile'] = StepFreqUser
    G2.createStepFreqUserFile()

if G2.cfg['$EnableRngProc'] == 'y':
    G2.createRngProcCmdFile(RngCmd)      # create rngcom config file
    if os.system('rngcom '+RngCmd) != 0: # run rngcom program
        ErrorExit('in range compression program')

# reset output params for later modules
G2.wrk['$Delay0'] = float(G2.cfg['$RadarDelayToStartSample'])+\
                     float(G2.cfg['$StartRngBinToProc'])/\
                     float(G2.cfg['$RadarA2DFreq'])
G2.wrk['$PRF'] = float(G2.cfg['$RadarPRF'])/\
                  float(G2.cfg['$ProcPresumRatio'])

#------------------------------
# Stepped Freq Processing Stage

if G2.cfg['$EnableStepFreqProc'] == 'y':
    if G2.cfg['$StepFreqMode'] == 'normal':
        G2.wrk['$StepFreqUserFile'] = 'null'  # reset for stepf cmd file
    G2.createStepFreqProcCmdFile(StepFreqCmd) # create stepf config file
    if os.system('stepf '+StepFreqCmd) != 0:  # run stepf program
        ErrorExit('in step frequency processing program')

if G2.cfg['$EnableStepFreqProc'] != 'off':
    # reset output params for ater modules
    G2.wrk['$PRF'] = float(float(G2.wrk['$PRF'])/\
                          int(G2.cfg['$NumberOfFreqSteps']))
    G2.parseParam(StepFreqLog,'Log File Version','$StepFreqLogVer','cfg')
    if G2.cfg['$StepFreqLogVer'] != '0.1':
        print ('WARNING - incorrect step freq log version - expecting 0.1 !!!')
    if G2.parseParam(StepFreqLog,'Number of range bins','$RngSamples','wrk') != 0:
        ErrorExit('in parsing range samples from step freq log')
    if G2.parseParam(StepFreqLog,'Output A2D frequency','$A2DFreq','wrk') != 0:
        ErrorExit('in parsing output A2D freq from step freq log')
    if G2.parseParam(StepFreqLog,'Number of range lines','$AzSamples','wrk') != 0:
        ErrorExit('in parsing az samples from step freq log')
    cornerInFile = StepFreqFile
else: # if StepFreqProc set to 'off'
    cornerInFile = RncFile
       
#------------
# Corner Turn
if G2.cfg['$EnableCornerTurn'] == 'y':
    AzSamplesLessOne = int(G2.wrk['$AzSamples'])-1
    cmd = 'corner '+cornerInFile+' '+CorFile+' '+G2.wrk['$RngSamples']+' 0 '+\
          AzSamplesLessOne+' 8 mem=50'
    if os.system(cmd) != 0:  # corner turn the data on disk
        ErrorExit('in corner turn program')     

#-------------------------
# Azimuth Processing Stage

if G2.calcAzProcParams() != 0:
    ErrorExit('in calculating az proc params')

if G2.cfg['$EnableAzProc'] == 'y':
    G2.createAzProcCmdFile(AzCmd)      # create azcom config file
    if os.system('azcom '+AzCmd) != 0: # run azcom program
        ErrorExit('in azimuth compression program')

if G2.cfg['$EnableAzProc']!= 'off':
    if G2.parseParam(AzLog,'Output azimuth samples','$AzSamples','wrk') != 0:
        ErrorExit('in parsing output azimuth samples from azcom log')
    
#------------------
# TIFF Output Stage
if G2.cfg['$EnableFloat2Tiff'] == 'y':
    if G2.cfg['$DetectMethod'] == 'cmplx':  # first run iq2mag
        cmd = 'iq2mag '+AzcFile+' '+MagFile+' '+G2.wrk['$RngSamples']+' '+\
              G2.wrk['$AzSamples']+' pow float' 
        if os.system(cmd) != 0: # convert IQ values to power
            ErrorExit('in IQ to mag conversion program')
        flt2byteInFile = MagFile
    else:                            # data not complex, so skip iq2mag       
        flt2byteInFile = AzcFile

    cmd = 'flt2byte '+flt2byteInFile+' '+ByteFile+' '+G2.wrk['$RngSamples']+\
          ' '+G2.wrk['$AzSamples']+' math='+G2.cfg['$Float2ByteMath']+\
          ' scale='+G2.cfg['$Float2ByteScale']
    if os.system(cmd) != 0:     # convert floats to unsigned chars (bytes) 
        ErrorExit('in float to byte conversion program')

    cmd = 'b2tif '+ByteFile+' '+TiffFile+' '+G2.wrk['$RngSamples']+' '+\
           G2.wrk['$AzSamples']+' orient=3'
    if os.system(cmd) != 0:     # make a TIFF format file from byte file
        ErrorExit('in byte to tiff conversion program')

#---------------
# Swap endian order for big endian output. Note that the bytes per value
# for the swapend prog is half the bytes per value used in the corner turn
# if the data is complex. The endian swap is performed in-place using the
# AzcFile.
if G2.cfg['$EnableEndianSwap'] == 'y' and G2.cfg['$OutputEndian'] == 'big':
    if G2.cfg['$DetectMethod'] == 'cmplx':
        ValsPerLine = 2*int(G2.wrk['$AzSamples'])
        SwapBytesPerVal = int(int(G2.cfg['$OutBytesPerVal'])/2)
    else:
        ValsPerLine = int(G2.wrk['$AzSamples'])
        SwapBytesPerVal = int(G2.cfg['$OutBytesPerVal'])
    cmd = 'swapend '+AzcFile+' '+AzcFile+' '+G2.wrk['$RngSamples']+\
          ' '+ValsPerLine+' '+SwapBytesPerVal
    if os.system(cmd) != 0:  # run swapend program
        ErrorExit('in endian swap of azimuth compression output')

#---------------
# Transpose for rngline output. After the corner turn, the AzcFile is 
# overwritten.
if G2.cfg['$EnableOrient'] == 'y' and G2.cfg['$OutputOrient'] == 'rngline':
    RngBinsToProcLessOne = int(G2.wrk['$RngSamples'])-1
    cmd = 'corner '+AzcFile+' '+CorFile+' '+G2.wrk['$AzSamples']+' 0 '+\
          RngBinsToProcLessOne+' '+G2.cfg['$OutBytesPerVal']+' mem=50'
    if os.system(cmd) != 0:  # run corner turn program
        ErrorExit('in corner turn after azimuth compression')
    print ('\nMoving '+CorFile+' to '+AzcFile)
    cmd = 'mv '+CorFile+' '+AzcFile
    if os.system(cmd) != 0:  # copy corner turned file over .azc file 
        ErrorExit('in copying over .azc file after orient corner turn')

#---------------
if G2.cfg['$EnableImageLog'] == 'y':
    print ('\n------------------')
    print ('Image Log Creation')
    if G2.writeImageLog(ImLog,MocLog,RngLog,AzLog) != 0: # write image log
        ErrorExit('in writing the image log file')
    
#---------------
# Clean Up Stage. Leaves the output image binary and TIFF files,
# the output image log file and the input config files.
if G2.cfg['$EnableCleanUp'] == 'y':
    print ('\n----------------------')
    print ('Cleaning of Temp Files')
    cmd = 'rm -f '+UnpkIMU+' '+UnpkDGPS+' '+MocMerge+' '+MocLog+' '+MocFile+\
          ' '+MocTxtFile+' '+MocPsFile+' '+DCFile+' '+RngCmd+' '+RngLog+\
          ' '+StepFreqCmd+' '+StepFreqLog+' '+StepFreqFile+\
          ' '+AzCmd+' '+AzLog+' '+RncFile+' '+CorFile+' '+ByteFile+' '+MagFile 
    if os.system(cmd) != 0: # remove temporary files
        ErrorExit('in cleaning up temporary files')

#---------------
# Done!!
print ('\nG2 Processor done!!')

    
