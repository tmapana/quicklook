"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2tools.py
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
The G2 SAR processor python tools. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-07-31
COMMENTS:
Initial python version.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-04
COMMENTS:
Calc FFT sizes automagically in C progs.

VERSION/AUTHOR/DATE : 0.3 / Jasper Horrell / 1999-09-14
COMMENTS:
Allow specification of interference suppression parameters.
RadarA2dFreq renamed RadarA2DFreq. Add UTC Hrs,Mins,Secs. 
1999-10-01 - allow notch filter in interference suppress.

VERSION/AUTHOR/DATE : 0.4 / Jasper Horrell / 1999-10-01
COMMENTS:
Allow notch interference suppression. Allow output endian and
orientation specification. Replace space with underscores in
image log file parameters.

VERSION/AUTHOR/DATE : 0.5 / Jasper Horrell / 2000-02-16
COMMENTS:
Add version control function to allow use of older proc config
files with updated processor version. Add checkForOnly function.

VERSION/AUTHOR/DATE : 0.6 / Jasper Horrell / 2000-03-24
COMMENTS:
Changes to allow for stepped freq processing. To simplify this
and future development, add a new dictionary to the G2Config
class called wrk (work). The strategy is to allow the "wrk"
parameters to change with time, but the "cfg" parameters remain
fixed (e.g. the step freq code changes the number of range bins
to be used in subsequent steps). Using this approach and updating
the "wrk" dictionary after any processing steps which change the
parameters, and setting up params for the next step as late as
possible avoids having to calculate every possible scenario
up front. The orginal params (such as radar A2D freq) are also
still available in the "cfg" dictionary (for image log file, etc).

VERSION/AUTHOR/DATE : 0.7 / Jasper Horrell / 2000-03-11
COMMENTS:
Add rng compress phase sign to step freq config file.
In rng compress config file, remove 'Name' from some fields
e.g. LogFileName -> LogFile.


=========================================================
"""

import string, sys
              

class G2Config:

    def __init__(self, ConfigDic={}):
        self.cfg = ConfigDic  # mostly for parameters static with time
        self.wrk = {}         # for time-varying parameters (e.g. RngSamples) 

    # Parse a config file and add keys (params) and
    # values to the config object's cfg dictionary
    def addCfgFile(self,FileName):
        f = open(FileName,'r')
        while 1:
            str1 = f.readline()
            if not str1:
                break
            if str1.find('=>') != -1:   # look for lines with '=>' 
                field = str1.split()
                self.cfg[field[0]] = field[len(field)-1]
        f.close()
        return 0

    def parseParam(self,FileName,SearchStr,KeyStr,dict):
        f = open(FileName,'r')
        SearchStrFound = 0
        while 1:    
            str1 = f.readline()
            if not str1:
                break
            if str1.find(SearchStr) != -1: # find line with SearchStr  
                field = str1.split()
                for i in range(len(field)):
                    if field[i] == ':':  # find the colon
                        if dict == 'cfg':
                            self.cfg[KeyStr] = field[i+1]
                        elif dict == 'wrk':
                            self.wrk[KeyStr] = field[i+1]
                        else:
                            return -1
                        SearchStrFound = 1
                        break
        f.close()
        if SearchStrFound != 1:
            return -1
        else:
            return 0

    def convSecsToTimeStr(self,inSecs,KeyStr):
        NumHrs = int(inSecs/3600.0)
        NumMins = int((inSecs-NumHrs*3600)/60.0)
        NumSecs = float(inSecs-NumHrs*3600-NumMins*60.0)
        if NumHrs < 10:
            HrsStr = '0'+NumHrs
        else:
            HrsStr = NumHrs
        if NumMins < 10:
            MinsStr = '0'+NumMins
        else:
            MinsStr = NumMins    
        if NumSecs < 10:
            SecsStr = '0'+NumSecs
        else:
            SecsStr = NumSecs       
        self.cfg[KeyStr] = HrsStr+'h'+MinsStr+':'+SecsStr
        return 0

    def convDegToDegMinSec(self,inDeg,KeyStr):
        inSecs = 3600.0*inDeg
        sign = 1
        if inSecs < 0.0:
            inSecs = -1.0*inSecs
            sign = -1
        NumDegs = int(inSecs/3600.0)
        NumMins = int((inSecs-NumDegs*3600)/60.0)
        NumSecs = float(inSecs-NumDegs*3600-NumMins*60.0)            
        self.cfg[KeyStr] = NumDegs*sign+':'+NumMins+':'+NumSecs+' deg:min:sec'
        return 0


    def displayCfgContents(self):
        for i in range(len(self.cfg)):
            print(i, self.cfg.keys()[i], self.cfg.values()[i])

    def displayWrkContents(self):
        for i in range(len(self.wrk)):
            print( i, self.wrk.keys()[i], self.wrk.values()[i]   )

    def versionControl(self):
        if float(self.cfg['$RadarConfigVersion']) < 0.4:
            self.cfg['$RawDataType'] = 'byte'
            self.cfg['$DCOffsetI'] = 'null'
            self.cfg['$DCOffsetQ'] = 'null'
            self.cfg['AveGroundSpeed'] = 'null'
            print( 'Version control - setting RawDataType to byte')
            self.cfg['$StepFreqMode'] = 'no'
            print( 'Version control - setting StepFreqMode to no')
            self.cfg['$NumberOfFreqSteps'] = 1
            self.cfg['$FirstStepCentreFreq'] = self.cfg['$RadarCarrierFreq']
            self.cfg['$StepFreqStepSize'] = self.cfg['$RadarChirpBandwidth']
            self.cfg['$StepFreqUserFile'] = 'null'
            
        if float(self.cfg['$ProcConfigVersion']) < 0.6:
            self.cfg['$Float2ByteMath'] = 'none'
            print( 'Version control - setting Float2ByteMath to none')

        if float(self.cfg['$ProcConfigVersion']) < 0.7:
            self.cfg['$EnableStepFreqProc'] = 'n'
            print( 'Version control - setting EnableStepFreqProc to no')
            self.cfg['$StepFreqWinConstTime'] = 0.08           
        return 0

#-----------------
# Cater for possibility that off is allowed for all modules even though
# only makes sense at this stage if core modules not in off mode.
    def disableAll(self):
        if self.cfg['$EnableUnpackIMU'] != 'off':
            self.cfg['$EnableUnpackIMU'] = 'n'
        if self.cfg['$EnableUnpackDGPS'] != 'off':
            self.cfg['$EnableUnpackDGPS'] = 'n'
        if self.cfg['$EnableMergeMocData'] != 'off':
            self.cfg['$EnableMergeMocData'] = 'n'
        if self.cfg['$EnableMocompCalc'] != 'off':
            self.cfg['$EnableMocompCalc'] = 'n'
        if self.cfg['$EnablePlotMotionError'] != 'off':
            self.cfg['$EnablePlotMotionError'] = 'n'
        if self.cfg['$EnableSniffDC'] != 'off':
            self.cfg['$EnableSniffDC'] = 'n'
        if self.cfg['$EnableRngProc'] != 'off':
            self.cfg['$EnableRngProc'] = 'n'
        if self.cfg['$EnableStepFreqProc'] != 'off':
            self.cfg['$EnableStepFreqProc'] = 'n'
        if self.cfg['$EnableCornerTurn'] != 'off':
            self.cfg['$EnableCornerTurn'] = 'n'
        if self.cfg['$EnableAzProc'] != 'off':      
            self.cfg['$EnableAzProc'] = 'n'
        if self.cfg['$EnableFloat2Tiff'] != 'off':
            self.cfg['$EnableFloat2Tiff'] = 'n'
        if self.cfg['$EnableEndianSwap'] != 'off':
            self.cfg['$EnableEndianSwap'] = 'n'
        if self.cfg['$EnableOrient'] != 'off':
            self.cfg['$EnableOrient'] = 'n'
        if self.cfg['$EnableImageLog'] != 'off':
            self.cfg['$EnableImageLog'] = 'n'
        if self.cfg['$EnableCleanUp'] != 'off':
            self.cfg['$EnableCleanUp'] = 'n'        
        return 0

#-------------------
    def checkForOnly(self):
        if  self.cfg['$EnableUnpackIMU'] == 'only':
            self.disableAll()
            self.cfg['$EnableUnpackIMU'] = 'y'
        elif self.cfg['$EnableUnpackDGPS'] == 'only':
            self.disableAll()
            self.cfg['$EnableUnpackDGPS'] = 'y'
        elif self.cfg['$EnableMergeMocData'] == 'only':
            self.disableAll()
            self.cfg['$EnableMergeMocData'] = 'y'
        elif self.cfg['$EnableMocompCalc'] == 'only':
            self.disableAll()
            self.cfg['$EnableMocompCalc'] = 'y'
        elif self.cfg['$EnablePlotMotionError'] == 'only':
            self.disableAll()
            self.cfg['$EnablePlotMotionError'] = 'y'
        elif self.cfg['$EnableSniffDC'] == 'only':
            self.disableAll()
            self.cfg['$EnableSniffDC'] = 'y'
        elif self.cfg['$EnableRngProc'] == 'only':
            self.disableAll()
            self.cfg['$EnableRngProc'] = 'y'
        elif self.cfg['$EnableStepFreqProc'] == 'only':
            self.disableAll()
            self.cfg['$EnableStepFreqProc'] = 'y'
        elif self.cfg['$EnableCornerTurn'] == 'only':
            self.disableAll()
            self.cfg['$EnableCornerTurn'] = 'y'
        elif self.cfg['$EnableAzProc'] == 'only':
            self.disableAll()
            self.cfg['$EnableAzProc'] = 'y'
        elif self.cfg['$EnableFloat2Tiff'] == 'only':
            self.disableAll()
            self.cfg['$EnableFloat2Tiff'] = 'y'
        elif self.cfg['$EnableEndianSwap'] == 'only':
            self.disableAll()
            self.cfg['$EnableEndianSwap'] = 'y'
        elif self.cfg['$EnableOrient'] == 'only':
            self.disableAll()
            self.cfg['$EnableOrient'] = 'y'
        elif self.cfg['$EnableImageLog'] == 'only':
            self.disableAll()
            self.cfg['$EnableImageLog'] = 'y'
        elif self.cfg['$EnableCleanUp'] == 'only':
            self.disableAll()
            self.cfg['$EnableCleanUp'] = 'y'
        return 0

#-----------------
    def initChecks(self): # called after inferParams()
        if ( int(self.cfg['$StartG2PRI'])+ int(self.cfg['$InputPRIsToProc']) >\
             int(self.cfg['$RadarAzSamples'])):
             print( 'ERROR - insufficient azimuth samples in raw file!')
             return 1
             
        # checks for stepped freq processing
        if self.cfg['$StepFreqMode'] == 'no' or \
           self.cfg['$NumberOfFreqSteps'] == '1':
            self.cfg['$EnableStepFreqProc'] = 'off' # don't calc new A2D freq, etc

        if self.cfg['$EnableStepFreqProc'] != 'off':
            if int(self.cfg['$StartG2PRI']) % int(self.cfg['$NumberOfFreqSteps']) != 0:
                print( 'ERROR - StartG2PRI must be multiple of NumberOfFreqSteps!')
                return 1
            if float(self.cfg['$RngComWinConstTime'])!=1.0:
                print( 'WARNING - range compress time domain window const should')
                print( '          be set to 1.0 for step freq processing!')
            
        # checks for module 'off' dependencies
        if self.cfg['$EnableUnpackIMU'] == 'off' and \
           self.cfg['$EnableUnpackDGPS'] != 'off':
            print( 'ERROR - UnpackIMU off while UnpackDGPS not off')
            return 1
        if self.cfg['$EnableMocompCalc'] != 'off' and \
           (self.cfg['$EnableUnpackIMU'] == 'off' or \
            self.cfg['$EnableUnpackDGPS'] == 'off' or \
            self.cfg['$EnableMergeMocData'] == 'off'):
            print( 'ERROR - IMU/DGPS/MergeMoc off while MocompCalc not off')
            return 1
        if self.cfg['$EnableMocompCalc'] == 'off' and \
           self.cfg['$EnablePlotMotionError'] != 'off':
            print( 'ERROR - MocompCalc off while PlotMotionError not off')
            return 1
        if self.cfg['$EnableSniffDC'] != 'off' and \
           self.cfg['$RawDataType'] != 'byte':
            print( 'ERROR - SniffDC must be off if raw data not byte')
            return 1
        return 0
        
       
#-------------------        
    def inferParams(self): # called near start

        # set up some working file names
        self.wrk['$RawDataFile']=self.cfg['$InputPath']+self.cfg['$RawDataFile']
        self.wrk['$LBRFile']=self.cfg['$InputPath']+self.cfg['$LBRFile']
        self.wrk['$DGPSFile']=self.cfg['$InputPath']+self.cfg['$DGPSFile']
        self.wrk['$StepFreqUserFile']=self.cfg['$InputPath']+\
                                      self.cfg['$StepFreqUserFile']        

        # general calcs
        PRIIncrLBR = 2*int(self.cfg['$RadarRngBins'])/512
        self.cfg['$PRIIncrLBR'] = PRIIncrLBR
        MocompLBRPRIOffset = int(float(self.cfg['$RadarMocTimeOffset'])*\
                             float(self.cfg['$RadarPRF'])+0.5)*PRIIncrLBR
        self.cfg['$MocompLBRPRIOffset'] = MocompLBRPRIOffset                      
               
        ProcStartLP = int(self.cfg['$DataStartLP'])+\
                      int(self.cfg['$StartG2PRI'])*PRIIncrLBR
        self.cfg['$ProcStartLP'] = ProcStartLP
        ProcEndLP = ProcStartLP+(int(self.cfg['$InputPRIsToProc'])-1)*PRIIncrLBR
        self.cfg['$ProcEndLP'] = ProcEndLP
        EndRngBinToProc = int(int(self.cfg['$StartRngBinToProc'])+\
                              int(self.cfg['$RngBinsToProc'])-1)
        self.cfg['$EndRngBinToProc'] = EndRngBinToProc                               
        MocompRefRngBin = int(int(self.cfg['$StartRngBinToProc'])+\
                              int(self.cfg['$RngBinsToProc'])/2)
        self.cfg['$MocompRefRngBin'] = MocompRefRngBin                      
 
        if self.cfg['$InterferenceSuppress'] == 'lms':
            self.cfg['$LmsFlg'] = 'Y'
            self.cfg['$NotchFlg'] = 'N'
        elif self.cfg['$InterferenceSuppress'] == 'notch':
            self.cfg['$LmsFlg'] = 'N'
            self.cfg['$NotchFlg'] = 'Y'
        elif self.cfg['$InterferenceSuppress'] == 'none':
            self.cfg['$LmsFlg'] = 'N'
            self.cfg['$NotchFlg'] = 'N'
        else:
            print( 'ERROR - unknown interference suppression method '+\
              self.cfg['$InterferenceSuppress']+'!')
            return 1
        return 0

#-----------------------------
# rngcom module requires a step freq user file for normal step freq mode.
    def createStepFreqUserFile(self):
        FreqSteps = int(self.cfg['$NumberOfFreqSteps'])
        UserFile = open(self.wrk['$StepFreqUserFile'],'w')
        UserFile.write(FreqSteps+'\n')
        CentreFreq = float(self.cfg['$FirstStepCentreFreq'])
        for step in range(FreqSteps):
            UserFile.write(CentreFreq+'\n')
            CentreFreq = CentreFreq + float(self.cfg['$StepFreqStepSize'])
        UserFile.close()
        return 0

#------------------
    def calcAzProcParams(self): # set up params for azimuth processing  
        if self.cfg['$DetectMethod'] != 'cmplx':
            self.cfg['$OutBytesPerVal'] = 4
        else:
            self.cfg['$OutBytesPerVal'] = 8

        if self.cfg['$DetectMethod'] == 'cmplx':
            self.cfg['$AzProcDetectMethod'] = 0
        elif self.cfg['$DetectMethod'] == 'mag':
            self.cfg['$AzProcDetectMethod'] = 1
        elif self.cfg['$DetectMethod'] == 'pow':
            self.cfg['$AzProcDetectMethod'] = 2
        elif self.cfg['$DetectMethod'] == 'powdB':
            self.cfg['$AzProcDetectMethod'] = 3
        else:
            print( 'ERROR - unknown detect method '+self.cfg['$DetectMethod']+'!')
            return 1
        return 0

#--------------------
    def writeImageLog(self,ILogName,MLogName,RLogName,ALogName):       
        
        # Parse relevant params
        self.parseParam(MLogName,'Antenna direction','$AntennaDirn','cfg')
        self.parseParam(MLogName,'Proc start DGPS UTC','$ImMocStartUTC','cfg')
        self.parseParam(MLogName,'Proc end DGPS UTC','$ImMocEndUTC','cfg')
        self.convSecsToTimeStr(float(self.cfg['$ImMocStartUTC']),'$ImMocStartUTCTime')
        self.convSecsToTimeStr(float(self.cfg['$ImMocEndUTC']),'$ImMocEndUTCTime')
        self.parseParam(MLogName,'Proc average altitude','$ImAveAlt','cfg')
        self.parseParam(MLogName,'Image near swath start latitude','$ImNearStartLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImNearStartLat']),'$ImNearStartLatDMS')
        self.parseParam(MLogName,'Image near swath start longitude','$ImNearStartLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImNearStartLong']),'$ImNearStartLongDMS')        
        self.parseParam(MLogName,'Image near swath start altitude','$ImNearStartAlt','cfg')
        self.parseParam(MLogName,'Image far swath start latitude','$ImFarStartLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImFarStartLat']),'$ImFarStartLatDMS')
        self.parseParam(MLogName,'Image far swath start longitude','$ImFarStartLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImFarStartLong']),'$ImFarStartLongDMS')
        self.parseParam(MLogName,'Image far swath start altitude','$ImFarStartAlt','cfg')
        self.parseParam(MLogName,'Image near swath end latitude','$ImNearEndLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImNearEndLat']),'$ImNearEndLatDMS')
        self.parseParam(MLogName,'Image near swath end longitude','$ImNearEndLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImNearEndLong']),'$ImNearEndLongDMS')
        self.parseParam(MLogName,'Image near swath end altitude','$ImNearEndAlt','cfg')
        self.parseParam(MLogName,'Image far swath end latitude','$ImFarEndLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImFarEndLat']),'$ImFarEndLatDMS')
        self.parseParam(MLogName,'Image far swath end longitude','$ImFarEndLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImFarEndLong']),'$ImFarEndLongDMS')
        self.parseParam(MLogName,'Image far swath end altitude','$ImFarEndAlt','cfg')
        self.parseParam(MLogName,'Image platform ref track start lat','$ImPlatRefStartLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImPlatRefStartLat']),'$ImPlatRefStartLatDMS')        
        self.parseParam(MLogName,'Image platform ref track start long','$ImPlatRefStartLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImPlatRefStartLong']),'$ImPlatRefStartLongDMS')
        self.parseParam(MLogName,'Image platform ref track start alt','$ImPlatRefStartAlt','cfg')
        self.parseParam(MLogName,'Image platform ref track end lat','$ImPlatRefEndLat','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImPlatRefEndLat']),'$ImPlatRefEndLatDMS')
        self.parseParam(MLogName,'Image platform ref track end long','$ImPlatRefEndLong','cfg')
        self.convDegToDegMinSec(float(self.cfg['$ImPlatRefEndLong']),'$ImPlatRefEndLongDMS')
        self.parseParam(MLogName,'Image platform ref track end alt','$ImPlatRefEndAlt','cfg')

        self.parseParam(RLogName,'Range sample spacing','$ImRngSampleSpacing','cfg')
        self.parseParam(RLogName,'Nominal range resolution','$ImNomRngRes','cfg')

        self.parseParam(ALogName,'Output azimuth sample spacing','$ImAzSampleSpacing','cfg')
        self.parseParam(ALogName,'Nominal azimuth resolution','$ImNomAzRes','cfg')
        self.parseParam(ALogName,'Output azimuth samples','$ImAzSamples','cfg')
        self.parseParam(ALogName,'Range samples to process','$ImRngSamples','cfg')
        self.parseParam(ALogName,'Output data type','$ImDataType','cfg')
        self.parseParam(ALogName,'Doppler bandwidth processed','$ImDopBW','cfg')
        self.parseParam(ALogName,'Processed near slant range','$ImNearRng','cfg')
        self.parseParam(ALogName,'Processed far slant range','$ImFarRng','cfg')
        
        # write image log
        ILog = open(ILogName,'w')
        ILog.write(\
        'G2 Processor (Ver 0.9) Image Log File\n'+\
        'Code: J.M. Horrell - UCT Radar Remote Sensing Group 1999-2000\n'\
        'Image_log_file_version              => 0.2\n\n'+\

        'Image_azimuth_samples               => '+self.cfg['$ImAzSamples']+'\n'+\
        'Image_range_samples                 => '+self.cfg['$ImRngSamples']+'\n'+\
        'Image_azimuth_sample_spacing        => '+self.cfg['$ImAzSampleSpacing']+' m\n'+\
        'Image_range_sample_spacing          => '+self.cfg['$ImRngSampleSpacing']+' m\n'+\
        'Image_Doppler_bandwidth_processed   => '+self.cfg['$ImDopBW']+' Hz\n'+\
        'Image_nominal_azimuth_resolution    => '+self.cfg['$AzComNomAzRes']+' m\n'+\
        'Image_nominal_range_resolution      => '+self.cfg['$ImNomRngRes']+' m\n'+\
        'Image_range_window_constant         => '+self.cfg['$RngComWinConstTime']+'\n'+\
        'Image_azimuth_window_constant_time  => '+self.cfg['$AzComWinConstTime']+'\n'+\
        'Image_azimuth_window_constant_freq  => '+self.cfg['$AzComWinConstFreq']+'\n'+\
        'Image_azimuth_looks                 => '+self.cfg['$NumAzLooks']+'\n'+\
        'Image_azimuth_look_overlap          => '+self.cfg['$AzLookOverlapFrac']+'\n'+\
        'Image_detection_method              => '+self.cfg['$DetectMethod']+'\n'+\
        'Image_data_type                     => '+self.cfg['$ImDataType']+'\n'+\
        'Image_data_endian_type              => '+self.cfg['$OutputEndian']+'\n'+\
        'Image_projection                    => slant_range\n'+\
        'Image_orientation                   => '+self.cfg['$OutputOrient']+'\n'+\
        'Image_near_slant_range              => '+self.cfg['$ImNearRng']+' m\n'+\
        'Image_far_slant_range               => '+self.cfg['$ImFarRng']+' m\n'+\
        'Image_motion_compensation           => '+self.cfg['$MoComp']+'\n'+\
        'Image_interference_suppression      => '+self.cfg['$InterferenceSuppress']+'\n'+\
        'Image_terrain_altitude              => '+self.cfg['$TerrainAlt']+' m (AMSL)\n'+\
        'Image_platform_average_altitude     => '+self.cfg['$ImAveAlt']+' m (AMSL)\n'+\
        'Image_near_swath_start_latitude     => '+self.cfg['$ImNearStartLat']+' deg ('\
                                                 +self.cfg['$ImNearStartLatDMS']+')\n')
        ILog.write(\
        'Image_near_swath_start_longitude    => '+self.cfg['$ImNearStartLong']+' deg ('\
                                                 +self.cfg['$ImNearStartLongDMS']+')\n'+\
        'Image_near_swath_start_altitude     => '+self.cfg['$ImNearStartAlt']+' m (AMSL)\n'+\
        'Image_far_swath_start_latitude      => '+self.cfg['$ImFarStartLat']+' deg ('\
                                                 +self.cfg['$ImFarStartLatDMS']+')\n'+\
        'Image_far_swath_start_longitude     => '+self.cfg['$ImFarStartLong']+' deg ('\
                                                 +self.cfg['$ImFarStartLongDMS']+')\n'+\
        'Image_far_swath_start_altitude      => '+self.cfg['$ImFarStartAlt']+' m (AMSL)\n'+\
        'Image_near_swath_end_latitude       => '+self.cfg['$ImNearEndLat']+' deg ('\
                                                 +self.cfg['$ImNearEndLatDMS']+')\n'+\
        'Image_near_swath_end_longitude      => '+self.cfg['$ImNearEndLong']+' deg ('\
                                                 +self.cfg['$ImNearEndLongDMS']+')\n'+\
        'Image_near_swath_end_altitude       => '+self.cfg['$ImNearEndAlt']+' m (AMSL)\n'+\
        'Image_far_swath_end_latitude        => '+self.cfg['$ImFarEndLat']+' deg ('\
                                                 +self.cfg['$ImFarEndLatDMS']+')\n'+\
        'Image_far_swath_end_longitude       => '+self.cfg['$ImFarEndLong']+' deg ('\
                                                 +self.cfg['$ImFarEndLongDMS']+')\n'+\
        'Image_far_swath_end_altitude        => '+self.cfg['$ImFarEndAlt']+' m (AMSL)\n'+\
        'Image_platform_ref_track_start_lat  => '+self.cfg['$ImPlatRefStartLat']+' deg ('\
                                                 +self.cfg['$ImPlatRefStartLatDMS']+')\n'+\
        'Image_platform_ref_track_start_long => '+self.cfg['$ImPlatRefStartLong']+' deg ('\
                                                 +self.cfg['$ImPlatRefStartLongDMS']+')\n'+\
        'Image_platform_ref_track_start_alt  => '+self.cfg['$ImPlatRefStartAlt']+' m (AMSL)\n'+\
        'Image_platform_ref_track_end_lat    => '+self.cfg['$ImPlatRefEndLat']+' deg ('\
                                                 +self.cfg['$ImPlatRefEndLatDMS']+')\n'+\
        'Image_platform_ref_track_end_long   => '+self.cfg['$ImPlatRefEndLong']+' deg ('\
                                                 +self.cfg['$ImPlatRefEndLongDMS']+')\n'+\
        'Image_platform_ref_track_end_alt    => '+self.cfg['$ImPlatRefEndAlt']+' m (AMSL)\n\n'+\

        'Raw_data_ID                         => '+self.cfg['$DataID']+'\n'+\
        'Raw_data_file                       => '+self.cfg['$RawDataFile']+'\n'+\
        'LBR_file                            => '+self.cfg['$LBRFile']+'\n'+\
        'DGPS_file                           => '+self.cfg['$DGPSFile']+'\n'+\
        'Radar_config_file                   => '+self.cfg['$RadarCfgFile']+'\n'+\
        'Step_freq_user_file                 => '+self.cfg['$StepFreqUserFile']+'\n'+\
        'Stepped_freq_mode                   => '+self.cfg['$StepFreqMode']+'\n'+\
        'Antenna_direction                   => '+self.cfg['$AntennaDirn']+'\n'+\
        'Radar_carrier_frequency             => '+self.cfg['$RadarCarrierFreq']+' Hz\n'+\
        'Radar_pulse_length                  => '+self.cfg['$RadarPulseLength']+' sec\n'+\
        'Radar_chirp_bandwidth               => '+self.cfg['$RadarChirpBandwidth']+' Hz\n'+\
        'Delay_to_raw_data_start_sample      => '+self.cfg['$RadarDelayToStartSample']+'\n'+\
        'Mocomp_start_DGPS_UTC               => '+self.cfg['$ImMocStartUTC']+' sec ('\
                                                 +self.cfg['$ImMocStartUTCTime']+')\n'+\
        'Mocomp_end_DGPS_UTC                 => '+self.cfg['$ImMocEndUTC']+' sec ('\
                                                 +self.cfg['$ImMocEndUTCTime']+')\n'+\
        'Raw_data_start_process_range_sample => '+self.cfg['$StartRngBinToProc']+'\n'+\
        'Raw_data_range_samples_to_process   => '+self.cfg['$RngBinsToProc']+'\n'+\
        'Raw_data_start_azimuth_sample       => '+self.cfg['$StartG2PRI']+'\n'+\
        'Raw_data_azimuth_samples_to_process => '+self.cfg['$InputPRIsToProc']+'\n'+\
        'Raw_data_start_LBR_PRI              => '+self.cfg['$DataStartLP']+'\n'+\
        'Raw_data_PRF                        => '+self.cfg['$RadarPRF']+' Hz\n'+\
        'Proc_presum_ratio                   => '+self.cfg['$ProcPresumRatio']+'\n'+\
        'Proc_average_ground_speed           => '+self.cfg['$ProcAveGrndSpeed']+' m/s\n'+\
        'Rng_compress_ref_phase_sign         => '+self.cfg['$RngComRefPhase']+'\n'+\
        'Az_compress_ref_phase_sign          => '+self.cfg['$AzComRefPhase']+'\n'+\
        'Mocomp_phase_sign                   => '+self.cfg['$MoCompPhaseSign']+'\n'+\
        'Mocomp_range_shift_sign             => '+self.cfg['$MoCompRngShiftSign']+'\n'\
        )
        ILog.close()
        return 0
     
#---------------------------        
    def createRngProcCmdFile(self,FileName):
        f = open(FileName,'w')
        f.write(\
        'G2-created rngcom command file\n'+\
        '$ProgramVersion (jmh)    => 1.3\n'+\
        '--------------------------------------\n\n'+\
        '--General (required)--\n'+\
        '$ScreenUpdateRate        => 20\n'+\
        '$LogFile                 => '+self.cfg['$RngLog']+'\n'+\
        '$InputFile               => '+self.wrk['$RawDataFile']+'\n'+\
        '$OutputFile              => '+self.cfg['$RncFile']+'\n'+\
        '$InputDataType           => '+self.wrk['$InputDataType']+'\n'+\
        '$OutputDataType          => 3\n'+\
        '$StartProcessPRI         => '+self.cfg['$StartG2PRI']+'\n'+\
        '$PreSumRatio             => '+self.cfg['$ProcPresumRatio']+'\n'+\
        '$PreSummedPulsesToUse    => '+self.wrk['$AzSamples']+'\n'+\
        '$InputFileRngBins        => '+self.cfg['$RadarRngBins']+'\n'+\
        '$StartRngBin             => '+self.cfg['$StartRngBinToProc']+'\n'+\
        '$RngBinsToProcess        => '+self.wrk['$RngSamples']+'\n'+\
        '$HeaderBytes             => 0\n'+\
        '$FooterBytes             => 0\n'+\
        '$InputDCOffsetI          => '+self.cfg['$DCOffsetI']+'\n'+\
        '$InputDCOffsetQ          => '+self.cfg['$DCOffsetQ']+'\n'+\
        '$InputIQRatio            => '+self.cfg['$IQRatio']+'\n\n'+\

        '--Misc (required)--\n'+\
        '$CarrierFreq [Hz - SRC,RW,MoC] => '+self.cfg['$RadarCarrierFreq']+'\n'+\
        '$StepFreqUserFile [note]       => '+self.wrk['$StepFreqUserFile']+'\n'+\
        '$A2DFreq [Hz]                  => '+self.wrk['$A2DFreq']+'\n'+\
        '$RngShiftInterpSize [RW,MoC]   => '+self.cfg['$RngShiftInterpSize']+'\n'+\
        '$Scale                         => '+self.cfg['$RngComScale']+'\n\n'+\

        '--Range compression specific (RC)--\n'+\
        '$RngComFlg [Y/N - SRC]         => '\
                     +self.cfg['$RngCompress'].upper()+'\n'+\
        '$RngComRefFuncPhaseSign [+-1]  => '+self.cfg['$RngComRefPhase']+'\n'+\
        '$RngComChirpBandwidth [Hz]     => '+self.cfg['$RadarChirpBandwidth']+'\n'+\
        '$RngComPulseLen [sec]          => '+self.cfg['$RadarPulseLength']+'\n'+\
        '$RngComWinConstTime            => '+self.cfg['$RngComWinConstTime']+'\n\n'+\

        '--Motion compensation specific (MoC)--\n'+\
        '$MoCompFlg [Y/N]                => '\
                      +self.cfg['$MoComp'].upper()+'\n'+\
        '$MoCompFileName                 => '+self.cfg['$MocFile']+'\n'+\
        '$MoCompRngShiftFlg [Y/N]        => '\
                      +self.cfg['$MoCompRngShiftFlg'].upper()+'\n'+\
        '$MoCompRngShiftSign [+-1]       => '+self.cfg['$MoCompRngShiftSign']+'\n'+\
        '$MoCompRngShiftIndex [note]     => 0\n'+\
        '$MoCompPhaseSign [+-1]          => '+self.cfg['$MoCompPhaseSign']+'\n'+\
        '$MoCompRngUpdates [note]        => 1\n'+\
        '$MoCompStartPRI [note]          => '+self.cfg['$StartG2PRI']+'\n\n'+\

        '--LMS interference suppression (LMS)--\n'+\
        '$LmsFlg [Y/N]                   => '+self.cfg['$LmsFlg']+'\n'+\
        '$LmsUpdateRate                  => '+self.cfg['$LmsUpdateRate']+'\n'+\
        '$LmsNumWeights                  => '+self.cfg['$LmsNumWeights']+'\n'+\
        '$LmsSidelobeOrder [note]        => '+self.cfg['$LmsSidelobeOrder']+'\n\n'+\

        '--Notch interference suppression (Notch)--\n'+\
        '$NotchFlg [Y/N]                 => '+self.cfg['$NotchFlg']+'\n'+\
        '$NotchUpdateRate                => '+self.cfg['$NotchUpdateRate']+'\n'+\
        '$NotchNumFFTLines [note]        => '+self.cfg['$NotchNumFFTLines']+'\n'+\
        '$NotchCutoff [dB - note]        => '+self.cfg['$NotchCutoff']+'\n'+\
        '$NotchMedianKernLen [note]      => '+self.cfg['$NotchMedianKernLen']+'\n\n'+\

        '--STC specific (STC)--\n'+\
        '$STCFlg [Y/N]                   => N\n'+\
        '$STCFileName [note]             => tmp.stc\n\n'+\

        '--SRC, Doppler centroid and range walk specific (SRC,DOPC,RW)--\n'+\
        '$SRCFlg [Y/N]                   => N\n'+\
        '$DopCentroid [Hz - note]        => 0.0\n'+\
        '$RngWalkRngShiftFlg [Y/N]       => N\n'+\
        '$RngWalkPhaseShiftFlg [Y/N]     => N\n'+\
        '$RngWalkAzBeamwidth [deg - RW]  => 60.0\n'+\
        '$SRCFocusRng [m - SRC]          => 0.0\n'+\
        '$NomGroundSpeed [m/s - SRC,RW]  => '+self.cfg['$ProcAveGrndSpeed']+'\n'+\
        '$SquintAngle [deg - SRC,RW]     => 0.0\n'+\
        '$InputPRF [Hz - RW,DOPC]        => '+self.wrk['$PRF']+'\n\n'+\
        
        'Notes - notes are not included in G2-created command files\n'+\
        '        (create template command file for module to view notes).\n'\
        )
        f.close()
        return 0


#--------------------------------
    def createStepFreqProcCmdFile(self,FileName):

        # Perform relevant calcs #

        f = open(FileName,'w')
        f.write(\
        'G2-created stepf command file (Stepped Freq Processing)\n'+\
        '$CmdFileVersion (rtl)    => 0.1\n'+\
        '-------------------------------------------------------\n\n'+\
      
        '--General--\n'+\
        '$LogFileName  (null for none)    => '+self.cfg['$StepFreqLog']+'\n'+\
        '$InputFileName                   => '+self.cfg['$RncFile']+'\n'+\
        '$OutputFileName                  => '+self.cfg['$StepFreqFile']+'\n'+\
        '$InputDataType [note]            => 3\n'+\
        '$PreSummedPulsesToUse [note]     => '+self.wrk['$AzSamples']+'\n'+\
        '$RngBinsToProcess                => '+self.wrk['$RngSamples']+'\n'+\
        '$InputDCOffsetI                  => 0.0\n'+\
        '$InputDCOffsetQ                  => 0.0\n'+\
        '$NarrowChirpBandwidth [Hz]       => '+self.cfg['$RadarChirpBandwidth']+'\n'+\
        '$NarrowPulseLen [sec]            => '+self.cfg['$RadarPulseLength']+'\n'+\
        '$InputA2DFreq [Hz]               => '+self.wrk['$A2DFreq']+'\n'+\
        '$InputStartSampleDelay [sec]     => '+self.wrk['$Delay0']+'\n'+\
        '$RngComWinConstTime              => '+self.cfg['$StepFreqWinConstTime']+'\n'+\
        '$RngComRefFuncPhaseSign [+-1]    => '+self.cfg['$RngComRefPhase']+'\n'+\
        '$StepFreqProcMode (normal/user)  => '+self.cfg['$StepFreqMode']+'\n\n'+\
        '--Stepped-Frequency (normal)--\n'+\
        '$StartCentreFrequency [Hz]       => '+self.cfg['$FirstStepCentreFreq']+'\n'+\
        '$StepSize [Hz]                   => '+self.cfg['$StepFreqStepSize']+'\n'+\
        '$NumberOfFreqSteps               => '+self.cfg['$NumberOfFreqSteps']+'\n\n'+\

        '--Stepped-Frequency (user)--\n'+\
        '$StepFreqUserFile [note]         => '+self.wrk['$StepFreqUserFile']+'\n\n'+\
        
        'Notes - notes are not included in G2-created command files\n'+\
        '        (create template command file for module to view notes).\n'\
        )
        f.close()
        return 0


#--------------------------              
    def createAzProcCmdFile(self,FileName):
        f = open(FileName,'w')
        f.write(\
        'G2-created azcom command file (SAR Azimuth Compression)\n'+\
        '$ProgramVersion (jmh)    => 1.1\n'+\
        '-------------------------------------------------------\n\n'+\
        '$ScreenUpdateRate             => 20\n'+\
        '$LogFileName                  => '+self.cfg['$AzLog']+'\n'+\
        '$InputStartSampleDelay        => '+self.wrk['$Delay0']+'\n'+\
        '$CarrierFreq [Hz]             => '+self.cfg['$RadarCarrierFreq']+'\n'+\
        '$InputPRF [Hz]                => '+self.wrk['$PRF']+'\n'+\
        '$NomGroundSpeed [m/s]         => '+self.cfg['$ProcAveGrndSpeed']+'\n'+\
        '$InputFileAzPts               => '+self.wrk['$AzSamples']+'\n'+\
        '$StartProcessAzPt             => 0\n'+\
        '$AzPtsToProcess               => '+self.wrk['$AzSamples']+'\n'+\
        '$InputFileRngBins             => '+self.wrk['$RngSamples']+'\n'+\
        '$StartProcessRngBin           => 0\n'+\
        '$RngBinsToProcess             => '+self.wrk['$RngSamples']+'\n'+\
        '$InputDCOffsetI               => 0.0\n'+\
        '$InputDCOffsetQ               => 0.0\n'+\
        '$InvFFTSizeReduc [pow of 2]   => '+self.cfg['$AzComInvFFTSizeReduc']+'\n'+\
        '$InputFileName                => '+self.cfg['$CorFile']+'\n'+\
        '$OutputFileName               => '+self.cfg['$AzcFile']+'\n'+\
        '$AppendExistOutFileFlg [Y/N]  => N\n'+\
        '$RngFocSegments               => '+self.cfg['$RngFocSegments']+'\n'+\
        '$RefFuncSign [+-1]            => '+self.cfg['$AzComRefPhase']+'\n'+\
        '$A2DFreq [Hz]                 => '+self.wrk['$A2DFreq']+'\n'+\
        '$NomAzRes [m - see note]      => '+self.cfg['$AzComNomAzRes']+'\n'+\
        '$WinConstTime [0.0-1.0]       => '+self.cfg['$AzComWinConstTime']+'\n'+\
        '$NumLooks                     => '+self.cfg['$NumAzLooks']+'\n'+\
        '$LookOverlapFrac [0.0-1.0]    => '+self.cfg['$AzLookOverlapFrac']+'\n'+\
        '$WinConstFreq [0.0-1.0]       => '+self.cfg['$AzComWinConstFreq']+'\n'+\
        '$RngCurvInterpSize            => '+self.cfg['$AzComRngCurvInterpSize']+'\n'+\
        '$RngCurvBatchSize             => '+self.cfg['$AzComRngCurvBatchSize']+'\n'+\
        '$PostSumRatio                 => 1\n'+\
        '$DetectMethod                 => '+self.cfg['$AzProcDetectMethod']+'\n'+\
        '$InputDataType                => 3\n'+\
        '$OutputDataType               => 3\n'+\
        '$Scale                        => '+self.cfg['$AzComScale']+'\n'+\
        '$ReportMax [1/0]              => 1\n\n'\

        'Notes - notes are not included in G2-created command files.\n'+\
        '        (create template command file for module to view notes).\n'\
        )
        f.close()
        return 0   
