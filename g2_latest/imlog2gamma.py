#!/usr/bin/env python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group (1999)
FILE NAME: imlog2gamma.py
CODE CONTROLLER: Jasper Horrell (jasper@eng.uct.ac.za)
DESCRIPTION:
The G2 SAR processor top-level controlling module.

VERSION/AUTHOR/DATE : 0.1 /
                      Jasper Horrell and Andrew Wilkinson /
                      1999-10-04
COMMENTS:
Python program to write out a header file for the gamma processor
using the info supplied in the G2 processor's image log file.

VERSION/AUTHOR/DATE :  0.2 / Richard Lord / 16 March 2001
COMMENTS: Due to stepped-frequency upgrade of G2 processor, the imlog
          file has changed slightly.

VERSION/AUTHOR/DATE :  /  /
COMMENTS:

e.g.

python  imlog2gamma.py  mtar-000-vv-3.imlog test1.par

=========================================================
"""

import sys, string , navtools
from math import *
from Numeric import *

#---Misc functions---#

def ErrorExit(cause):
    print '\nERROR - '+cause+'!!'
    sys.exit(-1)



#---Class definition---#

class LogUtil:

    def __init__(self, ConfigDic={}):
        self.cfg = ConfigDic

    # Parse an image log file and add keys (params) and
    # values to the log object's dictionary.
    def addLogFile(self,FileName):
        f = open(FileName,'r')
        while 1:
            str1 = f.readline()
            if not str1:
                break
            if string.find(str1,'=>') != -1:   # look for lines with '=>'
                field = string.split(str1)
                for i in range(len(field)):
                    if field[i] == '=>':  # find position of '=>' in line
                        self.cfg[field[0]] = field[i+1]
                        break
        f.close()
        return 0

    def displayContents(self):
        for i in range(len(self.cfg)):
            print i, self.cfg.keys()[i], ' : ',self.cfg.values()[i]   

    def writeGammaHeader(self,GamFileName):       
        
        # write header file
        HeadF = open(GamFileName,'w')
        HeadF.write(\
                    'Gamma Interferometric SAR Processor - Image Data Parameter File\n'\
                    'Created by utility imlog2gamma.py - JH & AJW 8.10.1999\n'\
                    'This file contains SASAR parameters extracted from an imlog file\n\n'\
        )
        HeadF.write('title:    '+self.cfg['Raw_data_ID']+'\n')
        HeadF.write('sensor:    SASAR\n')

        # Extract date from file ID
        ID = self.cfg['Raw_data_ID']
        year = ID[0:4]
        month = ID[4:6]
        day = ID[6:8]
        HeadF.write('date:      '+year+' '+month+' '+day+'\n')   # FIX

        c = 2.997E8
        azimuth_lines = int(self.cfg['Image_azimuth_samples'])
        range_pixel_spacing = float(self.cfg['Image_range_sample_spacing'])
        adc_sample_rate = c/(2.0*range_pixel_spacing)                                 # Currently calculated but should be given in imlog file
        azimuth_pixel_spacing = float(self.cfg['Image_azimuth_sample_spacing'])
        speed = float(self.cfg['Proc_average_ground_speed'])                # is this correct????
        prf = speed/azimuth_pixel_spacing                                   # must change later
        start_time = float(self.cfg['Mocomp_start_DGPS_UTC'])
        center_time = start_time + (azimuth_lines/2.0)/prf
        end_time = start_time + azimuth_lines/prf

        HeadF.write('start_time:                  '+`start_time`+' s - currently using Mocomp_start_DGPS_UTC\n')                     # FIX in log file
        HeadF.write('center_time:                 '+`center_time`+' s\n')
        HeadF.write('end_time:                    '+`end_time`+' s\n')
        HeadF.write('line_header_size:            0\n')
        HeadF.write('range_samples:               '+self.cfg['Image_range_samples']+'\n')
        HeadF.write('azimuth_lines:               '+`azimuth_lines`+'\n')
        HeadF.write('range_looks:                       1\n')  # This is an SLC file
        HeadF.write('azimuth_looks:                     1\n')
        HeadF.write('range_scale_factor:        1.0000000\n')
        HeadF.write('azimuth_scale_factor:      1.0000000\n')

        # Compute mid swath lat/long
        Image_near_swath_start_latitude = float( self.cfg['Image_near_swath_start_latitude'] )
        Image_near_swath_start_longitude = float( self.cfg['Image_near_swath_start_longitude'] )
        Image_far_swath_end_latitude = float( self.cfg['Image_far_swath_end_latitude'] )
        Image_far_swath_end_longitude = float( self.cfg['Image_far_swath_end_longitude'] )
        center_latitude = (Image_near_swath_start_latitude + Image_far_swath_end_latitude)/2.0
        center_longitude = (Image_near_swath_start_longitude + Image_far_swath_end_longitude)/2.0

        HeadF.write('center_latitude:             '+`center_latitude`+'  degrees\n')       #  check
        HeadF.write('center_longitude:            '+`center_longitude`+'  degrees\n')      # check
        HeadF.write('heading:                            0.0     degrees\n')                                # FIX
        HeadF.write('range_pixel_spacing:          '+`range_pixel_spacing`+' m\n')
        HeadF.write('azimuth_pixel_spacing:        '+`azimuth_pixel_spacing`+' m\n')

        # Obtain near,mid & far ranges
        near_range_slc = float(self.cfg['Image_near_slant_range'])
        far_range_slc = float(self.cfg['Image_far_slant_range']) 
        center_range_slc = (near_range_slc + far_range_slc)/2.0
     
        HeadF.write('near_range_slc:               '+`near_range_slc`+' m\n')
        HeadF.write('center_range_slc:          '+`center_range_slc`+' m\n')
        HeadF.write('far_range_slc:             '+`far_range_slc`+' m\n')      

        # Calculate incidence angle at mid swath
        Image_platform_average_altitude = float(self.cfg['Image_platform_average_altitude'])
        Image_terrain_altitude = float(self.cfg['Image_terrain_altitude'])
        height = Image_platform_average_altitude - Image_terrain_altitude
        incidence_angle = acos(height / center_range_slc) / 3.14159 * 180.0

        HeadF.write('incidence_angle:           '+`incidence_angle`+'  degrees\n')          # Computed at mid swath!

        HeadF.write('azimuth_deskew:          ON\n')                               # Check this parameter!!!
        HeadF.write('azimuth_angle:               90.0000   degrees\n')            # ZERO DOPPLER PROCESSING
        HeadF.write('radar_frequency:         '+self.cfg['Radar_carrier_frequency']+' Hz\n')
        HeadF.write('adc_sampling_rate:        '+`adc_sample_rate`+' Hz\n') #  24.0e+06   Hz\n')               # Calculated from range sampe spacing in m
        HeadF.write('chirp_bandwidth:           '+self.cfg['Radar_chirp_bandwidth']+' Hz\n')
        HeadF.write('prf:                       '+`prf`+' Hz\n')                                                 # !!! Obtained by velocity/az pixel spacing\n')
        HeadF.write('azimuth_proc_bandwidth:     '+self.cfg['Image_Doppler_bandwidth_processed']+' Hz\n')             # Bandwidth of azimuth data in Hz
        HeadF.write('doppler_polynomial:        00e+00  0.00000e+00  0.00000e+00  0.00000e+00  Hz Hz/m Hz/m^2 Hz/m^3\n')         # FIX
        HeadF.write('receiver_gain:                0.0000   dB\n')
        HeadF.write('calibration_gain:             0.0000   dB\n')



        # Calculate state vector positions and velocities at start and end of MOCOMP line
        
        P1_lat = float(self.cfg['Image_platform_ref_track_start_lat'])/180*pi           # convert to radians
        P1_long = float(self.cfg['Image_platform_ref_track_start_long'])/180*pi
        P1_alt = float(self.cfg['Image_platform_ref_track_start_alt'])
        P1_geod = array([P1_lat,P1_long,P1_alt])

        P3_lat = float(self.cfg['Image_platform_ref_track_end_lat'])/180*pi
        P3_long = float(self.cfg['Image_platform_ref_track_end_long'])/180*pi
        P3_alt = float(self.cfg['Image_platform_ref_track_end_alt'])
        P3_geod = array([P3_lat,P3_long,P3_alt])

        P1 = array(navtools.ToCartes(P1_geod))                 # Convert to XYZ coordinates
        P3 = array(navtools.ToCartes(P3_geod))

        # Calculate mid point vector i.e. P2
        P2 = P1 + (P3-P1)/2

        t1 = start_time     # i.e. use same times as start range line and end range lines
        t2 = (start_time+end_time)/2
        t3 = end_time

        # Vel = (P3-P1)/(t3-t1)
        V1 = (P3-P1)/(t3-t1)
        V2 = V1
        V3 = V1


        # Calculate sensor and scene radii from earth centre

        sar_to_earth_center = sqrt( dot(P1,P1) )    # use height at start of image

        L = navtools.ToCartes([P1_lat,P1_long,Image_terrain_altitude])                # Convert to XYZ coordinat

        earth_radius_below_sensor =  sqrt( dot(L,L) )

        HeadF.write('sar_to_earth_center:             '+`sar_to_earth_center` +'  m\n')
        HeadF.write('earth_radius_below_sensor:       '+`earth_radius_below_sensor`+'  m\n')
        HeadF.write('earth_semi_major_axis:           6378137.0000   m\n')
        HeadF.write('earth_semi_minor_axis:           6356752.3141   m\n')

        # I found by trial and error that GAMMA required a minimum of 3 state vectors for the baseline 
        # so I generated a third one at the mid point!

        HeadF.write('number_of_state_vectors:                    3\n')                 # Start and end coords of MOCOMP line
        HeadF.write('time_of_first_state_vector:          '+`t1`+'   s\n')
        HeadF.write('state_vector_interval:               '+`t2-t1`+'  s\n')
        HeadF.write('state_vector_position_1:   '+`P1[0]`+'  '+`P1[1]`+'  '+`P1[2]`+'  m   m   m\n')
        HeadF.write('state_vector_velocity_1:   '+`V1[0]`+'  '+`V1[1]`+'  '+`V1[2]`+'  m/s m/s m/s\n')
        HeadF.write('state_vector_position_2:   '+`P2[0]`+'  '+`P2[1]`+'  '+`P2[2]`+'  m   m   m\n')
        HeadF.write('state_vector_velocity_2:   '+`V2[0]`+'  '+`V2[1]`+'  '+`V2[2]`+'  m/s m/s m/s\n')
        HeadF.write('state_vector_position_3:   '+`P3[0]`+'  '+`P3[1]`+'  '+`P3[2]`+'  m   m   m\n')
        HeadF.write('state_vector_velocity_3:   '+`V3[0]`+'  '+`V3[1]`+'  '+`V3[2]`+'  m/s m/s m/s\n')
        HeadF.write('\n')
        HeadF.write('*************** END OF SLC DATA PARAMETERS ******************\n')

        Vel = sqrt( dot(V1,V1) )
        HeadF.write('Calculated Velocity1 (i.e. from State vectors):  '+`Vel` +'  m/s\n')
        d = P2-P1
        HeadF.write('Displacement vector P2-P1:  '+'d'+'  m\n')
        dist  = sqrt( dot(d,d) )
        HeadF.write('Displacement magnitude  |P2-P1|:  '+`dist` +'  m\n')
        HeadF.write('Height above ground:  '+`height`+'  m\n')
        HeadF.close()
        return 0


# azsamples = int(self.cfg['$ImAzSamples'])
# azsamples = azsamples + 1
# self.cfg['$ImAzSamples'] = `azsamples`


#------------#
# MAIN PROGRAM


print '\n------------------------------------------------------------------------'
print 'Prog: imlog2gamma.py - Ver 0.2 (jmh and ajw)\n'
print 'This utility converts the SASAR .imlog file to a GAMMA processor .par file'
print '\n------------------------------------------------------------------------'
if len(sys.argv) < 3:
    print 'Creates header file for Gamma processor based on SASAR image log.'
    print 'USAGE: saslog2gamma.py [G2_Image_Log_File] [Gamma_Header_File]'
    sys.exit(-1)

ImLogName = sys.argv[1]
GamHeaderName = sys.argv[2]

print '\nReading File: '+ImLogName + '\n'

# Create new log file utility object and parse image log file.
L = LogUtil()    # create a new log file utility object
if L.addLogFile(ImLogName) != 0: # add contents of proc config file to object
    ErrorExit('in parsing image log file '+ImLogName)

# L.displayContents()

print 'Creating File: '+GamHeaderName + '\n'

L.writeGammaHeader(GamHeaderName)

print 'Done!\n'