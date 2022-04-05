#!/usr/bin/python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2mocfilt.py
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
Mocomp preprocessing to combine IMU inertial data with DGPS
data. Does this by:
1) Reading in the unpacked DGPS and IMU files.
2) Interpolate the two data sets to same grid and subtract mean
   value from each.
3) Find difference curve and fit poly of order 2 (quadratic) 
   to this diff curve (effectively gives smooth estimate of difference).
4) Add difference poly back to the IMU data and also add the DGPS
   mean value to this (for correct image coords).  
5) Write to output file.

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-08-06
COMMENTS:
Initial version, modified from g2filt_dgps_imu.py (ver 0.1) in order
to smooth the DGPS data and perhaps incorporate inert data as well.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-13
COMMENTS:
This version based on g2analmoc.py (ver 0.2). First full working
version. Ignores the erroneous IMU height values.

VERSION/AUTHOR/DATE : 0.3 / Jasper Horrell / 1999-08-19
COMMENTS:
Add in motion comp (LBR) PRI offset. This is implemented here by
subtracting a LBR PRI offset from the LBR PRI ID which is written
to the output file.

VERSION/AUTHOR/DATE : 0.4 / Jasper Horrell / 2000-03-07
COMMENTS:
Made changes to the plots for more meaningful display (EUSAR2000
paper). No changes to main workings.

=========================================================

"""

import sys, string, time, Polynomial, Gnuplot, LeastSquares, Interpolation, grnn
from Numeric import *
from umath import *
from FFT import fft

### Functions ###
def polyEval(params,x):
    y = 0.0
    for i in range(len(params)):
        y = y + params[i]*pow(x,i)
    return y    

#################
### Main Prog ###
print '\n--------------------'
print 'Prog: g2mocfilt.py - Ver 0.4 (jmh)\n'  
if len(sys.argv) < 5:
    print 'Merges unpacked DGPS and unpacked IMU data for mocomp.'
    print 'USAGE: g2mocfilt.py [DGPS_unpkFile] [IMU_unpkFile] [OutFile] [MocOffset]'
    sys.exit()

MocPRIOffset = int(sys.argv[4]) 

## READ FROM UNPACKED DGPS FILE ### (unpacked in range of IMU records)
DGPSFileName = sys.argv[1]     
DGPSFile = open(DGPSFileName, 'r')
print 'Opened input unpacked DGPS file: ' + DGPSFileName

lines = 0L
field = string.split(DGPSFile.readline())
if not field:
    print 'ERROR - no lines read of DGPS file!!'
    sys.exit()

PRI, UTC, Lat, Long, Alt = [], [], [], [], []   # PRI of the LBR variety
while 1:  # Repeat for all lines in the DGPS file
    if not field:
        break
    lines = lines + 1L
    if len(field) != 5:   # check that expected length
        break
    PRI.append(float(field[0]))	
    UTC.append(float(field[1]))
    Lat.append(float(field[2]))
    Long.append(float(field[3]))
    Alt.append(float(field[4]))
 
    # Read next line of file
    field = string.split(DGPSFile.readline())

# tidy up
DGPSFile.close()	
print 'DGPS: records read: ' + `lines`
print 'DGPS: start and end PRI: '+`PRI[0]`+' / '+`PRI[len(PRI)-1]`


### READ FROM UNPACKED IMU FILE ###
IMUFileName = sys.argv[2]
IMUFile = open(IMUFileName, 'r')
print 'Opened input unpacked IMU file: ' + IMUFileName

IMULines = 0L
IMUFile.readline()  # skip header line

IMUPRI, IMULat, IMULong = [], [], []      # IMUPRI of the LBR variety
while 1:  # Repeat for all lines in the IMU file (ignore IMU alt values)
    field = string.split(IMUFile.readline())
    if not field:
        break
    IMUPRI.append(float(field[0]))	
    IMULat.append(float(field[2]))
    IMULong.append(float(field[3]))
    IMULines = IMULines + 1L

# tidy up
IMUFile.close()	
print 'IMU: records read: ' + `IMULines`
print 'IMU: start and end PRI: '+`IMUPRI[0]`+' / '+`IMUPRI[len(IMUPRI)-1]`



### START OF MAIN PROCESSING ### (now that all data read into lists)

OutFileName = sys.argv[3]
OutFile = open(OutFileName, 'w')
print 'Opened output file: ' + OutFileName

# select a subregion (useful for testing)
start, end = 0, len(PRI)    # DGPS
PRI = PRI[start:end]
UTC = UTC[start:end]
Lat = Lat[start:end]
Long = Long[start:end]
Alt = Alt[start:end]
IMUstart, IMUend = 0, len(IMUPRI)  # IMU
IMUPRI = IMUPRI[IMUstart:IMUend]
IMULat = IMULat[IMUstart:IMUend]
IMULong = IMULong[IMUstart:IMUend]

# create arrays (fixed size) for the interpolation 
ArrPRI = zeros((len(PRI)),Float)  # create array
ArrUTC = zeros((len(PRI)),Float)
ArrLat = zeros((len(PRI)),Float)
ArrLong = zeros((len(PRI)),Float)
ArrAlt = zeros((len(PRI)),Float)
ArrIMUPRI = zeros((len(IMUPRI)),Float)  # create array
ArrIMULat = zeros((len(IMUPRI)),Float)
ArrIMULong = zeros((len(IMUPRI)),Float)

# Populate arrays for interp. Keep actual LBR PRI's at this stage
# and remove mean values from data sets. Does not use the IMU height
# values as these are incorrect.
sumUTCVal = 0.0
sumLatVal = 0.0
sumLongVal = 0.0
sumAltVal = 0.0
for i in range(len(PRI)):               # DGPS
    ArrPRI[i] = PRI[i]
    ArrUTC[i] = UTC[i]
    ArrLat[i] = Lat[i]
    ArrLong[i] = Long[i]
    ArrAlt[i] = Alt[i]
    sumUTCVal = sumUTCVal + ArrUTC[i]
    sumLatVal = sumLatVal + ArrLat[i]
    sumLongVal = sumLongVal + ArrLong[i]
    sumAltVal = sumAltVal + ArrAlt[i]
UTCMean = sumUTCVal/len(PRI)    
LatMean = sumLatVal/len(PRI)
LongMean = sumLongVal/len(PRI)
AltMean = sumAltVal/len(PRI)
ArrUTC = ArrUTC - UTCMean
ArrLat = ArrLat - LatMean  # subtract mean values for analysis
ArrLong = ArrLong - LongMean
ArrAlt = ArrAlt - AltMean

sumIMULatVal = 0.0
sumIMULongVal = 0.0
for i in range(len(IMUPRI)):            # IMU
    ArrIMUPRI[i] = IMUPRI[i]
    ArrIMULat[i] = IMULat[i]
    ArrIMULong[i] = IMULong[i]
    sumIMULatVal = sumIMULatVal + ArrIMULat[i] 
    sumIMULongVal = sumIMULongVal + ArrIMULong[i]
LatIMUMean = sumIMULatVal/len(IMUPRI)
LongIMUMean = sumIMULongVal/len(IMUPRI)    
ArrIMULat = ArrIMULat - LatIMUMean  # subtract mean value for analysis
ArrIMULong = ArrIMULong - LongIMUMean

# Find interp function (IF) for DGPS and IMU data (requires arrays as args)    
IFUTC = Interpolation.InterpolatingFunction((ArrPRI,),ArrUTC)
IFLat = Interpolation.InterpolatingFunction((ArrPRI,),ArrLat)
IFIMULat = Interpolation.InterpolatingFunction((ArrIMUPRI,),ArrIMULat)
IFLong = Interpolation.InterpolatingFunction((ArrPRI,),ArrLong)
IFIMULong = Interpolation.InterpolatingFunction((ArrIMUPRI,),ArrIMULong)
IFAlt = Interpolation.InterpolatingFunction((ArrPRI,),ArrAlt)

# Evaluate interp function on regular grid. The same grid is used for both 
# data sets with use made of the DGPS start and end PRI's as these lie within
# the IMU values as a result of the unpacking program. During this step, the
# start PRI is removed from the InterpIndex list.
minInterpPRI = PRI[0]
maxInterpPRI = PRI[end-1]
UTCInterpEval = []
LatInterpEval, IMULatInterpEval = [],[]
LongInterpEval, IMULongInterpEval = [],[]
AltInterpEval, InterpIndex = [],[]
TrackTime = []

step = 640   # 160 is 10 G2 PRI's for 4096 rng bins
for index in range(minInterpPRI,maxInterpPRI,step):
    UTCInterpEval.append(IFUTC(index))
    LatInterpEval.append(IFLat(index))
    IMULatInterpEval.append(IFIMULat(index))
    LongInterpEval.append(IFLong(index))
    IMULongInterpEval.append(IFIMULong(index))
    AltInterpEval.append(IFAlt(index))    
    InterpIndex.append(index-minInterpPRI)  # add start PRI back later
    TrackTime.append((index-minInterpPRI)/(16.0*136.3636364))

print 'NB: track time for plots assumes LBR PRI incr 16, PRF 136.364 Hz!'    

"""
# plot interp data
gLat = Gnuplot.Gnuplot()
gLatInterpData = Gnuplot.Data(TrackTime,LatInterpEval,\
                 title='Interp DGPS latitude', with='lines')
gIMULatInterpData = Gnuplot.Data(TrackTime,IMULatInterpEval,\
                    title='Interp inertial latitude', with='lines')
gLat.ylabel('Latitude (deg with zero mean)')
gLat.xlabel('Track Time (secs)')                    
gLat.plot(gLatInterpData,gIMULatInterpData)
gLat.hardcopy(filename='lat_IMU_DGPS.ps',color=1,largefont=1)
time.sleep(5.0)

gLong = Gnuplot.Gnuplot()
gLongInterpData = Gnuplot.Data(TrackTime,LongInterpEval,\
                 title='Interp DGPS longitude', with='lines')
gIMULongInterpData = Gnuplot.Data(TrackTime,IMULongInterpEval,\
                    title='Interp inertial longitude', with='lines')
gLong.ylabel('Longitude (deg with zero mean)')                    
gLong.xlabel('Track Time (secs)')
gLong.plot(gLongInterpData,gIMULongInterpData)
gLong.hardcopy(filename='long_IMU_DGPS.ps',color=1,largefont=1)
time.sleep(5.0)
"""

# Calculate the  difference curves
LatDiff,LongDiff = [],[]
for i in range(len(InterpIndex)):
    LatDiff.append(LatInterpEval[i] - IMULatInterpEval[i])
    LongDiff.append(LongInterpEval[i] - IMULongInterpEval[i])


# Fit quadratic to lat and long diff curves and plot
order = 2
ballprk = (0.,0.,0.)
LatDiff2,LongDiff2 = [],[]
fitskip = 2  # speed up for poly fitting
fitpts = int(len(InterpIndex)/fitskip)
for i in range(fitpts):  # create tuples for poly fitting
    k = fitskip*i
    LatDiff2.append((InterpIndex[k],LatDiff[k]))
    LongDiff2.append((InterpIndex[k],LongDiff[k]))
LatDiffPoly = LeastSquares.polynomialLeastSquaresFit(order,ballprk,LatDiff2)
LongDiffPoly = LeastSquares.polynomialLeastSquaresFit(order,ballprk,LongDiff2)
LatDiffPolyEval = []
LongDiffPolyEval = []
for i in range(len(InterpIndex)): # evaluate polys at the interp points
    LatDiffPolyEval.append(polyEval(LatDiffPoly[0],InterpIndex[i]))
    LongDiffPolyEval.append(polyEval(LongDiffPoly[0],InterpIndex[i]))

# plot difference curves
gLatDiff = Gnuplot.Gnuplot()
gLatDiffData = Gnuplot.Data(TrackTime,LatDiff,\
               title='Difference latitude', with='lines')
gLatDiffPolyData = Gnuplot.Data(TrackTime,LatDiffPolyEval,\
                   title='Difference latitude poly', with='lines')
gLatDiff.ylabel('Latitude Difference (deg)')
gLatDiff.xlabel('Track Time (secs)') 
gLatDiff.plot(gLatDiffData,gLatDiffPolyData)
#gLatDiff.hardcopy(filename='lat_diff.ps',color=1,largefont=1)

gLongDiff = Gnuplot.Gnuplot()
gLongDiffData = Gnuplot.Data(TrackTime,LongDiff,\
                title='Difference longitude', with='lines')
gLongDiffPolyData = Gnuplot.Data(TrackTime,LongDiffPolyEval,\
                    title='Difference longitude poly', with='lines')
gLongDiff.ylabel('Latitude Difference (deg)')
gLongDiff.xlabel('Track Time (secs)') 
gLongDiff.plot(gLongDiffData,gLongDiffPolyData)
#gLongDiff.hardcopy(filename='long_diff.ps',color=1,largefont=1)

# set up alt smoothing
SmoothAlt = grnn.GenRegNN(InterpIndex,AltInterpEval)
SmoothAlt.stddev = 5000
#SmoothAlt.display(100,'sleep')

# Form output corrected data (incl PRI offset) and write out to ASCII file
OutPRI,OutUTC,OutLat,OutLong,OutAlt,OutAltPlot = [],[],[],[],[],[]
for i in range(len(InterpIndex)):
    OutPRI.append(InterpIndex[i]+minInterpPRI-MocPRIOffset)
    OutUTC.append(UTCInterpEval[i]+UTCMean)
    OutLat.append(IMULatInterpEval[i]+LatDiffPolyEval[i]+LatMean)
    OutLong.append(IMULongInterpEval[i]+LongDiffPolyEval[i]+LongMean)   
#    OutAlt.append(AltInterpEval[i]+AltMean)  # original interpolated alt
    OutAlt.append(SmoothAlt.eval(InterpIndex[i])+AltMean) # smoothed alt
    OutAltPlot.append(SmoothAlt.eval(InterpIndex[i])) # used for plot  
    OutFile.write(`int(OutPRI[i])`+' '+`OutUTC[i]`+' '+`OutLat[i]`+' '+\
               `OutLong[i]`+' '+`OutAlt[i]`+'\n')

# plot the smoothed altitude
gAlt = Gnuplot.Gnuplot()
gAltInterpData = Gnuplot.Data(TrackTime,AltInterpEval,\
                 title='Interp DGPS altitude', with='lines')
gAltSmoothData = Gnuplot.Data(TrackTime,OutAltPlot,\
                 title='Smoothed DGPS altitude', with='lines')
gAlt.ylabel('DGPS Altitude (m with zero mean)')
gAlt.xlabel('Track Time (secs)')
gAlt.plot(gAltInterpData,gAltSmoothData)
#gAlt.hardcopy(filename='alt_smooth.ps',color=1,largefont=1)
time.sleep(5.0)
    
    
# Tidy up
OutFile.close()

"""
gOutLatData = Gnuplot.Data(InterpIndex,OutLat,\   # remove mean 
                 title='Output lat data', with='lines')
gLatFinal = Gnuplot.Gnuplot()
gLatFinal.plot(gLatInterpData,gOutLatData)
gOutLongData = Gnuplot.Data(InterpIndex,OutLong,\ # remove mean
                 title='Output long data', with='lines') 
gLongFinal = Gnuplot.Gnuplot()
gLongFinal.plot(gLongInterpData,gOutLongData)
"""

#print 'Press return to continue...'
#sys.stdin.readline()
		
sys.exit(0)  # return zero on success 	
