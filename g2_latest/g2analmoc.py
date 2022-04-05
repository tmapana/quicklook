#!/usr/bin/python

"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2analmoc.py
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
Renamed from g2filt_dgps_imu.py. Used to plot motion spectra etc.
Part of prototyping stage with eventual goal of developing a mocomp
prefilter. Uses both DGPS and IMU data. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-08-05
COMMENTS:
Initial version - only DGPS data

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-11
COMMENTS:
Incl IMU data and redo substantially. Incorporate interpolating function.

=========================================================
"""

import sys, string, Polynomial, Gnuplot, LeastSquares, Interpolation
from Numeric import *
from umath import *
from FFT import fft

print '\n--------------------'
print 'Prog: g2analmoc.py - Ver 0.1 (jmh)\n'  
if len(sys.argv) != 3:
    print 'USAGE: g2analmoc.py [DGPS_unpkFile] [IMU_unpkFile]'
    sys.exit()

### PROCESS UNPACKED DGPS FILE ### (unpacked in range of IMU records)
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
print 'DGPS: records processed: ' + `lines`
print 'DGPS: start and end PRI: '+`PRI[0]`+' / '+`PRI[len(PRI)-1]`


### PROCESS UNPACKED IMU FILE ###
IMUFileName = sys.argv[2]
IMUFile = open(IMUFileName, 'r')
print 'Opened input unpacked IMU file: ' + IMUFileName

IMULines = 0L
IMUFile.readline()  # skip header line

IMUPRI, IMULat, IMULong = [], [], []      # IMUPRI of the LBR variety
while 1:  # Repeat for all lines in the IMU file
    field = string.split(IMUFile.readline())
    if not field:
        break
    IMUPRI.append(float(field[0]))	
    IMULat.append(float(field[2]))
    IMULong.append(float(field[3]))
    IMULines = IMULines + 1L

# tidy up
IMUFile.close()	
print 'IMU: records processed: ' + `IMULines`
print 'IMU: start and end PRI: '+`IMUPRI[0]`+' / '+`IMUPRI[len(IMUPRI)-1]`


### START OF MAIN PROCESSING ### (now that all data read into lists)

start, end = 0, len(PRI)    # DGPS
PRI = PRI[start:end]
Lat = Lat[start:end]
Long = Long[start:end]
Alt = Alt[start:end]

IMUstart, IMUend = 0, len(IMUPRI)  # IMU
IMUPRI = IMUPRI[IMUstart:IMUend]
IMULat = IMULat[IMUstart:IMUend]
IMULong = IMULong[IMUstart:IMUend]


# create arrays (fixed size) for the interpolation 
ArrPRI = zeros((len(PRI)),Float)  # create array
ArrLat = zeros((len(PRI)),Float)

ArrIMUPRI = zeros((len(IMUPRI)),Float)  # create array
ArrIMULat = zeros((len(IMUPRI)),Float)

# populate arrays for interp
sumVal = 0.0
for i in range(len(PRI)):               # DGPS
    ArrPRI[i] = PRI[i]
    ArrLat[i] = Lat[i]
    sumVal = sumVal + ArrLat[i] 
ArrLat = ArrLat - sumVal/len(PRI)  # subtract mean value for analysis

sumIMUVal = 0.0
for i in range(len(IMUPRI)):            # IMU
    ArrIMUPRI[i] = IMUPRI[i]
    ArrIMULat[i] = IMULat[i]
    sumIMUVal = sumIMUVal + ArrIMULat[i] 
ArrIMULat = ArrIMULat - sumIMUVal/len(IMUPRI)  # subtract mean value for analysis

# find interp function (IF) for DGPS and IMU data     
IFLat = Interpolation.InterpolatingFunction((ArrPRI,),ArrLat)  # needs arrays as args
IFIMULat = Interpolation.InterpolatingFunction((ArrIMUPRI,),ArrIMULat)

# evaluate interp function on regular grid (same grid for both data sets)
minInterpPRI = PRI[0]
maxInterpPRI = PRI[end-1]
LatInterpEval = []
IMULatInterpEval = []

# interp both data sets on same grid
step = 320
for index in range(minInterpPRI,maxInterpPRI,step):  # DGPS
    LatInterpEval.append(IFLat(index))
    IMULatInterpEval.append(IFIMULat(index))
index = arrayrange(minInterpPRI,maxInterpPRI,step)

# plot orig and interp
gLat = Gnuplot.Gnuplot()
#gLatData = Gnuplot.Data(PRI,Lat,title='Orig DGPS lat data', with='lines')
gLatInterpData = Gnuplot.Data(index,LatInterpEval,\
                 title='Interp DGPS lat data', with='lines')
#gIMULatData = Gnuplot.Data(IMUPRI,IMULat,title='Orig IMU lat data', with='lines')
gIMULatInterpData = Gnuplot.Data(index,IMULatInterpEval,\
                    title='Interp IMU lat data', with='lines')
#gLat.plot(gIMULatInterpData,gLatInterpData)
gLat.plot(gLatInterpData,gIMULatInterpData)

# plot diff curv
gdiff = Gnuplot.Gnuplot()
diff = []
for i in range(len(LatInterpEval)):
    diff.append(LatInterpEval[i] - IMULatInterpEval[i])

gDiffData = Gnuplot.Data(index,diff,title='Diff', with='lines')
gdiff.plot(gDiffData)


"""
# apply window, zero pad, find fft and plot
DGPSlatpad = []
IMUlatpad = []
for i in range(len(LatInterpEval)):
    winc = 0.08+0.92*math.pow( math.sin( (math.pi*i)/(len(LatInterpEval)-1.0 )),2.0)
    DGPSlatpad.append(LatInterpEval[i]*winc)
    IMUlatpad.append(IMULatInterpEval[i]*winc)

gwin = Gnuplot.Gnuplot()
gwindata = Gnuplot.Data(index,DGPSlatpad,\
           title='DGPS lat win data',with='linespoints')
gwinIMUdata = Gnuplot.Data(index,IMUlatpad,\
           title='IMU lat win data',with='linespoints')           
gwin.plot(gwindata,gwinIMUdata)


padpts = 10000
indx2 = arrayrange(0,padpts+len(index),1)

for i in range(padpts):
    DGPSlatpad.append(0.0)
    IMUlatpad.append(0.0)

minx = 0
maxx = int(len(indx2)/200)
gfft = Gnuplot.Gnuplot()
gfftData = Gnuplot.Data(indx2[minx:maxx],abs(fft(DGPSlatpad))[minx:maxx],\
           title='Freq domain DGPS data (padded)',with='linespoints')
gfftIMUData = Gnuplot.Data(indx2[minx:maxx],abs(fft(IMUlatpad))[minx:maxx],\
              title='Freq domain IMU data (padded)',with='linespoints')
gfft.plot(gfftIMUData,gfftData)



# plot fft without padding
minx = 0
maxx = int(len(index)/200)
gfft = Gnuplot.Gnuplot()
gfftData = Gnuplot.Data(index[minx:maxx],abs(fft(LatInterpEval))[minx:maxx],\
           title='Freq domain DGPS data',with='linespoints')
gfftIMUData = Gnuplot.Data(index[minx:maxx],abs(fft(IMULatInterpEval))[minx:maxx],\
              title='Freq domain IMU data',with='linespoints')
gfft.plot(gfftIMUData,gfftData)
"""

"""
###### Poly fitting stuff (for best results, remove start values)

def polyEval(params,x):
    y = 0.0
    for i in range(len(params)):
        y = y + params[i]*pow(x,i)
    return y    

# set up empty lists for poly fitting
Lat2, long2, Alt2 = [], [], []

# populate arrays for interp
for i in range(len(PRI)):               # DGPS
    Lat2.append((PRI[i]-PRI[0],Lat[i]-Lat[0]))   # create list of tuples for poly fits
    Long2.append((PRI[i]-PRI[0],Long[i]-Long[0]))
    Alt2.append((PRI[i]-PRI[0],Alt[i]-Alt[0]))

order2 = 4
ballprk = (0.,0.,0.,0.,0.)

LatPoly2 = LeastSquares.polynomialLeastSquaresFit(order2,ballprk,Lat2)
LongPoly2 = LeastSquares.polynomialLeastSquaresFit(order2,ballprk,Long2)
AltPoly2 = LeastSquares.polynomialLeastSquaresFit(order2,ballprk,Alt2)

LatPolyEval = []
LatPoly2Eval = []
LongPolyEval = []
LongPoly2Eval = []
AltPolyEval = []
AltPoly2Eval = []
for i in range(len(PRI)):
    LatPoly2Eval.append(polyEval(LatPoly2[0],PRI[i]))
    LongPoly2Eval.append(polyEval(LongPoly2[0],PRI[i]))
    AltPoly2Eval.append(polyEval(AltPoly2[0],PRI[i]))

gLat = Gnuplot.Gnuplot()
gLatData = Gnuplot.Data(PRI,Lat,title='Orig lat data', with='linespoints')
gLatPoly2Data = Gnuplot.Data(PRI,LatPoly2Eval,title='Poly2 lat data',with='linespoints')
gLong = Gnuplot.Gnuplot()
gLongData = Gnuplot.Data(PRI,Long,title='Orig long data', with='lines')
gLongPoly2Data = Gnuplot.Data(PRI,LongPoly2Eval,title='Poly2 long data',with='linespoints')
gAlt = Gnuplot.Gnuplot()
gAltData = Gnuplot.Data(PRI,Alt,title='Orig alt data', with='linespoints')
gAltPoly2Data = Gnuplot.Data(PRI,AltPoly2Eval,title='Poly2 alt data', with='linespoints')

gLat.plot(gLatData,gLatPoly2Data)
gLong.plot(gLongData,gLongPoly2Data)
gAlt.plot(gAltData,gAltPoly2Data)

"""

print 'Press return to continue...'
sys.stdin.readline()
		
	
