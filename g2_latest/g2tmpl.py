"""
=========================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2tmpl.py
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
The G2 SAR processor template writing functions. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-07-30
COMMENTS:
Initial python version.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-08-04
COMMENTS:
Calc FFT sizes automagically. Radar cfg file name moved to
proc cfg file.

VERSION/AUTHOR/DATE : 0.3 / Jasper Horrell / 1999-08-14
COMMENTS:
Remove mocomp prelims enable and add in separate mocomp 
enables. Also include mocomp IMU/DGPS merge option.

VERSION/AUTHOR/DATE : 0.4 / Jasper Horrell / 1999-08-20
COMMENTS:
Add in mocomp data time offset in radar config file.
Rearrange config files. 

VERSION/AUTHOR/DATE : 0.5 / Jasper Horrell / 1999-10-01
COMMENTS:
Allow notch interference suppression. Allow output endian and
orientation specification.

VERSION/AUTHOR/DATE : 0.6 / Jasper Horrell / 2000-02-16
COMMENTS:
Allow Flt2ByteMath field.

VERSION/AUTHOR/DATE : 0.7 / Jasper Horrell / 2000-03-24
COMMENTS:
Add fields for stepped frequency processing. Add "off"
option for certain modules in Proc Control section of
proc config. This (with the "n" option) allows greater
flexibility when restarting a run half-way through and
parameters have to be calculated as if the module had
been run. Add printHelp() function.

VERSION/AUTHOR/DATE : 0.8 / Jasper Horrell / 2000-08-04
COMMENTS:
Write out User Manual as HTML file.

VERSION/AUTHOR/DATE : 0.9 / Jasper Horrell / 2001-02-13
COMMENTS:
Minor changes to HTML user manual.


=========================================================
"""

import string, sys


#----------------------------------
# radar configuration file template

def createRadarCfgTmpl(fileName):
    outf = open(fileName,'w')
    outf.write("""SASAR VHF configuration file for G2 Processor
$RadarConfigVersion => 0.4    
=============================================
$DataID                                   => 19990721-tznn009-vv
$RawDataFile (no path)                    => tznn-h3.009
$LBRFile (no path)                        => tznn.009
$DGPSFile (no path)                       => tzaneen.gps
$RawDataType (byte/float - see note)      => byte
$DataStartLP (LBR PRI for raw data start) => 2415152
$RadarDelayToStartSample (secs)           => 10.6667e-06
$RadarPRF (as in raw data file - Hz)      => 136.363636
$RadarAzSamples (in raw data file)        => 44562
$RadarRngBins (in raw data file)          => 4096
$RadarA2DFreq (Hz)                        => 24.0e+06
$RadarCarrierFreq (nominal - Hz)          => 141.0e+06
$RadarPulseLength (sec)                   => 6.66667e-06
$RadarChirpBandwidth (Hz - zero for mono) => 12.0e+06
$RadarMocTimeOffset (sec)                 => 0.5
$TerrainAlt (m)                           => 1300

--Stepped Freq Setup (see note) --
$StepFreqMode (no/normal/user)            => no
$NumberOfFreqSteps (normal)               => 1
$FirstStepCentreFreq (normal - Hz)        => 141.0e+06
$StepFreqStepSize (normal - Hz)           => 12.0e+06
$StepFreqUserFile (user - no path)        => null

--Optional Params (else 'null' - see note)--
$DCOffsetI                                => null
$DCOffsetQ                                => null
$IQRatio                                  => null
$AveGroundSpeed (m/s)                     => null

Notes
=====

RawDataType - byte (8-bit I, 8-bit Q, unsigned char)
              float (32-bit IEEE float I and Q, litte endian)

Stepped Freq Setup:
    For stepped freq operation, the "RadarPulseLen" and
    "RadarChirpBandwidth" parameters are taken to be those of each
    narrow band pulse (assumed constant for a run). The step freq
    user file first line should be an integer which is the number
    of freq steps with the centre frequency of each step (in Hz)
    on subsequent lines.

Optional Params:
    These values are only used if the processor module where they
    are normally calculated is set to 'off'. For example the DC
    offsets and IQ ratio are usually calculated in the SniffDC
    module. The average ground speed is calculated in the motion
    compensation calculation module.""")
    outf.close()
    return 0


#--------------------------------------
# processor configuration file template

def createProcCfgTmpl(fileName):
    outf = open(fileName,'w')
    outf.write("""G2 Processor config for SASAR VHF
$ProcConfigVersion => 0.7
=================================

---General Parameters---
$RunID (root of output file names)        => ./run1-vv-1
$InputPath (location of raw data files)   => /scratch/sasar/
$RadarCfgFile (path optional)             => /scratch/sasar/run1-vv.cfg
$StartG2PRI (0 is start of raw data)      => 0
$InputPRIsToProc                          => 44560
$ProcPresumRatio                          => 5
$StartRngBinToProc                        => 0
$RngBinsToProc                            => 4096

---Proc Control [see note] ---
$EnableUnpackIMU (y/n/only/off)           => y
$EnableUnpackDGPS (y/n/only/off)          => y
$EnableMergeMocData (y/n/only/off)        => y
$EnableMocompCalc (y/n/only/off)          => y
$EnablePlotMotionError (y/n/only/off)     => y
$EnableSniffDC (y/n/only/off)             => y
$EnableRngProc (y/n/only)                 => y
$EnableStepFreqProc (y/n/only/off)        => y
$EnableCornerTurn (y/n/only)              => y
$EnableAzProc (y/n/only)                  => y
$EnableFloat2Tiff (y/n/only/off)          => y
$EnableEndianSwap (y/n/only/off)          => y
$EnableOrient (y/n/only/off)              => y
$EnableImageLog (y/n/only/off)            => y
$EnableCleanUp (y/n/only/off)             => y

---RngProcStageParams---
$RngCompress (y/n)                        => y
$RngComRefPhase (+-1)                     => -1
$RngComWinConstTime [see note]            => 0.08
$RngComScale                              => 1.0
$MoComp (y/n)                             => y
$MoCompRngShiftFlg (y/n)                  => y
$RngShiftInterpSize                       => 8
$MoCompRngShiftSign [+-1]                 => -1
$MoCompPhaseSign [+-1]                    => -1
$InterferenceSuppress (lms/notch/none)    => notch
$LmsUpdateRate                            => 1
$LmsNumWeights                            => 256
$LmsSidelobeOrder                         => 0
$NotchUpdateRate                          => 1000
$NotchNumFFTLines                         => 1000
$NotchCutoff (dB)                         => 2
$NotchMedianKernLen                       => 65

---StepFreqProcStageParams---
$StepFreqWinConstTime                     => 0.08

---AzProcStageParams---
$RngFocSegments                           => 256
$AzComRefPhase [+-1]                      => -1
$AzComNomAzRes [m]                        => 12.0
$AzComWinConstTime                        => 1.0
$AzComWinConstFreq                        => 0.08
$NumAzLooks                               => 4
$AzLookOverlapFrac [0.0-1.0]              => 0.5
$AzComRngCurvInterpSize                   => 4
$AzComInvFFTSizeReduc (power of 2)        => 2
$AzComRngCurvBatchSize                    => 256
$DetectMethod (cmplx/mag/pow/powdB)       => pow
$AzComScale                               => 1.0

---FloatToTiffStageParams---
$Float2ByteMath (none/powx/log/logpowx)   => none
$Float2ByteScale                          => 1.0e-5

---EndianSwapStageParams---
$OutputEndian (little/big)                => little

---OutputOrientStageParams---
$OutputOrient (azline/rngline)            => azline


NOTES:
=====
Proc Control:          
    y    - Run this module, performing relevant parameter
           calculations or reads from log files, and write
           log file, if any.
    n    - Do not run this module, but infer parameters and
           read from log files as if it had been run.
    only - Run only this module. Dynamic parameters are
           calculated or read from intermediate log files
           of the other modules unless these are set to "off"
    off  - Do not run this module and do not infer any
           parameters or read from its log file (i.e. as if
           module absent altogrether). Note that the 'off'
           state is only permitted for certain modules. Checks
           are performed for modules not in 'off' state which
           depend on a module which is in the 'off state.
    These options allow flexible processor configuration. For
    example, a run may be restarted from half-way, by setting
    the previously completed portions to "n". The correct
    parameters for the later modules will still be calculated
    or read from the intermediate log files.

    The main modules (available as standalone modules, unless
                      contra-indicated):
    
    UnpackIMU -
        Unpack the IMU records from the SASAR LBR file and write
        to an ASCII file (C executable).
    UnpackDGPS -
        Unpack the DGPS records to a more readable ASCII format
        and sync with the IMU data (Python code).
    MergeMocData -
        Merge the IMU and DGPS records to form a single ASCII
        file with the LBR PRI, Latitude, Longitude, etc. This
        does all the smoothing, interpolation, etc. (Python code
        which also requires the Scientific Python modules as
        freely available on the web).
    MocompCalc -
        Calculate the range shifts required for each range line
        from the merged motion data. Also creates the geocoding
        information. (C executable).
    PlotMotionError -
        Plot the motion compensation range shift as calculated
        by the MocompCalc module. (Python code - not standalone)
    SniffDC -
        Calculate the DC offsets and average I to Q value ratio
        from an analysis of part of the raw data (C executable)
    RngProc -
        Range compression, interference suppression and motion
        compensation correction implementation (C executable).
    StepFreqProc -
        Step frequency processing (C executable).
    CornerTurn -
        Corner turn the range compressed file (C executable).
    AzProc -
        Range curvature correction, azimuth compression and
        multilook (C executable).
    Float2Tiff -
        Convert floating point output from azimuth compression
        to TIFF file for easy viewing (various C executables).
    EndianSwap -
        Swap endian format (C executable).
    Orient -
        Covert from azimuth line format to range line format
        (uses corner turn C executable)
    ImageLog -
        Create image log file with geocoding info, etc. (Python
        code - not standalone)
    CleanUp -
        Remove all temporary data and log files (Python code -
        not standalone).

RngComWinConstTime - set to 1.0 if step freq processing.        
    """)      
    outf.close()
    return 0


#---------------------
# general help message

def printHelp():

    print (
'Quick help screen:\n',
'------------------\n\n'

'RUN PROCESSOR : (perform SAR processing - see "create template files" below)\n',
'usage: g2.py [processor_config_file_name]\n\n',

'CREATE TEMPLATE FILES: (create template processor- and radar config files)\n',
'usage: g2.py --templates\n\n',

'USER MANUAL : (install HTML user manual in current dir)\n',
'usage: g2.py --user-manual\n\n',

'QUICK HELP : (this screen)\n',
'usage: g2.py --help\n')

 #-------------------------------------------------------
# write out User Manual as HTML file in current directory

def createHTMLUserManual(fileName):
    outf = open(fileName,'w')
    outf.write("""
<!doctype html public "-//w3c//dtd html 4.0 transitional//en">
<html>
<head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <meta name="GENERATOR" content="Mozilla/4.72 [en] (X11; U; Linux 2.2.12-20 i586) [Netscape]">
   <title>G2 SAR Processor User Manual</title>
</head>
<body lang="en">

<h1>
<font color="#000099">G2 SAR Processor (Ver. 1.0) User Manual (Rev. B)</font></h1>

<h3>
Author: J.M. Horrell - Radar Remote Sensing Group, University of Cape Town</h3>
Note: this HTML file has been automatically generated by the G2 processor.
To update this document, simply copy any HTML source changes into the "createHTMLUserManual()"
function in the file "g2tmpl.py". This will allow users of the processor
to automatically generate the latest version of this document in their
current directory.
<p><a href="#overview">Overview</a>
<br><a href="#quick_start">Quick Start</a>
<br><a href="#design_philos">Design Philosophy</a>
<br><a href="#integrated_config">Integrated Processor Configuration</a>
<br><a href="#standalone_config">Standalone Module Configuration</a>
<br><a href="#queries">Queries</a>
<br>
<hr width="100%">

<a NAME="overview"></a>
<h3>Overview</h3>

<p>The G2 sythetic aperture radar (SAR) processor is based on the range-Doppler
algorithm and has been designed for the processing of airborne SAR data.
The processor is modular and flexible and can handle a wide range of SAR
processing tasks. This version of the processor includes a Python program,
"g2.py", which integrates the various modules and has been designed for
the processing of SASAR VHF data (South African SAR VHF sensor). In many
cases, it is possible to use the integrated version for processing of data
from systems other than the SASAR VHF system by simply changing input parameters.
Where the integrated processor is found not to be suitable for a particular
processing job, it might still be possible to configure the various modules
manually to get the job done (see below).&nbsp; In addition, if semi-automated
processing is required for a new SAR sensor, a new version of the overall
glue program in Python could be created which used the same core modules
(such as range compression, azimuth compression, etc.).

<a NAME="quick_start"></a>
<p><b>For a quick start to the integrated processor</b>, type from the
command line (no quotes): "g2.py --help". This will show the syntax required
for operation, how to generate templates of the required processor- and
radar configuration files and also how to generate this document in the
current directory.
<p>The G2 processor in its integrated form relies heavily on two ASCII
configuration files being properly set up, the processor configuration
file and the radar configuration file. These contain parameter names and
values and must be configured prior to a processing run. The processor
will automatically generate templates of these two configuration files
in the correct format and which include notes on parameter usage by typing
(no quotes) "g2.py --templates".
<ul>
<li>
<b>processor configuration file</b> - This is configured prior to processing
by the processor operator and contains all those parameters which the processor
operator usually would select (like the required azimuth resolution). For
normal operation, the processor config file would be the ONLY file with
which the processor operator needs to be concerned.</li>

<li>
<b>radar configuration file</b> - This is ideally written by the radar
control software and would not normally be altered by the processor operator.
Parameters in the radar config file include radar centre frequency and
pulse length, for example. The radar config file also contains names (without
paths) of the raw data file, LBR file, DGPS file, etc. The paths to those
files are set up in the processor configuration file.</li>
</ul>

<a NAME="design_philos"></a>
<h3>Design Philosophy and Flexibility</h3>

<p>The processor has been written to be as flexible as possible (required
for experimental systems) whilst retaining a simple operator interface.
A very modular approach has been taken with modules usually reading from
some input file and writing to an output file. Log files are also generated
by certain of the modules. These are read by later modules to set parameters
in the integrated processor.
<p>The operator can, through the processor configuration file, decide whether
to keep all the intermediate temporary files and temporary log files or
delete them at the end. For debugging a processing run, it is often very
useful having all the temporary log- and config (and even data) files available,
but these can take up a lot of disk space.
<p>Most of the modules are written in C (for speed) with certain modules
in the Python language. The overall glue for the processor is written in
Python. The C modules have been compiled to be standalone executables and
the Python modules to Python byte code which are called by the main Python
code. The Python code automatically configures the processing for the particular
module.
<p>For example, the azimuth processing module requires a configuration
ASCII text file which the Python code automatically generates before calling
the operating system to excute azimuth processing. As mentioned, this approach
also allows the flexibility of running the C executables independently
of the rest of the processor.

<a NAME="integrated_config"></a>
<h3>Integrated Processor Configuration</h3>

For basic operation, see the <a href="#quick_start">Quick Start</a> section.

In the integrated processor ("g2.py"), modules may be switched in and out
(y/n/only/off states are possible - see the notes in the processor configuration
file for details) which allows for very flexible operation. Using this
approach, it is possible to proceed half-way through a processing run,
stop and examine the intermediate files, and then continue where the processor
left off, for example. Another use might be to omit modules altogether
which are not relevant for the particular radar.
<p>Described in the section below are the standalone modules called by
the g2.py program. For more detail on the integrated processor operation,
a quick scan through the g2.py file&nbsp; will give the user a good idea
of the top-level operation. The g2.py file is well commented and should
be fairly readable even to those unfamiliar with the Python programming
language (note that in Python, logical grouping is done by indentation,
not brackets).

<a NAME="standalone_config"></a>
<h3>Standalone Module Configuration</h3>

As an alternative to the integrated processor approach, the major modules
also be run from the command line independently of the rest of the processor.
This is possible as these modules are in fact standalone executable programs,
As for the integrated processor (g2.py program), some of the major modules
also depend on their own ASCII configuration files being set up. The ability
to run the major modules separately from the command line significantly adds
to the flexibility available to the user, extending the processor's
usefulness well beyond the SASAR VHF system.

<p>The modules which may be run in standalone mode are (in normal order
of SASAR processing execution):
<ul>
<li>
<b>imu_unpack </b>- (C executable - SASAR specific) Unpacks the binary
format LBR motion records to an ASCII file. Type "imu_unpack" for usage.</li>

<li>
<b>g2unpk_dgps.py</b> - (Python script - SASAR specific) Parse DGPS file
from OmniStar system and extract relevant part to ASCII file for processing.
Type "g2unpk_dgps.py" or "python g2unpk_dgps.pyc" for usage.</li>

<li>
<b>g2mocfilt.py</b> (Python script - slightly SASAR specific - rewriting
this in C will speed things up) Merge the data streams from the IMU and
DGPS (including smoothing, interpolation, etc.). This as input&nbsp; the
two output ASCII files from the previous two steps. The output is a single
ASCII file with the merged latitude, longitude and altitude data. Type
"g2mocfilt.py" or "python g2mocfilt.pyc" for usage.</li>

<li>
<b>mocomp</b> - (C executable - slightly SASAR specific). Calculate the
range shift required for each pulse of the processing run. This takes as
input the output merged data file from the previous step and the output
is in the form of a binary file containing pulse numbers (4-byte integers)
and range shifts (4-byte floats). A log file may be generated which contains
the average ground speed calculated from the motion data for the scene
to be processed. Type "mocomp" for usage.</li>

<li>
<b>sniffdc </b>- (C executable) Find the DC offsets and average I/Q channel
ratio for the scene. The values calculated may be written to a log file.
Type "sniffdc" for usage.</li>

<li>
<b>rngcom</b> - (C executable) Motion compensation correction, interference
suppression and range compression. Type "rngcom" for usage. The automatically
generated template rngcom configuration file contains notes on setting
up parameters.</li>

<li>
<b>stepf </b>- (C executable) Perform stepped frequency processing. Combines
the returns from a stepped frequency burst of pulses into a single range
line (with combined bandwidth, upsampled, etc). Type "stepf" for usage.</li>

<li>
<b>corner</b> - (C executable) Corner turn (transpose) data. Reads from
a data file, corner turns and writes out again to a new file. Type "corner"
for usage.</li>

<li>
<b>azcom</b> - (C executable) Range curvature correction, azimuth compression
and multilook processing. Reads the corner turned file and writes out the
azimuth compressed file (image). Type "azcom" for usage.</li>

<li>
<b>iq2mag</b> - (C executable) Convert complex IQ data to magnitude/power.
Reads and writes to disk. Type "iq2mag" for usage.</li>

<li>
<b>flt2byte</b> - (C executable) Convert floating point binary data to
unsigned char binary data (scaled). Type "flt2byte" for usage.</li>

<li>
<b>b2tif</b> - (C executable) Write unsigned char binary data (image) to
TIFF format. Type "b2tif" for usage.</li>

<li>
<b>swapend</b> - (C executable) Swap endian format of binary data on disk.
This useful if moving data from i386 (little-endian) to Sun (big-endian),
for example. Type "swapend" for usage.</li>
</ul>

<a NAME="queries"></a>
<h3>Queries</h3>

<p>For developers, more detail on the G2 processor may be found in the
G2 processor design document. Queries or requests may be directed to:
<p>Jasper Horrell (email: jasper@eng.uct.ac.za) or Richard Lord (email:
rlord@rrsg.ee.uct.ac.za)
<br>UCT Radar Remote Sensing Group (web: http://rrsg.ee.uct.ac.za)
</body>
</html>

""")
