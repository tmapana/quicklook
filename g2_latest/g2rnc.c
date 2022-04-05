/*==========================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2rnc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Unauthorized use prohibited!!
Part of G2 SAR processor. Main range compression procedure. Also can perform
presumming, STC removal, motion comp, Doppler centroid removal, range walk
correction.

VERSION/AUTHOR/DATE : rngcom / Jasper Horrell / 1994
COMMENTS: 
Initial versions.

VERSION/AUTHOR/DATE : 1995-05-06 / Jasper Horrell / 1995-05-06
COMMENTS:
 Motion comp: 06-05-95 03:43pm
   If mcomflg='Y',motion compensation phase correction performed in time
   domain after range compression. If in addition, mcomrflg='Y', motion comp
   range shifts performed in time domain. Note: not possible to perform the
   range shifts without also performing the phase corrections, but is possible
   to perform the phase shifts without the range shifts. 
 Doppler centroid now may be removed. 
 Range compression may be omitted to still use other functions.
   This is not very efficient. 06-05-95 10:46am 
 Now Performs Secondary Range Compression (SRC) 06-05-95 10:52am 
 Now performs range walk range shifts in time domain if RngWalkRngShiftFlg='Y'
   (can be combined with motion comp range shifts) . Also performs corresponding
   range walk phase shifts if RngWalkPhaseShiftFlg='Y' (can be
   combined with motcomp phase shifts). The two range walk operations
   (shift and phase) are independent, unlike the motcomp operations. 

VERSION/AUTHOR/DATE : g2rng1a / Jasper Horrell / 1997-10-06
COMMENTS:
rngcom8.c modified for G2 SASAR ground processsor

VERSION/AUTHOR/DATE : g2rnc1b / Jasper Horrell / 1997-12-05
COMMENTS:
Allow compilation as a function or standalone. Add error checks in parsing, 
check program version in command file 

VERSION/AUTHOR/DATE : 1998-07-20 / Jasper Horrell / 1998-07-20
COMMENTS:
Add in option to read in and write out float data. Change max message for 
float output, allow specify direction of the motion comp range shift. 
Open spec file as text to improve DOS/Win32 compatibility  

VERSION/AUTHOR/DATE : 1999-01-20 / Jasper Horrell / 1999-01-20
COMMENTS:
Pass structure (as defined in g2rnc.h header file) to prog instead of the 
command file name when calling this code as a function. Allow separate I 
and Q input DC offset parameters (make these doubles). 

VERSION/AUTHOR/DATE : 1999-01-27 / Jasper Horrell / 1999-01-22
COMMENTS:
Add a few messages, rename some variables, re-organise spec file and 
RncFileStruct, add RngCom flag. Fixed bug where DC Offset removed in 
RngBinsToProcess (not X2) for zero Doppler centroid case. Changed some 
default parameter values for template command file. Move FUNC to g2func.h 
and rename to FUNC_RNGCOM. Fixed bug giving mostly zero output when 
DopCentoid zero. Shift time domain ref function to be origin centred (also 
remove the 'move' parameter). 

VERSION/AUTHOR/DATE : 1999-01-28 / Richard Lord, Jasper Horrell / 1999-01-28
COMMENTS:
Changes to incorporate Richard's interference suppression routines. Changed
num_weights to LmsNumWeights and altered intereference suppression error
checking (1999-02-02). Added check that mocomp file large enough and more
mocomp message output (1999-02-08). 
Incorporated LmsWeights from g2lms.c into g2rnc.c (1999-02-08).
Included LmsUpdateRate in g2lms.c in order to implement the standard 
LMS adaptive filter without using transfer functions. (1999-02-22).
Added LmsUpdateRate-note at end of file.  (1999-02-22).

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 1999-08-04
COMMENTS:
Calc range FFT size automagically. Divide scale by FFT size.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 1999-09-14
COMMENTS:
Changed interference suppression updates rates to Int4B type.

VERSION/AUTHOR/DATE : 1.2 / Jasper Horrell / 1999-09-16
COMMENTS:
Specify log file in command file (may be "null" in which case info 
written to stdout). Check input file size big enough. Include
'InputIQRatio' parameter.
19991011 - close log file.

VERSION/AUTHOR/DATE : 1.3 / Jasper Horrell / 2000-03-28
COMMENTS:
Allow for step freq processing motion compensation. This requires
that the mocomp phase correction is correct for each centre frequency
transmitted. Also, removed ability to be a C function (getting long
and messy).

VERSION/AUTHOR/DATE : 1.4 / Jasper Horrell / 2000-04-11
COMMENTS:
Correct bug for zero center shift of reference function when reference
function length is odd (was incorrectly swapping I and Q). Remove "Name"
from config file e.g. LogFileName -> LogFile. 

TO DO (maybe): 
--------------
For step freq mode, allow presumming, SRC, range walk
phase shift, and Doppler centroid removal. Also for step freq mode, 
allow interference suppression to average over multiple PRIs. The
original interference suppression assumes the same part of the spectrum
is available with each PRI.

Complete rewrite to make more modular (e.g. presum, interference suppress,
mocomp, etc.) Step freq certainly is adding complications.

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "g2func.h"
#include "g2rnc.h"
#include "g2rline.h"                          /* Richard */
#include "g2lms.h"                            /* Richard */
#include "g2notch.h"                          /* Richard */


/* Misc defns limited to this file */
#define PROG_VERSION "1.3"
#define LIMIT_UNCHAR 255


/* Function prototypes */
Int2B WriteRngComTmplCmdFile(char CmdFileName[]);

/*********************************************/
/* MAIN RNG COMPRESSION PROG (function */

Int2B main (int argc,char *argv[])
{
  FILE *outfile,*cmdfp,*infile,*stcfile,*mcomfile=NULL,*msg=NULL,
       *logfile=NULL,*stepfuserfile=NULL;
  char InputFileName[STRING_SPACE],OutputFileName[STRING_SPACE],
       STCFileName[STRING_SPACE],MoCompFileName[STRING_SPACE],
       ProgVersion[STRING_SPACE],LogFileName[STRING_SPACE]="null",
       StepFUserFileName[STRING_SPACE]="null";
  Int4B bin,         /* Loop counter */
        binb,bine,   /* Rng bin indices for motion comp */
        closebin,    /* Nearest bin from current for rng interp */
        Errors = 0,
	    FooterBytes, /* No. of bytes at end of each range line */
        FreqStep = 0, /* a counter for stepped freq mode freq steps */
        HeaderBytes, /* No. of bytes in headers */
        i,j,         /* Loop counter */
        indx,indx2,indxi,indxq,indxi2,indxq2, /* Array indices */
        inp,         /* Input pulse index for Doppler centroid removal*/
	    inskip,      /* Initial no. of bytes to skip to start of data */
        InPtSize=0,    /* Byte size of input complex pt (2 - chars,8 - float)*/
	    InputDataType, /* 0 - unsigned char (1 byte),
			     3 - float (4 bytes) */ 
	    InputFileRngBins,
	    InputFileSize,
	    InputFileSizeRequired,
        LmsUpdateRate,   /* # of range lines after which transfer function is updated */ /* Richard */
        NotchUpdateRate, /* # of range lines after which transfer function is updated */ /* Richard */
        MaxValueOutputPRI=-99,
        MaxValueOutputRngBin=-99,
        MoCompRngShiftIndex, /* which mocomp update to use for rng
				   shift */
        MoCompRngShiftSign, /* controls dirn of range shift (not phase shift)*/
        MoCompStartPRI,   /* Index of startpri in motion file */
        MoCompRngUpdates,
	    OverFlows=0,      /* Byte overflows after range compression */
        OutPtSize=0,    /* Byte size of output complex pt (2 - chars,8 - float)*/
	    OutputDataType,   /* 0 - unsigned char (1 byte),
			        3 - float (4 bytes) */ 
	    PreSumRatio,      /* Azimuth presum ratio */
	    pri0,        /* First PRI number in motion file*/
	    pulse,       /* Loop counter */
	    PreSummedPulsesToUse, 
        rnglen,      /* Length of range reference function (in rng bins)*/
        MoCompPhaseSign,     /* Sign of phase runout in motion comp */
        RefShiftX2 = 0, /* shift of reference function to origin */
        RngComRefFuncPhaseSign,
	    RngBinsToProcess,  
        RngBinsToProcessX2,     /* Twice RngBinsToProcess */
	    RngFFTSize,     /* Range FFT size (power of 2)*/
        RngFFTSizeX2,    /* Twice no of FFT pts */
        RngShiftInterpSize, 
	    skip,        /* Bytes to skip to required bin of next pulse*/
	    StartRngBin,  
        StartProcessPRI,
        StepFSteps=1, /* Number of step frequency steps */
        truncshift,  /* Range shift in bins rounded down */
        update,      /* Index of motion comp rng update */
        Zeros=0,       /* Zero counter */
        LmsNumWeights,    /* # of weights of LMS adaptive filter */           /* Richard */
        LmsNumWeightsX2,                                                      /* Richard */
        LmsSidelobeOrder, /* level of sidelobe suppression in LMS filter */   /* Richard */
        NotchNumFFTLines, /* # of rng spectra to average for notch filter */  /* Richard */
        NotchMedianKernLen;  /* length of median filter for notch filter */   /* Richard */

  char MoCompFlg,             /* (Y/N) */
       MoCompRngShiftFlg,     /* (Y/N) */
       RngComFlg,             /* Range Compression Flag (Y/N) */       
       RngWalkPhaseShiftFlg,  /* (Y/N) */
       SRCFlg,       /* Secondary Range Compression flag (Y/N) */
       STCFlg,       /* STC flag (Y/N) */
       StepFreqFlg = 'N',
       RngWalkRngShiftFlg,      /* (Y/N) */
       LmsFlg,                  /* (Y/N) */     /* Richard */
       NotchFlg;                /* (Y/N) */     /* Richard */
  double A2DFreq, 
         ang_m,      /* Half beamwidth minus squint angle (rad) */
         ang_p,      /* Half beamwidth plus squint angle (rad) */
         arg,        /* Argument of trig function */
         RngSampleSpacing,    /* Range bin size (m) */
	     CarrierFreq, /* Hz */ 
         DopCentroid, /* Dop centroid (Hz) to be removed (not for SRC calc) */
         dopcenSRC,  /* Dop centroid (Hz) for SRC calc */
         doprateSRC, /* Doppler rate (Hz/sec) for SRC calc */
         fracshift,  /* Fractional bin shift for rng interp */
 	     icor,       /* Motion comp I correction */
         ipart,      /* Stores dummy integer */
         InputDCOffsetI,
         InputDCOffsetQ,
         InputIQRatio=1.0,     	     
         InputPRF,
         lambda,     /* Carrrier freq wavelength (check use with stepped freq) */
         MaxValue=0.0,
         MaxAbsValue=-99.0,         
         NomGroundSpeed,  /* (m/sec) for SRC calc */
         NomRngRes,
         NotchCutoff,  /* cutoff level in dB above neighbouring signal level */   /* Richard */      
         oldchrpc=0.0, /* Chirp const before SRC */
         phcor,      /* Motion comp phase correction */
         qcor,       /* Motion comp Q correction */
         RngWalkAzBeamwidth,  /* (deg) for rng walk correction */
	     RngComChirpBandwidth, /* (Hz)*/
         RngComChirpConst=0.0,  /* Chirp constant */
	     RngComPulseLen,   /* (sec) */
         RngComWinConstTime, /* Window constant */	     
         rngshift,   /* Range shift for pulse (motion + walk) */
         Scale,      /* Reduction scale factor after range compression */
         ScreenUpdateRate,         
         SquintAngle, /* (deg) for SRC and range walk calc */
         SRCFocusRng, /* (m) */
	     t,          /* Time in rng reference */
         tmp,        /* Temporary double */
         tmpi, tmpq, /* Temp I and Q values */
	     walk_rshift=0.0, /* Rng walk lin range shift (m per presummed pulse)*/
	     walk_pshift=0.0, /* Rng walk lin phase shift (rad per pulse)*/
         win;        /* Window factor */
  float  *rcor=NULL;   /* Range corrections for motion comp */

  time_t timeS,timeE; /* Start and end times */

  unsigned char *IQInChar=NULL,  /* Pointer to array for data */
                *IQOutChar=NULL;
  float *IQInFloat=NULL,  /* Pointer to array for data */
        *IQOutFloat=NULL;
  double *comiq=NULL,         /* Pointer to range compressed data array */
         *dummyaddress=NULL,  /* Needed for FFT routine */
         *StepFCentreFreq=NULL,   /* Array of centre freqs for step freq mode */
	     *rngrefiq=NULL,      /* Range reference function I and Q values */
	     *rngrefiq2=NULL,     /* Copy of range reference function for RFI suppression */  /* Richard */
         *rfi=NULL,           /* Pointer to RFI suppression transfer function */       /* Richard */
         *LmsWeights = NULL,  /* Weights of LMS adaptive filter */    
         *sinc=NULL,          /* SINC function values for rng interp */
         *stc=NULL,           /* Pointer to STC array */
         *sumiq=NULL,         /* Array used for presumming */
         *tmpiq=NULL,         /* Temp storage array for rng interpolation */
         *tmp2iq=NULL;        /* Temp array used in ref func origin shift */

  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n------------\n");
  fprintf(msg,"Prog: RNGCOM (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (Copyright UCT 1994-2000)\n");

  if (argc != 2) { 
    fprintf(msg,"Range compression for SAR data.\n\n");
    fprintf(msg,"USAGE: rngcom [cmd file]\n");
    fprintf(msg,"To see command file structure, type `rngcom -tmpl'\n"); 
    fprintf(msg,"(generates a template command file `rngcom.tmp')\n\n");  
    exit(0); 
  }

  /* Write template command file, if requested */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0) {
     if (WriteRngComTmplCmdFile("rngcom.tmp")==0)
       fprintf(msg,"Template command file `rngcom.tmp' written!\n");
     else
       fprintf(msg,"ERROR - in writing template command file!\n");
     exit(0);
  }

  fprintf(msg,"\n");

  /* Open spec file */
  if ( (cmdfp = fopen(argv[1],"r")) == NULL) {
    fprintf(msg,"ERROR - command file %s not found/opened!\n",argv[1]);
    exit(1);
  }  


  /* READ PROCESSING SPECS */

  if (SeekP(cmdfp,"$ProgramVersion","=>",-1,0)) { 
    fprintf(msg,
      "WARNING - in parsing command file (prog version not found)!\n\n"); 
  } 
  else {
    ReadStr(cmdfp,ProgVersion,STRING_SPACE); 
  }

  Errors = 0;

  if (SeekP(cmdfp,"$ScreenUpdateRate","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&ScreenUpdateRate);
  if (SeekP(cmdfp,"$LogFile","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,LogFileName,STRING_SPACE);  
  if (SeekP(cmdfp,"$InputFile","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,InputFileName,STRING_SPACE);  
  if (SeekP(cmdfp,"$OutputFile","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,OutputFileName,STRING_SPACE);
  if (SeekP(cmdfp,"$InputDataType","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &InputDataType);
  if (SeekP(cmdfp,"$OutputDataType","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &OutputDataType);
  if (SeekP(cmdfp,"$StartProcessPRI","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&StartProcessPRI);
  if (SeekP(cmdfp,"$PreSumRatio","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&PreSumRatio);
  if (SeekP(cmdfp,"$PreSummedPulsesToUse","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&PreSummedPulsesToUse);
  if (SeekP(cmdfp,"$InputFileRngBins","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&InputFileRngBins);
  if (SeekP(cmdfp,"$StartRngBin","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&StartRngBin);
  if (SeekP(cmdfp,"$RngBinsToProcess","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&RngBinsToProcess);
  if (SeekP(cmdfp,"$HeaderBytes","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&HeaderBytes);
  if (SeekP(cmdfp,"$FooterBytes","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&FooterBytes);
  if (SeekP(cmdfp,"$InputDCOffsetI","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &InputDCOffsetI);
  if (SeekP(cmdfp,"$InputDCOffsetQ","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &InputDCOffsetQ);
  if (SeekP(cmdfp,"$InputIQRatio","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &InputIQRatio);

  if (SeekP(cmdfp,"$CarrierFreq","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &CarrierFreq);  
  if (SeekP(cmdfp,"$StepFreqUserFile","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,StepFUserFileName,STRING_SPACE);
  if (SeekP(cmdfp,"$A2DFreq","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &A2DFreq);
  if (SeekP(cmdfp,"$RngShiftInterpSize","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &RngShiftInterpSize);
  if (SeekP(cmdfp,"$Scale","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &Scale);

  if (SeekP(cmdfp,"$RngComFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&RngComFlg);
  if (SeekP(cmdfp,"$RngComRefFuncPhaseSign","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &RngComRefFuncPhaseSign);
  if (SeekP(cmdfp,"$RngComChirpBandwidth","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &RngComChirpBandwidth);
  if (SeekP(cmdfp,"$RngComPulseLen","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &RngComPulseLen); 
  if (SeekP(cmdfp,"$RngComWinConstTime","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &RngComWinConstTime);

  if (SeekP(cmdfp,"$MoCompFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&MoCompFlg);
  if (SeekP(cmdfp,"$MoCompFileName","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%s",MoCompFileName);
  if (SeekP(cmdfp,"$MoCompRngShiftFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&MoCompRngShiftFlg);
  if (SeekP(cmdfp,"$MoCompRngShiftSign","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&MoCompRngShiftSign);
  if (SeekP(cmdfp,"$MoCompRngShiftIndex","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &MoCompRngShiftIndex);
  if (SeekP(cmdfp,"$MoCompPhaseSign","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &MoCompPhaseSign);
  if (SeekP(cmdfp,"$MoCompRngUpdates","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &MoCompRngUpdates);
  if (SeekP(cmdfp,"$MoCompStartPRI","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &MoCompStartPRI);

  /* Richard */
  if (SeekP(cmdfp,"$LmsFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&LmsFlg);
  if (SeekP(cmdfp,"$LmsUpdateRate","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&LmsUpdateRate);
  if (SeekP(cmdfp,"$LmsNumWeights","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&LmsNumWeights);
  if (SeekP(cmdfp,"$LmsSidelobeOrder","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&LmsSidelobeOrder);

  /* Richard */
  if (SeekP(cmdfp,"$NotchFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&NotchFlg);
  if (SeekP(cmdfp,"$NotchUpdateRate","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&NotchUpdateRate);
  if (SeekP(cmdfp,"$NotchNumFFTLines","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&NotchNumFFTLines);
  if (SeekP(cmdfp,"$NotchCutoff","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&NotchCutoff);
  if (SeekP(cmdfp,"$NotchMedianKernLen","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&NotchMedianKernLen);

  if (SeekP(cmdfp,"$STCFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&STCFlg);
  if (SeekP(cmdfp,"$STCFileName","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,STCFileName,STRING_SPACE);

  if (SeekP(cmdfp,"$SRCFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&SRCFlg);
  if (SeekP(cmdfp,"$DopCentroid","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&DopCentroid);
  if (SeekP(cmdfp,"$RngWalkRngShiftFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&RngWalkRngShiftFlg);
  if (SeekP(cmdfp,"$RngWalkPhaseShiftFlg","=>",-1,0)) Errors++; 
  else ReadChar(cmdfp,&RngWalkPhaseShiftFlg);
  if (SeekP(cmdfp,"$RngWalkAzBeamwidth","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &RngWalkAzBeamwidth);
  if (SeekP(cmdfp,"$SRCFocusRng","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &SRCFocusRng);
  if (SeekP(cmdfp,"$NomGroundSpeed","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &NomGroundSpeed);
  if (SeekP(cmdfp,"$SquintAngle","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%lf", &SquintAngle);
  if (SeekP(cmdfp,"$InputPRF","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&InputPRF);

  fclose (cmdfp);

  if (Errors!=0) {
    fprintf(msg,"ERROR - %ld errors in parsing command file!\n\n",Errors);
    exit(1);
  }

  /* open log file, if specified */
  if (strcmp(LogFileName,"null")!=0) {
    if ( (logfile = fopen(LogFileName,"wt")) == NULL) {
      fprintf(msg,"ERROR - log file %s not opened!\n",LogFileName);
      exit(1);
    }    
    msg = logfile;  /* assign messages to the log file */
    fprintf(msg,"Prog: RNGCOM (Ver. %s) Log File\n",PROG_VERSION);
    fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 2000\n\n");
  }

  /* Check version IDs match */
  if (strcmp(ProgVersion,PROG_VERSION)!=0) {
    fprintf(msg,
      "WARNING - command file version (%s) not same as program (%s)!\n\n",
      ProgVersion,PROG_VERSION);
  }

  /* Check for lower case flags */
  if (RngComFlg=='y') RngComFlg='Y';
  if (RngComFlg=='n') RngComFlg='N';
  if (STCFlg=='y') STCFlg='Y';
  if (STCFlg=='n') STCFlg='N';
  if (MoCompFlg=='y') MoCompFlg='Y';
  if (MoCompFlg=='n') MoCompFlg='N';
  if (MoCompRngShiftFlg=='y') MoCompRngShiftFlg='Y';
  if (MoCompRngShiftFlg=='n') MoCompRngShiftFlg='N';
  if (SRCFlg=='y') SRCFlg='Y';
  if (SRCFlg=='n') SRCFlg='N';
  if (LmsFlg=='y') LmsFlg='Y';              /* Richard */
  if (LmsFlg=='n') LmsFlg='N';              /* Richard */
  if (NotchFlg=='y') NotchFlg='Y';          /* Richard */
  if (NotchFlg=='n') NotchFlg='N';          /* Richard */
  if (RngWalkRngShiftFlg=='y') RngWalkRngShiftFlg='Y';
  if (RngWalkRngShiftFlg=='n') RngWalkRngShiftFlg='N';
  if (RngWalkPhaseShiftFlg=='y') RngWalkPhaseShiftFlg='Y';
  if (RngWalkPhaseShiftFlg=='n') RngWalkPhaseShiftFlg='N';

  if (MoCompFlg=='N') MoCompRngShiftFlg='N';


  /* STEP FREQ STUFF - If step freq user file name not 'null', open step freq
   * user file and read number of freq steps and centre frequencies. Also close
   * step freq file. */
  if (strcmp(StepFUserFileName,"null")!=0) {
    if ( (stepfuserfile = fopen(StepFUserFileName,"rt")) == NULL) {
      fprintf(msg,"ERROR - step freq user file %s not opened!\n",StepFUserFileName);
      exit(1);
    }
  StepFreqFlg = 'Y';     
  fscanf(stepfuserfile,"%ld\n",&StepFSteps);
  }

  StepFCentreFreq = (double *)malloc(sizeof(double)*StepFSteps);
  if (StepFCentreFreq == NULL) {
    fprintf(msg,"ERROR - in step freq array mem allocation!\n"); exit(1);
  }
  
  if (StepFreqFlg =='Y') {
    for (i=0; i<StepFSteps; i++)
      fscanf(stepfuserfile,"%lf\n",&StepFCentreFreq[i]);
    fclose(stepfuserfile);
    /* check that SRC, presum, Doppler centroid removal, and range walk phase
     * shift not selected */
    if (SRCFlg=='Y' || PreSumRatio!=1 || DopCentroid!=0.0 ||
        RngWalkPhaseShiftFlg=='Y') {
      fprintf(msg,"ERROR - SRC, pre-sum, Doppler centroid removal, and range walk\n");
      fprintf(msg,"        phase shift not supported with step freq mode!\n");
      exit(1); 
    }
    if ( (LmsFlg=='Y' && LmsUpdateRate!=1) || (NotchFlg=='Y') ) {
      fprintf(msg,"ERROR - interference suppression with stepped freq mode is\n");
      fprintf(msg,"        only supported for update rate of unity with LMS filter.\n");
      exit(1); 
    }
  }
  else {  /* if not step freq mode */
    StepFCentreFreq[0] = CarrierFreq;
  }


  /* Check that processing parameters are valid */
  if ( (StartRngBin+RngBinsToProcess) > InputFileRngBins ) {
    fprintf(msg,"ERROR - Attempt to process beyond valid range data!\n");
    fprintf(msg,"%ld bins would be needed",(StartRngBin+RngBinsToProcess) );
    exit(1);
  }

  if ( MoCompRngShiftFlg=='Y' && MoCompFlg=='N' ) {
    fprintf(msg,"ERROR - Motion comp rng shift also requires phase shift!\n");
    exit(1);
  }

  if (MoCompRngShiftFlg=='Y' || RngWalkRngShiftFlg=='Y') {
    if ( (RngShiftInterpSize % 2 != 0 && RngShiftInterpSize != 1)
	     || (RngShiftInterpSize <= 0) ) {
      fprintf(msg,"ERROR - Rng shift interp size must be 1 or even!\n");
      exit(1);
    }
  }

  if ( MoCompFlg=='Y' && MoCompRngShiftFlg=='Y' )
    if (MoCompRngShiftIndex >= MoCompRngUpdates) {
      fprintf(msg,
	 "ERROR - MoCompRngShiftIndex >= MoCompRngUpdates!\n");
      exit(1);
    }

 /* Richard */
  if ( LmsFlg == 'Y' && NotchFlg == 'Y') {
    fprintf(msg,"ERROR - LmsFlg and NotchFlg may not both be set to Y!\n");
    exit(1);
  }

  /* Richard */
  if (LmsFlg == 'Y') {
    if ((LmsNumWeights < 2) || (LmsNumWeights > RngBinsToProcess)) {
      fprintf(msg,
        "ERROR - invalid value for LmsNumWeights %ld!\n",LmsNumWeights);
      exit(1);
    }
    if ((LmsSidelobeOrder < 0) || (LmsSidelobeOrder > 3)) {
      fprintf(msg,
        "ERROR - invalid value for LmsSidelobeOrder %ld!\n",LmsSidelobeOrder);
      exit(1);
    }
    if (LmsUpdateRate < 1) {
      fprintf(msg,
        "ERROR - invalid value for LmsUpdateRate %ld!\n",LmsUpdateRate);
      exit(1);
    }
  }  /* end if LmsFlg == Y */

  /* Richard */
  if (NotchFlg == 'Y') {
    if ((NotchNumFFTLines < 1) || (NotchNumFFTLines > PreSummedPulsesToUse)) {
      fprintf(msg,
        "ERROR - invalid value for NotchNumFFTLines %ld!\n",NotchNumFFTLines);
      exit(1);
    }
    if ((NotchMedianKernLen < 3) || (NotchMedianKernLen > RngBinsToProcess)) {
      fprintf(msg,
        "ERROR - invalid value for NotchMedianKernLen %ld!\n",NotchMedianKernLen);
      exit(1);
    }
    if (NotchUpdateRate < 1) {
      fprintf(msg,
        "ERROR - invalid value for NotchUpdateRate %ld!\n",NotchUpdateRate);
      exit(1);
    }
  }  /* end if NotchFlg == Y */

  if ( SRCFlg == 'Y' && RngComFlg == 'N') {
    fprintf(msg,"ERROR - no SRC without rng compression!\n");
    exit(1);
  }

  if (InputDataType != 0 && InputDataType != 3) {
    fprintf(msg,"ERROR - invalid input data type %ld!\n",InputDataType);
    exit(1);
  }

  if (OutputDataType != 0 && OutputDataType != 3) {
    fprintf(msg,"ERROR - invalid input data type %ld!\n",OutputDataType);
    exit(1);
  }

  /* Calculate length of range reference function (in rng bins) */
  rnglen=0;
  if (RngComFlg == 'Y') {
    rnglen = (Int4B)(RngComPulseLen*A2DFreq);
    if (rnglen > RngBinsToProcess/2)
      fprintf(msg,
	"\n WARNING: less than half output bins will be valid!\n\n");
  }

  /* Calc required FFT size */
  RngFFTSize = nextPowerOfTwo(RngBinsToProcess+rnglen);

  /* Get start time */
  timeS = time(NULL);

  /* Miscellaneous */
  RngBinsToProcessX2 = 2*RngBinsToProcess;
  RngFFTSizeX2 = 2*RngFFTSize;
  Scale = Scale/(double)RngFFTSize;  /* to remove FFT scaling effect */
  LmsNumWeightsX2 = 2*LmsNumWeights;                 /* Richard */
  RngSampleSpacing = C/(2.0*A2DFreq);
  lambda = C/CarrierFreq;
  inp = 0;
  if (RngComPulseLen != 0.0) {
    RngComChirpConst = RngComChirpBandwidth/RngComPulseLen;
  }
  if (RngComChirpBandwidth==0.0) {
    NomRngRes = C*RngComPulseLen/2.0;  /* used only for log */
  }
  else {
    NomRngRes = C/(2.0*RngComChirpBandwidth); /* only for log */
  }
  if (InputDataType == 0)
    { InPtSize = 2*sizeof(unsigned char); }
  else if (InputDataType == 3) 
    { InPtSize = 2*sizeof(float); }
  if (OutputDataType == 0) 
    { OutPtSize = 2*sizeof(unsigned char); }
  else if (OutputDataType == 3) 
    { OutPtSize = 2*sizeof(float); }  


  /* Calc SRC effective chirp */
  if (SRCFlg == 'Y') {
    dopcenSRC = 2.0*NomGroundSpeed*sin(SquintAngle*PI/180.0)/lambda;
    doprateSRC = -2.0*NomGroundSpeed*NomGroundSpeed/(SRCFocusRng*lambda);
    oldchrpc = RngComChirpConst;
    RngComChirpConst = oldchrpc*pow(1-(oldchrpc/doprateSRC)*
                 pow(dopcenSRC/CarrierFreq,2.0),-1.0);
  }

  /* Calc rng walk linear range shift */
  if (RngWalkRngShiftFlg=='Y' || RngWalkPhaseShiftFlg=='Y') {
    ang_m = (RngWalkAzBeamwidth/2.0 - SquintAngle)*PI/180.0;
    ang_p = (RngWalkAzBeamwidth/2.0 + SquintAngle)*PI/180.0;
    walk_rshift = (NomGroundSpeed*PreSumRatio/InputPRF)*
                  ( cos(ang_m)-cos(ang_p) )
                  / (sin(ang_p)*cos(ang_m)+sin(ang_m)*cos(ang_p));
  }

   /* Calc rng walk linear phase shift */
   if (RngWalkPhaseShiftFlg == 'Y')
     walk_pshift = (double)MoCompPhaseSign*4.0*PI*walk_rshift/lambda;


   /* MESSAGES */
   fprintf(msg,"MESSAGES:\n");
   fprintf(msg,"Input file                     : %s\n",InputFileName);
   fprintf(msg,"Output file                    : %s\n",OutputFileName);
   switch(InputDataType) {
     case 0:
       fprintf(msg,
	  "Input data type                : unsigned char (0 - 255)\n");
       break;
     case 3:
       fprintf(msg,
	  "Input data type                : float (+-3.4x10^(+-28))\n");
       break;
   }  /* end switch InputDataType */
   switch(OutputDataType) {
     case 0:
       fprintf(msg,
	  "Output data type               : unsigned char (0 - 255)\n");
       break;
     case 3:
       fprintf(msg,
	  "Output data type               : float (+-3.4x10^(+-28))\n");
       break;
   }  /* end switch OutputDataType */
   fprintf(msg,"Start PRI to process           : %ld\n",
	  StartProcessPRI);
   fprintf(msg,"Rng compress presum ratio      : %ld\n",
	     PreSumRatio);
   fprintf(msg,"Pre-summed pulses to use       : %ld\n",
	  PreSummedPulsesToUse);
   fprintf(msg,
      "Range bins to process          : %ld (%ld - %ld from a possible %ld)\n",
	  RngBinsToProcess,StartRngBin,
                StartRngBin+RngBinsToProcess-1,InputFileRngBins);
   fprintf(msg,"Header / footer bytes (input)  : %ld / %ld\n",
          HeaderBytes,FooterBytes);
   fprintf(msg,"Input DC offset I              : %f\n",InputDCOffsetI);
   fprintf(msg,"Input DC offset Q              : %f\n",InputDCOffsetQ);
   fprintf(msg,"Input I/Q ratio                : %f\n",InputIQRatio);
   fprintf(msg,"Carrier frequency              : %e Hz\n",
          CarrierFreq);
   if (StepFreqFlg == 'Y') {
     fprintf(msg,"Step freq user file            : %s\n", StepFUserFileName); 
     fprintf(msg,"Step freq steps                : %ld\n", StepFSteps);
   }         
   fprintf(msg,"A2D sampling frequency         : %e Hz\n",
          A2DFreq);
   fprintf(msg,"Range sample spacing           : %.3f m\n",
          RngSampleSpacing);
   fprintf(msg,"Range FFT size                 : %ld\n",
          RngFFTSize);
   fprintf(msg,"Multiplicative scale factor    : %f\n",
          Scale);                     

   if (RngComFlg == 'N') {
     fprintf(msg,"Rng compression                : no\n");
     fprintf(msg,"Nominal range resolution       : unknown\n");
   }
   else {
     fprintf(msg,"Range compression              : yes\n");
     fprintf(msg,"Reference func phase sign      : %ld\n",
             RngComRefFuncPhaseSign);
     fprintf(msg,"Rng compression chirp BW       : %.4f MHz\n",
             RngComChirpBandwidth*1.0e-6);
     fprintf(msg,"Pulse length                   : %.4e sec\n",
             RngComPulseLen);
     fprintf(msg,"Nominal range resolution       : %.4f m\n",NomRngRes);
     fprintf(msg,"Window time domain constant    : %.4f\n",
             RngComWinConstTime);
     fprintf(msg,"Range ref function length      : %ld rng bins\n",rnglen);
     fprintf(msg,"Num valid output bins          : %ld of %ld\n",
           RngBinsToProcess-rnglen,RngBinsToProcess);
   }

   if (MoCompFlg == 'N') {
     fprintf(msg,"Motion compensation            : no\n");
   }
   if (MoCompFlg=='Y') {
     fprintf(msg,"Motion comp. file              : %s\n",MoCompFileName);
     fprintf(msg,"Motion comp. start PRI to use  : %ld\n",MoCompStartPRI);
     fprintf(msg,"Motion comp. range updates     : %ld\n",MoCompRngUpdates);

     if (MoCompRngShiftFlg == 'Y') {
       fprintf(msg,"Motion compensation            : phase and rng shift\n");
       fprintf(msg,"Motion comp. range shift sign  : %ld\n",MoCompRngShiftSign);
       fprintf(msg,"Motion com. range shift index  : %ld\n",MoCompRngShiftIndex);
     }
     if (MoCompRngShiftFlg == 'N') {
       fprintf(msg,"Motion compensation            : phase shift only\n");
     }
     fprintf(msg,"Motion comp. phase sign        : %ld\n",MoCompPhaseSign);
   }  /* end if MoCompFlg == Y */  


   if (RngWalkRngShiftFlg == 'Y' && RngWalkPhaseShiftFlg == 'Y')
     fprintf(msg,"Range walk correction          : shift and phase\n");
   if (RngWalkRngShiftFlg == 'Y' && RngWalkPhaseShiftFlg == 'N')
     fprintf(msg,"Range walk correction          : shift only\n");

   if (MoCompRngShiftFlg == 'Y' || RngWalkRngShiftFlg == 'Y') {
     if (RngShiftInterpSize == 1) {
       fprintf(msg,"Range shifts                   : nearest neighbour\n");
     }
     else {
       fprintf(msg,"Range shifts                   : %ld pt interp\n",
	       RngShiftInterpSize);
     }
   } /* end if mocomp rng shift or rng walk rng shift */

   /* Richard */
   if (LmsFlg == 'N') {
     fprintf(msg,"LMS interference suppression   : no\n");
   }  
   else {
     fprintf(msg,"LMS interference suppression   : yes\n");
     fprintf(msg,"LMS update rate                : %ld\n",LmsUpdateRate);
     fprintf(msg,"LMS number of weights          : %ld\n",LmsNumWeights);
     fprintf(msg,"LMS sidelobe order             : %ld\n",LmsSidelobeOrder);
   }

   /* Richard */
   if (NotchFlg == 'N') {
     fprintf(msg,"Notch interference suppression : no\n");
   }  
   else {
     fprintf(msg,"Notch interference suppression : yes\n");
     fprintf(msg,"Notch update rate              : %ld\n",NotchUpdateRate);
     fprintf(msg,"Notch number of FFT lines      : %ld\n",NotchNumFFTLines);
     fprintf(msg,"Notch cutoff                   : %f\n",NotchCutoff);
     fprintf(msg,"Notch median length            : %ld\n",NotchMedianKernLen);
   }

   if (STCFlg == 'Y') {
     fprintf(msg,"STC removal during processing  : yes\n");
     fprintf(msg,"STC file                       : %s\n",STCFileName);
   }
   else {
     fprintf(msg,"STC removal during processing  : no\n");
   }
     
   if (SRCFlg == 'Y')
     fprintf(msg,"SRC old/new chirps             : %.4e/%.4e (Hz/sec)\n",
             oldchrpc,RngComChirpConst);
   if (DopCentroid != 0.0)
     fprintf(msg,"Doppler centroid               : removal %.4f Hz\n",
	     DopCentroid);

   fprintf(msg,"\n");

   
  /*****************************
   * OPEN FILES FOR PROCESSING *
   *****************************/

  /* Open input file and check big enough */
  if ( (infile = fopen (InputFileName, "rb") ) == NULL ) { 
    fprintf(msg,"ERROR - Input file %s not opened!\n",InputFileName); 
    exit(1); 
  }

  fseek(infile,0,SEEK_END);
  InputFileSize = ftell(infile);
  fseek(infile,0,SEEK_SET);  
  InputFileSizeRequired = PreSummedPulsesToUse*PreSumRatio*
        (InputFileRngBins*InPtSize+HeaderBytes+FooterBytes);
  if (InputFileSize<InputFileSizeRequired) {
    fprintf(msg,
      "ERROR - Input file size of %ld bytes is too small (%ld required)!\n",
      InputFileSize,InputFileSizeRequired);
    if (msg != stdout) {
      fprintf(stdout,
      "ERROR - Input file size of %ld bytes is too small (%ld required)!\n",
      InputFileSize,InputFileSizeRequired);
    }
    exit(1);
  }

  /* Open output file */
  if ( (outfile = fopen (OutputFileName, "wb") ) == NULL )
    { fprintf(msg,"ERROR - Output file not opened!\n"); exit(1); }

  /* Open motion comp data file, get to start of relevant data */
  if (MoCompFlg == 'Y') {
     if ( (mcomfile = fopen (MoCompFileName, "rb") ) == NULL )
       { fprintf(msg,"ERROR - Input motion file not found/opened!\n"); 
         exit(1); }

     fread(&pri0,sizeof(Int4B),1,mcomfile);
     fprintf(msg,"1st PRI in motion comp data file: %ld\n",pri0);
     
     /* Check mocomp file size large enough */
     fseek(mcomfile,0,SEEK_END);
     if ( ftell(mcomfile) < (4+4*MoCompRngUpdates)*
         (MoCompStartPRI-pri0+PreSummedPulsesToUse*PreSumRatio) ) {
       fprintf(msg,"ERROR - mocomp file size %ld too small!\n",ftell(mcomfile));
       exit(1);
     }
     else {  /* move file pointer back posn after read of first PRI*/
       fseek(mcomfile,sizeof(Int4B),SEEK_SET);  
     }  

     if (MoCompStartPRI >= pri0) 
       { fseek(mcomfile,(4+4*MoCompRngUpdates)*(MoCompStartPRI-pri0)+4,0); }
     else
       { fprintf(msg,"ERROR - attempt to seek beyond start of motion file!\n");
         fprintf(msg,"PRI0: %ld\n", pri0); exit(1); }

  }  /* end if MoCompFlg */ 



  /*********************************************
   * SKIP INPUT FILE TO START OF RELEVANT DATA *
   *********************************************/

   /* Get to start of required data (inskip in bytes)*/
   inskip = StartProcessPRI*(InPtSize*InputFileRngBins+HeaderBytes+FooterBytes)
            + HeaderBytes + InPtSize*StartRngBin; 

   fseek(infile, inskip, 0);

   /* Calc bytes to skip to required bin of next pulse after extraction */
   skip = InPtSize*(InputFileRngBins-StartRngBin-RngBinsToProcess)
          + FooterBytes + HeaderBytes + InPtSize*StartRngBin;

  /*****************************
   * ALLOCATE SPACE FOR ARRAYS *
   *****************************/

   /* General arrays */
   Errors = 0;

   sumiq    = (double *)malloc(sizeof(double)*RngFFTSizeX2);
   stc      = (double *)malloc(sizeof(double)*RngBinsToProcess);
   comiq    = (double *)malloc(sizeof(double)*RngFFTSizeX2);
   rfi       = (double *)malloc(sizeof(double)*RngFFTSizeX2);       /* Richard */
   rngrefiq = (double *)malloc(sizeof(double)*RngFFTSizeX2);

   if ( sumiq==NULL || stc==NULL || comiq==NULL || rfi==NULL 
       || rngrefiq==NULL ) Errors++;  /* Richard */

   if (InputDataType == 0) { 
     IQInChar  = (unsigned char *)malloc(sizeof(unsigned char)*RngBinsToProcessX2);
     if (IQInChar==NULL) Errors++;
   }
   else if (InputDataType == 3) {
     IQInFloat  = (float *)malloc(sizeof(float)*RngBinsToProcessX2);
     if (IQInFloat==NULL) Errors++;
   }
   if (OutputDataType == 0) { 
     IQOutChar  = (unsigned char *)malloc(sizeof(unsigned char)*RngBinsToProcessX2);
     if (IQOutChar==NULL) Errors++;
   }
   else if (OutputDataType == 3) {
     IQOutFloat  = (float *)malloc(sizeof(float)*RngBinsToProcessX2);
     if (IQOutFloat==NULL) Errors++;
   }   

   if (RngComFlg == 'Y') {
     tmp2iq   = (double *)malloc(sizeof(double)*RngFFTSizeX2);
     rngrefiq2 = (double *)malloc(sizeof(double)*RngFFTSizeX2);       /* Richard */
     if ( rngrefiq2 == NULL || tmp2iq ==NULL) Errors++;
   } 

   if (LmsFlg == 'Y') {
     LmsWeights  = (double *)malloc(sizeof(double)*LmsNumWeightsX2);  /* Richard */
     if ( LmsWeights == NULL ) Errors++; 

     /* Initialise weight array */
     for (i=0; i<LmsNumWeightsX2; i++) LmsWeights[i] = 0.0;             /* Richard */
   }

   if (Errors != 0)
     { fprintf(msg,"ERROR - in array memory allocation!\n"); exit(1); }


   /* Motion comp array */
   if (MoCompFlg=='Y') {
     rcor = (float *)malloc(sizeof(float)*MoCompRngUpdates);
     if (rcor==NULL) {
       fprintf(msg,"ERROR in rcor array memory allocation!\n");
       exit(1);
     }
   }

   /* Range interpolation arrays */
   if (MoCompRngShiftFlg=='Y' || RngWalkRngShiftFlg=='Y')
     if (RngShiftInterpSize != 0 && RngShiftInterpSize != 1) {
       tmpiq = (double *)malloc(sizeof(double)*
                 (RngBinsToProcess+RngShiftInterpSize-1)*2);
       sinc  = (double *)malloc(sizeof(double)*RngShiftInterpSize);
       if (tmpiq==NULL || sinc==NULL) {
         fprintf(msg,"ERROR in interpolation arrays' memory allocation!\n");
         exit(1);
       }
     }

  /**********************
   * READ IN STC VALUES *
   **********************/

  /* Open STC file */
  if (STCFlg=='Y') {
    if ( (stcfile = fopen (STCFileName, "rt") ) == NULL )
      { fprintf(msg,"ERROR - STC file not found/opened!\n"); exit(1); }

    /* Skip STC file until start of data */
    for (bin=0; bin<StartRngBin; bin++)
      fscanf (stcfile, "%lf", &tmp);

    /* Read in STC values */
    for (bin=0; bin<RngBinsToProcess; bin++)
      fscanf (stcfile, "%lf", &stc[bin]);

    /* Close STC file */
    fclose(stcfile);
  }
  else {  /* for no STC removal, set STC values to unity */
    for (bin=0; bin<RngBinsToProcess;bin++) 
      stc[bin] = 1.0; 
  }

  if (RngComFlg == 'Y') {

    /********************************************
    * GENERATE RANGE REFERENCE FUNCTION AND FFT *
    ********************************************/

    /* Initialize array */
    for (i=0; i<RngFFTSizeX2; i++) {
      tmp2iq[i] = 0.0;
    }

    /*  Generate reference function (time-reversed complex conjugate)
        (Zero padded at end)   */
    indxi = 0;
    indxq = 1;
    for (bin=0; bin<rnglen; bin++) {
      t = ( (double)bin*(double)rnglen/((double)rnglen-1.0)
             - (double)rnglen/2.0 )/A2DFreq;
      arg = (double)RngComRefFuncPhaseSign*PI*RngComChirpConst*t*t;
      win = RngComWinConstTime+(1.0-RngComWinConstTime)*pow(sin(PI*(Int4B)bin/
	        (Int4B)(rnglen-1.0)),2.0);
      tmp2iq[indxi++] = cos( arg )*win;
      tmp2iq[indxq++] = -sin( arg )*win;
      indxi++;
      indxq++;
    }  /* end for bin loop */

    /* shift time domain func to origin centered */
    RefShiftX2 = -2*(Int4B)(rnglen/2.0); /* -2 moves 2nd cmplx value to first etc */
    for (i=0; i<RngFFTSizeX2;i++) {
       rngrefiq[i] = tmp2iq[(i-RefShiftX2)%RngFFTSizeX2];       
    }  
 
    /* Fourier transform reference function */
    dummyaddress = rngrefiq - 1;
    CFFT(dummyaddress,RngFFTSize,-1);  /* Richard */

    /* Generate copy of reference function used for RFI suppression */    /* Richard */
    for (i=0; i<RngFFTSizeX2; i++) rngrefiq2[i] = rngrefiq[i];

  } /* end if RngComFlg == Y */

  /***************************************
  * REPEAT FOR EACH PRESUMMED RANGE LINE *
  ***************************************/

  fprintf(msg,"Processing the %ld output range lines....\n",
	 PreSummedPulsesToUse);

  Errors = 0;
  for (pulse=0; pulse<PreSummedPulsesToUse; pulse++) {

    /**********************************************************************
    * Read in range line, remove Doppler centroid and perform presumming *
    *                                                                    *
    * Rescale for presum and unSTC data                                  *
    * (end of array is zero padded to FFT size)                          *
    **********************************************************************/

     G2RncReadLine(
          infile,
          RngFFTSize,
          RngBinsToProcess,
          PreSumRatio,
          InputDataType,
          skip,
          pulse,
          &inp,
          IQInChar,
          IQInFloat,
          DopCentroid,
          InputPRF,
          InputDCOffsetI,
          InputDCOffsetQ,
          InputIQRatio,
          stc,
          sumiq                    /* Range line */
          );


     /******************************************************
      * Implement LMS Adaptive Filter                      *
      * and multiply with Reference Function if necessary  *
      ******************************************************/   /* Richard */

     if (LmsFlg == 'Y') {

       if (modf((double)pulse / (double)LmsUpdateRate,&ipart)==0) {
         lms(
             sumiq,               /* Noisy range line */
             RngFFTSize,
             RngBinsToProcess,
             (double)LmsUpdateRate,
             LmsNumWeights,
             LmsSidelobeOrder,
             LmsWeights,          /* Weights of adaptive filter */
             rfi                  /* Transfer function */
            );

        /* Multiply with reference function if necessary */   

         if (LmsUpdateRate>1) {
           if (RngComFlg=='N') {
             for (i=0; i<RngFFTSizeX2; i++) rngrefiq[i] = rfi[i];
           }
           else {
             indxi = 0;  
             indxq = 1;
             for (bin=0; bin<RngFFTSize; bin++) {
               rngrefiq[indxi] = (rngrefiq2[indxi] * rfi[indxi]) -
                                 (rngrefiq2[indxq] * rfi[indxq]);
               rngrefiq[indxq] = (rngrefiq2[indxi] * rfi[indxq]) + 
                                 (rngrefiq2[indxq] * rfi[indxi]);
               indxi++;
               indxi++;
               indxq++;
               indxq++;
             } /* end for loop */
           } /* end else */
         } /* end if LmsUpdateRate */
       }  /* end if modf */
     }  /* end if LmsFlg */


     /******************************************************
      * Update Notch Filter Transfer Function              *
      * and multiply with Reference Function               *
      ******************************************************/  /* Richard */

     if (NotchFlg == 'Y') {

       if (modf((double)pulse / (double)NotchUpdateRate,&ipart)==0) {
	 
         notch(
               NotchNumFFTLines,
               NotchMedianKernLen,
               NotchCutoff,
               infile,
               RngFFTSize,
               RngBinsToProcess,
               PreSummedPulsesToUse,
               PreSumRatio,
               InputDataType,
               skip,
               pulse,
               &inp,            
               IQInChar,
               IQInFloat,
               DopCentroid,
               InputPRF,
               InputDCOffsetI,
               InputDCOffsetQ,
               InputIQRatio,
               stc,
               sumiq,               
               rfi                 
              );

         /* Multiply with reference function */   
   
         if (RngComFlg=='N') {
           for (i=0; i<RngFFTSizeX2; i++) rngrefiq[i] = rfi[i];
         }
         else {
           indxi = 0;  
           indxq = 1;
           for (bin=0; bin<RngFFTSize; bin++) {
             rngrefiq[indxi] = (rngrefiq2[indxi] * rfi[indxi]) -
                               (rngrefiq2[indxq] * rfi[indxq]);
             rngrefiq[indxq] = (rngrefiq2[indxi] * rfi[indxq]) + 
                               (rngrefiq2[indxq] * rfi[indxi]);
             indxi++;
             indxi++;
             indxq++;
             indxq++;
           } /* end for loop */
         } /* end else */
       } /* end if modf */
     } /* end if NotchFlg */





     /******************************************************
      *                                                    *
      * Perform range compression, if required             *
      *                                                    *
      ******************************************************/ 

     if (RngComFlg=='Y' || (LmsFlg=='Y' && LmsUpdateRate>1) || NotchFlg=='Y')  {  /* Richard */


       /*************************
       * Fourier transform data *
       *************************/
       dummyaddress = sumiq - 1;
       CFFT(dummyaddress,RngFFTSize,-1);  /* Richard */


       /***************************************************************
       * Perform compression in freq domain (Multiply in freq domain) *
       ***************************************************************/
       indxi = 0;
       indxq = 1;
       for (bin=0; bin<RngFFTSize; bin++) {
         comiq[indxi] = sumiq[indxi]*rngrefiq[indxi]
                        - sumiq[indxq]*rngrefiq[indxq];
         comiq[indxq] = sumiq[indxi]*rngrefiq[indxq]
                        + sumiq[indxq]*rngrefiq[indxi];
         indxi++;
         indxi++;
         indxq++;
         indxq++;
       }

       /**************
       * Inverse FFT *
       **************/
       dummyaddress = comiq - 1;
       CFFT(dummyaddress,RngFFTSize,1);   /* Richard */

     } /* end if RngComFlg == Y or LmsFlg==Y or NotchFlg==Y */
     else {  /* Skip range compression */ 
        for (j=0;j<RngFFTSizeX2;j++)
           comiq[j] = sumiq[j];
     }


     /*****************************************
      * Read in motion comp data, if required *
      *****************************************/
     if (MoCompFlg=='Y') {
        fread (rcor,sizeof(float),MoCompRngUpdates,mcomfile);
        /* Move pointer in mot comp file to next posn */
        if (pulse != (PreSummedPulsesToUse-1))
          fseek(mcomfile,(4+4*MoCompRngUpdates)*(PreSumRatio-1)+4,1);
     }

     /*********************************************************
      * PHASE correction for rng walk in time domain or else  *
      * defer to combine with motion comp phase correction    *
      *********************************************************/
      if ( MoCompFlg=='N' && RngWalkPhaseShiftFlg=='Y' ) {
        icor = cos(walk_pshift*(double)pulse);
        qcor = sin(walk_pshift*(double)pulse);
        indxi = 0;
        indxq = indxi+1;
        for (bin=0; bin<RngBinsToProcess; bin++) {
          tmpi = comiq[indxi];
          tmpq = comiq[indxq];
          comiq[indxi] = tmpi*icor - tmpq*qcor;
          comiq[indxq] = tmpi*qcor + tmpq*icor;
	      indxi++;
	      indxi++;
	      indxq++;
	      indxq++;
        }  /* end for bin */
      }  /* end if MoCompFlg ... */

     /****************************************************
      * PHASE correction for motion comp. in time domain *
      * (with rng walk phase correction, if applicable)  *
      ****************************************************/
      if ( MoCompFlg == 'Y' ) {

    	/* Calc start and end rng bin indices for first mcomp segment*/
        binb = 0;
        bine = (Int4B)((double)RngBinsToProcess/(double)MoCompRngUpdates)-1;

        /* Repeat for each mocomp rng update */
        for (update=0;update<MoCompRngUpdates;update++) {

          phcor = (double)MoCompPhaseSign*4.0*PI*(double)rcor[update]*
                  StepFCentreFreq[FreqStep]/C;
          if (RngWalkPhaseShiftFlg=='Y')   /* Add rng walk phase correction*/
            phcor = phcor+walk_pshift*(double)pulse;
          icor = cos(phcor);
          qcor = sin(phcor);

          /* Perform phase correction */
          indxi = 2*binb;
          indxq = 2*binb+1;
          for (bin=binb; bin<(bine+1); bin++) {
            tmpi = comiq[indxi];
            tmpq = comiq[indxq];
            comiq[indxi] = tmpi*icor - tmpq*qcor;
            comiq[indxq] = tmpi*qcor + tmpq*icor;
	        indxi++;
	        indxi++;
	        indxq++;
	        indxq++;
          }  /* end for bin */

          /* Update start and end bins */
          binb = bine+1;
          bine += (Int4B)((double)RngBinsToProcess/(double)MoCompRngUpdates);

        }  /* End for update (mcomp rng update) loop */

        /* Increment frequency step counter and check in range */
        FreqStep++;
        if (FreqStep == StepFSteps) FreqStep = 0;

      }   /* End if MoCompFlg cond */


     /***************************************************************
     * Calc range shift, if required (+ve range shift for move away)*
     ***************************************************************/
     rngshift = 0.0;
     if ( MoCompRngShiftFlg=='Y' )
       rngshift = (double)MoCompRngShiftSign*(double)rcor[MoCompRngShiftIndex];
     if ( RngWalkRngShiftFlg=='Y' )
       rngshift = rngshift + (walk_rshift*pulse-walk_rshift*
			      PreSummedPulsesToUse/2.0);


     /*********************************************************
     * RNG SHIFTS in time domain (interpolation), if required *
     **********************************************************/
     if ( (MoCompRngShiftFlg=='Y'|| RngWalkRngShiftFlg=='Y')
	      && rngshift != 0.0 ) {

       /* NEAREST NEIGHBOUR correction */
       if (RngShiftInterpSize == 1) {

         if (rngshift < 0.0) {   /* For a shift closer */
           closebin = (Int4B)(rngshift/RngSampleSpacing - 0.5); 
           indxi = 0;
           indxq = indxi+1;
           for (bin=0;bin<RngBinsToProcess+closebin;bin++) { /* closebin -ve */
             comiq[indxi] = comiq[indxi-2*closebin];
             comiq[indxq] = comiq[indxq-2*closebin];
             indxi++; indxi++;
             indxq++; indxq++;
           }
           for (bin=0;bin>closebin; bin--) { /* Zero pad at end */
             comiq[indxi] = 0.0;
             comiq[indxq] = 0.0;
             indxi++; indxi++;
             indxq++; indxq++;
           }
         }                 /* end shift closer */

         if (rngshift >= 0.0) {  /* For a shift away */
           closebin = (Int4B)(rngshift/RngSampleSpacing + 0.5);
           indxi = 2*(RngBinsToProcess-1);
           indxq = indxi+1;
           for (bin=0;bin<RngBinsToProcess-closebin;bin++) {
             comiq[indxi] = comiq[indxi-2*closebin];
             comiq[indxq] = comiq[indxq-2*closebin];
             indxi--; indxi--;
             indxq--; indxq--;
           }
           for (bin=0;bin<closebin; bin++) { /* Zero pad at begin */
             comiq[indxi] = 0.0;
             comiq[indxq] = 0.0;
             indxi--; indxi--;
             indxq--; indxq--;
           }
         }  /* end if closebin >= 0 */

       }  /* end nearest neighbour (if RngShiftInterpSize == 1) */

       /* INTERPOLATION shifting (matched for rounding down) */
       else {

         if (rngshift<0.0) {  /* For a shift CLOSER */

           truncshift = (Int4B)(rngshift/RngSampleSpacing);
           fracshift = (double)truncshift-rngshift/RngSampleSpacing; /* note +ve */

           /* Zero array */
           for (i=0;i<(2*(RngBinsToProcess+RngShiftInterpSize-1));i++)
	         tmpiq[i]=0.0;

           /* Read comiq into tmpiq array */
           indx=RngShiftInterpSize - 2;
           indx2 = 0;
           for (j=0;j<2*RngBinsToProcess;j++)
             tmpiq[indx++] = comiq[indx2++];

           /* Calc sinc terms */
           for (i=0;i<RngShiftInterpSize;i++) { 
             arg = PI*(fracshift + (double)(RngShiftInterpSize/2) - 1.0
                   - (double)i);
             if (arg<1.0e-5 && arg>-1.0e-5) { 
               sinc[i]=1.0; 
             }
             else { 
               sinc[i] = (sin(arg))/arg; 
             }
           }

           indxi = 0;
           indxq = indxi+1;
           for (bin=0;bin<RngBinsToProcess;bin++) {
             tmpi = 0.0;
             tmpq = 0.0;
             indxi2 = 2*(bin - truncshift);  /* note truncshift -ve */
             indxq2 = indxi2+1;

             /* check full interp still in bounds of tmpiq array */
             if (indxi2 > 2*RngBinsToProcess - 2) {
               comiq[indxi++] = 0.0;   /* zero pad if out of bounds */
               comiq[indxq++] = 0.0;
               indxi++;
               indxq++;
              }
              else {  /* if within bounds, do interpolation */ 
                for (i=0;i<RngShiftInterpSize;i++) { 
                  tmpi += tmpiq[indxi2++]*sinc[i];
                  tmpq += tmpiq[indxq2++]*sinc[i];
                  indxi2++;
                  indxq2++;
                }

                comiq[indxi++] = tmpi;
                comiq[indxq++] = tmpq;
                indxi++;
                indxq++;
              }
            } /* end for bin loop */

          }                     /* end shift (CLOSER) */


          if (rngshift > 0.0) { /* For a shift AWAY */
            truncshift = (Int4B)(rngshift/RngSampleSpacing);
            fracshift = rngshift/RngSampleSpacing-(double)truncshift; /* note +ve */

            /* Zero array */
            for (i=0;i<(2*(RngBinsToProcess+RngShiftInterpSize-1));i++)
	          tmpiq[i]=0.0;

            /* Read comiq into tmpiq array */
            indx=RngShiftInterpSize;         /* note differs from -ve shift */
            indx2 = 0;
            for (j=0;j<2*RngBinsToProcess;j++)
              tmpiq[indx++] = comiq[indx2++];

            /* Calc sinc terms */
            for (i=0;i<RngShiftInterpSize;i++) { 
               arg = PI*(fracshift + (double)(RngShiftInterpSize/2) - 1.0
                      - (double)i);
               if (arg<1.0e-5 && arg>-1.0e-5) { 
                 sinc[i]=1.0; 
               }
               else { 
                 sinc[i] = (sin(arg))/arg; 
               }
            }

           indxi = 2*(RngBinsToProcess-1);
           indxq = indxi+1;
           for (bin=0;bin<RngBinsToProcess;bin++) {
             tmpi = 0.0;
             tmpq = 0.0;
             indxi2 = 2*(RngBinsToProcess+RngShiftInterpSize-2-bin
			          - truncshift);
             indxq2 = indxi2+1;

             /* check full interp still in bounds of tmpiq array */
             if ( indxi2 < 2*RngShiftInterpSize - 2) {
               comiq[indxi--] = 0.0;   /* zero pad if out of bounds */
               comiq[indxq--] = 0.0;
               indxi--;
               indxq--;
             }
             else {  /* if within bounds, do interpolation */
               for (i=0;i<RngShiftInterpSize;i++) { 
                 tmpi += tmpiq[indxi2--]*sinc[i];
                 tmpq += tmpiq[indxq2--]*sinc[i];
                 indxi2--;
                 indxq2--;
               }

               comiq[indxi--] = tmpi;
               comiq[indxq--] = tmpq;
               indxi--;
               indxq--;
             }  /* end else */

           }  /* end bin loop */

         }  /* end shift (AWAY) */

       }  /* end INTERPOLATION SHIFTING */

     }  /* end RNG SHIFTS */


     /*************************************************************
     * ReSTC range-compressed data, scale and write to output file                                            *
     *************************************************************/
     indxi = 0;
     indxq = 1;
     for (bin=0; bin<RngBinsToProcess; bin++) {
	   tmpi = comiq[indxi]*stc[bin]*Scale;
       if (fabs(tmpi) > MaxAbsValue) {
         MaxValue = tmpi;
         MaxAbsValue = fabs(tmpi);
         MaxValueOutputPRI = pulse;
         MaxValueOutputRngBin = bin;
       }  
	   tmpq = comiq[indxq]*stc[bin]*Scale;
       if (fabs(tmpq) > MaxAbsValue) {
         MaxValue = tmpq;
         MaxAbsValue = fabs(tmpq);
         MaxValueOutputPRI = pulse;
         MaxValueOutputRngBin = bin;
       }  	  

       if (OutputDataType == 0) {
         if (tmpi > 128.0) {  /* Check for I overflows */
	       tmpi = 128.0; 
	       OverFlows++; 
	     }
	     else if (tmpi < -127.0) { 
	       tmpi = -127.0; 
	       OverFlows++; 
	     }
	     if (tmpq > 128.0) {  /* Check for Q overflows */
	       tmpq = 128.0; 
	       OverFlows++; 
	     }
	     else if (tmpq < -127.0) { 
	       tmpq = -127.0; 
	       OverFlows++; 
	     }
	     IQOutChar[indxi] = (unsigned char)(tmpi+127.5);
	     IQOutChar[indxq] = (unsigned char)(tmpq+127.5);
         if (IQOutChar[indxi]==127 && IQOutChar[indxq]==127) 
           Zeros++; /* check for zeros */
       } /* end if OutputDataType == 0 */
       else if (OutputDataType == 3) {
         IQOutFloat[indxi] = (float)tmpi; 
         IQOutFloat[indxq] = (float)tmpq;
       }

	   indxi++; indxi++;
	   indxq++; indxq++;
	 }  /* End bin loop */

     /* Write range line to output file */
     if (OutputDataType == 0) {
       fwrite(IQOutChar,sizeof(unsigned char),RngBinsToProcessX2,outfile);
     }
     else if (OutputDataType == 3) {
       fwrite(IQOutFloat,sizeof(float),RngBinsToProcessX2,outfile);
     }

     /* Report progress to stdout */
     if (modf((double)pulse/ScreenUpdateRate,&ipart)==0)
       fprintf(stdout,"PRI : %ld / rng shift : %.4f m\r",pulse,rngshift);


  } /* End pulse loop */

  /* Report on overflows, zeros and max value - to assist scaling byte output*/
  if (OutputDataType == 0)
    fprintf(msg,"Percentage overflows/zeros : %.4e / %.4e\n", 
           (double)OverFlows*100.0/((double)PreSummedPulsesToUse*
                  (double)RngBinsToProcess),
           (double)Zeros*100.0/((double)PreSummedPulsesToUse*
                  (double)RngBinsToProcess) );     

  if (OutputDataType == 0) {
    fprintf (msg,"Max value after processing : %8.6e (incl. DC offset 127)\n",
             MaxValue+127.5);
  }
  else if (OutputDataType == 3) {
    fprintf (msg,"Max value after processing : %8.6e (float output - no DC offset)\n",
             MaxValue);
  }

  fprintf (msg,"Max processed : output PRI %ld / rng bin %ld\n",
           MaxValueOutputPRI,MaxValueOutputRngBin);

  /* Get end time */
  timeE = time(NULL);
  fprintf(msg,"\nRNG COMPRESSION DONE - in %ld secs (%.2f min)\n",
	   timeE-timeS,(double)(timeE-timeS)/60.0);
  if (msg!=stdout) {
    fprintf(stdout,"\nRng Compression Done!\n");
  }


  /* Close log file */
  if (strcmp(LogFileName,"null")!=0) {
    fclose(logfile);
  }    

  /****************************************
  * Free memory used in range compression *
  *****************************************/
  free (stc);
  free (comiq);
  free (rngrefiq);
  free (StepFCentreFreq);
  free (rfi);    /* Richard */
  if (RngComFlg == 'Y') {
    free (rngrefiq2);
    free (tmp2iq);
  }

  /***********************************
  * FREE ARRAY SPACE AND CLOSE FILES *
  ************************************/
  if (LmsFlg=='Y')        /* Richard */
    free (LmsWeights);

  if (MoCompFlg=='Y') 
    free (rcor);

  free (sumiq);
  if (InputDataType == 0) {
    free (IQInChar); 
  }
  else if (InputDataType == 3) {
    free (IQInFloat); 
  }
  if (OutputDataType == 0) {
    free (IQOutChar); 
  }
  else if (OutputDataType == 3) {
    free (IQOutFloat); 
  }

  fclose (infile);
  fclose (outfile);

  return(0);

}  /* End main program */



/************************************************************************/
/* FUNCTION : WriteRngComTmplCmdFile() */
Int2B WriteRngComTmplCmdFile(char CmdFileName[])
{
  FILE *Out; 
  if ( (Out = fopen (CmdFileName, "w") ) == NULL ) return(1);

  fprintf(Out,"Command file for Rngcom (SAR range compression)\n");
  fprintf(Out,"$ProgramVersion (jmh) => %s\n",PROG_VERSION);
  fprintf(Out,"-----------------------------------------------\n\n");
  fprintf(Out,"--General (required)--\n");
  fprintf(Out,"$ScreenUpdateRate                      => 10\n");
  fprintf(Out,"$LogFile ('null' for none)             => null\n");
  fprintf(Out,"$InputFile                             => tmpg2.raw\n");
  fprintf(Out,"$OutputFile                            => tmpg2.rnc\n");
  fprintf(Out,"$InputDataType [see note]              => 0\n");
  fprintf(Out,"$OutputDataType [see note]             => 3\n");
  fprintf(Out,"$StartProcessPRI                       => 0\n");
  fprintf(Out,"$PreSumRatio                           => 1\n");
  fprintf(Out,"$PreSummedPulsesToUse                  => 1001\n");
  fprintf(Out,"$InputFileRngBins                      => 2048\n");
  fprintf(Out,"$StartRngBin                           => 0\n");
  fprintf(Out,"$RngBinsToProcess                      => 2048\n");
  fprintf(Out,"$HeaderBytes                           => 0\n");
  fprintf(Out,"$FooterBytes                           => 0\n");
  fprintf(Out,"$InputDCOffsetI                        => 127.0\n");
  fprintf(Out,"$InputDCOffsetQ                        => 127.0\n");
  fprintf(Out,"$InputIQRatio                          => 1.0\n\n");

  fprintf(Out,"--Misc (required)--\n");
  fprintf(Out,"$CarrierFreq [Hz - SRC,RW,MoC]         => 141.0e+06\n");
  fprintf(Out,"$StepFreqUserFile [note]               => null\n");
  fprintf(Out,"$A2DFreq [Hz]                          => 12.0e+06\n");
  fprintf(Out,"$RngShiftInterpSize [RW,MoC - note]    => 8\n");
  fprintf(Out,"$Scale                                 => 1.0e-6\n\n");

  fprintf(Out,"--Range compression specific (RC)--\n");
  fprintf(Out,"$RngComFlg [Y/N - SRC]                 => N\n");
  fprintf(Out,"$RngComRefFuncPhaseSign [+-1]          => -1\n");
  fprintf(Out,"$RngComChirpBandwidth [Hz]             => 0.0\n");
  fprintf(Out,"$RngComPulseLen [sec]                  => 83.333e-09\n");
  fprintf(Out,"$RngComWinConstTime                    => 0.08\n\n");

  fprintf(Out,"--Motion compensation specific (MoC)--\n");
  fprintf(Out,"$MoCompFlg [Y/N]                       => N\n");
  fprintf(Out,"$MoCompFileName                        => tmpg2.moc\n");
  fprintf(Out,"$MoCompRngShiftFlg [Y/N]               => N\n");
  fprintf(Out,"$MoCompRngShiftSign [+-1 - note]       => 1\n");
  fprintf(Out,"$MoCompRngShiftIndex [note]            => 0\n");
  fprintf(Out,"$MoCompPhaseSign [+-1]                 => -1\n");
  fprintf(Out,"$MoCompRngUpdates [note]               => 1\n");
  fprintf(Out,"$MoCompStartPRI [note]                 => 0\n\n");

  fprintf(Out,"--LMS interference suppression (LMS)--\n");     /* Richard */
  fprintf(Out,"$LmsFlg [Y/N]                          => N\n");
  fprintf(Out,"$LmsUpdateRate [note]                  => 1\n");
  fprintf(Out,"$LmsNumWeights                         => 512\n");
  fprintf(Out,"$LmsSidelobeOrder [note]               => 0\n\n");

  fprintf(Out,"--Notch interference suppression (Notch)--\n"); /* Richard */
  fprintf(Out,"$NotchFlg [Y/N]                        => N\n");
  fprintf(Out,"$NotchUpdateRate                       => 100\n");
  fprintf(Out,"$NotchNumFFTLines [note]               => 100\n");
  fprintf(Out,"$NotchCutoff [dB - note]               => 3\n");
  fprintf(Out,"$NotchMedianKernLen [note]             => 33\n\n");
 
  fprintf(Out,"--STC specific (STC)--\n");
  fprintf(Out,"$STCFlg [Y/N]                          => N\n");
  fprintf(Out,"$STCFileName [note]                    => tmpg2.stc\n\n");

  fprintf(Out,"--SRC, Doppler centroid and range walk specific (SRC,DOPC,RW)--\n");
  fprintf(Out,"$SRCFlg [Y/N]                          => N\n");
  fprintf(Out,"$DopCentroid [Hz - note]               => 0.0\n");
  fprintf(Out,"$RngWalkRngShiftFlg [Y/N]              => N\n");
  fprintf(Out,"$RngWalkPhaseShiftFlg [Y/N]            => N\n"); 
  fprintf(Out,"$RngWalkAzBeamwidth [deg - RW]         => 60.0\n");
  fprintf(Out,"$SRCFocusRng [m - SRC]                 => 0.0\n");
  fprintf(Out,"$NomGroundSpeed [m/s - SRC,RW]         => 72.0\n");
  fprintf(Out,"$SquintAngle [deg - SRC,RW]            => 0.0\n");
  fprintf(Out,"$InputPRF [Hz - RW,DOPC]               => 136.364\n");

  fprintf(Out,"\nNotes:\n");
  fprintf(Out,"------\n\n");
  fprintf(Out,"InputDataType : 0 - unsigned char IQ (2*1 bytes per point)\n");
  fprintf(Out,"              : 3 - float IQ (2*4 bytes per point)\n\n");
  fprintf(Out,"OutputDataType : 0 - unsigned char IQ\n");
  fprintf(Out,"                     (2*1 bytes per point - DC offset 127)\n");
  fprintf(Out,"               : 3 - float IQ (2*4 bytes per point)\n\n");
  fprintf(Out,
    "StepFreqUserFile : 'null' for no stepped freq mode, else ASCII file with number\n");
  fprintf(Out,"          of freq steps on first line followed by centre frequency of\n");
  fprintf(Out,"          each step, each on a new line, in the transmit order.\n\n");
  fprintf(Out,"RngShiftInterpSize : 0 - none\n");
  fprintf(Out,"                   : 1 - nearest neighbour\n");
  fprintf(Out,"                   : else even\n\n");
  fprintf(Out,
    "MoCompRngShiftSign  : Controls direction of range shift for motion\n");
  fprintf(Out,
    "                      comp. and is independent of phase shift.\n\n");
  fprintf(Out,
    "MoCompRngShiftIndex : The motion comp. range update to use for the range\n");
  fprintf(Out,
    "                      shift. While multiple phase shifts are possible per\n");
  fprintf(Out,
    "                      line, only a single range shift is possible. Must be\n");
  fprintf(Out,
    "                      less than MoCompRngUpdates.\n\n");
  fprintf(Out,
    "MoCompRngUpdates : Ensure precalculated motion comp rng shifts fit in with\n");
  fprintf(Out,
    "                   selected portion to process. If rng bins is not a multiple\n");
  fprintf(Out,
    "                   of num of mocomp rng updates, a few rng bins at far swath\n");
  fprintf(Out,
    "                   may be zeroed due to truncation (not serious).\n\n");
  fprintf(Out,
    "MoCompStartPRI : Must be updated with StartProcessPRI.\n\n");  

  /* Richard */
  fprintf(Out,"LmsUpdateRate    : 1     - implement standard LMS adaptive filter\n");
  fprintf(Out,"                           (i.e. without using transfer functions)\n\n");
  fprintf(Out,"LmsSidelobeOrder : 0     - no sidelobe supression\n");
  fprintf(Out,"                 : 1,2,3 - increasing level of sidelobe supression\n\n");
  fprintf(Out,
    "NotchNumFFTLines : Number of range spectra to average,\n");
  fprintf(Out,
    "                   usually equal to NotchUpdateRate.\n\n");
  fprintf(Out,
    "NotchCutoff   : Number of dB above average signal strength, where notches are\n");
  fprintf(Out,
    "                inserted. [Number of dB = 10 * log((abs(range spectrum))^2) ]\n\n");
  fprintf(Out,
    "NotchMedianKernLen : Length of median filter, used to find average signal\n");
  fprintf(Out,
    "                     strength. If interference peaks are lumped very close\n");
  fprintf(Out,
    "                     together, try increasing this value.\n\n");
  fprintf(Out,
    "STCFileName : ASCII file containing STC curve as space- or line-delimited\n");
  fprintf(Out,
    "              floats. Best to remove STC for rng compress or interference.\n\n");
  fprintf(Out,"DopCentroid : Removed in time domain.\n");

  fclose(Out);
  return(0);

} /* End WriteRngComTmplCmdFile function */

