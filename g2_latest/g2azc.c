/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1994-1999
FILE NAME: g2azc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
Unauthorized use prohibited!!
Part of G2 SAR processor. Range-Doppler azimuth compression stage. 
Also includes range curvature correction, multilook etc. 

VERSION/AUTHOR/DATE : Gp2az7d / Jasper Horrell / pre 1997-11-28
COMMENTS:
Version of azcom7d.c modified for SASAR GP2 ground processor. Different 
header file, cleaned up %f specifiers, Int2B and Int4B types.

VERSION/AUTHOR/DATE : g2azc1a / Jasper Horrell / pre 1997-11-28
COMMENTS:
Change into a function. Add option to append output file.

VERSION/AUTHOR/DATE : g2azc1b / Jasper Horrell / 1997-12-02
COMMENTS:
(1997-11-28) Allow compilation as a function or standalone.
(1997-12-02) Tidy, display az coord of max, rename variable AzresLook to 
  NomAzRes (also changes in cmd file), add parse error checks, more help 
  in template command file.

VERSION/AUTHOR/DATE : G2Azc 1997-12-05 / Jasper Horrell / 1997-12-05
COMMENTS:
Allow power output in dB. Check program version in command file.

VERSION/AUTHOR/DATE : 1997-12-08 / Jasper Horrell / 1997-12-08
COMMENTS:
Rename RngFocUpdates to RngFocSegments, ensure num focus updates does 
not change the num of output rng bins, allow RngfocSegments to be < 1 
(for max updates), fix bugs in dB output option, fix dB bug for zero

VERSION/AUTHOR/DATE : 1998-11-24 / Jasper Horrell / 1998-01-21 - 1998-11-24
COMMENTS:
Fixed range focus segment bug.
Reorder inverse FFT to fix phase hassles.
Change "max at output az/rng" message.
Add in option to read in float data
Open spec file as text to improve DOS/Win32 compatibility.
Add note to command file about delay to first range bin.
Allow specify separate DC offsets for input I and Q.
Minor changes to notes in template cmd file.

VERSION/AUTHOR/DATE : 1999-01-20 / Jasper Horrell / 1999-01-20
COMMENTS:
Pass structure (as defined in g2azc.h header file) to prog instead of 
the command file name when calling this code as a function.

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Move FUNC to g2func.h and rename to FUNC_AZCOM. Update with this 
header info. Changed direction of FFT in CFFT calls so that first 
forward FFT then inverse FFT.
1999-08-03 - cosmetic message change. Also remove repetitive message about
FFT size too small. 

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 1999-08-04
COMMENTS: 
Calcs FFT size automagically. 

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 1999-09-16
COMMENTS:
Specify log file in command file (may be "null" in which case info 
written to stdout). Check input file size big enough.
19991011 - close log file.

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "g2func.h"
#include "g2azc.h"


/* SPECIAL DEFINES */
#define REORDER_IFFT 1


/* Misc defns limited to this file */
#define PROG_VERSION "1.1"
#define DC_OFFSET_UNCHAR 127
#define DC_OFFSET_UNSHORT 32767
#define LIMIT_UNCHAR 255
#define LIMIT_UNSHORT 65535
#define DB_VALUE_FOR_ZERO -135

/* Function prototypes */
Int2B WriteAzComTmplCmdFile(char CmdFileName[]);
Int2B RngCurvStraighten(double A2DFreq,
		        float  **CorrectedBatchFreq,
		        Int4B  FFTSize,
		        Int4B  MasterRefFreqPtsD2,
		        Int4B  *NextRngBinToStraighten,
		        double RngBinSize,
		        Int4B  RngCurvBatchSize,
		        double RngCurvCalc,
		        float  **RngCurvedBatchFreq,
		        Int4B  RngCurvInterpSize,
		        double *RngCurvRefRng,
		        Int4B  *StraightenedRngBins,
		        Int4B  *StraightenedRngBinsUsed
		        );


/*******************************************************************/
/* MAIN AZCOM PROGRAM/FUNCTION */    

#if FUNC_AZCOM
Int2B G2Azc (struct AzcFileStruct Cmd)
#else
Int2B main (int argc,char *argv[])
#endif
{
  FILE *OutputFile,*cmdfp=NULL,*InputFile,*msg=NULL,*logfile=NULL;
  char InputFileName[STRING_SPACE],OutputFileName[STRING_SPACE],
       AppendExistOutFileFlg,ProgVersion[STRING_SPACE]="",
       LogFileName[STRING_SPACE]="null";

  Int2B    SizeFactor=1;     /* used to double size of complex output arrays */
  register Int4B i,k,indx,indx2,indxi,indxq,indxi2,indxq2;
  Int4B    AzPtsToProcess,   /* No. of az pts to process */
           AzPtsToProcessX2, 
           AzPtsPrePostSum,  /* Az points immediately prior to post sum */
           AzRealPtsPrePostSum, /* Scaled by two for complex output */
           bin,              /* Range bin index */
           CurvArrayRngBins=0, /* Rng bins in range curv array */
           dBMinCounter=0,
           DetectMethod,  /* (0-none, 1-mag, 2-power,3-power dB) */
           Errors=0,
           FarRngPad=0,   /* No. of rng curv int far rng bins zero padded */
	       FFTSize,
           FFTSizeX2,        /* Twice the FFT pts */
           FirstLookFreqStartPt,
	       FocusSegment,      /* Focussing segment index */
	       InputDataType,    /* 0 - unsigned char (1 byte),
			        3 - float (4 bytes) */
           InPtSize,         /* 2 - unsigned char, 8 - float */
           InputFileAzPts,
           InputFileRngBin,   
           InputFileRngBins,  
	       InputFileSize,
	       InputFileSizeRequired,
           InvFFTSize,       /* (power of 2) */
           InvFFTSizeReduc,  /* (pow of 2) */
           InvFFTSizeX2,
           look,             /* Look index */          
           LookFreqStartPt, 
           LookFreqPts,     
#if REORDER_IFFT
	       LookFreqPtsX2,
#endif
	       loopstart,        /* Start indx in loop */
           MasterRefFreqPtsD2, /* Half total freq pts spanned by master ref fu. */
           MasterRefTimePts,   /* Master az ref func length (time domain pts)*/
           MaxValueInputRngBin=0, 
           MaxValueAzPixel=0,
           move,             /* Time domain shift at output */
           mvbytes=0,          /* No of bytes to move in curvblock array */
           NextInputRngBin=0,   /* Used in check for end of valid range data */
           NextRngBinToStraighten,
           NumLooks,         /* No. of looks */
           OutputAzPts,      /* No. of output azimuth points */
           OutputAzRealPts,  /* Scaled by 2, for complex output */
	       OutputDataType,   /* 0 - unsigned char (1 byte),
			            1 - unsigned short int (2 bytes),
			            2  - long int (4 bytes),
			            3 - float (4 bytes)
			            4 - double (8 bytes) */
           OutputDCOffset=0,
           OutputLimit=LIMIT_UNCHAR,
			       /* Max possible value for output data type */
           OverFlows=0,           /* Overflows for rng bin */
           PostSumRatio,   
           RefFuncSign,       
           ReportMax,
           RngBinEndFocSeg,    
           RngBinStartFocSeg,
           RngBinsToProcess,  
           RngCurvBatchSize,    /* Rng bins for RC straightening per batch */
           RngCurvInterpSize,   /* Length (pts) of interpolator */
           RngCurvInterpSizeD2, /* Half interp. size (speedup variable) */
           RngCurvSpace,        /* Bytes allocated to rng curved array*/
           RngFocSegments,       /* Number of updates */  
	       ScreenUpdateRate,
           SkipInputBytes,      /* Input file bytes to skip after reading each line */
           StartProcessAzPt,     
           StartProcessRngBin,   
           StraightenedRngBins,     /* No. of bins straightened bins (out of batch) */
           StraightenedRngBinsUsed, /* (out of RngCurvBatchSize)*/
           StraightenedSpace,   /* Bytes allocated to rng curv straightened array*/
           Zeros=0;               /* Zeros for bin */

  double A2DFreq,      
         AzResFull,      /* Az res divided by no. of effective looks */
         CarrierFreq,    
         EffLooks,   
         FarProcRng,      
         FocusRng,
         InputDCOffsetI,
         InputDCOffsetQ,
         InputPRF,
         InputStartSampleDelay, /* Time delay to input data first range bin	*/
         ipart,           /* Stores dummy integer */
         LookDopBw,          /* Bandwidth of each look (Hz)*/
         LookOverlapFrac, 
         magsum,          /* Sum of mags (postsumming) */
         MaxAbsValue=-1.0,
         MaxCurvRngBins,  /* Range curvature at far swath (bins) */
         MaxValue=0.0,    /* for helping to set the scale factor */
         MidCurvRngBins,  /* Range curvature at mid swath (bins)*/
         NearProcRng,
         NomAzRes,        /* assuming window brodening factor of unity */
	     NomGroundSpeed,
         RefFuncPhase,     
         RngBinSize,
         RngCurvCalc=0.0, /* Speed up calc */
         RngCurvRefRng,   /* Closest range in rng curv correction calc */
         Scale,           /* Multiplicative  scale factor for output */
         StartProcessRng,
    	 SumI,SumQ,
         TgtRng,          /* Range to focusing target */
         tmp,tmpi,tmpq,
         TotalDopBw,      /* Total Doppler bandwidth needed for looks (Hz)*/
         Wavelength,     
         WindowFactor,             /* Window scaling factor */
         WinConstTime,            /* Window constant - time domain*/
         WinConstFreq;           /* Window constant - freq domain */
  
  time_t StartTime,EndTime; /* Start- and end times */

  unsigned char *OutputUnChar=NULL,
                *IQCharInput=NULL;   /* Input IQ unsigned char data */  
  unsigned Int2B *OutputUnShort=NULL;
  Int4B *OutputLongInt=NULL;
  float **RngCurvedBatchFreq=NULL,  /* Rng curved data for batch - freq domain */
        **CorrectedBatchFreq=NULL;  /* Rng curv straightened data batch  - freq domain*/
  float *RngCurvedBlock=NULL,       /* Base array of 2-D rng curved data
         	       			  array for batch - floats to save space */
        *OutputFloat=NULL,
        *IQFloatInput=NULL,
        *CorrectedBatchBlock=NULL;  /* Base array of 2-D curv corrected array for batch */
  double **LookRefFunc=NULL;        /* ref func for each look - freq domain */
  double *InputIQ=NULL,          /* Input IQ data */
         *CorrectedBinFreq=NULL, /* Straightened data -  freq domain, single range bin */
         *MasterAzRef=NULL,
         *LookRefBlock=NULL,        /* Base array of look ref functions */ 
         *dummyaddress=NULL,
         *Window=NULL,  /* array of window scaling factors */
         *ShiftedMasterAzRef=NULL, /* zero shifted */
         *FocusedAzLine=NULL,
         *OutputDouble=NULL,
         *PostSummedAzLine=NULL,
         *DetectedAzLine=NULL; /* Array of output detected data */
#if REORDER_IFFT
  double *DechirpedLook=NULL;
#endif


#if FUNC_AZCOM

  /* Assign variables from the structure passed */
  msg = Cmd.msg;
  strcpy(ProgVersion,Cmd.Version);
  ScreenUpdateRate      = Cmd.ScreenUpdateRate;
  strcpy(LogFileName,Cmd.LogFileName);
  InputStartSampleDelay = Cmd.InputStartSampleDelay;
  CarrierFreq           = Cmd.CarrierFreq;
  InputPRF              = Cmd.InputPRF;
  NomGroundSpeed        = Cmd.NomGroundSpeed;
  InputFileAzPts        = Cmd.InputFileAzPts;
  StartProcessAzPt      = Cmd.StartProcessAzPt;
  AzPtsToProcess        = Cmd.AzPtsToProcess;
  InputFileRngBins      = Cmd.InputFileRngBins;
  StartProcessRngBin    = Cmd.StartProcessRngBin;
  RngBinsToProcess      = Cmd.RngBinsToProcess;
  InputDCOffsetI        = Cmd.InputDCOffsetI; /* new */
  InputDCOffsetQ        = Cmd.InputDCOffsetQ; /* new */
  InvFFTSizeReduc       = Cmd.InvFFTSizeReduc;
  strcpy(InputFileName,Cmd.InputFileName);
  strcpy(OutputFileName,Cmd.OutputFileName);
  AppendExistOutFileFlg = Cmd.AppendExistOutFileFlg;  
  RngFocSegments        = Cmd.RngFocSegments;
  RefFuncSign           = Cmd.RefFuncSign;
  A2DFreq               = Cmd.A2DFreq;
  NomAzRes              = Cmd.NomAzRes; 
  WinConstTime          = Cmd.WinConstTime;
  NumLooks              = Cmd.NumLooks;
  LookOverlapFrac       = Cmd.LookOverlapFrac;
  WinConstFreq          = Cmd.WinConstFreq;
  RngCurvInterpSize     = Cmd.RngCurvInterpSize;
  RngCurvBatchSize      = Cmd.RngCurvBatchSize;
  PostSumRatio          = Cmd.PostSumRatio;
  DetectMethod          = Cmd.DetectMethod;
  InputDataType         = Cmd.InputDataType; /* new */
  OutputDataType        = Cmd.OutputDataType; 
  Scale                 = Cmd.Scale;
  ReportMax             = Cmd.ReportMax;

  fprintf (msg,"\nAZIMUTH COMPRESSION STAGE...\n");

#else /* called from command line with cmd file */

  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n-----------\n");
  fprintf(msg,"Prog: AZCOM (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (Copyright UCT 1994-1999)\n");

  if (argc != 2) { 
    fprintf(msg,"Azimuth compression for SAR data.\n\n");
    fprintf(msg,"USAGE: azcom [cmd file]\n");
    fprintf(msg,"To see command file structure, type `azcom -tmpl'\n"); 
    fprintf(msg,"(generates a template command file `azcom.tmp')\n\n");  
    exit(1); 
  }

  /* Write template command file, if requested */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0) {
     if (WriteAzComTmplCmdFile("azcom.tmp")==0)
       fprintf(msg,"Template command file `azcom.tmp' written!\n");
     else
       fprintf(msg,"ERROR - in writing template command file!!\n");
     exit(1);
  }

  fprintf(msg,"\n");

  /* Open spec file */
  if ( (cmdfp = fopen(argv[1],"r")) == NULL) {
    fprintf(msg,"ERROR - command file %s not found/opened!\n",argv[1]);
    exit(1);
  }


  /* PARSE COMMAND FILE */
  
  /* Read version ID of command file. */
  if (SeekP(cmdfp,"$ProgramVersion","=>",-1,0)) { 
    fprintf(msg,
      "WARNING - in parsing command file (program version not found)!\n\n"); 
  } 
  else {
    ReadStr(cmdfp,ProgVersion,STRING_SPACE);
  }

  Errors=0;

  if (SeekP(cmdfp,"$ScreenUpdateRate","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld",&ScreenUpdateRate);
  if (SeekP(cmdfp,"$LogFileName","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,LogFileName,STRING_SPACE); 
  if (SeekP(cmdfp,"$InputStartSampleDelay","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&InputStartSampleDelay);
  if (SeekP(cmdfp,"$CarrierFreq","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&CarrierFreq);
  if (SeekP(cmdfp,"$InputPRF","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&InputPRF);
  if (SeekP(cmdfp,"$NomGroundSpeed","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&NomGroundSpeed);
  if (SeekP(cmdfp,"$InputFileAzPts","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &InputFileAzPts);
  if (SeekP(cmdfp,"$StartProcessAzPt","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &StartProcessAzPt);
  if (SeekP(cmdfp,"$AzPtsToProcess","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &AzPtsToProcess);
  if (SeekP(cmdfp,"$InputFileRngBins","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &InputFileRngBins);
  if (SeekP(cmdfp,"$StartProcessRngBin","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &StartProcessRngBin);
  if (SeekP(cmdfp,"$RngBinsToProcess","=>",-1,0)) Errors++;
  else fscanf(cmdfp,"%ld", &RngBinsToProcess);
  if (SeekP(cmdfp,"$InputDCOffsetI","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf", &InputDCOffsetI);
  if (SeekP(cmdfp,"$InputDCOffsetQ","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf", &InputDCOffsetQ);
  if (SeekP(cmdfp,"$InvFFTSizeReduc","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &InvFFTSizeReduc);
  if (SeekP(cmdfp,"$InputFileName","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,InputFileName,STRING_SPACE); 
  if (SeekP(cmdfp,"$OutputFileName","=>",-1,0)) Errors++; 
  else ReadStr(cmdfp,OutputFileName,STRING_SPACE);
  if (SeekP(cmdfp,"$AppendExistOutFileFlg","=>",-1,0)) Errors++;
  else ReadChar(cmdfp,&AppendExistOutFileFlg);
  if (SeekP(cmdfp,"$RngFocSegments","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &RngFocSegments);
  if (SeekP(cmdfp,"$RefFuncSign","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &RefFuncSign);
  if (SeekP(cmdfp,"$A2DFreq","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&A2DFreq);
  if (SeekP(cmdfp,"$NomAzRes","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&NomAzRes);
  if (SeekP(cmdfp,"$WinConstTime","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&WinConstTime);
  if (SeekP(cmdfp,"$NumLooks","=>",-1,0)) Errors++;
  else fscanf(cmdfp,"%ld",&NumLooks);
  if (SeekP(cmdfp,"$LookOverlapFrac","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&LookOverlapFrac);
  if (SeekP(cmdfp,"$WinConstFreq","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&WinConstFreq);
  if (SeekP(cmdfp,"$RngCurvInterpSize","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &RngCurvInterpSize);
  if (SeekP(cmdfp,"$RngCurvBatchSize","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &RngCurvBatchSize);
  if (SeekP(cmdfp,"$PostSumRatio","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &PostSumRatio);
  if (SeekP(cmdfp,"$DetectMethod","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &DetectMethod);
  if (SeekP(cmdfp,"$InputDataType","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &InputDataType);
  if (SeekP(cmdfp,"$OutputDataType","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &OutputDataType);
  if (SeekP(cmdfp,"$Scale","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%lf",&Scale);
  if (SeekP(cmdfp,"$ReportMax","=>",-1,0)) Errors++; 
  else fscanf(cmdfp,"%ld", &ReportMax);
  
  fclose(cmdfp);

  if (Errors!=0) {
    fprintf(msg,"ERROR - %ld errors in parsing command file!\n\n",Errors);
    exit(1);
  } 

#endif

  /* open log file, if specified */
  if (strcmp(LogFileName,"null")!=0) {
    if ( (logfile = fopen(LogFileName,"wt")) == NULL) {
      fprintf(msg,"ERROR - log file %s not opened!\n",LogFileName);
      exit(1);
    }    
    msg = logfile;  /* assign messages to the log file */
    fprintf(msg,"Prog: AZCOM (Ver. %s) Log File\n",PROG_VERSION);
    fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1999\n\n");
  }


  /* Check version IDs match */
  if (strcmp(ProgVersion,PROG_VERSION)!=0) {
    fprintf(msg,
      "WARNING - command version (%s) not same as program (%s)!\n\n",
       ProgVersion,PROG_VERSION);
  }

  /* Calc misc and speedup variables */
  if (InputDataType == 0)
    { InPtSize = 2*sizeof(unsigned char); }
  else if (InputDataType == 3)
    { InPtSize = 2*sizeof(float); }
  else
    { fprintf(msg,"ERROR - input data type %ld invalid!\n",InputDataType); exit(1); }
  SkipInputBytes = InPtSize*(InputFileAzPts-AzPtsToProcess); /* after each input line read */

  if (RngCurvInterpSize>1)
    {RngCurvInterpSizeD2 = RngCurvInterpSize/2;}
  else
    {RngCurvInterpSizeD2 = 0;}

  Wavelength = C/CarrierFreq;
  EffLooks = (double)NumLooks - ((double)NumLooks-1.0)*LookOverlapFrac;
  AzResFull = NomAzRes/EffLooks;
  AzPtsToProcessX2 = AzPtsToProcess*2;
  RngBinSize = C/(2.0*A2DFreq);
  StartProcessRng = C*InputStartSampleDelay/2.0 + StartProcessRngBin*
                       RngBinSize;

  /* Calculate range curvature at mid- and far processed swath (the max)*/
  NearProcRng = StartProcessRng + RngBinSize*
                (double)RngBinsToProcess/2.0;  /* mid processing swath */
  MidCurvRngBins = (NearProcRng/RngBinSize)*
    ( sqrt(1.0+ 1.0/( 16.0*AzResFull*AzResFull/
		      (Wavelength*Wavelength)-1.0 ))-1.0 );
  FarProcRng = StartProcessRng + RngBinSize*
		(double)RngBinsToProcess;      /* far processed swath */
  MaxCurvRngBins = (FarProcRng/RngBinSize)*
    ( sqrt(1.0+ 1.0/( 16.0*AzResFull*AzResFull/
		      (Wavelength*Wavelength)-1.0 ))-1.0 );

  /* Calc FFT size, based on far swath requirements*/
  MasterRefTimePts = (Int4B)((FarProcRng*InputPRF/NomGroundSpeed)*
		        pow(4.0*AzResFull*AzResFull/
		        (Wavelength*Wavelength)-0.25,-0.5)+0.5);
  FFTSize = nextPowerOfTwo(MasterRefTimePts+AzPtsToProcess);
  FFTSizeX2 = FFTSize*2;
  InvFFTSize = (Int4B)((double)FFTSize/(double)InvFFTSizeReduc);
  InvFFTSizeX2 = 2*InvFFTSize;
  OutputAzPts = (Int4B)((double)AzPtsToProcess/
                ((double)InvFFTSizeReduc*(double)PostSumRatio));

  if (DetectMethod==0)
    SizeFactor = 2;  /* to double size of complex output arrays */ 
  else /* for DetectMethod 1, 2 or 3 */
    SizeFactor = 1;

  OutputAzRealPts = OutputAzPts*SizeFactor;
  AzPtsPrePostSum = OutputAzPts*PostSumRatio; /* az pts prior to postsum */
  AzRealPtsPrePostSum = AzPtsPrePostSum*SizeFactor;

  if (DetectMethod == 2 || DetectMethod == 3)
    {  Scale = Scale/((double)FFTSize*(double)FFTSize); } /* note scale changed here */
  else 
    { Scale = Scale/(double)FFTSize; }
  InputFileRngBin = StartProcessRngBin;  

  /* Check and arrange output (DC offsets applied only to complex
       output case and only for unsigned char and unsigned short output) */

  if (OutputDataType!=0 && OutputDataType!=1 && OutputDataType!=2 &&
      OutputDataType!=3 && OutputDataType!=4) {
    fprintf(msg,"ERROR - output data type (%ld) invalid!!\n",OutputDataType);
    exit(1);
  }

  if (OutputDataType==0) OutputLimit = LIMIT_UNCHAR;
  else if (OutputDataType==1) OutputLimit = LIMIT_UNSHORT;
  
  switch(DetectMethod) {
    case 0:    /* complex output */
      if (NumLooks != 1) {
	    fprintf(msg,"ERROR - complex output only with single look!\n");
	    exit(1);
	  }
      switch(OutputDataType) { 
        case 0: /* unsigned char output */
	      OutputDCOffset = DC_OFFSET_UNCHAR; break;
        case 1: /* unsigned short */
	      OutputDCOffset = DC_OFFSET_UNSHORT; break;
	    case 2: /* Int4B output (process as for double) */
	    case 3: /* float (process as for double) */
	    case 4: /* double */
	      OutputDCOffset = 0; break; 
      }   /* end switch OutputDataType */ 
      break;
    case 1:    /* magnitude detected (process as for power) */
    case 2:    /* power detected */      
    case 3:    /* power detected dB */
      OutputDCOffset = 0; break;
    default:
      fprintf(msg,"ERROR - detection option (%ld) invalid!!\n",DetectMethod);
      exit(1);
    }  /* end switch DetectMethod */  

  /* Check FFTSize and InvFFTSizeReduc is power of 2 */ 
  tmp = (double)((Int4B)( log((double)FFTSize+0.1)/log(2.0) ) );
  if ( pow(2.0,tmp) != (double)FFTSize ) { 
    fprintf(msg,"ERROR - FFTSize %ld not power of two!\n",FFTSize); 
    exit(1); 
  }
  tmp = (double)((Int4B)( log((double)InvFFTSizeReduc+0.1)/log(2.0) ) );
  if ( pow(2.0,tmp) != (double)InvFFTSizeReduc ) { 
    fprintf(msg,"ERROR - InvFFTSizeReduc %ld not power of two!\n",
            InvFFTSizeReduc); 
    exit(1);
  }
  
  
  /* Check look size in freq domain */
  LookDopBw = NomGroundSpeed/NomAzRes;
  LookFreqPts = (Int4B)((LookDopBw*(double)FFTSize/InputPRF)+0.5);

#if REORDER_IFFT
  LookFreqPtsX2 = LookFreqPts*2;
#endif

  if (LookFreqPts > InvFFTSize) {
    fprintf(msg,"\nERROR : InvFFTSize (%ld) < LookFreqPts (%ld)\n",
                   InvFFTSize,LookFreqPts);
    exit(1);
  }

  /* Check total bandwidth used less than InputPRF */
  TotalDopBw = EffLooks*LookDopBw;
  MasterRefFreqPtsD2 = (Int4B)(TotalDopBw*(double)FFTSize/
			       (2.0*(double)InputPRF));
  if( TotalDopBw > InputPRF ) {
    fprintf(msg,
      "ERROR - Total Doppler bandwidth needed (%.3f Hz) is > PRF (%.3f Hz)\n",
             TotalDopBw, InputPRF);
    exit(1);
  }

  /* Check processing block in limits */
  if (StartProcessRngBin+RngBinsToProcess > InputFileRngBins) {
    fprintf(msg,"ERROR - attempt to process beyond valid range data!\n");
    exit(1);
  }

  if (StartProcessAzPt+AzPtsToProcess > InputFileAzPts) {
    fprintf(msg,"ERROR - attempt to process beyond valid az data!\n");
    exit(1);
  }

  /* Check focus updates in range */
  if (RngFocSegments > RngBinsToProcess) { 
    fprintf(msg,
      "ERROR - Rng focus segments %ld out of range (max rng bins %ld)!\n",
            RngFocSegments,RngBinsToProcess); 
    exit(1);
  }
  else if (RngFocSegments < 1) { 
    RngFocSegments  = RngBinsToProcess;
  }

  /* Messages */
   fprintf(msg,"MESSAGES:\n");
   fprintf(msg,"Input file                     : %s\n",InputFileName);
   fprintf(msg,"Output file                    : %s\n",OutputFileName);
   fprintf(msg,"Appending output file          : %c\n",
           AppendExistOutFileFlg);
   fprintf(msg,"Nominal azimuth resolution     : %.3f m\n",NomAzRes);
   fprintf(msg,"Num looks                      : %ld\n",NumLooks);
   fprintf(msg,"Look overlap fraction          : %.3f\n",LookOverlapFrac);
   fprintf(msg,"Doppler bandwidth processed    : %.3f Hz\n",TotalDopBw);
   fprintf(msg,"Master ref func (time domain)  : %ld pts (far swath)\n",
	   MasterRefTimePts);
   fprintf(msg,"FFT size                       : %ld pts\n",FFTSize);
   fprintf(msg,"Look size in freq domain       : %ld pts\n",LookFreqPts);
   fprintf(msg,"Input azimuth samples          : %ld\n",AzPtsToProcess);
   fprintf(msg,"Output azimuth samples         : %ld\n",OutputAzPts);
   fprintf(msg,"Input azimuth sample spacing   : %.4f m\n",
	           NomGroundSpeed/InputPRF);
   fprintf(msg,"Output azimuth sample spacing  : %.4f m\n",
               NomGroundSpeed*InvFFTSizeReduc*PostSumRatio/InputPRF);
   fprintf(msg,
      "Range samples to process       : %ld from %ld (samples %ld - %ld)\n",
	 RngBinsToProcess,InputFileRngBins,StartProcessRngBin,
		StartProcessRngBin+RngBinsToProcess-1);
   fprintf(msg,
      "Processed near slant range     : %.3f m\n",StartProcessRng);
   fprintf(msg,
      "Processed far slant range      : %.3f m\n",FarProcRng); 			
   fprintf(msg,"Range bins per focus segment   : %ld\n",
           (Int4B)(RngBinsToProcess/RngFocSegments));
   fprintf(msg,"Range curv at mid / far swath  : %.2f / %.2f bins\n",
	   MidCurvRngBins,MaxCurvRngBins);
   if (RngCurvInterpSize == 0)
     {fprintf(msg,"Range curvature correction     : none\n");}
   else if (RngCurvInterpSize == 1)
     {fprintf(msg,"Range curvature correction     : nearest neighbour \n");}
   else
     {fprintf(msg,"Range curvature correction     : %ld pt interpolation\n",
	     RngCurvInterpSize);}
   
   switch(DetectMethod)
     {
     case 0:
       fprintf(msg,
        "Detection method               : none (complex output)\n"); 
       break;
     case 1:
       fprintf(msg,"Detection method               : magnitude\n");
       break;
     case 2:
       fprintf(msg,"Detection method               : power\n");
       break;
     case 3:
       fprintf(msg,"Detection method               : power in dB\n");
     }  /* end switch DetectMethod */

   switch(InputDataType)
     {
     case 0:
       fprintf(msg,
          "Input data type                : unsigned char (1 byte I, 1 byte Q)\n");
       break;
     case 3:  
      fprintf(msg,
          "Input data type                : float (4 bytes I, 4 bytes Q)\n");
      break;
     }  /* end switch InputDataType */

   switch(OutputDataType)
     {
     case 0:
       fprintf(msg,
	  "Output data type               : unsigned char (0 - 255)\n");
       break;
     case 1:
       fprintf(msg,
      "Output data type               : unsigned short int (0 - 65535)\n");
       break;
     case 2:
       fprintf(msg,
	  "Output data type               : long int (+- 2 147 483 647)\n");
       break;
     case 3:
       fprintf(msg,
	  "Output data type               : float (+-3.4x10^(+-28))\n");
       break;
     case 4:
       fprintf(msg,
	  "Output data type               : double (+- 1.7*10^(+-308))\n");
     }  /* end switch OutputDataType */

   fprintf (msg,
      "DC offset for output           : %ld\n",(Int4B)OutputDCOffset);
   
   if (Scale==0)
     {
     fprintf(msg,
	"Auto scaling not yet implemented - scale factor set to unity!\n");
     if (DetectMethod==2 || DetectMethod==3) 
       { Scale = 1.0/((double)FFTSize*(double)FFTSize); }
     else 
       { Scale = 1.0/(double)FFTSize; }
     }
   else if (Scale==1)
     { fprintf(msg,"No scaling of output (scale factor = 1).\n"); }
   else
     {
     if (DetectMethod == 2 || DetectMethod == 3) 
       { fprintf(msg,"Multiplicative scale factor    : %6.4e\n",
                 Scale*(double)FFTSize*(double)FFTSize);}
     else 
       { fprintf(msg,"Multiplicative scale factor    : %6.4e\n",
                 Scale*(double)FFTSize);}    
     }

   
  /******************************************************************
  * Open output and input files and move input file pointer to start*
  ******************************************************************/

  if (AppendExistOutFileFlg=='Y' || AppendExistOutFileFlg=='y')
    {
    if ( (OutputFile = fopen (OutputFileName, "ab") ) == NULL )
       { fprintf(msg,"ERROR: Output file not opened!\n");exit(1); }
    }
  else
    {
    if ( (OutputFile = fopen (OutputFileName, "wb") ) == NULL )
      { fprintf(msg,"ERROR: Output file not opened!\n");exit(1); }
    }

  /* Open input file and check big enough */
  if ( (InputFile = fopen (InputFileName, "rb") ) == NULL ) { 
    printf ("ERROR: Input file not found/opened!\n");
    exit(1); 
  }

  fseek(InputFile,0,SEEK_END);
  InputFileSize = ftell(InputFile);
  fseek(InputFile,0,SEEK_SET);  
  InputFileSizeRequired = InputFileAzPts*InputFileRngBins*InPtSize;
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
  

  /* For no rng curv correction or too close to edge*/
  if ( StartProcessRngBin<RngCurvInterpSizeD2 || RngCurvInterpSizeD2==0)
    {
    fseek (InputFile,StartProcessRngBin*InputFileAzPts*InPtSize +
	     StartProcessAzPt*InPtSize,0);
    }
  else  /* where the extra near-rng rng curv interp bins can be used */
    {
    fseek(InputFile,(StartProcessRngBin-RngCurvInterpSizeD2)*
	           InputFileAzPts*InPtSize + StartProcessAzPt*InPtSize,0);
    }

  /********************************
  * ALLOCATE SPACE FOR ARRAYS *
  *********************************/
   Errors = 0;

   InputIQ = (double *)malloc(sizeof(double)*FFTSizeX2);
   CorrectedBinFreq = (double *)malloc(sizeof(double)*FFTSizeX2);
   FocusedAzLine = (double *)malloc(sizeof(double)*InvFFTSizeX2);
   if (PostSumRatio!=1)
     DetectedAzLine = (double*)malloc(sizeof(double)*AzRealPtsPrePostSum);
   PostSummedAzLine = (double*)malloc(sizeof(double)*OutputAzRealPts);
   MasterAzRef=(double *)malloc(sizeof(double)*FFTSizeX2);
   ShiftedMasterAzRef=(double *)malloc(sizeof(double)*FFTSizeX2);
   LookRefBlock=(double *)malloc(sizeof(double)*NumLooks*LookFreqPts*2);
   LookRefFunc=(double **)malloc(sizeof(double *)*NumLooks);

   for (i=0;i<NumLooks;i++)
      LookRefFunc[i]= &LookRefBlock[i*LookFreqPts*2];

   if ( InputIQ==NULL || CorrectedBinFreq==NULL || FocusedAzLine==NULL
        || MasterAzRef==NULL
	|| ShiftedMasterAzRef==NULL || LookRefBlock==NULL || LookRefFunc==NULL
        || PostSummedAzLine==NULL )
     { Errors++; }

   if (InputDataType == 0)
     { 
     IQCharInput=(unsigned char *)malloc(sizeof(unsigned char)*AzPtsToProcessX2);
     if (IQCharInput == NULL) { Errors++; }
     }
   else if (InputDataType == 3)
     {
     IQFloatInput=(float *)malloc(sizeof(float)*AzPtsToProcessX2);
     if (IQFloatInput == NULL) { Errors++; }
     }

   if (Errors != 0)
    {
     fprintf(msg,"ERROR - in memory allocation of processing arrays\n");
     exit(1);
     }

   if (PostSumRatio!=1)
     if(DetectedAzLine==NULL)
       {
       fprintf(msg,"ERROR - in memory allocation of DetectedAzLine array\n");
       exit(1);
       }
   
   Errors = 0;
   switch(OutputDataType)
     {
     case 0:
       OutputUnChar=(unsigned char *)malloc(sizeof(unsigned char)
					    *OutputAzRealPts);
       if (OutputUnChar==NULL) Errors++; break;
     case 1:
       OutputUnShort=(unsigned Int2B *)malloc(sizeof(unsigned Int2B)
						  *OutputAzRealPts);
       if (OutputUnShort==NULL) Errors++; break;
     case 2:
       OutputLongInt=(Int4B *)malloc(sizeof(Int4B)*OutputAzRealPts);
       if(OutputLongInt==NULL) Errors++; break;
     case 3:
       OutputFloat=(float *)malloc(sizeof(float)*OutputAzRealPts);
       if(OutputFloat==NULL) Errors++; break;
     case 4:
       OutputDouble=(double *)malloc(sizeof(double)*OutputAzRealPts);
       if(OutputDouble==NULL) Errors++;
     }  /* end switch OutputDataType */

   if (Errors != 0)
     {fprintf(msg,"ERROR - in memory allocation of output arrays\n");exit(1);}

#if REORDER_IFFT
   DechirpedLook = (double *)malloc(sizeof(double)*LookFreqPts*2);
   if(DechirpedLook==NULL)
     {
     fprintf(msg,"ERROR - in mem allocation of DechirpedLook array\n");
     exit(1);
     }
#endif

   /***********************************************
   * PERFORM IF RANGE CURVATURE CORRECTION *
   *************************************************/

   if (RngCurvInterpSize != 0)
    {

     /* Allocate 2-D array spaces. Note RngCurvedBatchFreq contains extra
	range bins to allow for interpolation at edges and curvature while
	CorrectedBatchFreq contains only RngCurvBatchSize straightened bins */

     if (RngCurvInterpSize == 1)  /* Nearest neighbour */
       { CurvArrayRngBins = (Int4B)MaxCurvRngBins+1+RngCurvBatchSize;
	 mvbytes = sizeof(float)*FFTSizeX2*((Int4B)MaxCurvRngBins+1); }
     else                /* For interpolation */
       {
       CurvArrayRngBins = RngCurvInterpSize+(Int4B)MaxCurvRngBins+
			    1+RngCurvBatchSize;
       mvbytes = sizeof(float)*FFTSizeX2*(RngCurvInterpSize+
			    (Int4B)MaxCurvRngBins+1);
       }

     RngCurvSpace = sizeof(float)*CurvArrayRngBins*FFTSizeX2;
     StraightenedSpace = sizeof(float)*RngCurvBatchSize*FFTSizeX2;

     RngCurvedBlock=(float *)malloc(RngCurvSpace);
     RngCurvedBatchFreq=(float **)malloc(sizeof(float *)*CurvArrayRngBins);
     CorrectedBatchBlock=(float *)malloc(StraightenedSpace);
     CorrectedBatchFreq=(float **)malloc(sizeof(float *)*RngCurvBatchSize);

     for (i=0;i<CurvArrayRngBins;i++)
       RngCurvedBatchFreq[i]= &RngCurvedBlock[i*FFTSizeX2];
     for (i=0;i<RngCurvBatchSize;i++)
       CorrectedBatchFreq[i]= &CorrectedBatchBlock[i*FFTSizeX2];

     if (RngCurvedBlock==NULL || RngCurvedBatchFreq==NULL ||
	 CorrectedBatchBlock==NULL || CorrectedBatchFreq==NULL )
      {
      fprintf (msg,"ERROR - in rng curv array memory allocation\n");
      exit(1);
      }

     fprintf(msg,
	"Range curv array allocated     : %ld bytes (%ld rng bin(s))\n",
	    RngCurvSpace+StraightenedSpace,RngCurvBatchSize);


     /* FILL RngCurvedBatchFreq (First RngCurvInterpSize/2 bins kept
            for interpolation i.e. RngCurvedBatchFreq[RngCurvInterpSizeD2][i]
	   corresponds to the StartProcessRngBin ) */

     loopstart = 0;
     /* Zero pad near edge, if necessary for interpolation */
     if (StartProcessRngBin<RngCurvInterpSizeD2)
       {
       loopstart = RngCurvInterpSizeD2;
       for (k=0;k<RngCurvInterpSizeD2;k++)
	 for(i=0;i<FFTSizeX2;i++) RngCurvedBatchFreq[k][i]=0.0;
       }

     FarRngPad = 0;
     NextInputRngBin = StartProcessRngBin;

     /* Read in data, FFT, zero centre and store in rng curved array */
     for (k=loopstart;k<CurvArrayRngBins;k++)
        {

         /* Initialize array for data */
         for (i=0; i<FFTSizeX2; i++) InputIQ[i]=0.0;

         /* Read in az line (if exists) and remove DC offset.
            (If does not exist, zero pad in rng occurs) */
         if (NextInputRngBin < InputFileRngBins) { 
	   if (InputDataType == 0) { 
             fread (IQCharInput,sizeof(unsigned char),AzPtsToProcessX2,InputFile);
             indx = 0;
             for (i=0;i<AzPtsToProcess;i++) {
               InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetI;
               indx++;
               InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetQ;
               indx++;
             } /* end for loop */
           }
           else if (InputDataType == 3) {
             fread (IQFloatInput,sizeof(float),AzPtsToProcessX2,InputFile);
             indx = 0;
             for (i=0;i<AzPtsToProcess;i++) {
               InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetI;
               indx++;
               InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetQ;
               indx++;      
             } /* end for loop */
           }
         }
         else { 
           FarRngPad++;
         }

         /* Increment file pointer, if valid */
         if (NextInputRngBin < InputFileRngBins-1)
	   fseek(InputFile,SkipInputBytes,1);

         /* Increment rng check */
         NextInputRngBin++;

         /* FFT data */
         dummyaddress = InputIQ - 1;
         CFFT(dummyaddress,FFTSize,-1);

         /* Shift to zero centered and store in RngCurvedBatchFreq */
         for (i=0; i<FFTSize;i++)
            RngCurvedBatchFreq[k][i] = (float)InputIQ[i+FFTSize];
         for (i=FFTSize; i<FFTSizeX2;i++)
	    RngCurvedBatchFreq[k][i] = (float)InputIQ[i-FFTSize];

        } /* End k loop */

     /* Misc for initial rng curv straightening */
     RngCurvCalc = Wavelength*InputPRF/(2.0*NomGroundSpeed*FFTSize);
     RngCurvCalc *= RngCurvCalc;        /* Square RngCurvCalc */
     NextRngBinToStraighten = StartProcessRngBin;
     RngCurvRefRng = StartProcessRng;

     /* STRAIGHTEN RngCurvBatchSize BINS WITH INTERPOLATION
            (straighten only the bandwidth used) */
     if (RngCurvStraighten(A2DFreq,
		      CorrectedBatchFreq,
		      FFTSize,
		      MasterRefFreqPtsD2,
		      &NextRngBinToStraighten,
		      RngBinSize,
		      RngCurvBatchSize,
		      RngCurvCalc,
		      RngCurvedBatchFreq,
		      RngCurvInterpSize,
		      &RngCurvRefRng,
		      &StraightenedRngBins,
		      &StraightenedRngBinsUsed
		      )!=0) {
       fprintf(msg,"ERROR - in RngCurvStraighten routine!!\n");
       exit(1);
     }

   }   /* End if rng curv correction */

   /* Get start time */
   StartTime = time(NULL);
 
   /**********************************
   * REPEAT FOR EACH FOCUSING RANGE *
   **********************************/

   /* Initialize start and end bins for focusing segment. Note that these
      are indexed from the start of processing, not from the first input
      file range bin. */
   RngBinStartFocSeg = 0;
   RngBinEndFocSeg   = (Int4B)(RngBinsToProcess/RngFocSegments) - 1;

   for (FocusSegment=0;FocusSegment<RngFocSegments;FocusSegment++)
     {

     /***********************************************
      * GENERATE AZIMUTH REFERENCE FUNCTION *
      ***********************************************/

      /* Initialize master reference function */
      for (k=0;k<FFTSizeX2;k++) MasterAzRef[k] = 0.0;

      /* Calc focussing range and ref func length for segment */
      FocusRng = StartProcessRng + RngBinSize*
		   (double)(RngBinStartFocSeg+RngBinEndFocSeg)/2.0;
      MasterRefTimePts = (Int4B)((FocusRng*InputPRF/NomGroundSpeed)*
			    pow(4.0*AzResFull*AzResFull/
		            (Wavelength*Wavelength)-0.25,-0.5)+0.5);      
      move = 2*(Int4B)(MasterRefTimePts/(2.0*(double)InvFFTSizeReduc));
      
      indxi = 0;
      indxq = 1;
      for (k=0; k<MasterRefTimePts; k++)
         {
          WindowFactor = WinConstTime+(1.0-WinConstTime)*pow(sin(PI*(double)k/
					((double)MasterRefTimePts-1.0)),2.0);
          TgtRng=sqrt(FocusRng*FocusRng+pow(((double)k-(double)
		  (MasterRefTimePts-1)/2.0)*NomGroundSpeed/InputPRF,2.0));
          RefFuncPhase=(double)RefFuncSign*4.0*PI*CarrierFreq*TgtRng/C;
          MasterAzRef[indxi++] = cos(RefFuncPhase)*WindowFactor;
          MasterAzRef[indxq++] = -sin(RefFuncPhase)*WindowFactor;/*Cmplx conj */
          indxi++;
          indxq++;
         }

       /* Fourier transform master reference function */
      dummyaddress = MasterAzRef - 1;
      CFFT(dummyaddress,FFTSize,-1);

      /* Shift MasterAzRef to zero centered */
      for (i=0; i<FFTSize;i++)
         ShiftedMasterAzRef[i] = MasterAzRef[i+FFTSize];
      for (i=FFTSize; i<FFTSizeX2;i++)
         ShiftedMasterAzRef[i] = MasterAzRef[i-FFTSize];

      /* Allocate array space for freq window function and calc window */
      if ( (Window = (double *)malloc(sizeof(double)*LookFreqPts)) == NULL )
       {
       fprintf (msg,"ERROR - in Window array memory allocation\n");
       exit(1);
       }

      for (i=0;i<LookFreqPts;i++)
        Window[i] = WinConstFreq+(1.0-WinConstFreq)*pow(sin(PI*(double)i/
		 ((double)LookFreqPts-1.0)),2.0);


      /* CALCULATE REFERENCE FUNCTION FOR EACH LOOK */

      /* Calc the freq posn of start of first look */
      FirstLookFreqStartPt =(Int4B)((FFTSize - EffLooks*LookFreqPts)/2);

      for (look=0; look<NumLooks; look++)
         {
          LookFreqStartPt = (Int4B)(FirstLookFreqStartPt +
			     look*(1.0-LookOverlapFrac)*LookFreqPts);
          indxi=2*LookFreqStartPt;
          indxq=indxi+1;
	  indx = 0;
          for(i=0;i<LookFreqPts;i++)
	     {
              LookRefFunc[look][indx++] = ShiftedMasterAzRef[indxi++]
		                           *Window[i];
              LookRefFunc[look][indx++] = ShiftedMasterAzRef[indxq++]
		                           *Window[i];
              indxi++;
              indxq++;
             }

         }  /* End look loop */


   /*********************************************
     * REPEAT PROCESSING FOR EACH RANGE BIN *
     **********************************************/

     for (bin = RngBinStartFocSeg; bin < (RngBinEndFocSeg+1); bin++)
     {   
   
      /* Initialize array for output*/
      if (PostSumRatio != 1)
	for (i=0;i<AzRealPtsPrePostSum;i++) DetectedAzLine[i] = 0.0;
      else /* if no post summing */
	for (i=0;i<OutputAzRealPts;i++) PostSummedAzLine[i] = 0.0; 
      
      /* PERFORM IF NO RNG CURV CORRECTION */
      if (RngCurvInterpSize == 0)
        {
         /* Initialize array for data */
         for (i=0; i<FFTSizeX2; i++) InputIQ[i]=0.0;

         /* Read in azimuth line, remove DC offset,  and seek to next posn*/
         if (InputDataType == 0) { 
           fread (IQCharInput,sizeof(unsigned char),AzPtsToProcessX2,InputFile);
           indx = 0; 
           for (i=0;i<AzPtsToProcess;i++) {
             InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetI;
             indx++;
             InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetQ;
             indx++;
           } /* end for loop */
         }
         else if (InputDataType == 3) {
           fread (IQFloatInput,sizeof(float),AzPtsToProcessX2,InputFile);
           indx = 0;
           for (i=0;i<AzPtsToProcess;i++) {
             InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetI;
             indx++;
             InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetQ;
             indx++;
           } /* end for loop */
         }

         if ( FocusSegment!=RngFocSegments-1 || bin!=RngBinEndFocSeg )
	   fseek(InputFile,SkipInputBytes,1);

	 /* Fourier transform azimuth line */
         dummyaddress = InputIQ - 1;
         CFFT(dummyaddress,FFTSize,-1);

         /* Shift InputIQ to zero centered */
         for (i=0; i<FFTSize;i++)
            CorrectedBinFreq[i] = InputIQ[i+FFTSize];
	 for (i=FFTSize; i<FFTSizeX2;i++)
            CorrectedBinFreq[i] = InputIQ[i-FFTSize];
        }  /* end if no rng cuv correction */

      else  /* PERFORM FOR RNG CURV CORRECTION */
        {
         /* Get rng curv straightened az line */
         if ( StraightenedRngBinsUsed < RngCurvBatchSize )
           {
	   for (i=0; i<FFTSizeX2;i++)
             CorrectedBinFreq[i] =
		 (double)CorrectedBatchFreq[StraightenedRngBinsUsed][i];
	   StraightenedRngBinsUsed++;
           }  /* end if StraightenedBinsUsed */
         else /* If no more straightened bins are available */
           {

            /* Move existing data in RngCurvedBatchFreq */
            memmove(RngCurvedBlock,&RngCurvedBlock[FFTSizeX2*
		       RngCurvBatchSize],mvbytes);
 
	    /* Read in RngCurvBatchSize bins from input file */
            for (k=0;k<RngCurvBatchSize;k++)
               {

	       /* Initialize array for data */
               for (i=0; i<FFTSizeX2; i++) InputIQ[i]=0.0;

               /* Read in az line (if exists) and remove DC offset.
                              (If does not exist, zero pad in rng occurs) */
	       if (NextInputRngBin < InputFileRngBins) {
		 if (InputDataType == 0) { 
                   fread (IQCharInput,sizeof(unsigned char),AzPtsToProcessX2,InputFile);
                   indx = 0 ;
                   for (i=0;i<AzPtsToProcess;i++) {
                     InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetI;
                     indx++;
                     InputIQ[indx] = (double)IQCharInput[indx]-InputDCOffsetQ;
                     indx++;
                   } /* end for loop */
                 }
                 else if (InputDataType == 3) {
                   fread (IQFloatInput,sizeof(float),AzPtsToProcessX2,InputFile);
                   indx = 0;
                   for (i=0;i<AzPtsToProcess;i++) {
                     InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetI;
                     indx++;
                     InputIQ[indx] = (double)IQFloatInput[indx]-InputDCOffsetQ;      
                     indx++;
                   } /* end for loop */
                 }	
	       }
               else { 
                 FarRngPad++;
               }

               /* Increment file pointer if valid */
               if (NextInputRngBin < InputFileRngBins-1)
		 fseek(InputFile,SkipInputBytes,1);

               /* Increment rng check */
	       NextInputRngBin++;

               /* FFT data */
               dummyaddress = InputIQ - 1;
               CFFT(dummyaddress,FFTSize,-1);

               /* Shift to zero centered and store in RngCurvedBatchFreq */
               indx = CurvArrayRngBins - RngCurvBatchSize + k;
               for (i=0; i<FFTSize;i++)
                  RngCurvedBatchFreq[indx][i] = (float)InputIQ[i+FFTSize];
               for (i=FFTSize; i<FFTSizeX2;i++)
                  RngCurvedBatchFreq[indx][i] = (float)InputIQ[i-FFTSize];

               }  /* End k loop */

            /* STRAIGHTEN RngCurvBatchSize BINS WITH INTERPOLATION
                         (NextRngBinToStraighten and RngCurvRefRng correct from
	           last interpolation)*/
            if (RngCurvStraighten(A2DFreq,
		          CorrectedBatchFreq,
		          FFTSize,
		          MasterRefFreqPtsD2,
		          &NextRngBinToStraighten,
		          RngBinSize,
		          RngCurvBatchSize,
		          RngCurvCalc,
		          RngCurvedBatchFreq,
		          RngCurvInterpSize,
		          &RngCurvRefRng,
		          &StraightenedRngBins,
		          &StraightenedRngBinsUsed)!=0) {
              fprintf(msg,"ERROR - in RngCurvStraighten routine!!\n");
              exit(1);
		    }

            /* Get rng curv straightened az line */
	    for (i=0; i<FFTSizeX2;i++)
              CorrectedBinFreq[i] =
		 (double)CorrectedBatchFreq[StraightenedRngBinsUsed][i];
            StraightenedRngBinsUsed++;
           }  /* end else no more straightened bins available */

        }  /* end else perform rng curv correction */


      /* REPEAT PROCESSING FOR EACH LOOK */
      for (look=0; look<NumLooks; look++)
       {

        /* Multiply in frequency domain */
        LookFreqStartPt = (Int4B)(FirstLookFreqStartPt +
			  look*(1.0-LookOverlapFrac)*LookFreqPts);
	indxi = 2*LookFreqStartPt;
        indxq = indxi+1;
	indxi2 = 0;
	indxq2 = 1;
        for (i=0; i<LookFreqPts; i++)
	   {
#if REORDER_IFFT
	    DechirpedLook[indxi2] =  CorrectedBinFreq[indxi]*
			       LookRefFunc[look][indxi2]
			     -CorrectedBinFreq[indxq]*
			       LookRefFunc[look][indxq2];
	    DechirpedLook[indxq2] =  CorrectedBinFreq[indxi]*
			       LookRefFunc[look][indxq2]
			     +CorrectedBinFreq[indxq]*
			       LookRefFunc[look][indxi2];
#else
	    FocusedAzLine[indxi2] =  CorrectedBinFreq[indxi]*
			       LookRefFunc[look][indxi2]
			     -CorrectedBinFreq[indxq]*
			       LookRefFunc[look][indxq2];
	    FocusedAzLine[indxq2] =  CorrectedBinFreq[indxi]*
			       LookRefFunc[look][indxq2]
			     +CorrectedBinFreq[indxq]*
			       LookRefFunc[look][indxi2];
#endif
	    indxi++; indxi++;
            indxq++; indxq++;
	    indxi2++; indxi2++;
	    indxq2++; indxq2++;
           }

#if REORDER_IFFT
	/* zero pad array */
	for (i=0;i<InvFFTSizeX2;i++) FocusedAzLine[i]=0.0;

	/* Re-order array for IFFT */
	indx = 2*(InvFFTSize - (FFTSize/2 - LookFreqStartPt));
	while (indx >= InvFFTSizeX2)
	  indx = indx - InvFFTSizeX2;
	for (i=0; i<LookFreqPtsX2;i++)
	  {
	  FocusedAzLine[indx++] = DechirpedLook[i];
	  if (indx == InvFFTSizeX2) indx = 0;
	  }

#else
	/* zero pad rest of focused az line in freq domain */
	for (i=2*LookFreqPts;i<InvFFTSizeX2;i++) FocusedAzLine[i]=0.0;
#endif

	/* Inverse FFT back to time domain */
	dummyaddress = FocusedAzLine - 1;
	CFFT(dummyaddress,InvFFTSize,1);

	/* Az shift, scale and add to Detected array. Scaling here is
	   slightly inefficient, if lots of looks. Data added directly
	   to PostSummedAzLine if no post summing */
	indxi = move;
	indxq = indxi+1;
	indx = 0; /* initialise, in case complex output */
	for (i=0;i<AzPtsPrePostSum;i++)
	  {
	  if (DetectMethod == 2 || DetectMethod==3)     /* power output */
	    {
	    tmpi = FocusedAzLine[indxi++];
	    tmpq = FocusedAzLine[indxq++];
	    if (PostSumRatio==1)
	      {PostSummedAzLine[i] += (tmpi*tmpi+tmpq*tmpq)*Scale;}
	    else
	      {DetectedAzLine[i] += (tmpi*tmpi+tmpq*tmpq)*Scale;}     
	    }
          else if (DetectMethod == 1) /* mag output */
            {
            tmpi = FocusedAzLine[indxi++];
            tmpq = FocusedAzLine[indxq++];
            if (PostSumRatio==1)
	      {PostSummedAzLine[i] += sqrt(tmpi*tmpi+tmpq*tmpq)*Scale;}
	    else
	      {DetectedAzLine[i] += sqrt(tmpi*tmpi+tmpq*tmpq)*Scale;}
	    }
	  else                        /* complex output (single look) */
	    {
            if (PostSumRatio==1)
	      {
	      PostSummedAzLine[indx++] = FocusedAzLine[indxi++]*Scale;
	      PostSummedAzLine[indx++] = FocusedAzLine[indxq++]*Scale;
	      }
            else
	      {
	      DetectedAzLine[indx++] = FocusedAzLine[indxi++]*Scale;
	      DetectedAzLine[indx++] = FocusedAzLine[indxq++]*Scale;
              }  
	    }  /* end else complex output */
          indxi++;
	  indxq++;
          }  /* end for i loop */
 
       } /* end look loop */

 
       /* PERFORM POST SUM  */
       if (PostSumRatio!=1)
         {
         if (DetectMethod==0)  /* complex post sum */
	   {
           SumI = 0.0;
	   SumQ = 0.0;
	   indx = 0;
           indx2 = 0;
           for (i = 0; i < AzPtsPrePostSum; i++)
             {
             SumI += DetectedAzLine[indx++];
	     SumQ += DetectedAzLine[indx++];
             if ( ((i+1)%PostSumRatio) == 0 )
               {
               PostSummedAzLine[indx2++] = SumI*Scale;
	       PostSummedAzLine[indx2++] = SumQ*Scale;
	       SumI = 0.0;
	       SumQ = 0.0;
               }
             } /* end for i loop */
	   }   /* end if complex post sum */
	 else  /* mag or power output post sum */
	   {
	   magsum = 0.0;
           indx = 0;
           for (i = 0; i < AzPtsPrePostSum; i++)
             {
             magsum += DetectedAzLine[i];
             if ( ((i+1)%PostSumRatio) == 0 )
               {
               PostSummedAzLine[indx++] = magsum*Scale;
               magsum = 0.0;
               }
             } /* end for i loop */
	   }  /* end else for mag or power output post sum */
       	 }  /* end if post sum */

 

       /* CHECK FOR MAX ABS VALUE, if required */
       if (ReportMax!=0)
	 for (i=0; i<OutputAzRealPts; i++)
            if (fabs(PostSummedAzLine[i]) > MaxAbsValue)
              {
              if (DetectMethod != 3)
		{ MaxValue = PostSummedAzLine[i]; } /* no DC offset added */
              else
                { MaxValue = 10.0*log10(PostSummedAzLine[i]); }
              MaxAbsValue = fabs(PostSummedAzLine[i]);
              MaxValueInputRngBin = InputFileRngBin;
              MaxValueAzPixel = (Int4B)((double)i/(double)SizeFactor); 
             }			

      /* CONVERT TO dB IF REQUIRED (needs to be after max value check
         as small values map to large negative dB values)  */
       if (DetectMethod == 3)
         {
         for (i=0; i<OutputAzRealPts; i++)
           {
           if (PostSummedAzLine[i] != 0.0)
             { PostSummedAzLine[i] = 10.0*log10(PostSummedAzLine[i]);}
           else
             { 
             PostSummedAzLine[i] = DB_VALUE_FOR_ZERO; 
             dBMinCounter++;
             }
	   }
         }

       
       /* CHECK FOR OVERFLOWS AND ADD DC OFFSET, if applicable. DC
	  offset only added for OutputDataTypes 0 and 1 and only if
	  data complex. Overflows and zeros checked only for
	  OutputDataTypes 0 and 1 (checked for all detection methods).
	  Overflows and zeros for complex data are counted for both I and Q
	  components (thus may be higher than the no. of az pts) */
       if (OutputDataType < 2)
         {
	 if (DetectMethod == 0) {indxi=0; indxq=1;}
	 for (i=0;i<OutputAzPts;i++)
	   {

	   if (DetectMethod==0)  /* complex output */
	     {
             tmpi = PostSummedAzLine[indxi] + OutputDCOffset + 0.5;
             tmpq = PostSummedAzLine[indxq] + OutputDCOffset + 0.5;

	     if (tmpi > (double)OutputLimit)  /* check for I channel */
	       { tmpi = OutputLimit; OverFlows++;}
	     else if (tmpi < 1.0)
	       { tmpi = 0.0; OverFlows++;}
	     else if (tmpi < (OutputDCOffset+1) && tmpi >= OutputDCOffset )
	       { Zeros++; }   

	     if (tmpq > (double)OutputLimit)  /* check for Q channel */
	       { tmpq = OutputLimit; OverFlows++;}
	     else if (tmpq < 1.0)
	       { tmpq = 0.0; OverFlows++;}
	     else if (tmpq < (OutputDCOffset+1) && tmpq >= OutputDCOffset )
	       { Zeros++; }   

             if (OutputDataType==0) /* unsigned char output */
               {
	       OutputUnChar[indxi++] = (unsigned char)tmpi;
               OutputUnChar[indxq++] = (unsigned char)tmpq;
               }   
             else  /* must be unsigned short */
               {
	       OutputUnShort[indxi++] = (unsigned Int2B)tmpi;
               OutputUnShort[indxq++] = (unsigned Int2B)tmpq;
	       }
	     indxi++;
	     indxq++;
	     }   /* end if complex output */
	   else  /* mag, power or power dB output */
	     {
             tmp = PostSummedAzLine[i]+0.5;
	     if (tmp > (double)OutputLimit)
	       { tmp = OutputLimit; OverFlows++; }
	     else if (tmp < 1.0)
	       { Zeros++; }

	     if (OutputDataType==0) /* unsigned char output */
               OutputUnChar[i] = (unsigned char)tmp;   
             else  /* must be unsigned short */
               OutputUnShort[i] = (unsigned Int2B)tmp;
	     } /* end else mag or power ouput */

	   }  /* end for i < OutputAzPts loop */

	 }  /* end if OutputDataType < 2 (CHECK FOR OVERFLOWS)*/

       
      /* ASSIGN OUTPUT ARRAYS FOR HIGHER DATA TYPES. No DC offsets
	 added for these types. The output arrays for the lower types
	 have already been assigned in the previous section. */
      if (OutputDataType >= 2)  
	{
	if (OutputDataType==2)
	  for (i=0;i<OutputAzRealPts;i++)
	    OutputLongInt[i]=(Int4B )(PostSummedAzLine[i]+0.5);   
        else if (OutputDataType==3)	    
          for (i=0;i<OutputAzRealPts;i++)
	    OutputFloat[i] = (float)PostSummedAzLine[i];
	else  /* must be of type 4 - double */
          for (i=0;i<OutputAzRealPts;i++)
	    OutputDouble[i] = (double)PostSummedAzLine[i];
	}  /* end if OutputDatatype >= 2 */
       
      /* WRITE DATA TO OUTPUT FILE */
      if (OutputDataType==0)
        fwrite(OutputUnChar,sizeof(unsigned char),
	       OutputAzRealPts,OutputFile);
      else if (OutputDataType==1)
        fwrite(OutputUnShort,sizeof(unsigned Int2B),
	       OutputAzRealPts,OutputFile);
      else if (OutputDataType==2)
        fwrite(OutputLongInt,sizeof(Int4B),
	       OutputAzRealPts,OutputFile);
      else if (OutputDataType==3)
        fwrite(OutputFloat,sizeof(float),
	       OutputAzRealPts,OutputFile);
      else   /* must be of type 4 - double */
        fwrite(OutputDouble,sizeof(double),
	       OutputAzRealPts,OutputFile);	

      /* Report progress to stdout */
      if (modf((double)bin/ScreenUpdateRate,&ipart)==0) {  
         fprintf(stdout,
		"Focus rng / output rng bin     : %.4f m / %ld\r",
		    FocusRng,bin); 
      }
      
      /* Update InputFileRngBin for next loop interation */
      InputFileRngBin++; 
 
    }   /* End bin within focus segment loop */

     
     /* Update focusing segment. Enlarge last focus segment in case
	rng focus updates not a factor of the rng bins to process
        (ensures all rng bins processed) */
     RngBinStartFocSeg += (Int4B)(RngBinsToProcess/RngFocSegments);
     RngBinEndFocSeg   += (Int4B)(RngBinsToProcess/RngFocSegments);
     if (FocusSegment == RngFocSegments-2)
       RngBinEndFocSeg = RngBinsToProcess - 1;    /* fixed 1998-01-21 */

     /* Free window function */
     free(Window);

  }   /* END FocusSegment LOOP */


  if (msg==stdout)
    { fprintf(msg,"%60s\n"," "); }
   
  if (RngCurvInterpSize!=0)
    fprintf(msg,
       "Extra zeros used at far rng    : %ld (by rng curv correction)\n",
	   FarRngPad);

  /* Report on overflows, zeros and max value - to assist scaling */
  if (OutputDataType < 2) 
    {
      fprintf(msg,"Percentage overflows/zeros     : %.4e / %.4e\n", 
           (double)OverFlows*100.0/((double)OutputAzRealPts*
                  (double)RngBinsToProcess),
           (double)Zeros*100.0/((double)OutputAzRealPts*
                  (double)RngBinsToProcess) );     
    }


  if (ReportMax!=0)
    {
    if (OutputDCOffset != 0.0)
      { 
      fprintf (msg,
        "Max value after processing     : %8.6e (incl. DC offset %ld)\n",
		 MaxValue+(double)OutputDCOffset,OutputDCOffset); 
      }
    else 
      {
      fprintf (msg,
        "Max value after processing     : %8.6e (no DC offset)\n",
          MaxValue); 
      }
    fprintf (msg,"Max at output az / rng pixel   : %ld / %ld\n",
             MaxValueAzPixel,MaxValueInputRngBin-StartProcessRngBin);  
    }

  if (DetectMethod == 3)
    fprintf(msg,
        "Zeros set to min dB val. (%d): %ld\n",
         DB_VALUE_FOR_ZERO,dBMinCounter); 
  
  /* Get end time */
  EndTime = time(NULL);
  fprintf (msg,"\nAZ COMPRESSION DONE - in %ld secs (%.2f min)\n",
	   EndTime-StartTime,(double)(EndTime-StartTime)/60.0);
  if (msg!=stdout) {
    fprintf(stdout,"\nAz Compression Done!\n");
  }


  
  /******************************
  * Free memory and close files *
  ******************************/

  if (strcmp(LogFileName,"null")!=0) {
    fclose(logfile);    /* first close log file so can see what happened */
  }
 
  if (RngCurvInterpSize != 0)
    {
     free (CorrectedBatchBlock);
     free (CorrectedBatchFreq);
     free (RngCurvedBlock);
     free (RngCurvedBatchFreq);
    }

  free (InputIQ);
  free (MasterAzRef);
  free (LookRefFunc);
  free (LookRefBlock);
  free (FocusedAzLine);
  free (ShiftedMasterAzRef);
  if (PostSumRatio!=1) free (DetectedAzLine);
  free (PostSummedAzLine);

  if (InputDataType == 0) free(IQCharInput);
  else if (InputDataType == 3) free(IQFloatInput);
 
  if   (OutputDataType==0) free(OutputUnChar);
  else if (OutputDataType==1) free(OutputUnShort);
  else if (OutputDataType==2) free(OutputLongInt);
  else if (OutputDataType==3) free(OutputFloat);
  else free(OutputDouble);

#if REORDER_IFFT
  free (DechirpedLook);
#endif


  fclose(InputFile);
  fclose(OutputFile);

  return(0);  /* for successful completion */
}  /* End function g2azc */


/************************************************************************/
/* FUNCTION : WriteAzComTmplCmdFile() */
Int2B WriteAzComTmplCmdFile(char CmdFileName[])
{
  FILE *OutFile; 
  if ( (OutFile = fopen (CmdFileName, "w") ) == NULL ) return(0);

  fprintf(OutFile,"Command file for Azcom (SAR Azimuth Compression)\n");
  fprintf(OutFile,"$ProgramVersion (jmh) => %s\n",PROG_VERSION);
  fprintf(OutFile,"------------------------------------------------\n\n");
  fprintf(OutFile,"$ScreenUpdateRate                    => 10\n");
  fprintf(OutFile,"$LogFileName ('null' for none)       => null\n");
  fprintf(OutFile,"$InputStartSampleDelay [see note]    => 9.0e-05\n");
  fprintf(OutFile,"$CarrierFreq [Hz]                    => 1.41e+08\n");
  fprintf(OutFile,"$InputPRF [Hz]                       => 500.0\n");
  fprintf(OutFile,"$NomGroundSpeed [m/s]                => 250.0\n");
  fprintf(OutFile,"$InputFileAzPts                      => 1001\n");
  fprintf(OutFile,"$StartProcessAzPt                    => 0\n");
  fprintf(OutFile,"$AzPtsToProcess                      => 1001\n");
  fprintf(OutFile,"$InputFileRngBins                    => 2048\n");
  fprintf(OutFile,"$StartProcessRngBin                  => 0\n");
  fprintf(OutFile,"$RngBinsToProcess [see note]         => 2048\n");
  fprintf(OutFile,"$InputDCOffsetI                      => 127.0\n");
  fprintf(OutFile,"$InputDCOffsetQ                      => 127.0\n");
  fprintf(OutFile,"$InvFFTSizeReduc [power of 2]        => 1\n");
  fprintf(OutFile,"$InputFileName                       => tmpg2.cor\n");
  fprintf(OutFile,"$OutputFileName                      => tmpg2.azc\n");
  fprintf(OutFile,"$AppendExistOutFileFlg [Y/N]         => N\n");
  fprintf(OutFile,"$RngFocSegments [see note]           => 128\n");
  fprintf(OutFile,"$RefFuncSign [+-1]                   => -1\n");
  fprintf(OutFile,"$A2DFreq [Hz]                        => 30.0e+06\n");
  fprintf(OutFile,"$NomAzRes [m - see note]             => 20.0\n");
  fprintf(OutFile,"$WinConstTime [0.0-1.0 - see note]   => 0.08\n");
  fprintf(OutFile,"$NumLooks                            => 1\n");
  fprintf(OutFile,"$LookOverlapFrac [0.0-1.0]           => 0.5\n");
  fprintf(OutFile,"$WinConstFreq [0.0-1.0 - see note]   => 1.0\n");
  fprintf(OutFile,"$RngCurvInterpSize [see note]        => 8\n");
  fprintf(OutFile,"$RngCurvBatchSize [see note]         => 256\n");
  fprintf(OutFile,"$PostSumRatio                        => 1\n");
  fprintf(OutFile,"$DetectMethod [see note]             => 2\n");
  fprintf(OutFile,"$InputDataType [see note]            => 0\n");
  fprintf(OutFile,"$OutputDataType [see note]           => 0\n");
  fprintf(OutFile,"$Scale [see note]                    => 2.8e-07\n");
  fprintf(OutFile,"$ReportMax [1/0]                     => 1\n");

  fprintf(OutFile,"\nNotes:\n");
  fprintf(OutFile,"------\n\n");
  fprintf(OutFile,
    "InputStartSampleDelay : [in secs]\n"); 
  fprintf(OutFile,
    "     This value is coresponds to the start of the input file and is\n");
  fprintf(OutFile,
    "     independent of which range bin is selected here for processing.\n");
  fprintf(OutFile,
    "     Note the value may not be the same as in the range compression\n");
  fprintf(OutFile,
    "     program. There are two reasons for this:\n");;
  fprintf(OutFile,
    "       1) Early range compression version introduce a delay of half pulse\n");
  fprintf(OutFile,
    "          length (subtract here - only for rngcom version < 1999-01-27).\n");
  fprintf(OutFile,
    "       2) A subsection of range bins may have been selected in the\n");
  fprintf(OutFile,
    "          compression/corner turn.\n\n");      
  fprintf(OutFile,
    "RngBinsToProcess : If required for the range curvature correction,\n");
  fprintf(OutFile,
    "                   the program will automatically make use of more\n"); 
  fprintf(OutFile,
    "                   range bins than specified here, if available.\n\n"); 
  fprintf(OutFile,
    "RngFocSegments : -1 - max number of az. ref. function updates\n");
  fprintf(OutFile,
    "               : else less than or equal to RngBinsToProcess\n\n");

  fprintf(OutFile,
    "NomAzRes : The nominal azimuth resolution in metres assuming a\n");
  fprintf(OutFile,
    "           window broadening factor of unity.\n\n");
  fprintf(OutFile,
    "WinConstTime : The constant for the window to be applied over the\n");
  fprintf(OutFile,
    "               length of the time-domain reference function.\n");
  fprintf(OutFile,
    "             : 1.0 - rectangular window\n");
  fprintf(OutFile,
    "             : 0.08 - Hamming window\n\n");
  fprintf(OutFile,
    "WinConstFreq : The constant for the window to be applied over the\n");
  fprintf(OutFile,
    "               length of each look in the azimuth frequency domain.\n");
  fprintf(OutFile,
    "             : 1.0 - rectangular window\n");
  fprintf(OutFile,
    "             : 0.08 - Hamming window\n\n");
  fprintf(OutFile,"RngCurvInterpSize : 0 - none\n");
  fprintf(OutFile,"                  : 1 - nearest neighbour\n");
  fprintf(OutFile,"                  : even (8 pt or more recommended)\n\n");
  fprintf(OutFile,
    "RngCurvBatchSize : The number of range bins to undergo range\n");
  fprintf(OutFile,
    "                   curvature correction at once. A smaller number\n");
  fprintf(OutFile,
    "                   will reduce memory requirements, but run slower.\n\n"); 
  fprintf(OutFile,"DetectMethod : 0 - none (complex output)\n");
  fprintf(OutFile,"             : 1 - magnitude\n");
  fprintf(OutFile,"             : 2 - power\n");
  fprintf(OutFile,"             : 3 - power in dB\n\n");
  fprintf(OutFile,"OutputDataType : 0 - unsigned char\n");
  fprintf(OutFile,"               : 1 - unsigned short int (2 bytes)\n");
  fprintf(OutFile,"               : 2 - long int (4 bytes)\n");
  fprintf(OutFile,"               : 3 - float (4 bytes)\n");
  fprintf(OutFile,"               : 4 - double (8 bytes)\n\n");
  fprintf(OutFile,"InputDataType : 0 - unsigned char IQ (2*1 bytes per point)\n");
  fprintf(OutFile,"              : 3 - float IQ (2*4 bytes per point)\n\n");
  fprintf(OutFile,"Scale : 0 - auto (not implemented yet)\n");
  fprintf(OutFile,"      : 1 - none\n");
  fprintf(OutFile,"      : other (double)\n\n");

  fclose(OutFile);
  return(0);
} 


/*****************************************************************/
/* FUNCTION RngCurvStraighten */

Int2B RngCurvStraighten(double A2DFreq,
		        float  **CorrectedBatchFreq,
		        Int4B  FFTSize,
		        Int4B  MasterRefFreqPtsD2,
		        Int4B  *NextRngBinToStraighten,
		        double RngBinSize,
		        Int4B  RngCurvBatchSize,
		        double RngCurvCalc,
		        float  **RngCurvedBatchFreq,
		        Int4B  RngCurvInterpSize,
		        double *RngCurvRefRng,
		        Int4B  *StraightenedRngBins,
		        Int4B  *StraightenedRngBinsUsed
		        )

{
  Int4B k,i,indxi,indxq,indxi2,indxq2,
           closebin,         /* Closest bin to rng curved data */
           FFTSizeX2,        /* Twice FFT Size */
           intbin,           /* Interpolation bin */
           RngBinForCorrection,
           RngCurvInterpSizeD2=0, 
           tmpindx;
  
  float    SincFactor,
           tmpifloat,tmpqfloat,tmpifloat2,tmpqfloat2;
  
  double   arg,
           InterpPosnDesired,
           RngCurvFreq,     /* Range curvature in freq domain */
           RngCurvTime2Way;

  /* Misc */
  FFTSizeX2 = 2*FFTSize;
  if (RngCurvInterpSize>1) {RngCurvInterpSizeD2 = RngCurvInterpSize/2;}
  
  /* Initialize straightened array */
     for (k=0;k<RngCurvBatchSize;k++)
       for (i=0;i<FFTSizeX2;i++) CorrectedBatchFreq[k][i] = 0.0;
    
     *StraightenedRngBins = 0;
     for (RngBinForCorrection = *NextRngBinToStraighten;
	  RngBinForCorrection<(*NextRngBinToStraighten+RngCurvBatchSize);
	  RngBinForCorrection++)
        {
         for (k=0; k<MasterRefFreqPtsD2; k++) /* Correct freq pairs */
            {

             /* Calc RngCurvFreq */
             RngCurvFreq = *RngCurvRefRng*(1.0/sqrt(1.0-k*k*RngCurvCalc)
			       - 1.0);
             RngCurvTime2Way = 2*RngCurvFreq/C;
             closebin = (Int4B)(RngCurvFreq/RngBinSize + 0.5);
             intbin = closebin+*StraightenedRngBins;

             indxi = 2*k+FFTSize;
             indxq = indxi+1;
             indxi2 = FFTSize-(k+1)*2;
             indxq2 = indxi2+1;

             if (RngCurvInterpSize == 1)  /* Nearest neighbour */
               {
                CorrectedBatchFreq[*StraightenedRngBins][indxi]=
		  RngCurvedBatchFreq[intbin][indxi];
                CorrectedBatchFreq[*StraightenedRngBins][indxq]=
		  RngCurvedBatchFreq[intbin][indxq];
                CorrectedBatchFreq[*StraightenedRngBins][indxi2]=
		  RngCurvedBatchFreq[intbin][indxi2];
                CorrectedBatchFreq[*StraightenedRngBins][indxq2]=
		  RngCurvedBatchFreq[intbin][indxq2];
	       }
             else /* Perform interpolation (+ve and -ve freqs simul) */
               {
                InterpPosnDesired = A2DFreq*RngCurvTime2Way +
		                     (double)RngCurvInterpSizeD2;
                tmpindx = closebin;
                tmpifloat = 0.0;
                tmpqfloat = 0.0;
                tmpifloat2 = 0.0;
                tmpqfloat2 = 0.0;
                for (i=0;i<RngCurvInterpSize;i++)
                  {
                   arg = PI*( InterpPosnDesired-(double)(tmpindx++) );
                   if (arg<1.0e-5 && arg>-1.0e-5)
		     {SincFactor = 1.0;}
                   else
		     {SincFactor = (float)(sin(arg)/arg);}
                   tmpifloat += RngCurvedBatchFreq[intbin][indxi]
		                  *SincFactor;
                   tmpqfloat += RngCurvedBatchFreq[intbin][indxq]
		                  *SincFactor;
                   tmpifloat2 += RngCurvedBatchFreq[intbin][indxi2]
		                  *SincFactor;
                   tmpqfloat2 += RngCurvedBatchFreq[intbin++][indxq2]
		                  *SincFactor;
                  }
                CorrectedBatchFreq[*StraightenedRngBins][indxi]
		  = tmpifloat;
                CorrectedBatchFreq[*StraightenedRngBins][indxq]
		  = tmpqfloat;
                CorrectedBatchFreq[*StraightenedRngBins][indxi2]
		  = tmpifloat2;
                CorrectedBatchFreq[*StraightenedRngBins][indxq2]
		  = tmpqfloat2;
               }   /* End else interpolate */
 
            } /* End k (freq pt) loop */

         /* Increment RngCurvRefRng and StraightenedRngBins for
	    straightened array */
         *RngCurvRefRng += RngBinSize;
         *StraightenedRngBins += 1;   /* don't use ++ here, unless brackets*/
	 
        }  /* End RngBinForCorrection loop */

     /* Incr NextRngBinToStraighten
        (bins up to (not incl) NextRngBinToStraighten 
         have been straightened) */
     *NextRngBinToStraighten += RngCurvBatchSize;

     /* Initialize StraightenedRngBins_used to indicate
            CorrectedBatchFreq array primed */
     *StraightenedRngBinsUsed = 0;

     return(0);
}  /* end RngCurvStraighten function */

