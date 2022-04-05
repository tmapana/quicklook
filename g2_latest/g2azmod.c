/* Program to perform azimuth compression (multiple range bins)
   using FFT approach. */
/* Author: J.M. Horrell  (Copyright 1994-1997) */
/* UNAUTHORIZED USE PROHIBITED !! */

/* CHANGES: */
/* Ver. Gp2az7d - version of azcom7d.c modified for SASAR GP2 ground
    processor. Different header file, cleaned up %f specifiers, Int2B
    and Int4B types. */
/* Ver g2azc1a - change into a function. Add option to append output file. */
/* Ver g2azc1b - (1997-11-28) allow compilation as a function or standalone */
/*             - (1997-12-02) tidy, display az coord of max, rename variable
                 AzresLook to NomAzRes (also changes in cmd file), add parse
                 error checks, more help in template command file. */
/* Ver G2Azc (1997-12-05) - allow power output in dB, check program
                 version in command file */
/* Ver. (1997-12-08) - rename RngFocUpdates to RngFocSegments, ensure num 
        focus updates does not change the num of output rng bins, allow 
        RngfocSegments to be < 1 (for max updates), fix bugs in dB output 
        option, fix dB bug for zero */             

/* 1998-01-21 - fixed range focus segment bug. */
/* g2azc_mod.c (1998-02-26) - modified version to give float complex output after 
           the  az FFT and also after the range curvature correction interpolation
           (in range-Doppler domain). Files written specified below in define. 
           The size of the two output R-D files should be the same if all the 
           input file range bins are processed. Else the corrected file will 
           contain only the range bins processed while the R-D uncorrected file
           contains extra rbins for the range curv straightening. If crash on 
           running, then be skeptical as to correct R-D output. Try adjusting 
           processing parameters so as in reasonable range. */
/* 1998-03-05 - added in R-D spread az res compensation function.
           Only compensates for one range */
/* 1998-03-12 - reorder IFFT */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "g2func.h"

/* SET AS FUNCTION OR STANDALONE */
#define FUNC 0


/* SET TO PRODUCE RANGE-DOPPLER DOMAIN OUTPUT OR ENABLE COMPENSATION
OF AZ RES FOR RANGE-DOPPLER SPREADING (NOT BOTH!!!!)*/
#define R_D_OUTPUT 1
#define R_D_COMP 0

/* SET TO REORDER IFFT */
#define REORDER_IFFT 1


#define R_D_FILE_NAME "rd.dat"
#define R_D_CORR_FILE_NAME "rd_corr.dat"
#define R_D_COMP_FILE_NAME "rdcompfile.dat"
#define R_D_COMP_OUTPUT_FILE_NAME "rdcompoutfile.dat"
#define R_D_COMP_FILE_DATA_POINTS 32768

/* Misc defns limited to this file */
#define PROG_VERSION "1997-12-08"
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

#if FUNC
Int2B G2Azc (FILE *msg, char SpecFileName[])
#else
Int2B main (int argc,char *argv[])
#endif
{
  FILE *OutputFile,*cmdfp,*InputFile;
#if !FUNC 
  FILE *msg=NULL;
#endif
#if R_D_OUTPUT
  FILE *RDFile, *RDCorrFile;
#elif R_D_COMP
  FILE *RDCompFile;
  /*  FILE *RDCompOutputFile;*/
#endif


  char InputFileName[STRING_SPACE],OutputFileName[STRING_SPACE],
       AppendExistOutFileFlg,ProgVersion[STRING_SPACE];

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
           InputDCOffset,
           FarRngPad=0,   /* No. of rng curv int far rng bins zero padded */
           FFTSize,       
           FFTSizeX2,        /* Twice the FFT pts */
           FirstLookFreqStartPt,
	   FocusSegment,      /* Focussing segment index */
           InvFFTSize,       /* (power of 2) */
           InvFFTSizeX2,
           InputFileAzPts,
           InputFileRngBin,   
           InputFileRngBins,  
           InvFFTSizeReduc,  /* (pow of 2) */
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
         FocusRng,
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

#if R_D_OUTPUT
  float *RDArray=NULL,
        *RDCorrArray=NULL;
#elif R_D_COMP
  float *RDCompArray=NULL;
#endif
 

#if FUNC
  fprintf (msg,"\nAZIMUTH COMPRESSION STAGE...\n");

  /* Open spec file as binary */
  if ( (cmdfp = fopen (SpecFileName, "rb") ) == NULL )
    {
    fprintf(msg,"ERROR: Command file %s not found/opened\n",SpecFileName );
    exit(1);
    }
#else

  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\nProg: AZCOM (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (Copyright UCT 1997)\n");

  if (argc != 2)
    { 
    fprintf(msg,"Azimuth compression for SAR data.\n\n");
    fprintf(msg,"USAGE: azcom [cmd file]\n");
    fprintf(msg,"To see command file structure, type `azcom -tmpl'\n"); 
    fprintf(msg,"(generates a template command file `azcom.tmp')\n\n");  
    exit(1); 
    }

  /* Write template command file, if requested */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0)
    {
     if (WriteAzComTmplCmdFile("azcom.tmp"))
       fprintf(msg,"Template command file `azcom.tmp' written!\n");
     exit(1);
    }

  fprintf(msg,"\n");

  /* Open spec file */
  if ( (cmdfp = fopen(argv[1],"rb")) == NULL)
    {
    fprintf(msg,"ERROR - command file %s not found/opened!\n",argv[1]);
    exit(1);
    }
  
#endif

  /*********************
  * Parse command file *
  **********************/
  
  /* Check version IDs match */
  if (!SeekP(cmdfp,"$ProgramVersion","=>",-1,0)) 
    { fprintf(msg,
      "WARNING - in parsing command file (program version not found)!\n\n"); } 
  else
    {
    ReadStr(cmdfp,ProgVersion,STRING_SPACE); 
    if (strcmp(ProgVersion,PROG_VERSION)!=0)
      {
      fprintf(msg,
        "WARNING - command file version (%s) not same as program (%s)!\n\n",
        ProgVersion,PROG_VERSION);
      }
    } /* end else */

  Errors=0;

  if (!SeekP(cmdfp,"$ScreenUpdateRate","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&ScreenUpdateRate);
  if (!SeekP(cmdfp,"$InputStartSampleDelay","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&InputStartSampleDelay);
  if (!SeekP(cmdfp,"$CarrierFreq","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&CarrierFreq);
  if (!SeekP(cmdfp,"$InputPRF","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&InputPRF);
  if (!SeekP(cmdfp,"$NomGroundSpeed","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&NomGroundSpeed);
  if (!SeekP(cmdfp,"$InputFileAzPts","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &InputFileAzPts);
  if (!SeekP(cmdfp,"$StartProcessAzPt","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &StartProcessAzPt);
  if (!SeekP(cmdfp,"$AzPtsToProcess","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &AzPtsToProcess);
  if (!SeekP(cmdfp,"$InputFileRngBins","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &InputFileRngBins);
  if (!SeekP(cmdfp,"$StartProcessRngBin","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &StartProcessRngBin);
  if (!SeekP(cmdfp,"$RngBinsToProcess","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &RngBinsToProcess);
  if (!SeekP(cmdfp,"$InputDCOffset","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &InputDCOffset);
  if (!SeekP(cmdfp,"$FFTSize","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &FFTSize);
  if (!SeekP(cmdfp,"$InvFFTSizeReduc","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &InvFFTSizeReduc);
  if (!SeekP(cmdfp,"$InputFileName","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,InputFileName,STRING_SPACE); 
  if (!SeekP(cmdfp,"$OutputFileName","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,OutputFileName,STRING_SPACE);
  if (!SeekP(cmdfp,"$AppendExistOutFileFlg","=>",-1,0)) Errors++;
  ReadChar(cmdfp,&AppendExistOutFileFlg);
  if (!SeekP(cmdfp,"$RngFocSegments","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &RngFocSegments);
  if (!SeekP(cmdfp,"$RefFuncSign","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &RefFuncSign);
  if (!SeekP(cmdfp,"$A2DFreq","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&A2DFreq);
  if (!SeekP(cmdfp,"$NomAzRes","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&NomAzRes);
  if (!SeekP(cmdfp,"$WinConstTime","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&WinConstTime);
  if (!SeekP(cmdfp,"$NumLooks","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &NumLooks);
  if (!SeekP(cmdfp,"$LookOverlapFrac","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&LookOverlapFrac);
  if (!SeekP(cmdfp,"$WinConstFreq","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&WinConstFreq);
  if (!SeekP(cmdfp,"$RngCurvInterpSize","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &RngCurvInterpSize);
  if (!SeekP(cmdfp,"$RngCurvBatchSize","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &RngCurvBatchSize);
  if (!SeekP(cmdfp,"$PostSumRatio","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &PostSumRatio);
  if (!SeekP(cmdfp,"$DetectMethod","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &DetectMethod);
  if (!SeekP(cmdfp,"$OutputDataType","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &OutputDataType);
  if (!SeekP(cmdfp,"$Scale","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&Scale);
  if (!SeekP(cmdfp,"$ReportMax","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld", &ReportMax);
  
  fclose(cmdfp);

  if (Errors!=0)
    {fprintf(msg,"WARNING - %ld errors in parsing command file!\n\n",Errors);} 

  /* Calc misc and speedup variables */
  if (RngCurvInterpSize>1)
    {RngCurvInterpSizeD2 = RngCurvInterpSize/2;}
  else
    {RngCurvInterpSizeD2 = 0;}
  Wavelength = C/CarrierFreq;
  EffLooks = (double)NumLooks - ((double)NumLooks-1.0)*LookOverlapFrac;
  AzResFull = NomAzRes/EffLooks;
  FFTSizeX2 = FFTSize*2;
  AzPtsToProcessX2 = AzPtsToProcess*2;
  InvFFTSize = (Int4B)((double)FFTSize/(double)InvFFTSizeReduc);
  InvFFTSizeX2 = 2*InvFFTSize;
  OutputAzPts = (Int4B)((double)AzPtsToProcess/
                ((double)InvFFTSizeReduc*(double)PostSumRatio));
  if (DetectMethod==0)
    SizeFactor = 2;  /* to double size of complex output arrays */ 
  else /* for DetectMethod 2, 3 or 4 */
    SizeFactor = 1;
  OutputAzRealPts = OutputAzPts*SizeFactor;
  AzPtsPrePostSum = OutputAzPts*PostSumRatio; /* az pts prior to postsum */
  AzRealPtsPrePostSum = AzPtsPrePostSum*SizeFactor;
  SkipInputBytes = 2*(InputFileAzPts-AzPtsToProcess); /* after each input line read */
  RngBinSize = C/(2.0*A2DFreq);
  StartProcessRng = C*InputStartSampleDelay/2.0 + StartProcessRngBin*
                       RngBinSize;
  if (DetectMethod == 2 || DetectMethod == 3)
    {  Scale = Scale/((double)FFTSize*(double)FFTSize); } /* note scale changed here */
  else 
    { Scale = Scale/(double)FFTSize; }
  InputFileRngBin = StartProcessRngBin;  

  /* Check and arrange output (DC offsets applied only to complex
       output case and only for unsigned char and unsigned short output) */

  if (OutputDataType!=0 && OutputDataType!=1 && OutputDataType!=2 &&
      OutputDataType!=3 && OutputDataType!=4)
    {
    fprintf(msg,"ERROR - output data type (%ld) invalid!!\n",OutputDataType);
    exit(1);
    }

  if (OutputDataType==0) OutputLimit = LIMIT_UNCHAR;
  else if (OutputDataType==1) OutputLimit = LIMIT_UNSHORT;
  
  switch(DetectMethod)
    {
    case 0:    /* complex output */
      if (NumLooks != 1)
	{
	fprintf(msg,"ERROR - complex output only with single look!\n");
	exit(1);
	}
      switch(OutputDataType)
	{ 
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
  if ( pow(2.0,tmp) != (double)FFTSize )
    { 
    fprintf(msg,"ERROR - FFTSize %ld not power of two!\n",FFTSize); 
    exit(1); 
    }
  tmp = (double)((Int4B)( log((double)InvFFTSizeReduc+0.1)/log(2.0) ) );
  if ( pow(2.0,tmp) != (double)InvFFTSizeReduc )
    { 
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

  if (LookFreqPts > InvFFTSize)
    {
    fprintf(msg,"\nERROR : InvFFTSize (%ld) < LookFreqPts (%ld)\n",
                   InvFFTSize,LookFreqPts);
    exit(1);
    }

  /* Check total bandwidth used less than InputPRF */
  TotalDopBw = EffLooks*LookDopBw;
  MasterRefFreqPtsD2 = (Int4B)(TotalDopBw*(double)FFTSize/
			       (2.0*(double)InputPRF));
  if( TotalDopBw > InputPRF )
    {
    fprintf(msg,"ERROR - Total Doppler bandwidth needed (%.3f Hz) is > PRF\n",
             TotalDopBw);
    exit(0);
    }

  /* Check processing block in limits */
  if (StartProcessRngBin+RngBinsToProcess > InputFileRngBins)
   {
   fprintf(msg,"ERROR - attempt to process beyond valid range data!\n");
   exit(0);
   }

  if (StartProcessAzPt+AzPtsToProcess > InputFileAzPts)
   {
   fprintf(msg,"ERROR - attempt to process beyond valid az data!\n");
   exit(0);
   }

  /* Calculate range curvature at mid- and far processed swath (the max)*/
  FocusRng = StartProcessRng + RngBinSize*
                (double)RngBinsToProcess/2.0;  /* mid processing swath */
  MidCurvRngBins = (FocusRng/RngBinSize)*
    ( sqrt(1.0+ 1.0/( 16.0*AzResFull*AzResFull/
		      (Wavelength*Wavelength)-1.0 ))-1.0 );
  FocusRng = StartProcessRng + RngBinSize*
                (double)RngBinsToProcess;      /* far processed swath */
 
#if R_D_OUTPUT
  /* Modify this for R_D_OUTPUT to allocate large enough array in range 
     for full az bandwidth correction. May need to add one rbin, if crashes*/
 MaxCurvRngBins = 2 + (FocusRng/RngBinSize)*
    ( sqrt(1.0+ 1.0/( 16.0*(NomGroundSpeed/InputPRF)*(NomGroundSpeed/InputPRF)/
		      (Wavelength*Wavelength)-1.0 ))-1.0 );
#else
 MaxCurvRngBins = (FocusRng/RngBinSize)*
    ( sqrt(1.0+ 1.0/( 16.0*AzResFull*AzResFull/
		      (Wavelength*Wavelength)-1.0 ))-1.0 );
#endif


  /* Check FFT size large enough at far swath*/

  MasterRefTimePts = (Int4B)((FocusRng*InputPRF/NomGroundSpeed)*
		        pow(4.0*AzResFull*AzResFull/
		        (Wavelength*Wavelength)-0.25,-0.5)+0.5);
  if (MasterRefTimePts+AzPtsToProcess>FFTSize)
    fprintf(msg,
     "\nWARNING: FFT size (%ld pts) < data+reference (%ld pts) at far swath\n",
        FFTSize,MasterRefTimePts+AzPtsToProcess);

  /* Check focus updates in range */
  if (RngFocSegments > RngBinsToProcess)
    { 
    fprintf(msg,
      "ERROR - Rng focus segments %ld out of range (max rng bins %ld)!\n",
            RngFocSegments,RngBinsToProcess); 
    exit(1);
    }
  else if (RngFocSegments < 1)
    { RngFocSegments  = RngBinsToProcess; }    

  /* Messages */
   fprintf(msg,"Input / output files           : %s / %s\n",
	  InputFileName,OutputFileName);
   fprintf(msg,"Appending output file          : %c\n",
           AppendExistOutFileFlg);
   fprintf(msg,"Az res / look(s) / bandwidth   : %.3f m / %ld / %.3f Hz\n",
           NomAzRes,NumLooks,TotalDopBw);
   fprintf(msg,"Master ref func (time domain)  : %ld pts (far swath)\n",
	   MasterRefTimePts);
   fprintf(msg,"FFT size                       : %ld pts\n",FFTSize);
   fprintf(msg,"Look size in freq domain       : %ld pts\n",LookFreqPts);
   fprintf(msg,"Input / output azimuth points  : %ld / %ld\n",
	  AzPtsToProcess,OutputAzPts);
   fprintf(msg,"Input / output az bin sizes    : %.4fm / %.4fm\n",
	  (NomGroundSpeed/InputPRF),
           (NomGroundSpeed*InvFFTSizeReduc*PostSumRatio/InputPRF));
   fprintf(msg,
      "Range bins to process          : %ld from %ld (bins %ld - %ld)\n",
	 RngBinsToProcess,InputFileRngBins,StartProcessRngBin,
                StartProcessRngBin+RngBinsToProcess-1);
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

#if R_D_OUTPUT
   fprintf(msg,"Producing range-Doppler domain output files.\n");
#endif



   
  /******************************************************************
  * Open output and input files and move input file pointer to start*
  ******************************************************************/

#if R_D_OUTPUT
   if ( (RDFile = fopen (R_D_FILE_NAME, "wb") ) == NULL )
      { fprintf(msg,"ERROR: Output R-D file not opened!\n");exit(1); }
   if ( (RDCorrFile = fopen (R_D_CORR_FILE_NAME, "wb") ) == NULL )
      { fprintf(msg,"ERROR: Output corrected R-D file not opened!\n");exit(1); }
#elif R_D_COMP
   if (R_D_COMP_FILE_DATA_POINTS > FFTSize)
     { fprintf(msg,"ERROR: RD_COMP_DATA_POINTS > FFTSize!\n");exit(1); }   

   if ( (RDCompFile = fopen (R_D_COMP_FILE_NAME,"rt") ) == NULL )
      { fprintf(msg,"ERROR: Input R-D comp file not opened!\n");exit(1); }   

   /*   if ( (RDCompOutputFile = fopen (R_D_COMP_OUTPUT_FILE_NAME,"wt") ) == NULL )
      { fprintf(msg,"ERROR: R-D comp output file not opened!\n");exit(1); }   */


   /* Allocate mem */
   RDCompArray = (float *)malloc(sizeof(float)*FFTSize);

   if (RDCompArray==NULL)
     { fprintf(msg,"ERROR: In RDCompArray allocation!\n");exit(1); }  

   /* Read in data, ensuring matched to FFTSize */
   indx=0;
   i = (Int4B)(FFTSize/(Int4B)(R_D_COMP_FILE_DATA_POINTS));
   while (indx<FFTSize)
     {
     fscanf(RDCompFile,"%f\n",&RDCompArray[indx++]);
     for (k=0;k<i-1;k++)  
       {
       RDCompArray[indx] = RDCompArray[indx-1]; /* Fill out */
       indx++;
       }
     }


    /* Convert data to inverse */
    for (i=0;i<FFTSize;i++)
      {
      if (RDCompArray[i] > 1.0e-12) 
        { RDCompArray[i] = 1.0/RDCompArray[i]; }
      else 
        { RDCompArray[i]=0.0; }  /* also removes any stray small -ves */
      }

    /* Clean up */ 
    fclose(RDCompFile);

#endif

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

  if ( (InputFile = fopen (InputFileName, "rb") ) == NULL )
    { printf ("ERROR: Input file not found/opened!\n"); exit(1); }

  /* For no rng curv correction or too close to edge*/
  if ( StartProcessRngBin<RngCurvInterpSizeD2 || RngCurvInterpSizeD2==0)
    {
    fseek (InputFile,StartProcessRngBin*InputFileAzPts*2 +
	     StartProcessAzPt*2,0);
    }
  else  /* where the extra near-rng rng curv interp bins can be used */
    {
    fseek(InputFile,(StartProcessRngBin-RngCurvInterpSizeD2)*
	           InputFileAzPts*2 + StartProcessAzPt*2,0);
    }

  /********************************
  * ALLOCATE SPACE FOR ARRAYS *
  *********************************/

   InputIQ = (double *)malloc(sizeof(double)*FFTSizeX2);
   CorrectedBinFreq = (double *)malloc(sizeof(double)*FFTSizeX2);
   FocusedAzLine = (double *)malloc(sizeof(double)*InvFFTSizeX2);
   if (PostSumRatio!=1)
     DetectedAzLine = (double*)malloc(sizeof(double)*AzRealPtsPrePostSum);
   PostSummedAzLine = (double*)malloc(sizeof(double)*OutputAzRealPts);
   IQCharInput=(unsigned char *)malloc(sizeof(unsigned char)*AzPtsToProcessX2);
   MasterAzRef=(double *)malloc(sizeof(double)*FFTSizeX2);
   ShiftedMasterAzRef=(double *)malloc(sizeof(double)*FFTSizeX2);
   LookRefBlock=(double *)malloc(sizeof(double)*NumLooks*LookFreqPts*2);
   LookRefFunc=(double **)malloc(sizeof(double *)*NumLooks);

#if R_D_OUTPUT
   RDCorrArray = (float *) malloc(sizeof(float)*FFTSizeX2);
   RDArray = (float *) malloc(sizeof(float)*FFTSizeX2);
   if (RDCorrArray==NULL || RDArray==NULL)
     {
     fprintf(msg,"ERROR - in memory allocation of R-D arrays\n");
     exit(1);
     }

#endif

   for (i=0;i<NumLooks;i++)
      LookRefFunc[i]= &LookRefBlock[i*LookFreqPts*2];

   if ( InputIQ==NULL || CorrectedBinFreq==NULL || FocusedAzLine==NULL
        || IQCharInput==NULL ||	MasterAzRef==NULL
	|| ShiftedMasterAzRef==NULL || LookRefBlock==NULL || LookRefFunc==NULL
        || PostSummedAzLine==NULL )
     {
     fprintf(msg,"ERROR - in memory allocation of processing arrays\n");
     exit(1);
     }

   if (PostSumRatio!=1)
     if(DetectedAzLine==NULL)
       {
       fprintf(msg,"ERROR - in memory allocation of DetectedAzLine array\n");
       exit(0);
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
     {fprintf(msg,"ERROR - in memory allocation of output arrays\n");exit(0);}  
#if REORDER_IFFT
   DechirpedLook = (double *)malloc(sizeof(double)*LookFreqPts*2);
   if(DechirpedLook==NULL)
     {
     fprintf(msg,"ERROR - in mem allocation of DechirpedLook array\n");
     exit(0);
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
         if (NextInputRngBin < InputFileRngBins)
          { fread (IQCharInput,1,AzPtsToProcessX2,InputFile);
            for (i=0;i<AzPtsToProcessX2;i++) InputIQ[i] =
		   (double)IQCharInput[i]-(double)InputDCOffset;}
         else
          { FarRngPad++;}

         /* Increment file pointer, if valid */
         if (NextInputRngBin < InputFileRngBins-1)
	   fseek(InputFile,SkipInputBytes,1);

         /* Increment rng check */
         NextInputRngBin++;

         /* FFT data */
         dummyaddress = InputIQ - 1;
         CFFT(dummyaddress,FFTSize,1);

         /* Shift to zero centered and store in RngCurvedBatchFreq */
         for (i=0; i<FFTSize;i++)
            RngCurvedBatchFreq[k][i] = (float)InputIQ[i+FFTSize];
         for (i=FFTSize; i<FFTSizeX2;i++)
            RngCurvedBatchFreq[k][i] = (float)InputIQ[i-FFTSize];

#if R_D_OUTPUT
         if (NextInputRngBin <= InputFileRngBins)   
           { 
           for (i=0; i<FFTSizeX2;i++)
               RDArray[i] = RngCurvedBatchFreq[k][i];
           fwrite(RDArray,sizeof(float),FFTSizeX2,RDFile);
           }
#endif

        } /* End k loop */

     /* Misc for initial rng curv straightening */
     RngCurvCalc = Wavelength*InputPRF/(2.0*NomGroundSpeed*FFTSize);
     RngCurvCalc *= RngCurvCalc;        /* Square RngCurvCalc */
     NextRngBinToStraighten = StartProcessRngBin;
     RngCurvRefRng = StartProcessRng;
     
     /* STRAIGHTEN RngCurvBatchSize BINS WITH INTERPOLATION
            (straighten only the bandwidth used) */
     RngCurvStraighten(A2DFreq,
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
		      );

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
      if (MasterRefTimePts+AzPtsToProcess>FFTSize)
	fprintf(msg,
	   "\nWARNING: FFT size (%ld pts) < data+reference (%ld pts)\n",
                  FFTSize,MasterRefTimePts+AzPtsToProcess);
      
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
      CFFT(dummyaddress,FFTSize,1);

      /* Shift MasterAzRef to zero centered */
      for (i=0; i<FFTSize;i++)
         ShiftedMasterAzRef[i] = MasterAzRef[i+FFTSize];
      for (i=FFTSize; i<FFTSizeX2;i++)
         ShiftedMasterAzRef[i] = MasterAzRef[i-FFTSize];

      /* Allocate array space for freq window function and calc window */
      if ( (Window = (double *)malloc(sizeof(double)*LookFreqPts)) == NULL )
       {
       fprintf (msg,"ERROR - in Window array memory allocation\n");
       exit(0);
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
         fread (IQCharInput,1,AzPtsToProcessX2,InputFile);
	 for (i=0;i<AzPtsToProcessX2;i++)
	   InputIQ[i] = (double)IQCharInput[i] - (double)InputDCOffset;
         if ( FocusSegment!=RngFocSegments-1 || bin!=RngBinEndFocSeg )
	   fseek(InputFile,SkipInputBytes,1);

	 /* Fourier transform azimuth line */
         dummyaddress = InputIQ - 1;
         CFFT(dummyaddress,FFTSize,1);

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

#if R_D_OUTPUT
            for (i=0; i<FFTSizeX2;i++)
               RDCorrArray[i] = (float)CorrectedBinFreq[i];
            fwrite(RDCorrArray,sizeof(float),FFTSizeX2,RDCorrFile);
#endif

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
               if (NextInputRngBin < InputFileRngBins)
                 {
		 fread (IQCharInput,1,AzPtsToProcessX2,InputFile);
                 for (i=0;i<AzPtsToProcessX2;i++)
                    InputIQ[i] = (double)IQCharInput[i]-
			             (double)InputDCOffset;
		 }
               else
                 { FarRngPad++;}

               /* Increment file pointer if valid */
               if (NextInputRngBin < InputFileRngBins-1)
		 fseek(InputFile,SkipInputBytes,1);

               /* Increment rng check */
               NextInputRngBin++;

               /* FFT data */
               dummyaddress = InputIQ - 1;
               CFFT(dummyaddress,FFTSize,1);

               /* Shift to zero centered and store in RngCurvedBatchFreq */
               indx = CurvArrayRngBins - RngCurvBatchSize + k;
               for (i=0; i<FFTSize;i++)
                  RngCurvedBatchFreq[indx][i] = (float)InputIQ[i+FFTSize];
               for (i=FFTSize; i<FFTSizeX2;i++)
                  RngCurvedBatchFreq[indx][i] = (float)InputIQ[i-FFTSize];

#if R_D_OUTPUT
         if (NextInputRngBin <= InputFileRngBins)   
           { 
           for (i=0; i<FFTSizeX2;i++)
               RDArray[i] = RngCurvedBatchFreq[k][i];
           fwrite(RDArray,sizeof(float),FFTSizeX2,RDFile);
           }
#endif


               }  /* End k loop */

            /* STRAIGHTEN RngCurvBatchSize BINS WITH INTERPOLATION
                         (NextRngBinToStraighten and RngCurvRefRng correct from
	           last interpolation)*/
            RngCurvStraighten(A2DFreq,
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
		      );

            /* Get rng curv straightened az line */
            for (i=0; i<FFTSizeX2;i++)
              CorrectedBinFreq[i] =
		 (double)CorrectedBatchFreq[StraightenedRngBinsUsed][i];
            StraightenedRngBinsUsed++;

#if R_D_OUTPUT
            for (i=0; i<FFTSizeX2;i++)
               RDCorrArray[i] = (float)CorrectedBinFreq[i];
            fwrite(RDCorrArray,sizeof(float),FFTSizeX2,RDCorrFile);
#endif

           }  /* end else no more straightened bins available */

        }  /* end else perform rng curv correction */


#if !R_D_OUTPUT

#if R_D_COMP
  indx = 0;
  for (i=0;i<FFTSize;i++)    /* multiply by R-D spread ampl. compensation */
    {
    tmp = (double)RDCompArray[i];
    CorrectedBinFreq[indx] = CorrectedBinFreq[indx]*tmp;
    indx++;
    CorrectedBinFreq[indx] = CorrectedBinFreq[indx]*tmp;
    /*    if (bin==22) 
      fprintf(RDCompOutputFile,"%15.8e\n",
             CorrectedBinFreq[indx]*CorrectedBinFreq[indx]+ 
             CorrectedBinFreq[indx-1]*CorrectedBinFreq[indx-1]); */
    indx++;
    }

#endif


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
        CFFT(dummyaddress,InvFFTSize,-1);

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

#endif

      /* Report progress to log */
      if (msg==stdout)
        { if (modf((double)bin/ScreenUpdateRate,&ipart)==0)   
            fprintf(msg,
		"Focus rng / output rng bin     : %.4f m / %ld\r",
		    FocusRng,bin); }
      
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
    fprintf (msg,"Max at output az /rng pixel    : %ld / %ld\n",
             MaxValueAzPixel,MaxValueInputRngBin);  
    }

  if (DetectMethod == 3)
    fprintf(msg,
        "Zeros set to min dB val. (%d): %ld\n",
         DB_VALUE_FOR_ZERO,dBMinCounter); 
  
  /* Get end time */
  EndTime = time(NULL);
  fprintf (msg,"AZ COMPRESSION DONE - in %ld secs (%.2f min)\n",
	   EndTime-StartTime,(double)(EndTime-StartTime)/60.0);

  
  /******************************
  * Free memory and close files *
  ******************************/
 
  if (RngCurvInterpSize != 0)
    {
     free (CorrectedBatchBlock);
     free (CorrectedBatchFreq);
     free (RngCurvedBlock);
     free (RngCurvedBatchFreq);
    }

  free (InputIQ);
  free (MasterAzRef);
  free (IQCharInput);
  free (LookRefFunc);
  free (LookRefBlock);
  free (FocusedAzLine);
  free (ShiftedMasterAzRef);
  if (PostSumRatio!=1) free (DetectedAzLine);
  free (PostSummedAzLine);

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

#if R_D_OUTPUT
  free(RDArray);
  free(RDCorrArray);
  fclose(RDFile);
  fclose(RDCorrFile);
#elif R_D_COMP
  free(RDCompArray);
  /* fclose(RDCompOutputFile);*/
#endif


  return (1);
}  /* End function g2azc */


/************************************************************************/
/* FUNCTION : WriteAzComTmplCmdFile() */
Int2B WriteAzComTmplCmdFile(char CmdFileName[])
{
  FILE *OutFile; 
  if ( (OutFile = fopen (CmdFileName, "wt") ) == NULL ) return(0);

  fprintf(OutFile,"Command file for Azcom (SAR Azimuth Compression)\n");
  fprintf(OutFile,"$ProgramVersion (jmh) => %s\n",PROG_VERSION);
  fprintf(OutFile,"------------------------------------------------\n\n");
  fprintf(OutFile,"$ScreenUpdateRate                    => 10\n");
  fprintf(OutFile,"$InputStartSampleDelay [sec]         => 9.0e-05\n");
  fprintf(OutFile,"$CarrierFreq [Hz]                    => 1.41e+08\n");
  fprintf(OutFile,"$InputPRF [Hz]                       => 500.0\n");
  fprintf(OutFile,"$NomGroundSpeed [m/s]                => 250.0\n");
  fprintf(OutFile,"$InputFileAzPts                      => 1001\n");
  fprintf(OutFile,"$StartProcessAzPt                    => 0\n");
  fprintf(OutFile,"$AzPtsToProcess                      => 1001\n");
  fprintf(OutFile,"$InputFileRngBins                    => 2048\n");
  fprintf(OutFile,"$StartProcessRngBin                  => 0\n");
  fprintf(OutFile,"$RngBinsToProcess [see note]         => 2048\n");
  fprintf(OutFile,"$InputDCOffset                       => 127\n");
  fprintf(OutFile,"$FFTSize [power of 2]                => 4096\n");
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
  fprintf(OutFile,"$OutputDataType [see note]           => 0\n");
  fprintf(OutFile,"$Scale [see note]                    => 2.8e-07\n");
  fprintf(OutFile,"$ReportMax [1/0]                     => 1\n");

  fprintf(OutFile,"\nNotes:\n");
  fprintf(OutFile,"------\n\n");
  fprintf(OutFile,"Input data must be unsigned char I and Q\n\n");
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
  fprintf(OutFile,"Scale : 0 - auto (not implemented yet)\n");
  fprintf(OutFile,"      : 1 - none\n");
  fprintf(OutFile,"      : other (double)\n\n");

  fclose(OutFile);
  return(1);
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


#if R_D_OUTPUT
         for (k=0; k<FFTSize/2; k++) /* Correct freq pairs for whole FFTSize*/
#else
	 for (k=0; k<MasterRefFreqPtsD2; k++) /* Correct freq pairs */
#endif
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

     return(1);
}  /* end RngCurvStraighten function */

