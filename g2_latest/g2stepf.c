/*===================================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2stepf.c
CODE CONTROLLER: Richard Lord

DESCRIPTION: 
Part of G2 SAR processor. Called after range compression stage. 
Processing of stepped-frequency waveforms, by reconstruction
of target reflectivity spectrum.

VERSION/AUTHOR/DATE : 0.1 / Richard Lord / 2000-03-31
COMMENTS: 
Initial version.

2000-03-30 (jmh) - renamed executable from "g2stepf" to "stepf", 
removed the !'s from the SeekP calls, and changed the main routine's
return value from 1 to 0 for successful completion. Added a printed
line above the program name, for neat display in G2 processor.

VERSION/AUTHOR/DATE :
COMMENTS:

=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"
#include "sf_procs.h"
#include "g2geth.h"


/* Misc defns limited to this file */
#define G2STEPF_VERSION "0.1"
#define G2STEPF_LOG_VERSION "0.1"

/* Function prototypes */
int WriteStepfTmplCmdFile(char CmdFileName[]);


/***************************************/
/* MAIN STEP-FREQUENCY PROCESSING PROG */
/***************************************/

int main (int argc,char *argv[])
{
  FILE *infile=NULL,
       *outfile=NULL,
       *cmdfp=NULL,
       *logfile=NULL,
       *userfile=NULL,
       *msg=NULL;
  
  char InputFileName[STRING_SPACE],
       OutputFileName[STRING_SPACE],
       CmdFileVersion[STRING_SPACE],
       LogFileName[STRING_SPACE],
       StepFreqUserFile[STRING_SPACE],
       StepFreqProcMode[STRING_SPACE];

  Int4B RngBinsToProcess,
        PreSummedPulsesToUse,
        RngBinsToProcessX2,
        NewRngBinsToProcess,
        NewRngBinsToProcessX2,
        RngFFTSize,
        RngFFTSizeX2,
        AzmLinesToProcess,
        newlength,
        NumberOfFreqSteps,
        Errors,
        RngComRefFuncPhaseSign,
        RangeCompressInputData,
        upsample,
        j,step,azmbin,                 /* Loop counters */
        InputDataType;                 /* 0 - unsigned char (1 byte),
                                          3 - float (4 bytes) */ 


  double *rline=NULL,
         *newrline=NULL,
         *Hsys=NULL,
         *fcenter_n=NULL,
         NarrowChirpBandwidth,
         TotalBandwidth,
         NarrowPulseLen,
         InputA2DFreq,
         OutputA2DFreq,
         InputStartSampleDelay,
         RngComWinConstTime,
         StartCentreFrequency,
         StepSize,
         fcenter,
         fcenter_min,
         fcenter_max,
         fshift,
         InputDCOffsetI,
         InputDCOffsetQ;


  /* Miscellaneous */
  RangeCompressInputData = 0;    /* 1=yes, 0=no */
  upsample = 8;
  Errors = 0;

  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n-----------\n");  
  fprintf(msg,"Prog: STEPF (Ver %s)\n",G2STEPF_VERSION);
  fprintf(msg,"Code: R.T. Lord (Copyright UCT 2000)\n");

  if (argc != 2) {
    fprintf(msg,"Stepped-Frequency Processing for SAR data.\n\n");
    fprintf(msg,"USAGE: stepf [cmd file]\n");
    fprintf(msg,"To see command file structure, type `stepf -tmpl'\n"); 
    fprintf(msg,"(generates a template command file `stepf.tmp')\n\n");  
    exit(1); 
  }

  fprintf(msg,"\n"); 

  /* Write template command file, if requested */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0) {
     if (WriteStepfTmplCmdFile("g2stepf.tmp")==0) {
       fprintf(msg,"Template command file `g2stepf.tmp' written!\n\n");
     }
     exit(1);
  }

  /* Open spec file */
  if ((cmdfp = fopen(argv[1],"r")) == NULL) {
    fprintf(msg,"ERROR - command file %s not found/opened!\n\n",argv[1]);
    exit(1);
  }  

  /* Read processing specs */
  if (SeekP(cmdfp,"$CmdFileVersion","=>",-1,0)) { 
    fprintf(msg,"WARNING - in parsing command file (prog version not found)!\n\n"); 
  } 
  else {
    ReadStr(cmdfp,CmdFileVersion,STRING_SPACE); 
  }


  if (SeekP(cmdfp,"$LogFileName","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,LogFileName,STRING_SPACE);  

  /* Create logfile */
  if (strcmp(LogFileName,"null")!=0) {
    if ( (logfile = fopen (LogFileName, "w") ) == NULL )
      fprintf(msg,"WARNING - Log file %s not created!\n\n",LogFileName); 
    else {
      msg = logfile;
      fprintf(msg,"G2 Stepped-Frequency Processor (Ver %s) Log File\n",G2STEPF_VERSION);
      fprintf(msg,"Code: R.T. Lord - UCT Radar Remote Sensing Group 2000\n"); 
      fprintf(msg,"Log File Version : %s\n\n",G2STEPF_LOG_VERSION); 
    }
  }


  fprintf(msg,"Input Parameters:\n");
  if (SeekP(cmdfp,"$InputFileName","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,InputFileName,STRING_SPACE);  
  fprintf(msg,"InputFilename                : %s\n", InputFileName);
  if (SeekP(cmdfp,"$OutputFileName","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,OutputFileName,STRING_SPACE);
  fprintf(msg,"OutputFileName               : %s\n", OutputFileName);
  if (SeekP(cmdfp,"$InputDataType","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%ld", &InputDataType);
  if (SeekP(cmdfp,"$PreSummedPulsesToUse","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&PreSummedPulsesToUse);
  fprintf(msg,"PreSummedPulsesToUse         : %ld\n", PreSummedPulsesToUse);
  if (SeekP(cmdfp,"$RngBinsToProcess","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&RngBinsToProcess);
  fprintf(msg,"RngBinsToProcess             : %ld\n", RngBinsToProcess);
  if (SeekP(cmdfp,"$InputDCOffsetI","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &InputDCOffsetI);
  if (SeekP(cmdfp,"$InputDCOffsetQ","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &InputDCOffsetQ);
  if (SeekP(cmdfp,"$NarrowChirpBandwidth","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &NarrowChirpBandwidth);
  fprintf(msg,"NarrowChirpBandwidth (Hz)    : %.5e\n", NarrowChirpBandwidth);
  if (SeekP(cmdfp,"$NarrowPulseLen","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &NarrowPulseLen); 
  fprintf(msg,"NarrowPulseLen (sec)         : %.5e\n", NarrowPulseLen);
  if (SeekP(cmdfp,"$InputA2DFreq","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &InputA2DFreq);
  fprintf(msg,"InputA2DFreq (Hz)            : %.5e\n", InputA2DFreq);
  if (SeekP(cmdfp,"$InputStartSampleDelay","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &InputStartSampleDelay);
  fprintf(msg,"InputStartSampleDelay (sec)  : %.5e\n", InputStartSampleDelay);
  if (SeekP(cmdfp,"$RngComWinConstTime","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &RngComWinConstTime);
  fprintf(msg,"RngComWinConstTime           : %.5e\n", RngComWinConstTime);
  if (SeekP(cmdfp,"$RngComRefFuncPhaseSign","=>",-1,0)) Errors++; 
  else fscanf(cmdfp, "%ld", &RngComRefFuncPhaseSign);
  fprintf(msg,"RngComRefFuncPhaseSign       : %ld\n", RngComRefFuncPhaseSign);
  if (SeekP(cmdfp,"$StepFreqProcMode","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,StepFreqProcMode,STRING_SPACE);  
  fprintf(msg,"StepFreqProcMode             : %s\n", StepFreqProcMode);
  if (SeekP(cmdfp,"$StartCentreFrequency","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &StartCentreFrequency);
  if (strcmp(StepFreqProcMode,"normal")==0)
    fprintf(msg,"StartCentreFrequency (Hz)    : %.5e\n", StartCentreFrequency);
  if (SeekP(cmdfp,"$StepSize","=>",-1,0)) Errors++; 
  fscanf(cmdfp, "%lf", &StepSize);
  if (strcmp(StepFreqProcMode,"normal")==0)
    fprintf(msg,"StepSize (Hz)                : %.5e\n", StepSize);
  if (SeekP(cmdfp,"$NumberOfFreqSteps","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&NumberOfFreqSteps);
  fprintf(msg,"NumberOfFreqSteps            : %ld\n\n", NumberOfFreqSteps);
  if (SeekP(cmdfp,"$StepFreqUserFile","=>",-1,0)) Errors++; 
  ReadStr(cmdfp,StepFreqUserFile,STRING_SPACE);  
  if (strcmp(StepFreqProcMode,"user")==0)
    fprintf(msg,"StepFreqUserFile             : %s\n\n", StepFreqUserFile);

  fclose (cmdfp);

  if (Errors!=0) {
    fprintf(msg,"ERROR - %ld errors in parsing command file!\n\n",Errors);
    exit(1); 
  }

  /* Check version IDs match */
  if (strcmp(CmdFileVersion,G2STEPF_VERSION)!=0) {
    fprintf(msg,"WARNING - command file version (%s) not same as program (%s)!\n\n",
      CmdFileVersion,G2STEPF_VERSION);
  }

  /* Error checking */
  if (PreSummedPulsesToUse < NumberOfFreqSteps) {
    fprintf(msg,"ERROR - Number of pulses to process (%ld) less than\n", PreSummedPulsesToUse);
    fprintf(msg,"number of frequency steps (%ld)!\n\n", NumberOfFreqSteps);
    exit(1); 
  }
  if (InputDataType != 0 && InputDataType != 3) {
    fprintf(msg,"ERROR - Invalid input data type %ld!\n\n",InputDataType);
    exit(1); 
  }
  if ((strcmp(StepFreqProcMode,"normal")!=0) && 
      (strcmp(StepFreqProcMode,"user")!=0)) {
    fprintf(msg,"ERROR - Invalid step-frequency processing mode: %s!\n\n", StepFreqProcMode);
    exit(1); 
  }


 /*****************************
  * OPEN FILES FOR PROCESSING *
  *****************************/

  /* Open input file */
  if ( (infile = fopen (InputFileName, "rb") ) == NULL ) {
    fprintf(msg,"ERROR - Input file %s not opened!\n\n",InputFileName); 
    exit(1); 
  }

  /* Open output file */
  if ( (outfile = fopen (OutputFileName, "wb") ) == NULL ) {
    fprintf(msg,"ERROR - Output file %s not opened!\n\n", OutputFileName);
    exit(1); 
  }

  /* Open user file */
  if (strcmp(StepFreqProcMode,"user")==0) {
    if ( (userfile = fopen (StepFreqUserFile, "r") ) == NULL ) {
      fprintf(msg,"ERROR - User file %s not opened!\n\n", StepFreqUserFile);
      exit(1); 
    }
  }


  /* Miscellaneous */
  fcenter_min = 1e15;
  fcenter_max = 0.0;
  RngBinsToProcessX2 = 2*RngBinsToProcess;

  /* Find number of frequency steps */
  if (strcmp(StepFreqProcMode,"user")==0) {
    fscanf(userfile, "%ld", &NumberOfFreqSteps);
  }
  if (NumberOfFreqSteps > 1000) {
    fprintf(msg,"ERROR - Number of frequency steps (%ld) greater than 1000.\n", NumberOfFreqSteps);
    exit(1); 
  }

  /* Find centre frequencies of narrowband pulses */
  fcenter_n = (double *)malloc(sizeof(double)*NumberOfFreqSteps);
  if ( fcenter_n==NULL ) {
    fprintf(msg,"ERROR - in array memory allocation!\n");
    exit (1); 
  }
  if (strcmp(StepFreqProcMode,"user")==0) {
    for (j=0; j<NumberOfFreqSteps; j++) {
      fscanf(userfile, "%lf", &fcenter_n[j]);
      if (fcenter_n[j] < fcenter_min) fcenter_min=fcenter_n[j];
      if (fcenter_n[j] > fcenter_max) fcenter_max=fcenter_n[j];
    }
  }
  else {
    for (j=0; j<NumberOfFreqSteps; j++) fcenter_n[j] = StartCentreFrequency + (j*StepSize);
    fcenter_min = fcenter_n[0];
    fcenter_max = fcenter_n[NumberOfFreqSteps-1];
  }

  /* Total radar bandwidth */
  TotalBandwidth = (fcenter_max+(NarrowChirpBandwidth/2.0))
    -(fcenter_min-(NarrowChirpBandwidth/2.0));
  fprintf(msg,"Total radar bandwidth after reconstruction (Hz): %.5e\n",TotalBandwidth);

  /* Centre frequency of reconstructed spectrum */
  fcenter = ((fcenter_max+(NarrowChirpBandwidth/2.0))+(fcenter_min-(NarrowChirpBandwidth/2.0)))/2.0;
  fprintf(msg,"Centre frequency of reconstructed spectrum (Hz): %.5e\n",fcenter);

  /* Number of wideband output range lines */
  AzmLinesToProcess = floor(PreSummedPulsesToUse / NumberOfFreqSteps);
  fprintf(msg,"Number of range lines of output file           : %ld\n",AzmLinesToProcess);

  /* Number of wideband output range bins */
  RngFFTSize =  pow(2.0,ceil(log((double)RngBinsToProcess) / log(2.0)));
  RngFFTSizeX2 = 2*RngFFTSize;
  newlength = ceil(((fcenter-fcenter_min)/InputA2DFreq)*RngFFTSize)
    +ceil(((fcenter_max-fcenter)/InputA2DFreq)*RngFFTSize) + RngFFTSize;
  NewRngBinsToProcess =  pow(2.0,ceil(log((double)newlength) / log(2.0)));
  fprintf(msg,"Number of range bins of output file            : %ld\n",NewRngBinsToProcess);
  NewRngBinsToProcessX2 = 2*NewRngBinsToProcess;

  /* New output A2D frequency */
  OutputA2DFreq = ((double)NewRngBinsToProcess/(double)RngFFTSize)*InputA2DFreq;
  fprintf(msg,"Output A2D frequency (Hz)                      : %.5e\n\n",OutputA2DFreq);


 /*****************************
  * ALLOCATE SPACE FOR ARRAYS *
  *****************************/

  rline = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  newrline = (double *)malloc(sizeof(double)*NewRngBinsToProcessX2);
  Hsys = (double *)malloc(sizeof(double)*NewRngBinsToProcessX2);
  if ( rline==NULL || newrline==NULL || Hsys==NULL ) Errors++;
  if (Errors != 0) {
    fprintf(msg,"ERROR - in array memory allocation!\n"); 
    exit (1);
  }


 /**************************
  * GET COMPRESSION FILTER *
  **************************/

  G2GetHsys(
       msg,
       infile,
       RngBinsToProcess,
       RngFFTSize,
       NewRngBinsToProcess,
       NumberOfFreqSteps,
       RngComRefFuncPhaseSign,
       upsample,
       fcenter,
       NarrowChirpBandwidth,
       InputA2DFreq,
       NarrowPulseLen,
       InputStartSampleDelay,
       RngComWinConstTime,
       fcenter_n,
       Hsys
       );


 /**********************************
  * FOR EVERY BURST OF FREQUENCIES *
  **********************************/

  for (azmbin=0; azmbin<AzmLinesToProcess; azmbin++)
  {

    /* Initialise new range line array */
    for (j=0; j<NewRngBinsToProcessX2; j++) newrline[j] = 0.0;

    for (step=0; step<NumberOfFreqSteps; step++)
    {
      fshift = (double)(((fcenter_n[step]-fcenter)/InputA2DFreq)*(double)RngFFTSize);


      /* Get range line from input file */
      G2GetRangeLine(
                msg,
                infile,
                RngBinsToProcess,
                RngFFTSize,
                InputDataType,
                InputDCOffsetI,
                InputDCOffsetQ,
                rline
                );
                

      /* Take FFT of range line, circular shift, phase compensation */
      G2FFTRangeLine(
                msg,
                RngFFTSize,
                InputA2DFreq,
                InputStartSampleDelay,
                rline
                );


      /* Range-compress range line if required */
      if (RangeCompressInputData == 1) {
        G2RangeCompress(
                msg,
                RngFFTSize,
                RngComRefFuncPhaseSign,
                NarrowChirpBandwidth,
                InputA2DFreq,
                NarrowPulseLen,
                rline
                );
      }


      /* Coherently add shifted subspectra */
      G2AddSubSpectrum(
                msg,
                RngFFTSize,
                NewRngBinsToProcess,
                step,
                fshift,
                rline,
                newrline
                );

    }


    /* Apply phase compensation */
    G2PhaseComp(
         msg,
         RngFFTSize,
         NewRngBinsToProcess,
         InputA2DFreq,
         InputStartSampleDelay,
         newrline
         );


    /* Multiply with compression filter and take IFFT */
    G2GetRangeProfile(
         msg,
         NewRngBinsToProcess,
         newrline,
         Hsys
         );


    /* Write range line to file */
    G2PutRangeLine( 
              msg,
              outfile,
              NewRngBinsToProcess,
              newrline
              );

  }


  /***********************************
  * FREE ARRAY SPACE AND CLOSE FILES *
  ************************************/

  free (fcenter_n);
  free (rline);
  free (newrline);
  free (Hsys);
  if (infile != NULL) fclose (infile);
  if (outfile != NULL) fclose (outfile);
  if (logfile != NULL) fclose (logfile);
  if (userfile != NULL) fclose (userfile);

  fprintf(stdout,"Stepped-frequency processing done!\n\n");

  return(0);

}  /* End G2stepf procedure */



/***********************************************************************/
/* FUNCTION : WriteStepfTmplCmdFile() */
int WriteStepfTmplCmdFile(char CmdFileName[])
{
  FILE *Out; 
  if ( (Out = fopen (CmdFileName, "w") ) == NULL ) return(1);

  fprintf(Out,"Command file for G2stepf (SAR stepped-frequency processing)\n");
  fprintf(Out,"$CmdFileVersion (rtl) => %s\n",G2STEPF_VERSION);
  fprintf(Out,"-----------------------------------------------------------\n\n");

  fprintf(Out,"--General --\n");
  fprintf(Out,"$LogFileName  ('null' for none)        => null\n");
  fprintf(Out,"$InputFileName                         => sim.bin\n");
  fprintf(Out,"$OutputFileName                        => test.bin\n");
  fprintf(Out,"$InputDataType [see note]              => 3\n");
  fprintf(Out,"$PreSummedPulsesToUse [see note]       => 4\n");
  fprintf(Out,"$RngBinsToProcess                      => 320\n");
  fprintf(Out,"$InputDCOffsetI                        => 0.0\n");
  fprintf(Out,"$InputDCOffsetQ                        => 0.0\n");
  fprintf(Out,"$NarrowChirpBandwidth [Hz]             => 12.0e+06\n");
  fprintf(Out,"$NarrowPulseLen [sec]                  => 10.0e-06\n");
  fprintf(Out,"$InputA2DFreq [Hz]                     => 24.0e+06\n");
  fprintf(Out,"$InputStartSampleDelay [sec]           => 2.669e-05\n");
  fprintf(Out,"$RngComRefFuncPhaseSign                => -1\n");
  fprintf(Out,"$RngComWinConstTime                    => 0.08\n");
  fprintf(Out,"$NumberOfFreqSteps                     => 4\n");
  fprintf(Out,"$StepFreqProcMode (normal/user)        => normal\n\n");

  fprintf(Out,"--Stepped-Frequency (normal)--\n");
  fprintf(Out,"$StartCentreFrequency [Hz]             => 124.8e+06\n");
  fprintf(Out,"$StepSize [Hz]                         => 10.8e+06\n\n");

  fprintf(Out,"--Stepped-Frequency (user)--\n");
  fprintf(Out,"$StepFreqUserFile [see note]           => user.txt\n\n\n");


  fprintf(Out,"Notes:\n");
  fprintf(Out,"------\n\n");

  fprintf(Out,"InputDataType : 0 - unsigned char IQ (2*1 bytes per point)\n");
  fprintf(Out,"              : 3 - float IQ (2*4 bytes per point)\n\n");

  fprintf(Out,"PreSummedPulsesToUse : These are the un-combined, narrow-bandwidth pulses\n\n");

  fprintf(Out,"StepFreqUserFile : ASCII file, with centre frequency of each step\n");
  fprintf(Out,"                   on a new line, in the order they were transmitted.\n\n");

  fclose(Out);

  return(0);

} /* End WriteStepfTmplCmdFile function */


