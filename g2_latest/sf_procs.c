/*===================================================================
COPYRIGHT: UCT Radar Remote Sensing Group 2000
FILE NAME: sf_procs.c
CODE CONTROLLER: Richard Lord

DESCRIPTION: 
Collection of procedures used by stepped-frequency processor.

List of procedures:
G2AddSubSpectrum
G2FFTRangeLine
G2GetRangeLine
G2GetPointSim
G2PhaseComp
G2PutRangeLine
G2RangeCompress
G2GetRangeProfile



VERSION/AUTHOR/DATE : 0.1 / Richard Lord / 2000-03-31
COMMENTS: 

=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"
#include "sf_procs.h"


/*********************************************************************
* Procedure: G2AddSubSpectrum
* 1) Shifts subspectrum to required location
* 2) Coherently adds subspectrum to reconstructed spectrum
* 3) Perform phase compensation
**********************************************************************/
int G2AddSubSpectrum(
              FILE *msg,
              Int4B RngFFTSize,
              Int4B NewRngBinsToProcess,
              Int4B step,
              double fshift,
              double rline[],
              double newrline[]
              )

{
  Int4B RngFFTSizeX2,
        NewRngBinsToProcessX2,
        firstbin, lastbin,
        i,indx,
        indxi, indxq,
        shift,
        Errors;

  double constant,
         tempval,
         *dummyaddress=NULL,
         *phase=NULL,
         *temp=NULL;


  /* Miscellaneous */
  Errors=0;
  RngFFTSizeX2 = 2*RngFFTSize;
  NewRngBinsToProcessX2 = 2*NewRngBinsToProcess;


  /* Allocate space for arrays */
  phase = (double *)malloc(sizeof(double)*NewRngBinsToProcess);
  temp = (double *)malloc(sizeof(double)*NewRngBinsToProcessX2);
  if ( phase==NULL || temp==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n");
     exit (1);
  }


  /* Copy rline into middle of temp */
  firstbin = NewRngBinsToProcess - RngFFTSize;
  lastbin = firstbin + RngFFTSizeX2;
  indx = 0;
  for (i=0; i<firstbin; i++) {
    temp[i] = 0.0;
  }
  for (i=firstbin; i<lastbin; i++) {
    temp[i] = rline[indx];
    indx++;
  }
  for (i=lastbin; i<NewRngBinsToProcessX2; i++) {
    temp[i] = 0.0;
  }


  /* Create phase ramp */
  constant = (-2.0 * PI) / ((double)NewRngBinsToProcess);
  shift = NewRngBinsToProcess / 2;
  for (i=0; i<NewRngBinsToProcess; i++) {
    phase[i] = (((double)((i+shift)%NewRngBinsToProcess)*constant) + PI) * fshift;
  }


  /* Take FFT */
  dummyaddress = temp - 1;
  CFFT(dummyaddress,NewRngBinsToProcess,-1);


  /* Multiply with phase ramp */
  indxi = 0;
  indxq = 1;
  for (i=0; i<NewRngBinsToProcess; i++) {
    tempval = temp[indxi];
    temp[indxi] = tempval*cos(phase[i]) - temp[indxq]*sin(phase[i]);
    temp[indxq] = tempval*sin(phase[i]) + temp[indxq]*cos(phase[i]);
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }


  /* Take IFFT */
  dummyaddress = temp - 1;
  CFFT(dummyaddress,NewRngBinsToProcess,1);


  /* Add shifted and scaled subspectrum coherently to newrline */
  for (i=0; i<NewRngBinsToProcessX2; i++) {
    newrline[i] += (temp[i] / NewRngBinsToProcess);
  }


  /* Free arrays */
  free (phase);
  free (temp);

  return(0);
}



/*********************************************************************
* Procedure: G2FFTRangeLine
* 1) Takes FFT of range line
* 2) Circular shift
* 3) Phase compensation
**********************************************************************/
int G2FFTRangeLine(
           FILE *msg,
           Int4B RngFFTSize,
           double InputA2DFreq,
           double InputStartSampleDelay,
           double rline[]
           )

{
  double tempval,
         phase_gradient,
         *tempiq=NULL,
         *phase=NULL,
         *dummyaddress=NULL;
         

  Int4B RngFFTSizeX2,
        shift,
        Errors,
        i, indxi, indxq;


  /* Initialise */
  Errors=0;
  RngFFTSizeX2 = 2*RngFFTSize;


  /* Allocate space for arrays */
  tempiq = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  phase = (double *)malloc(sizeof(double)*RngFFTSize);
  if ( tempiq==NULL || phase==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n"); 
     exit (1); 
  }


  /* Make copy of range line */
  for (i=0; i<RngFFTSizeX2; i++) tempiq[i] = rline[i];


  /* Take FFT */
  dummyaddress = tempiq - 1;
  CFFT(dummyaddress,RngFFTSize,-1);


  /* Circular shift FFT and rescale */
  shift = -(Int4B)RngFFTSize;
  for (i=0; i<RngFFTSizeX2;i++) {
    rline[i] = tempiq[(i-shift)%RngFFTSizeX2] / RngFFTSize;       
  }  


  /* Create phase compensation ramp */
  phase_gradient = -2.0*PI*(InputA2DFreq/(double)(RngFFTSize-1))*InputStartSampleDelay;
  for (i=0; i<RngFFTSize;i++) {
    phase[i] = ((double)i - ((double)RngFFTSize / 2.0)) * phase_gradient;
  }


  /* Multiply with phase ramp */
  indxi = 0;
  indxq = 1;
  for (i=0; i<RngFFTSize; i++) {
    tempval = rline[indxi];
    rline[indxi] = tempval*cos(phase[i]) - rline[indxq]*sin(phase[i]);
    rline[indxq] = tempval*sin(phase[i]) + rline[indxq]*cos(phase[i]);
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }


  /* Free Arrays */
  free (tempiq);
  free (phase);

  return(0);
}



/*********************************************************************
* Procedure: G2GetRangeLine
* 1) Read in range line
* 2) Output array is rline 
**********************************************************************/
int G2GetRangeLine( 
           FILE *msg,
           FILE *infile,
           Int4B RngBinsToProcess,
           Int4B RngFFTSize,
           Int4B InputDataType,
           double InputDCOffsetI,
           double InputDCOffsetQ,
           double rline[]
           )

{
  Int4B RngBinsToProcessX2,
        RngFFTSizeX2,
        Errors,
        j, indxi, indxq;

  unsigned char *IQInChar=NULL;

  float *IQInFloat=NULL;

  /* Initialise */
  Errors=0;
  RngBinsToProcessX2 = 2*RngBinsToProcess;
  RngFFTSizeX2 = 2*RngFFTSize;


  /* Allocate space for arrays */
  if (InputDataType == 0) { 
    IQInChar  = (unsigned char *)malloc(sizeof(unsigned char)*RngBinsToProcessX2);
    if (IQInChar==NULL) Errors++;
  }
  else if (InputDataType == 3) {
    IQInFloat  = (float *)malloc(sizeof(float)*RngBinsToProcessX2);
    if (IQInFloat==NULL) Errors++;
  }
  if (Errors != 0) {
    fprintf(msg,"ERROR - in array memory allocation!\n"); exit (1);
  }


  /* Read in required range bins */
  if (InputDataType == 0) { 
    if ( fread(IQInChar,sizeof(unsigned char),RngBinsToProcessX2,infile) 
          != RngBinsToProcessX2 ) { 
      fprintf(msg,"ERROR - in read from input file!\n"); exit(1); 
    }
  }
  else if (InputDataType == 3) {
    if ( fread(IQInFloat,sizeof(float),RngBinsToProcessX2,infile)
         != RngBinsToProcessX2) { 
      fprintf(msg,"ERROR - in read from input file!\n"); exit(1); 
    }
  }


  /* Initialise range line array */
  for (j=0; j<RngFFTSizeX2; j++) rline[j] = 0.0;


  /* Read range line */
  indxi = 0;
  indxq = 1;
  for (j=0; j<RngBinsToProcess;j++) {
    if (InputDataType == 0) { 
      rline[indxi] = (double)IQInChar[indxi]-InputDCOffsetI;
      rline[indxq] = (double)IQInChar[indxq]-InputDCOffsetQ;
    }
    else if (InputDataType == 3) { 
      rline[indxi] = (double)IQInFloat[indxi]-InputDCOffsetI;
      rline[indxq] = (double)IQInFloat[indxq]-InputDCOffsetQ;
    }
    indxi++; indxi++;
    indxq++; indxq++;
  }


  /* Free Arrays */
  if (InputDataType == 0) {
    free (IQInChar); 
  }
  else if (InputDataType == 3) {
    free (IQInFloat); 
  }

  return(0);
}



/*********************************************************************
* Procedure: G2GetPointSim
* 1) Gets point target simulation
**********************************************************************/
int G2GetPointSim(
       FILE *msg,
       Int4B RngFFTSize,
       Int4B upsample,
       double fcenter,
       double NarrowChirpBandwidth,
       double InputA2DFreq,
       double NarrowPulseLen,
       double InputStartSampleDelay,
       double RangeTarget,
       double rline[]
       )


{
  Int4B RngFFTSizeX2,
        startbin, endbin,
        numsamples,
        Errors=0,
        indxi, indxq,
        indx2i, indx2q,
        j;

  double constant1,
         constant2,
         RngFFTSizeUpsample,
         binlength,
         binspacing,
         timeshift,
         tstep,
         timeval,
         phase,
         gamma,
         *tempiq=NULL,
         *dummyaddress=NULL;


  /* Miscellaneous */
  RngFFTSizeX2 = 2*RngFFTSize;
  gamma = NarrowChirpBandwidth / NarrowPulseLen;
  binlength = C / (2.0 * InputA2DFreq * (double)upsample);
  constant1 = -4.0 * PI * fcenter * (RangeTarget / C);
  constant2 = PI * gamma;
  tstep = 1.0 / (InputA2DFreq * (double)upsample);
  timeshift = (2.0 * RangeTarget) / C;
  binspacing = InputA2DFreq / (double)RngFFTSize;
  numsamples = (Int4B)((floor(((NarrowChirpBandwidth / binspacing) / 2.0) + 0.5)) * 2);


  /* Allocate space for arrays */
  tempiq = (double *)malloc(sizeof(double)*RngFFTSizeX2*upsample);
  if ( tempiq==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n");
     exit (1);
  }


  /* Initialise arrays */
  for (j=0; j<RngFFTSizeX2; j++) rline[j] = 0.0;
  for (j=0; j<RngFFTSizeX2*upsample; j++) tempiq[j] = 0.0;


  /* Find point target return */
  RngFFTSizeUpsample = (double)(RngFFTSize*upsample);
  startbin = (Int4B)(floor(((RngFFTSizeUpsample - (NarrowPulseLen/tstep)) / 2.0) + 0.5) + 1);
  endbin = (Int4B)(floor((double)startbin + (NarrowPulseLen/tstep) + 0.5));
  for (j=startbin; j<endbin; j++) {
    timeval = (((double)j * tstep) + InputStartSampleDelay) - timeshift;
    phase = constant1 + (constant2 * (timeval * timeval));
    tempiq[j*2] = cos(phase);
    tempiq[(j*2)+1] = sin(phase);
  }


  /* Take FFT */
  dummyaddress = tempiq - 1;
  CFFT(dummyaddress,RngFFTSize*upsample,-1);


  /* Copy relevant portion of FFT to range line */
  indxi=0;
  indxq=1;
  for (j=0; j<=(numsamples/2); j++) {
    rline[indxi] = tempiq[indxi] / RngFFTSizeUpsample;
    rline[indxq] = tempiq[indxq] / RngFFTSizeUpsample;
    indxi++; indxi++;
    indxq++; indxq++;
  }

  indxi=RngFFTSizeX2-numsamples;
  indxq=indxi+1;
  indx2i=(RngFFTSizeX2*upsample)-numsamples;
  indx2q=indx2i+1;
  for (j=(numsamples/2); j<numsamples; j++) {
    rline[indxi] = tempiq[indx2i] / RngFFTSizeUpsample;
    rline[indxq] = tempiq[indx2q] / RngFFTSizeUpsample;
    indxi++; indxi++;
    indxq++; indxq++;
    indx2i++; indx2i++;
    indx2q++; indx2q++;
  }


  /* Take IFFT */
  dummyaddress = rline - 1;
  CFFT(dummyaddress,RngFFTSize,1);


  /* Free Arrays */
  free (tempiq);

  return(0);
}



/*********************************************************************
* Procedure: G2PhaseComp
* 1) Phase compensation
**********************************************************************/
int G2PhaseComp(
           FILE *msg,
           Int4B RngFFTSize,
           Int4B NewRngBinsToProcess,
           double InputA2DFreq,
           double tshift,
           double newrline[]
           )

{
  double tempval,
         phase_gradient,
         *phase=NULL;

  Int4B Errors,
        i, indxi, indxq;


  /* Initialise */
  Errors=0;

  /* Allocate space for arrays */
  phase = (double *)malloc(sizeof(double)*NewRngBinsToProcess);
  if ( phase==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n");
     exit (1);
  }


  /* Create phase compensation ramp */
  phase_gradient = 2.0 * PI * (InputA2DFreq / (double)(RngFFTSize - 1)) * tshift;
  for (i=0; i<NewRngBinsToProcess;i++) {
    phase[i] = ((double)i - (NewRngBinsToProcess / 2.0)) * phase_gradient;
  }


  /* Multiply with phase ramp */
  indxi = 0;
  indxq = 1;
  for (i=0; i<NewRngBinsToProcess; i++) {
    tempval = newrline[indxi];
    newrline[indxi] = tempval*cos(phase[i]) - newrline[indxq]*sin(phase[i]);
    newrline[indxq] = tempval*sin(phase[i]) + newrline[indxq]*cos(phase[i]);
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }


  /* Free Arrays */
  free (phase);

  return(0);
}



/*********************************************************************
* Procedure: G2PutRangeLine
* 1) Writes range line to file
**********************************************************************/
int G2PutRangeLine( 
           FILE *msg,
           FILE *outfile,
           Int4B RngBinsToProcess,
           double rline[]
           )

{
  Int4B RngBinsToProcessX2,
        Errors,
        j;

  float *IQOutFloat=NULL;


  /* Initialise */
  Errors=0;
  RngBinsToProcessX2 = 2*RngBinsToProcess;


  /* Allocate space for arrays */
  IQOutFloat  = (float *)malloc(sizeof(float)*RngBinsToProcessX2);
  if (IQOutFloat==NULL) Errors++;
  if (Errors != 0) {
    fprintf(msg,"ERROR - in array memory allocation!\n"); exit (1);
  }


  for (j=0; j<RngBinsToProcessX2; j++) IQOutFloat[j] = (float)rline[j];

  fwrite(IQOutFloat,sizeof(float),RngBinsToProcessX2,outfile);


  /* Free Arrays */
  free (IQOutFloat); 

  return(0);
}



/*********************************************************************
* Procedure: G2RangeCompress
* 1) Performs range compression
**********************************************************************/
int G2RangeCompress(
       FILE *msg,
       Int4B RngFFTSize,
       Int4B RngComRefFuncPhaseSign,
       double NarrowChirpBandwidth,
       double InputA2DFreq,
       double NarrowPulseLen,
       double rline[]
       )


{
  Int4B RngFFTSizeX2,
        RefShiftX2,
        Errors,
        shift,
        indxi, indxq,
        rnglen,
        j;

  double tstep,
         timeval,
         gamma,
         phase,
         tempval,
         *dummyaddress=NULL,
         *tempiq=NULL,
         *matchfil=NULL;


  /* Miscellaneous */
  Errors=0;
  RngFFTSizeX2 = 2*RngFFTSize;
  gamma = NarrowChirpBandwidth / NarrowPulseLen;
  tstep = 1.0 / InputA2DFreq;
  rnglen = (Int4B)(NarrowPulseLen*InputA2DFreq);


  /* Allocate space for arrays */
  matchfil = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  tempiq = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  if ( tempiq==NULL || matchfil==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n");
     exit (1);
  }


  /* Initialise array */
  for (j=0; j<RngFFTSizeX2; j++) matchfil[j] = 0.0;


  /* Generate reference function */
  indxi = 0;
  indxq = 1;
  for (j=0; j<rnglen; j++) {
    timeval = ( (double)j*(double)rnglen/((double)rnglen-1.0)
              - (double)rnglen/2.0 )/InputA2DFreq;
    phase = (double)RngComRefFuncPhaseSign*PI*gamma*timeval*timeval;
    matchfil[indxi++] = cos(phase);
    matchfil[indxq++] = -sin(phase);
    indxi++;
    indxq++;
  }


  /* Shift time domain func to origin centered */
  RefShiftX2 = -(Int4B)((floor(rnglen/2.0))*2);
  for (j=0; j<RngFFTSizeX2;j++) {
     tempiq[j] = matchfil[(j-RefShiftX2)%RngFFTSizeX2];       
  }  


  /* Take FFT of reference function */
  dummyaddress = tempiq - 1;
  CFFT(dummyaddress,RngFFTSize,-1);


  /* Circular shift FFT */
  shift = -(Int4B)RngFFTSize;
  for (j=0; j<RngFFTSizeX2;j++) {
    matchfil[j] = tempiq[(j-shift)%RngFFTSizeX2];       
  }  


  /* Multiply range line with matched filter */
  indxi = 0;
  indxq = 1;
  for (j=0; j<RngFFTSize;j++) {
    tempval = rline[indxi];
    rline[indxi] = rline[indxi]*matchfil[indxi] - rline[indxq]*matchfil[indxq];
    rline[indxq] = tempval*matchfil[indxq] + rline[indxq]*matchfil[indxi];
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }


  /* Free Arrays */
  free (tempiq);
  free (matchfil);

  return(0);
}



/*********************************************************************
* Procedure: G2GetRangeProfile
* 1) Multiplies reconstructed spectrum with compression filter Hsys
* 2) Performs cyclic shift
* 3) Applies IFFT
**********************************************************************/
int G2GetRangeProfile(
       FILE *msg,
       Int4B NewRngBinsToProcess,
       double newrline[],
       double Hsys[]
       )

{
  Int4B NewRngBinsToProcessX2,
        Errors,
        shift,
        j, indxi, indxq;

  double *tempiq=NULL,
         *dummyaddress=NULL;


  /* Initialise */
  Errors=0;
  NewRngBinsToProcessX2 = 2*NewRngBinsToProcess;


  /* Allocate space for arrays */
  tempiq = (double *)malloc(sizeof(double)*NewRngBinsToProcessX2);
  if ( tempiq==NULL ) Errors++;
  if (Errors != 0) {
     fprintf(msg,"ERROR - in array memory allocation!\n");
     exit (1);
  }


  /* Multiply spectrum with compression filter */
  indxi = 0;
  indxq = 1;
  for (j=0; j<NewRngBinsToProcess;j++) {
    tempiq[indxi] = newrline[indxi]*Hsys[indxi] - newrline[indxq]*Hsys[indxq];
    tempiq[indxq] = newrline[indxi]*Hsys[indxq] + newrline[indxq]*Hsys[indxi];
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }



  /* Circular shift FFT */
  shift = -NewRngBinsToProcess;
  for (j=0; j<NewRngBinsToProcessX2;j++) {
     newrline[j] = tempiq[(j-shift)%NewRngBinsToProcessX2];       
  }  


  /* Take IFFT */
  dummyaddress = newrline - 1;
  CFFT(dummyaddress,NewRngBinsToProcess,1);


  /* Free Array */
  free (tempiq);

  return(0);
}


