/*===================================================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2geth.c
CODE CONTROLLER: Richard Lord

DESCRIPTION: 
Part of G2 SAR processor. Called from G2 stepped-frequency processor
Performs following steps:
   1) Constructs compression filter Hsys

VERSION/AUTHOR/DATE : 0.1 / Richard Lord / 2000-03-31
COMMENTS: 

VERSION/AUTHOR/DATE :
COMMENTS:

=====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"
#include "sf_procs.h"


int G2GetHsys(
       FILE *msg,
       FILE *infile,
       Int4B RngBinsToProcess,
       Int4B RngFFTSize,
       Int4B NewRngBinsToProcess,
       Int4B NumberOfFreqSteps,
       Int4B RngComRefFuncPhaseSign,
       Int4B upsample,
       double fcenter,
       double NarrowChirpBandwidth,
       double InputA2DFreq,
       double NarrowPulseLen,
       double InputStartSampleDelay,
       double RngComWinConstTime,
       double fcenter_n[],
       double Hsys[]
       )

{
  Int4B RngBinsToProcessX2,
        RngFFTSizeX2,
        NewRngBinsToProcessX2,
        step,
        winlength,
        Errors=0,
        buffer,
        point1, point2,
        j, indxi, indxq;

  double fshift,
         binlength,
         RangeTarget,
         time_shift,
         win,
         fcenter_min,
         fcenter_max,
         *rline=NULL,
         *newrline=NULL,
         temp1i, temp1q,
         val1;


  /* Miscellaneous */
  fcenter_min = 1e15;
  fcenter_max = 0.0;
  for (j=0; j<NumberOfFreqSteps; j++) {
    if (fcenter_n[j] < fcenter_min) fcenter_min=fcenter_n[j];
    if (fcenter_n[j] > fcenter_max) fcenter_max=fcenter_n[j];
  }
  RngBinsToProcessX2 = 2*RngBinsToProcess;
  RngFFTSizeX2 = 2*RngFFTSize;
  NewRngBinsToProcessX2 = 2*NewRngBinsToProcess;
  binlength = C / (2.0 * InputA2DFreq);
  RangeTarget = ((RngFFTSize / 2.0) * binlength) + (InputStartSampleDelay * (C / 2.0));
  time_shift = ((2.0 * RangeTarget) / C);
 

  /* Allocate space for arrays */
  rline = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  newrline = (double *)malloc(sizeof(double)*NewRngBinsToProcessX2);
  if ( rline==NULL || newrline==NULL ) Errors++;
  if (Errors != 0) { 
    fprintf(msg,"ERROR - in array memory allocation!\n"); exit (1);
  }


  /* Initialise arrays */
  for (j=0; j<NewRngBinsToProcessX2; j++) Hsys[j] = 0.0;
  for (j=0; j<NewRngBinsToProcessX2; j++) newrline[j] = 0.0;

  for (step=0; step<NumberOfFreqSteps; step++)
    {

    fshift = (double)(((fcenter_n[step]-fcenter)/InputA2DFreq)*(double)RngFFTSize);


    /* Get point target return */
    G2GetPointSim(
         msg,
         RngFFTSize,
         upsample,
         fcenter_n[step],
         NarrowChirpBandwidth,
         InputA2DFreq,
         NarrowPulseLen,
         InputStartSampleDelay,
         RangeTarget,
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


    /* Range-compress range line */
    G2RangeCompress(
              msg,
              RngFFTSize,
              RngComRefFuncPhaseSign,
              NarrowChirpBandwidth,
              InputA2DFreq,
              NarrowPulseLen,
              rline
              );
  

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


  /* Phase compensation */
  G2PhaseComp(
       msg,
       RngFFTSize,
       NewRngBinsToProcess,
       InputA2DFreq,
       time_shift,
       newrline
       );



  /* Invert spectrum */

  buffer = 0;
  point1 = (Int4B)((NewRngBinsToProcess / 2) + floor(((fcenter_min - fcenter - (NarrowChirpBandwidth / 2.0)) / InputA2DFreq) * (double)RngFFTSize) + buffer - 0);
  point2 = (Int4B)((NewRngBinsToProcess / 2) + floor(((fcenter_max - fcenter + (NarrowChirpBandwidth / 2.0)) / InputA2DFreq) * (double)RngFFTSize) - buffer - 1);

  /* 
  point1 is the first nonzero Hsys,
  point2 is the first zero Hsys.
  */

  temp1i = newrline[(point1)*2];
  temp1q = newrline[((point1)*2)+1];
  val1 = sqrt((temp1i * temp1i) + (temp1q * temp1q));


  indxi = 2*point1;
  indxq = indxi+1;
  winlength = point2-point1;
  for (j=0; j<winlength;j++) 
  {
    win = RngComWinConstTime+(1.0-RngComWinConstTime)*pow(sin(PI*(Int4B)j/(Int4B)(winlength-1.0)),2.0);
    Hsys[indxi] = newrline[indxi] / ((newrline[indxi]*newrline[indxi]) + (newrline[indxq]*newrline[indxq])) * win * val1;
    Hsys[indxq] = -newrline[indxq] / ((newrline[indxi]*newrline[indxi]) + (newrline[indxq]*newrline[indxq])) * win * val1;
    indxi++;
    indxi++;
    indxq++;
    indxq++;
  }


  /* Free array space */
  free (rline);
  free (newrline);

  return(0);
}










