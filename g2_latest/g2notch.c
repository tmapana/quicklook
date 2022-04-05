/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2notch.c
CODE CONTROLLER: Richard Lord
DESCRIPTION: 
Part of G2 SAR processor. Called from range compression stage. 
Notch filter for interference suppression. 

VERSION/AUTHOR/DATE : 1999-01-28 / Richard Lord / 1999-01-28
COMMENTS: 
Initial version.
Median Filter Bug fixed on 1999-01-29

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"
#include "g2rline.h"

void notch(
           Int4B NotchNumFFTLines,
           Int4B NotchMedianKernLen,
           double NotchCutoff,
           FILE *infile,
           Int4B RngFFTSize,
           Int4B RngBinsToProcess,
           Int4B PreSummedPulsesToUse,
           Int4B PreSumRatio,
           Int4B InputDataType,
           Int4B skip,
           Int4B pulse,
           Int4B *inp,
           unsigned char IQInChar[],
           float IQInFloat[],
           double DopCentroid,
           double InputPRF,
           double InputDCOffsetI,
           double InputDCOffsetQ,
           double InputIQRatio,
           double stc[],
           double sumiq[],
           double rfi[]
          )


{
  FILE  *msg=NULL;
  
  long file_pos;

  Int4B RngFFTSizeX2,
        RngBinsToProcessX2,
        i, j,
        indxi, indxq,
        inp_temp,
        Errors = 0;

  double *dummyaddress = NULL,
         *sumiq2 = NULL,
         *avg_fft = NULL,
         *shifted_arr = NULL,
         *sort_arr = NULL,
         *median_arr = NULL,
         small_nonzero_number = 0.000001;

  RngFFTSizeX2 = RngFFTSize * 2;
  RngBinsToProcessX2 = RngBinsToProcess * 2;


  /*****************************
   * ALLOCATE SPACE FOR ARRAYS *
   *****************************/

  sumiq2 = (double *)malloc(sizeof(double)*RngFFTSizeX2);
  avg_fft = (double *)malloc(sizeof(double)*RngFFTSize);
  shifted_arr = (double *)malloc(sizeof(double)*RngFFTSize);
  sort_arr = (double *)malloc(sizeof(double)*NotchMedianKernLen);
  median_arr = (double *)malloc(sizeof(double)*RngFFTSize);

  if (sumiq2==NULL || avg_fft==NULL || shifted_arr==NULL ||
      sort_arr==NULL || median_arr==NULL) Errors++;

  if (Errors != 0) {
    fprintf(msg,"ERROR - in array memory allocation!\n"); 
    exit (1);
  }


  /***************************************************
   * OBTAIN AVERAGE RANGE SPECTRUM IN dB             *
   ***************************************************/

  /* Make copy of sumiq in sumiq2 */
  for (i=0; i<RngFFTSizeX2; i++) sumiq2[i] = sumiq[i];

  dummyaddress = sumiq2 - 1;
  CFFT(dummyaddress,RngFFTSize,-1);

  indxi = 0;
  indxq = 1;
  for (i=0; i<RngFFTSize; i++) {
    avg_fft[i] = (sumiq2[indxi] * sumiq2[indxi]) + (sumiq2[indxq] * sumiq2[indxq]);
    indxi++; indxi++;
    indxq++; indxq++;
  }
  
  if ((pulse + NotchNumFFTLines) > PreSummedPulsesToUse) {
     NotchNumFFTLines = PreSummedPulsesToUse - pulse;
  }

  if (NotchNumFFTLines > 1) {

    /* save parameters */
    file_pos = ftell(infile);
    inp_temp = *inp;

    for (j=1; j<NotchNumFFTLines; j++) {

      /* obtain next range line, contained in sumiq2 */
      G2RncReadLine(
                    infile,
                    RngFFTSize,
                    RngBinsToProcess,
                    PreSumRatio,
                    InputDataType,
                    skip,
                    pulse,
                    &inp_temp,
                    IQInChar,
                    IQInFloat,
                    DopCentroid,
                    InputPRF,
                    InputDCOffsetI,
                    InputDCOffsetQ,
                    InputIQRatio,
                    stc,
                    sumiq2
                   );

      dummyaddress = sumiq2 - 1;
      CFFT(dummyaddress,RngFFTSize,-1);

      /* add sumiq2 to avg_fft */
      indxi = 0;
      indxq = 1;
      for (i=0; i<RngFFTSize; i++) {
        avg_fft[i] += (sumiq2[indxi] * sumiq2[indxi]) + (sumiq2[indxq] * sumiq2[indxq]);
        indxi++; indxi++;
        indxq++; indxq++;
      }

    }

    /* restore parameters */
    fseek(infile, file_pos, SEEK_SET);
    *inp = inp_temp;
  }

  /* convert to dB */
  for (i=0; i<RngFFTSize; i++) {
    if (avg_fft[i] <= 0.0) avg_fft[i] = small_nonzero_number;
    avg_fft[i] = 10.0 * log10(avg_fft[i]);
  }


  /*****************************
   * APPLY MEDIAN FILTER       *
   *****************************/

  /* shift avg_fft */
  j = 0;
  for (i=((RngFFTSize/2)+1); i<RngFFTSize; i++) {
    shifted_arr[j] = avg_fft[i];
    j++;
  }
  j = (RngFFTSize/2)-1;
  for (i=0; i<=(RngFFTSize/2); i++) {
    shifted_arr[j] = avg_fft[i];
    j++;
  }

  /* apply median filter */
  for (i=0; i<((NotchMedianKernLen - 1) / 2); i++) 
    median_arr[i] = shifted_arr[i];

  for (i=(RngFFTSize-NotchMedianKernLen)+((NotchMedianKernLen-1)/2);i<RngFFTSize; i++) 
    median_arr[i] = shifted_arr[i];

  for (i=0; i<=(RngFFTSize - NotchMedianKernLen); i++) {
    for (j=0; j<NotchMedianKernLen; j++) sort_arr[j] = shifted_arr[i+j];    
    median_arr[i + ((NotchMedianKernLen - 1) / 2)] = median(sort_arr, NotchMedianKernLen);
  }

  /* re-shift median_arr to correspond with original avg_fft array */
  j = 0;
  for (i=((RngFFTSize/2)-1); i<RngFFTSize; i++) {
    shifted_arr[j] = median_arr[i];
    j++;
  }
  j = (RngFFTSize/2)+1;
  for (i=0; i<=((RngFFTSize/2)-2); i++) {
    shifted_arr[j] = median_arr[i];
    j++;
  }


  /*****************************
   * CREATE NOTCH FILTER       *
   *****************************/

  indxi = 0;  
  indxq = 1;
  for (i=0; i<RngFFTSize; i++) {
  
    if ((avg_fft[i] - shifted_arr[i]) > NotchCutoff) 
      rfi[indxi] = 0.0;
    else 
      rfi[indxi] = 1.0;

    rfi[indxq] = 0.0;
    indxi++; indxi++;
    indxq++; indxq++;
  }


  /*********************
   * TAPER OFF NOTCHES *
   *********************/

  if (rfi[0] == 0) {
    rfi[RngFFTSizeX2-2] *= 0.5;
    rfi[2] *= 0.5;
  }
  if (rfi[RngFFTSizeX2-2] == 0) {
    rfi[RngFFTSizeX2-4] *= 0.5;
    rfi[0] *= 0.5;
  }

  indxi = 2;  
  for (i=2; i<(RngFFTSizeX2-2); i+=2) {
    if (rfi[i] == 0) {
      rfi[i-2] *= 0.5;
      rfi[i+2] *= 0.5;
    }
  }


  /*********************
   * FREE ARRAY SPACE  *
   *********************/

  free (sumiq2);
  free (avg_fft);
  free (shifted_arr);
  free (sort_arr);
  free (median_arr);

}


