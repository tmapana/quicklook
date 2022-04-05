/*==========================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2rline.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Called from G2 range compression and notch filter.
Performs following steps:
   1) Read in range line
   2) remove Doppler centroid and perform presumming
   3) Rescale for presum and unSTC data (end of array is zero 
      padded to FFT size)
   4) Output array is sumiq 

VERSION/AUTHOR/DATE : 1999-01-28 / Jasper Horrell / 1999-01-28
COMMENTS: 
Initial version extracted from g2rnc.c (range compression prog).

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"

void G2RncReadLine( 
           FILE *infile,
           Int4B RngFFTSize,
           Int4B RngBinsToProcess,
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
           double sumiq[]
           )

{
  FILE  *msg=NULL;

  Int4B RngFFTSizeX2,
        RngBinsToProcessX2,
        i, j, bin,
        indxi, indxq;

  double cenphase,
         ceni, cenq,
         dati=0.0, datq=0.0;


  /* Initialize */
  RngFFTSizeX2 = RngFFTSize * 2;
  RngBinsToProcessX2 = RngBinsToProcess * 2;


  /* Initialize presum array*/
  for (i=0; i<RngFFTSizeX2; i++) sumiq[i] = 0.0;

  /* Perform presumming */
  for (i=0; i<PreSumRatio; i++) {

    /* Read in required range bins */
    if (InputDataType == 0) { 
      if ( fread(IQInChar,sizeof(unsigned char),RngBinsToProcessX2,infile) 
            != RngBinsToProcessX2 ) { 
        fprintf(msg,"ERROR - in read from input file!\n"); exit(1); 
      }
    } /* end if InputDataType */
    else if (InputDataType == 3) {
      if ( fread(IQInFloat,sizeof(float),RngBinsToProcessX2,infile)
           != RngBinsToProcessX2) { 
        fprintf(msg,"ERROR - in read from input file!\n"); exit(1); 
      }
    } /* end else if */

    cenphase = -2.0*PI*DopCentroid*(double)(*inp)/(double)InputPRF; /* Richard */
    *inp = *inp + 1;  /* Richard */

    /* Remove Doppler centroid and add to presum array */

    if (DopCentroid != 0.0) {
      ceni = cos(cenphase);
      cenq = sin(cenphase);
      indxi = 0;
      indxq = 1;
      for (j = 0; j < RngBinsToProcess; j++) {
         if (InputDataType == 0) {
           dati = ((double)IQInChar[indxi])/InputIQRatio-InputDCOffsetI;
	       datq = (double)IQInChar[indxq]-InputDCOffsetQ;  
         }
         else if (InputDataType == 3) {
 	       dati = ((double)IQInFloat[indxi])/InputIQRatio-InputDCOffsetI;
	       datq = (double)IQInFloat[indxq]-InputDCOffsetQ;
         }
         sumiq[indxi] += dati*ceni - datq*cenq;
         sumiq[indxq] += datq*ceni + dati*cenq;
         indxi++; indxi++;
         indxq++; indxq++;
       } /* end for loop */
    }
    else  { /* for zero doppler centroid */
      indxi = 0;
      indxq = 1;
      for (j=0; j<RngBinsToProcess;j++) {
        if (InputDataType == 0) { 
          sumiq[indxi] += ((double)IQInChar[indxi])/InputIQRatio-InputDCOffsetI;
          sumiq[indxq] += (double)IQInChar[indxq]-InputDCOffsetQ;
        }
        else if (InputDataType == 3) { 
          sumiq[indxi] += ((double)IQInFloat[indxi])/InputIQRatio-InputDCOffsetI;
          sumiq[indxq] += (double)IQInFloat[indxq]-InputDCOffsetQ;
        }
        indxi++; indxi++;
        indxq++; indxq++;
      }  /* end for j loop */  
    } /* end else for zero Doppler centroid */

    /* Skip data in input file until start of next data segment */
    fseek(infile,skip,1);

  }  /* End presum loop */


  /**********************************************
  * Rescale for presum and unSTC data          *
  * (end of array is zero padded to FFT size)  *
  **********************************************/

  j=0;
  for (bin=0; bin<RngBinsToProcess; bin++) {
    sumiq[j] =  sumiq[j]/((double)PreSumRatio*stc[bin]);
    j++;
    sumiq[j] =  sumiq[j]/((double)PreSumRatio*stc[bin]);
    j++;
  } /* end for bin loop */

} /* end G2RncReadLine function */
