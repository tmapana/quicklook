/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2lms.c
CODE CONTROLLER: Richard Lord
DESCRIPTION: 
Part of G2 SAR processor. Called from range compression stage. 
LMS adaptive filter for interference suppression. 

VERSION/AUTHOR/DATE : 1999-01-28 / Richard Lord / 1999-01-28
COMMENTS: 
Initial version.

VERSION/AUTHOR/DATE : 1999-02-22 / Richard Lord / 1999-02-22
COMMENTS:
Added standard implementation of LMS adaptive filter, i.e.
without using transfer functions.
Added LmsUpdateRate in parameter list.

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "g2func.h"

void lms(
        double sumiq[],
        Int4B RngFFTSize,
        Int4B RngBinsToProcess,
        double LmsUpdateRate,
        Int4B LmsNumWeights,
        Int4B LmsSidelobeOrder,
        double LmsWeights[],
        double rfi[]
        )

{
  FILE  *msg=NULL;

  Int4B RngFFTSizeX2,
        RngBinsToProcessX2,
        LmsNumWeightsX2,
        i, j, bin,
        indxi, indxq, indxi2, indxq2,
        constant,
        delay = 1,
        numtimes = 1;
  
  double *t_shift = NULL,
         *sumiq_pad = NULL,
         *dummyaddress = NULL,
         avg_power,
         max_mu,
         mu,
         muX2,
         first_div = 100.0,
         second_div = 1.0,
         e_i, e_q,
         total_i, total_q,
         temp,
         hr, hi;


  /* Initialise */
  LmsNumWeightsX2 = LmsNumWeights * 2;
  RngFFTSizeX2 = RngFFTSize * 2;
  RngBinsToProcessX2 = RngBinsToProcess * 2;
  constant = LmsNumWeights - 1 + delay;


  /* Find mu */
  avg_power = 0.0;
  indxi = 0;
  indxq = 1;
  for (bin=0; bin<RngBinsToProcess; bin++) {
    avg_power += (sumiq[indxi] * sumiq[indxi]) + (sumiq[indxq] * sumiq[indxq]);
    indxi++; indxi++;
    indxq++;indxq++;
  }
  avg_power /= (double)RngBinsToProcess;
  max_mu = 1.0 / (((double)LmsNumWeights + 1.0) * avg_power);

  mu = max_mu / first_div;  
  muX2 = mu * 2.0;


  /*******************************************
   * Implement standard LMS adaptive filter  * 
   *******************************************/

  if (LmsUpdateRate==1) { 

    /* Allocate space for arrays */
    sumiq_pad = (double *)malloc(sizeof(double)*((RngBinsToProcess + constant) * 2));
    if (sumiq_pad==NULL) {
      fprintf(msg,"ERROR - in array memory allocation!\n"); exit (1);
    }

    /* Create zero-padded array */
    for (i=0; i < (constant * 2); i++) sumiq_pad[i] = 0.0;
    for (i=0; i < RngBinsToProcess * 2; i++) sumiq_pad[i + (constant * 2)] = sumiq[i];

    for (bin=0; bin<RngBinsToProcess; bin++) {
      total_i = 0.0;  
      total_q = 0.0;
      indxi = bin * 2;
      indxq = (bin * 2) + 1;
      indxi2 = 0;
      indxq2 = 1;
      for (i=0; i<LmsNumWeights; i++) {
        total_i += (sumiq_pad[indxi] * LmsWeights[indxi2]) - (sumiq_pad[indxq] * LmsWeights[indxq2]);
        total_q += (sumiq_pad[indxi] * LmsWeights[indxq2]) + (sumiq_pad[indxq] * LmsWeights[indxi2]);
          
        indxi++;  indxi++;
        indxq++;  indxq++;
        indxi2++; indxi2++;
        indxq2++; indxq2++;
      }

      /* Compute Error */
      e_i = sumiq_pad[(bin + constant) * 2] - total_i;
      e_q = sumiq_pad[((bin + constant) * 2) + 1] - total_q;

      /* Update Weights */
      indxi = bin * 2;
      indxq = (bin * 2) + 1;
      for (i=0; i<LmsNumWeights; i++) {
        LmsWeights[i*2]     += muX2 * ((e_i * sumiq_pad[indxi]) + (e_q * sumiq_pad[indxq]));
        LmsWeights[(i*2)+1] += muX2 * ((e_q * sumiq_pad[indxi]) - (e_i * sumiq_pad[indxq]));
        
        indxi++; indxi++;
        indxq++; indxq++;
      }

      /* Update Output */
      sumiq[bin * 2]       = e_i;
      sumiq[(bin * 2) + 1] = e_q;
    }


  /* Free array space */
  free (sumiq_pad);


  } /* end if LmsUpdateRate */

  else {


    /**********************************************************
     * Implement LMS adaptive filter using transfer functions *
     **********************************************************/
  
    /* Allocate space for arrays */
    t_shift = (double *)malloc(sizeof(double)*RngFFTSizeX2);
    if (t_shift==NULL) {
      fprintf(msg,"ERROR - in array memory allocation!\n"); exit (1);
    }

    /* Find weights, iterating range line numtimes */

    for (j=0; j<numtimes; j++) {
      for (bin=0; bin<=(RngBinsToProcess - LmsNumWeights - delay); bin++) {
        total_i = 0.0;  
        total_q = 0.0;
        indxi = bin * 2;
        indxq = (bin * 2) + 1;
        indxi2 = 0;
        indxq2 = 1;
        for (i=0; i<LmsNumWeights; i++) {
          total_i += (sumiq[indxi] * LmsWeights[indxi2]) - (sumiq[indxq] * LmsWeights[indxq2]);
          total_q += (sumiq[indxi] * LmsWeights[indxq2]) + (sumiq[indxq] * LmsWeights[indxi2]);
          
          indxi++;  indxi++;
          indxq++;  indxq++;
          indxi2++; indxi2++;
          indxq2++; indxq2++;
        }

        /* Compute Error */
        e_i = sumiq[(bin + LmsNumWeights - 1 + delay) * 2] - total_i;
        e_q = sumiq[((bin + LmsNumWeights - 1 + delay) * 2) + 1] - total_q;

        /* Update Weights */
        indxi = bin * 2;
        indxq = (bin * 2) + 1;
        for (i=0; i<LmsNumWeights; i++) {
          LmsWeights[i*2]     += muX2 * ((e_i * sumiq[indxi]) + (e_q * sumiq[indxq]));
          LmsWeights[(i*2)+1] += muX2 * ((e_q * sumiq[indxi]) - (e_i * sumiq[indxq]));
        
          indxi++; indxi++;
          indxq++; indxq++;
        }
      }
      muX2 /= second_div;
    }    
  

    /*******************************
     * Calculate Transfer Function *
     *******************************/

    /* Initialise arrays */
    for (i=0; i<RngFFTSizeX2; i++) {
      t_shift[i] = 0.0;
      rfi[i] = 0.0;
    }

    indxi = 0;
    for (j=(RngFFTSize/2); j<RngFFTSize; j++) {
      temp = -((((double)j / RngFFTSize) * 2.0 * PI) - PI) * (double)delay;
      t_shift[indxi++] = cos(temp);
      t_shift[indxi++] = sin(temp);
    }

    indxi = RngFFTSize;
    for (j=0; j<(RngFFTSize/2); j++) {
      temp = -((((double)j / RngFFTSize) * 2.0 * PI) - PI) * (double)delay;
      t_shift[indxi++] = cos(temp);
      t_shift[indxi++] = sin(temp);
    }

    indxi = 0;
    indxi2 = LmsNumWeightsX2 - 2;
    indxq2 = LmsNumWeightsX2 - 1;
    for (j=0; j<LmsNumWeights; j++) {
      rfi[indxi++] = LmsWeights[indxi2--];
      rfi[indxi++] = LmsWeights[indxq2--];
    
      indxi2--; indxq2--;
    }

    dummyaddress = rfi - 1;
    CFFT(dummyaddress,RngFFTSize,-1);

    indxi = 0;
    indxq = 1;
    indxi2 = 0;
    indxq2 = 1;
    for (i=0; i<RngFFTSize; i++) {
      temp = rfi[indxi];
    
      rfi[indxi] = 1.0 - ((rfi[indxi] * t_shift[indxi2]) - (rfi[indxq] * t_shift[indxq2]));
      rfi[indxq] = -((temp * t_shift[indxq2]) + (rfi[indxq] * t_shift[indxi2]));

      indxi++;  indxi++;
      indxq++;  indxq++;
      indxi2++; indxi2++;
      indxq2++; indxq2++;
    }


    /********************
     * Remove Sidelobes *
     ********************/

    indxi = 0;  
    indxq = 1;

    switch(LmsSidelobeOrder) {

      case 1: for (bin=0; bin<RngFFTSize; bin++) {
                hr = rfi[indxi];
                hi = rfi[indxq];

                rfi[indxi] = 2*hr - hr*hr + hi*hi;
                rfi[indxq] = 2 * hi * (1 - hr);

                indxi++; indxi++;
                indxq++; indxq++;
              }
              break;

      case 2: for (bin=0; bin<RngFFTSize; bin++) {
                hr = rfi[indxi];
                hi = rfi[indxq];

                rfi[indxi] = 3*hr - 3*hr*hr + 3*hi*hi + hr*hr*hr - 3*hr*hi*hi;
                rfi[indxq] = hi * (3 - 6*hr + 3*hr*hr - hi*hi);

                indxi++; indxi++;
                indxq++; indxq++;
              }
              break;

      case 3: for (bin=0; bin<RngFFTSize; bin++) {
                hr = rfi[indxi];
                hi = rfi[indxq];

                rfi[indxi] = 4*hr - 6*hr*hr + 6*hi*hi + 4*hr*hr*hr - 12*hr*hi*hi -
                             hr*hr*hr*hr + 6*hr*hr*hi*hi - hi*hi*hi*hi;
                rfi[indxq] = -4*hi * (hr - 1 + hi) * (hr - 1 - hi) * (hr - 1);
              
                indxi++; indxi++;
                indxq++; indxq++;
              }

    } /* end switch */
  

    /* Free array space */
    free (t_shift);

  } /* end else */

}




