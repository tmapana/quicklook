/*===================================================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2geth.h
CODE CONTROLLER: Richard Lord

DESCRIPTION: 
Part of G2 SAR processor. Header file for g2geth.c 

VERSION/AUTHOR/DATE : 0.1 / Richard Lord / 2000-03-31
COMMENTS: 

VERSION/AUTHOR/DATE :
COMMENTS:

=====================================================================*/

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
       );


