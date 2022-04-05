/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2rline.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2rline.c 

VERSION/AUTHOR/DATE : 1999-01-28 / Jasper Horrell / 1999-01-28
COMMENTS: 
Initial version.

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

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
           );

