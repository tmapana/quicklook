/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2notch.h
CODE CONTROLLER: Richard Lord
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2notch.c (notch 
interference filter) 

VERSION/AUTHOR/DATE : 1999-01-28 / Richard Lord / 1999-01-28
COMMENTS: 
Initial version.

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

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
           );


