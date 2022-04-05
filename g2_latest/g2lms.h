/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2lms.h
CODE CONTROLLER: Richard Lord
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2lms.c 

VERSION/AUTHOR/DATE : 1999-01-28 / Richard Lord / 1999-01-28
COMMENTS: 
Initial version.
LmsWeights now defined in main procedure (1999-02-08).

VERSION/AUTHOR/DATE : 1999-02-22 / Richard Lord / 1999-02-22
COMMENTS:
Added LmsUpdateRate in parameter list.

=========================================*/

void lms(
        double sumiq[],
        Int4B RngFFTSize,
        Int4B RngBinsToProcess,
        double LmsUpdateRate,
        Int4B LmsNumWeights,
        Int4B LmsSidelobeOrder,
        double LmsWeights[],
        double rfi[]
        );

