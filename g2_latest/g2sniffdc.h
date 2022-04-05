/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2sniffdc.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2sniffdc.c. 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Add this header info.

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/



/* Header file for g2sniffdc program */

struct SniffDCStruct {
  FILE *msg;
  char  Version[STRING_SPACE];
  char  InputFileName[STRING_SPACE];
  Int4B Rows;
  Int4B Columns;
  Int4B StartRow;
  float AveI;   /* this returned */
  float AveQ;   /* this returned */
};

Int2B SniffDC(struct SniffDCStruct *Cmd);  /* pass pointer to allow 
                                              return values */