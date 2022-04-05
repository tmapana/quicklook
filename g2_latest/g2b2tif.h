/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2b2tif.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2b2tif.c 

VERSION/AUTHOR/DATE : 1999-01-27 / Jasper Horrell / 1999-01-27
COMMENTS: 
Initial version.

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

struct B2TifStruct {
  FILE  *msg;
  char  Version[STRING_SPACE];
  char  InputFileName[STRING_SPACE];
  char  OutputFileName[STRING_SPACE];
  Int4B Rows;
  Int4B Cols;
  char Endian[STRING_SPACE];
  Int2B Orient;
};

/* Prototypes for functions in g2b2tif.c file */
Int2B B2Tif (struct B2TifStruct Cmd);
