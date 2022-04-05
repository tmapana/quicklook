/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2cor.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2cor.c 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Updated with this header info.

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

struct CorFileStruct {
  FILE  *msg;
  char  Version[STRING_SPACE];
  char  InputFileName[STRING_SPACE];
  char  OutputFileName[STRING_SPACE];
  Int4B MaxInputCols;
  Int4B BytesPerValue;
  Int4B RowHeaderBytesToSkip;
  Int4B RowFooterBytesToSkip;  
  Int4B InputStartRow;
  Int4B InputEndRow;
  Int4B InputStartCol;
  Int4B InputEndCol;
  double MaxMem;
};

/* Prototypes for functions in g2cor.c file */
Int2B G2Cor (struct CorFileStruct Cmd);
