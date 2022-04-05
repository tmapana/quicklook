/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2moc.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2moc.c 

VERSION/AUTHOR/DATE : 1999-01-31 / Jasper Horrell / 1999-01-31
COMMENTS: 
Original.

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

struct MocCmdStruct {
  FILE  *msg;
  char  Version[STRING_SPACE];
  char  LBRFileName[STRING_SPACE];
  char  OutFileName[STRING_SPACE];
  char  OutTextFileName[STRING_SPACE];
  char  LogFileName[STRING_SPACE];
  Int4B RefRngBin;     /* set to DEFAULT_VAL_I, if not used */
  Int4B ProcStartPRI;  /* set to DEFAULT_VAL_I, if not used */
  Int4B ProcEndPRI;    /* set to DEFAULT_VAL_I, if not used */
  Int4B StartG2PRI;    /* set to DEFAULT_VAL_I, if not used */
  Int4B EndG2PRI;      /* set to DEFAULT_VAL_I, if not used */
  Int4B KernSize;      /* set to DEFAULT_VAL_I, if not used */
  char AntennaDirn;
  char  Radar;
  char  HBR;  
  double TerrainAlt;
  double NomAlt;
};

/* Prototypes for functions in g2cor.c file */
Int2B G2Moc (struct MocCmdStruct Cmd);
