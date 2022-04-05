/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2sniffdc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Program to sniff out DC value of I and Q channels 
where input are unsigned chars. May start at some offset in the file and 
specify pts to investigate. Uses slightly convoluted way of finding 
averages so as to avoid precision loss.

VERSION/AUTHOR/DATE : pre 1999-01-26 / Jasper Horrell / pre 1999-01-26
COMMENTS: 
(1998-10-02) Rewrote to take command line args.
(1998-10-14) Cosmetic change.
(1999-01-15) Minor changes to file open order etc.
(1999-01-25) Allow use as function or command line. Renamed to g2sniffdc.c
             and g2sniffdc.h (new header file) 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS:
Move FUNC to g1func.h and rename to FUNC_SNIFFDC.

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1999-08-02
COMMENTS:
Option to write output to log file for G2 processor

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-09-20
COMMENTS:
Check on gain imbalances as well as the DC offsets. Remove ability to
be compiled as a function. Add option to start and end at some arbitrary column.

=========================================*/

                 
#include "g2func.h"

#define PROG_VERSION "0.2"

Int2B main (int argc, char *argv[])
{
  FILE	*InFile=NULL,*msg=NULL,*LogFile=NULL;
  char InputFileName[STRING_SPACE],ProgVersion[STRING_SPACE],
       LogFileName[STRING_SPACE];
  unsigned char *InIQ=NULL;
  Int4B	i,indx,line,InRows,InCols,StartRow=0,StartCol=0,
        ColsToProc=DEFAULT_VAL_I,RowsToProc=DEFAULT_VAL_I,
        InPtSize=2; /* unsigned char complex input */
  double AveI=0.0,AveQ=0.0,AveDevI=0.0,AveDevQ=0.0,
         RowDevI,RowDevQ,RowI,RowQ;
  Int2B CreateLogFile=0;

  msg = stdout;
  fprintf(msg,"\n------------\n");
  fprintf(msg,"Prog SniffDC (Ver. %s) - Code J.M.Horrell\n",PROG_VERSION);

  if (argc < 4) { 
    fprintf(msg,"(Sniffs out DC values in unsigned char IQ input.)\n");
    fprintf(msg,"\nUSAGE: sniffdc [InFile] [InRows] [InCols] <optional params>\n");
    fprintf(msg,"\nwhere:\n"); 
    fprintf(msg,"      InRows  -  input file rows\n");
    fprintf(msg,"      InCols  -  size of each row in complex points\n");
    fprintf(msg,"Optional params (at end of line, no braces):\n");
    fprintf(msg,"      <StartRow=n> - row at which to start processing\n");
    fprintf(msg,"      <RowsToProc=n> - num rows to process (default end of file)\n");    
    fprintf(msg,"      <StartCol=n> - col at which to start proc (default 0)\n");
    fprintf(msg,"      <ColsToProc=n> - num columns to process (default end of row)\n");
    fprintf(msg,"      <LogF=s>     - name of log file\n");
    exit(1); 
  }

  /* Misc */
  strcpy(InputFileName,argv[1]);
  sscanf(argv[2], "%ld", &InRows);
  sscanf(argv[3], "%ld", &InCols);

  /* check for optional parameters */
  for (i=4; i<argc; i++) {

    if (strncmp(argv[i],"StartRow=",9)==0) {
      sscanf(argv[i],"StartRow=%ld",&StartRow);
    }
    else if (strncmp(argv[i],"RowsToProc=",11)==0) {
      sscanf(argv[i],"RowsToProc=%ld",&RowsToProc);
    }
    else if (strncmp(argv[i],"StartCol=",9)==0) {
      sscanf(argv[i],"StartCol=%ld",&StartCol);
    }
    else if (strncmp(argv[i],"ColsToProc=",11)==0) {
      sscanf(argv[i],"ColsToProc=%ld",&ColsToProc);
    }
    else if (strncmp(argv[i],"LogF=",5)==0) {
      sscanf(argv[i],"LogF=%s",LogFileName);
      CreateLogFile = 1;
      if ( (LogFile = fopen (LogFileName, "wt") ) == NULL ) {
         fprintf (msg,"ERROR: log file %s not opened!\n",LogFileName);
         exit(1);
      }
      msg = LogFile;
      fprintf(msg,"Prog SniffDC (Ver. %s) - Code J.M.Horrell\n",PROG_VERSION);     
    }  /* end else if */   
    else {
      fprintf(msg,"ERROR - parameter %s unknown!\n",argv[i]);
      exit(1);
    }     

  }  /* end for loop */

  /* Check that dimensions set correctly */
  if (ColsToProc==DEFAULT_VAL_I) {
    ColsToProc = InCols-StartCol;
  }
  if (RowsToProc==DEFAULT_VAL_I) {
    RowsToProc = InRows-StartRow;
  }
  if (StartCol+ColsToProc>InCols) {
    fprintf(msg,"ERROR - Not enough input cols (%ld) - %ld are required!\n",
    InCols,StartCol+ColsToProc);
  }
  if (StartRow+RowsToProc>InRows) {
    fprintf(msg,"ERROR - Not enough input rows (%ld) -%ld are required!\n",
    InRows,StartRow+RowsToProc);
  }

  /* Messages */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"Input file          : %s\n",InputFileName);
  fprintf(msg,"Input file rows     : %ld\n",InRows);
  fprintf(msg,"Input file columns  : %ld\n",InCols);
  fprintf(msg,"Start row           : %ld\n",StartRow);
  fprintf(msg,"Rows to process     : %ld\n",RowsToProc);
  fprintf(msg,"Start column        : %ld\n",StartCol);
  fprintf(msg,"Columns to process  : %ld\n",ColsToProc);

  /* Open input file */
  if ( (InFile = fopen (InputFileName, "rb") ) == NULL ) {
     fprintf (msg,"ERROR: input file %s not opened!\n",InputFileName);
     exit(1);
  }
  /* Allocate memory */
  InIQ = (unsigned char *)malloc(sizeof(unsigned char)*ColsToProc*InPtSize);
  if (InIQ==NULL) {
    fprintf(msg,"ERROR - in array memory allocation!\n");
    exit(1);
  }
 
  /*-----------------*/
  /* FIND DC OFFSETS */

  fseek(InFile,(StartRow*InCols+StartCol)*InPtSize,SEEK_SET); /* move to start proc */
  for (line = 0; line < RowsToProc; line++) {
    fread (InIQ, sizeof(unsigned char), ColsToProc*InPtSize, InFile); /* uchar only!*/
    indx = 0;
    RowI = 0.0;
    RowQ = 0.0;
    for (i = 0; i < ColsToProc; i++) {  /* find sum of I and sum of Q for row */     
      RowI += (float)InIQ[indx++];
      RowQ += (float)InIQ[indx++];
     }
    AveI += RowI/(float)ColsToProc; /* add average for row to overall sum */
    AveQ += RowQ/(float)ColsToProc;
    fseek(InFile,(InCols-ColsToProc)*InPtSize,SEEK_CUR); /* update file pointer */
  } /* end for line loop */
  AveI = AveI/(float)RowsToProc; /* average for data processed */
  AveQ = AveQ/(float)RowsToProc;
  fprintf(msg,"\nAverage I value     : %f\n", AveI);
  fprintf(msg,"Average Q value     : %f\n", AveQ);

  /*--------------------------*/
  /* INVESTIGATE GAIN BALANCE */
  
  fseek(InFile,(StartRow*InCols+StartCol)*InPtSize,SEEK_SET);
  for (line = 0; line < RowsToProc; line++) {
    indx = 0;
    RowDevI = 0.0;
    RowDevQ = 0.0;
    for (i = 0; i < ColsToProc; i++) { /* find sum of deviations for row */
      RowDevI += abs((float)InIQ[indx++]-AveI);
      RowDevQ += abs((float)InIQ[indx++]-AveQ);
    }
    AveDevI += RowDevI/(float)ColsToProc; /* add average for row to overall sum */
    AveDevQ += RowDevQ/(float)ColsToProc;
    fseek(InFile,(InCols-ColsToProc)*InPtSize,SEEK_CUR); /* update file pointer */
  } /* end for line loop */
  AveDevI = AveDevI/(float)RowsToProc; /* average for data processed */
  AveDevQ = AveDevQ/(float)RowsToProc;
  fprintf(msg,"Average I deviation : %f\n",AveDevI);
  fprintf(msg,"Average Q deviation : %f\n",AveDevQ);
  fprintf(msg,"Average I/Q ratio   : %f\n",AveDevI/AveDevQ);
  
  fprintf(stdout,"DONE!\n");
  
  /* Tidy up */
  free(InIQ);
  fclose(InFile);
  if (CreateLogFile)
    fclose(LogFile);

  return(0);  /* for successful completion */

}  /* end sniffdc code */
