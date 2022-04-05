/*==========================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2cor.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Program to corner turn a binary data set.

VERSION/AUTHOR/DATE : pre 1998-12-10 / Jasper Horrell / 1994 - 1998-12-10
COMMENTS: 
Version of corner.c modified for G2 SAR processor. Allow compilation as 
a function or standalone. 
(1997-12-03) Basically a complete rewrite. Remove presum option 
   (actually used to skip rather than presum), rename variables and
   change command file structure.
(1997-12-05) Check program version in command file.
(1998-07-20) Open spec file as text to improve DOS/Win32 compatibility
(1998-12-09) Allow to be run from command line with params.
(1998-12-10) Cosmetic changes to initial help.

VERSION/AUTHOR/DATE : 1999-01-18 / Jasper Horrell / 1999-01-18
COMMENTS:
Substantial changes to input parameters. If function, no longer cmd file 
option, but rather a structure containing the parameters is passed. If from 
command line, either use cmd file or run with params (some optional). 
Specify MaxMem instead of proc segment size, etc

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS:
Fixed bugs (MaxMem calc of columns and Int4B out of range.           
Move define FUNC to g2func.h (rename to FUNC_CORNER).
1999-08-03 - cosmetic message change
1999-09-16 - added else's to the statements to parse cmd file.
2000-03-04 - removed sizeof(unsigned char) from malloc of UnCharBlock
             (no change to functionality, but more logical code).

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

#include "g2func.h"
#include "g2cor.h"

/* Misc defns limited to this file */
#define PROG_VERSION "1999-01-26"
#define MB2BYTES 1048576.0
#define MEM_DEFAULT 100.0

/* Function prototypes */
Int2B WriteCornerTmplCmdFile(char CmdFileName[]);


/******************************************************************/
/* MAIN CORNER TURN PROGRAM/FUNCTION */

#if FUNC_CORNER
Int2B G2Cor (struct CorFileStruct Cmd) 
#else
Int2B main (int argc,char *argv[])
#endif
{
  FILE *OutFile,*cmdfp=NULL,*InFile,*msg=NULL;
  char ProgVersion[STRING_SPACE]="",InputFileName[STRING_SPACE],
       OutputFileName[STRING_SPACE];
  Int4B byte,
        BytesPerValue,
	    BytesToSkip,       /* to required column of next row */
	    BytesToSkipAtStart,
        column,
	    ColsPerProcSegment=1,    
        DataRowBytes=0,
        Errors = 0,
        i,    
        InputColsToUse,
        InputEndCol,
        InputEndRow,
        InputFileSize,
        InputRowsToUse,   
        InputStartCol,
        InputStartRow,
        LastProcSegmentCols,
        MaxInputCols,
        MaxInputRows,
        NumProcSegments,      
        row,
        RowHeaderBytesToSkip,
        RowFooterBytesToSkip,
        segment;
  double ActualMem,
         ColsPerProcSegmentCalc,
         MaxMem;
	
  time_t timeS,timeE; /* Start and end times */

  unsigned char **UnCharArray, /* Pointer to pointer to the array */
		*UnCharBlock,  /* Pointers to arrays for data */
		*InputRow,*InputCol;
#if FUNC_CORNER

  /* Assign variables from the structure passed */
  msg = Cmd.msg;
  strcpy(ProgVersion,Cmd.Version);
  strcpy(InputFileName,Cmd.InputFileName);
  strcpy(OutputFileName,Cmd.OutputFileName);
  MaxInputCols         = Cmd.MaxInputCols;
  BytesPerValue        = Cmd.BytesPerValue;
  RowHeaderBytesToSkip = Cmd.RowHeaderBytesToSkip;
  RowFooterBytesToSkip = Cmd.RowFooterBytesToSkip;  
  InputStartRow        = Cmd.InputStartRow;
  InputEndRow          = Cmd.InputEndRow;
  InputStartCol        = Cmd.InputStartCol;
  InputEndCol          = Cmd.InputEndCol;
  MaxMem               = Cmd.MaxMem;

  fprintf(msg,"\nCORNER TURN STAGE...\n");
  
  /* Check version IDs match */
  if (strcmp(ProgVersion,PROG_VERSION)!=0) {
    fprintf(msg,
      "WARNING - command version (%s) not same as program (%s)!\n\n",
      ProgVersion,PROG_VERSION);
  }

#else /* called from command line */
  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n------------\n");
  fprintf(msg,"Prog: CORNER (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1994-1999\n");

  /* If wrong number of arguments, print help and exit */
  if ( (argc != 2 && argc < 7) || argc > 12) { 
    fprintf(msg,"Corner-turn (transpose) binary data set.\n\n");
    fprintf(msg,"Choice of two USAGE methods:\n");
    fprintf(msg,
      "  1) To use command line parameters, syntax :\n");
    fprintf(msg,
      "       corner InFile OutFile MaxInputCols StartRow EndRow BytesPerValue\n");
    fprintf(msg,
      "     Optional params (use at end of line):\n");
    fprintf(msg,
      "       <startcol=n> - input start column (excl. header/footer)\n");
    fprintf(msg,
      "       <endcol=n>   - input end column (excl. header/footer)\n");
    fprintf(msg,
      "       <header=n>   - header bytes to skip for each row\n");
    fprintf(msg,
      "       <footer=n>   - footer bytes to skip for each row\n");
    fprintf(msg,
      "       <mem=n>      - max memory to use in MBytes\n"); 
    fprintf(msg,"     e.g. 'corner test.in test.cor 512 0 199 2 mem=100'\n\n");
    fprintf(msg,"  2) To use command file, syntax :\n");
    fprintf(msg,"       corner [cmd file]\n");
    fprintf(msg,"     Type 'corner -tmpl' to generate a commented template\n");
    fprintf(msg,"     command file (named 'corner.tmp').\n\n");
    exit(1); 
  }

  /* Write template command file, if requested, and exit */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0) {
     if (WriteCornerTmplCmdFile("corner.tmp")==0)
       fprintf(msg,"Template command file `corner.tmp' written!\n");
     else
       fprintf(msg,"ERROR - in writing template command file\n");     
     exit(1);
  }


  /* IF COMMAND LINE USAGE WITH SPEC FILE */
  if (argc==2) {
 
    if ( (cmdfp = fopen(argv[1],"r")) == NULL) {
      fprintf(msg,"ERROR - command %s not found/opened!\n",argv[1]);
      exit(1);
    }

    /* read processing specs */
    if (SeekP(cmdfp,"$ProgramVersion","=>",-1,0)) {
      fprintf(msg,
        "WARNING - in parsing command file (prog version not found)!\n\n"); 
    } 
    else {
      ReadStr(cmdfp,ProgVersion,STRING_SPACE); 
    }

    Errors = 0;
 
    if (SeekP(cmdfp,"$InputFileName","=>",-1,0)) Errors++; 
    else ReadStr(cmdfp,InputFileName,STRING_SPACE);  
    if (SeekP(cmdfp,"$OutputFileName","=>",-1,0)) Errors++; 
    else ReadStr(cmdfp,OutputFileName,STRING_SPACE); 
    if (SeekP(cmdfp,"$MaxInputCols","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&MaxInputCols);
    if (SeekP(cmdfp,"$BytesPerValue","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&BytesPerValue);
    if (SeekP(cmdfp,"$RowHeaderBytesToSkip","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&RowHeaderBytesToSkip);
    if (SeekP(cmdfp,"$RowFooterBytesToSkip","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&RowFooterBytesToSkip);
    if (SeekP(cmdfp,"$InputStartRow","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&InputStartRow);
    if (SeekP(cmdfp,"$InputEndRow","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&InputEndRow);
    if (SeekP(cmdfp,"$InputStartCol","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&InputStartCol);
    if (SeekP(cmdfp,"$InputEndCol","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%ld",&InputEndCol);         
    if (SeekP(cmdfp,"$MaxMem","=>",-1,0)) Errors++; 
    else fscanf(cmdfp,"%lf",&MaxMem);

    fclose (cmdfp);

    if (Errors!=0) {
      fprintf(msg,"ERROR - %ld errors in parsing command file!\n\n",Errors);
      exit(1);
    }

    /* Check version IDs match */
    if (strcmp(ProgVersion,PROG_VERSION)!=0) {
      fprintf(msg,
        "WARNING - command file version (%s) not same as program (%s)!\n\n",
         ProgVersion,PROG_VERSION);
    }

  } /* end if argc == 2 (i.e. command line with spec file)  */


  /* IF COMMAND LINE USAGE WITH PARAMETERS */  
  if (argc > 6) {  
    strcpy(InputFileName,argv[1]);
    strcpy(OutputFileName,argv[2]);
    sscanf(argv[3],"%ld",&MaxInputCols);
    sscanf(argv[4],"%ld",&InputStartRow);
    sscanf(argv[5],"%ld",&InputEndRow);
    sscanf(argv[6],"%ld",&BytesPerValue);
    
    /* set up other params to default values */
    InputStartCol = 0;
    InputEndCol = MaxInputCols - 1; 
    RowHeaderBytesToSkip = 0;
    RowFooterBytesToSkip = 0;
    MaxMem = MEM_DEFAULT;
  }  /* end if argc > 6 */

  /* Check for optional arguments (occurs if argc > 7) */
  for (i=7;i<argc;i++) {
    if (strncmp(argv[i],"startcol=",9)==0) { 
      sscanf(argv[i],"startcol=%ld",&InputStartCol); 
    }
    else if (strncmp(argv[i],"endcol=",7)==0) { 
      sscanf(argv[i],"endcol=%ld",&InputEndCol); 
    }
    if (strncmp(argv[i],"header=",7)==0) { 
      sscanf(argv[i],"header=%ld",&RowHeaderBytesToSkip); 
    }
    else if (strncmp(argv[i],"footer=",7)==0) { 
      sscanf(argv[i],"footer=%ld",&RowFooterBytesToSkip); 
    }
    else if (strncmp(argv[i],"mem=",4)==0) { 
      sscanf(argv[i],"mem=%lf",&MaxMem); 
    }
  }


#endif

  /* Get start time */
  timeS = time(NULL);
  

  /* Misc */
  InputRowsToUse = InputEndRow - InputStartRow + 1;
  InputColsToUse = InputEndCol - InputStartCol + 1;
  ColsPerProcSegmentCalc = (double)MaxMem * MB2BYTES / 
    ((double)InputRowsToUse * (double)BytesPerValue);

  if ( ColsPerProcSegmentCalc > (double)InputColsToUse ) {
    ColsPerProcSegment = InputColsToUse;
  }
  else if (ColsPerProcSegment < 1) {
    ColsPerProcSegment = 1;
    fprintf(msg,
     "\nWARNING - Cols per proc segment < 1 - reset to unity!\n\n");
  }
  else {
    ColsPerProcSegment = (Int4B)ColsPerProcSegmentCalc;
  }

  ActualMem = (double)(ColsPerProcSegment*InputRowsToUse*BytesPerValue)
                 / MB2BYTES;
  
  NumProcSegments = (Int4B)((double)InputColsToUse/
                            (double)ColsPerProcSegment);
  if (InputColsToUse % ColsPerProcSegment != 0) {
    NumProcSegments++; /* cater for non multiples */
    LastProcSegmentCols = InputColsToUse % ColsPerProcSegment;
  }
  else {
    LastProcSegmentCols = ColsPerProcSegment;
  }

  BytesToSkip = BytesPerValue*(MaxInputCols-ColsPerProcSegment) +
                RowHeaderBytesToSkip + RowFooterBytesToSkip;
  DataRowBytes = ColsPerProcSegment*BytesPerValue;

  /* Messages */
  fprintf(msg,"Input file                         : %s\n",InputFileName);
  fprintf(msg,"Output file                        : %s\n",OutputFileName);
  fprintf(msg,"Max input columns                  : %ld\n",MaxInputCols);
  fprintf(msg,"Bytes per data value               : %ld\n",BytesPerValue);
  fprintf(msg,"Row header / footer bytes to skip  : %ld / %ld\n",
          RowHeaderBytesToSkip,RowFooterBytesToSkip);   
  fprintf(msg,"Processing input rows              : %ld - %ld\n",
          InputStartRow,InputEndRow);
  fprintf(msg,"Processing input columns           : %ld - %ld\n",
          InputStartCol,InputEndCol);
  fprintf(msg,"Memory usage (max / actual)        : %.3f / %.3f MBytes\n",
         MaxMem,ActualMem);
  fprintf(msg,"Cols per proc segment              : %ld\n\n",
          ColsPerProcSegment);

  /* Open files */
  if ( (InFile = fopen (InputFileName, "rb") ) == NULL )
    { fprintf(msg,"ERROR: Input file not found/opened!\n"); exit(1); }
  if ( (OutFile = fopen (OutputFileName, "wb") ) == NULL )
    { fprintf(msg,"ERROR: Output file not opened!\n"); exit(1); }


  /* Check that processing is within confines of data set */
  fseek(InFile,0,SEEK_END);
  InputFileSize = ftell(InFile);
  fseek(InFile,0,SEEK_SET);
  
  MaxInputRows = (Int4B)((double)InputFileSize / 
                 (double)(MaxInputCols*BytesPerValue
                 + RowHeaderBytesToSkip+RowFooterBytesToSkip));

  if ( (InputEndCol-InputStartCol+1) > MaxInputCols ) {
     fprintf (msg,"ERROR - Attempt to process beyond valid data!\n");
     fprintf (msg,"%ld columns would be needed (%ld is the max)\n", 
                  InputEndCol+1,MaxInputCols );
     fprintf (msg,"File size : %ld bytes\n",InputFileSize);             
     exit(1);
  }

  if ( (InputEndRow-InputStartRow+1) > MaxInputRows ) {
     fprintf (msg,"ERROR - Attempt to process beyond valid data!\n");
     fprintf (msg,"%ld rows would be needed (%ld is the max)\n", 
                   InputEndRow+1,MaxInputRows );
     fprintf (msg,"File size : %ld bytes\n",InputFileSize);       
     exit(1);
  } 
  
  
 /****************************************************
   * ALLOCATE SPACE FOR ARRAYS AND ALLOCATE ADDRESSES *
   ****************************************************/

   UnCharArray = (unsigned char **)malloc(sizeof(unsigned char *)
                  *InputRowsToUse);
   UnCharBlock = (unsigned char *)malloc(InputRowsToUse*
                   ColsPerProcSegment*BytesPerValue); 
   InputRow = (unsigned char *)malloc(sizeof(unsigned char)
                  *BytesPerValue*ColsPerProcSegment);
   InputCol = (unsigned char *)malloc(sizeof(unsigned char)
                  *BytesPerValue*InputRowsToUse);
   
   if ( (UnCharArray==NULL)||(UnCharBlock==NULL)||
	(InputRow==NULL)||(InputCol==NULL) )
     {fprintf(msg,"ERROR : in array memory allocation!\n"); exit(1);}
     
   /* Allocate addresses */
   for (row=0; row<InputRowsToUse; row++)
     UnCharArray[row] = &UnCharBlock[row*ColsPerProcSegment*BytesPerValue];


  /**************************
  * REPEAT FOR EACH SEGMENT *
  ***************************/

  for (segment=0; segment<NumProcSegments; segment++)
  {
  fprintf(msg,"Segment %ld...\n",segment); 

  /* cater for straggler columns at end */
  if (segment == NumProcSegments - 1) {
    ColsPerProcSegment = LastProcSegmentCols;
    BytesToSkip = BytesPerValue*(MaxInputCols-ColsPerProcSegment) +
                  RowHeaderBytesToSkip + RowFooterBytesToSkip;
    DataRowBytes = ColsPerProcSegment*BytesPerValue;
  }

  /* Get to start of required data (BytesToSkipAtStart in bytes)*/
  BytesToSkipAtStart = InputStartRow*(BytesPerValue*MaxInputCols
                       + RowHeaderBytesToSkip + RowFooterBytesToSkip) 
                       + RowHeaderBytesToSkip + BytesPerValue*InputStartCol;
  fseek(InFile,BytesToSkipAtStart,0);
     
  /**********************
   * PERFORM EXTRACTION *
   **********************/

   /* Read in data for each row */
   for (row=0; row<InputRowsToUse; row++) {

     /* Read in data for each */
     fread(InputRow,1,DataRowBytes,InFile);
     for (byte=0; byte<DataRowBytes; byte++)
       UnCharArray[row][byte] = InputRow[byte];
      
     /* Skip input data file until start of next data segment */
     fseek(InFile,BytesToSkip,1);

     if (msg==stdout) { 
       if ((row%20)==0)
        fprintf(msg,"Row %ld read\r", row); 
     }

   }  /* End row loop */

   if (msg==stdout) printf("\n");

  /***********************************************
   * WRITE DATA TO OUTPUT FILE IN AZ LINE FORMAT *
   ***********************************************/

  /* Write each column */
  for (column=0; column<ColsPerProcSegment; column++)
    {
    byte = 0;   
    for (row=0; row<InputRowsToUse; row++)
      {
      for (i=0;i<BytesPerValue;i++)
        InputCol[byte++] = UnCharArray[row][column*BytesPerValue+i];
      }

    fwrite(InputCol,1,BytesPerValue*InputRowsToUse,OutFile);

    if (msg==stdout)
      { if((column%20)==0)
           fprintf(msg,"Column %ld written\r",column); }

    }  /* End column loop */

  if (msg==stdout) printf("\n");

   
  /*Increment start column for next segment*/  
  InputStartCol += ColsPerProcSegment;

   
  }  /* End segment (NumProcSegments) loop */

  /* Get end time */
  timeE = time(NULL);
  fprintf(msg,"CORNER TURN done - in %ld secs (%.2f min)\n",
	  timeE-timeS,(double)(timeE-timeS)/60.0);

  /* Tidy up */
  free (UnCharArray[0]);     /* Frees UnCharBlock block */
  free (UnCharArray);
  free (InputRow);
  free (InputCol);

  fclose (InFile);
  fclose (OutFile);

  return(0);  /* on success */
}  /* End g2cor function */


/************************************************************************/
/* FUNCTION : WriteCornerTmplCmdFile() */
Int2B WriteCornerTmplCmdFile(char CmdFileName[])
{
  FILE *OutFile; 
  if ( (OutFile = fopen (CmdFileName, "w") ) == NULL ) return(1);

  fprintf(OutFile,"Command file for Corner (Corner turn program)\n");
  fprintf(OutFile,"$ProgramVersion (jmh) => %s\n",PROG_VERSION);
  fprintf(OutFile,"---------------------------------------------\n\n");
  fprintf(OutFile,"$InputFileName               => tmpg2.rnc\n");
  fprintf(OutFile,"$OutputFileName              => tmpg2.cor\n");
  fprintf(OutFile,"$MaxInputCols [see note]     => 2048\n");
  fprintf(OutFile,"$BytesPerValue [see note]    => 2\n");
  fprintf(OutFile,"$RowHeaderBytesToSkip        => 0\n");
  fprintf(OutFile,"$RowFooterBytesToSkip        => 0\n");
  fprintf(OutFile,"$InputStartRow               => 0\n");
  fprintf(OutFile,"$InputEndRow                 => 1001\n");
  fprintf(OutFile,"$InputStartCol               => 0\n");
  fprintf(OutFile,"$InputEndCol                 => 2048\n");
  fprintf(OutFile,"$MaxMem [in MB - see note]   => 100.0\n");

  fprintf(OutFile,"\nNotes:\n");
  fprintf(OutFile,"------\n\n");
  fprintf(OutFile,
   "Assumes the input file is written row by row. The byte length of a row\n");
  fprintf(OutFile,
   "in the input file is : (BytesPerValue*MaxInputCols) +\n");
  fprintf(OutFile,
   "                       RowHeaderBytesToSkip + RowFooterBytesToSkip.\n\n");
  fprintf(OutFile,
   "MaxInputCols  : Excludes header and footer bytes.\n\n");
  fprintf(OutFile,
   "BytesPerValue : For complex data, this is the byte size of the full\n");
  fprintf(OutFile,
   "                complex number (I or Q together) and must thus be\n");
  fprintf(OutFile,
   "                even and at least 2.\n\n");
  fprintf(OutFile,
   "MaxMem : the maximum memory usage for the corner turn in MBytes. To\n");
  fprintf(OutFile,
   "         reduce mem usage, the input array may be processed in blocks\n");
  fprintf(OutFile,
   "         of input columns. The program respects the start and end rows\n");
  fprintf(OutFile,
   "         and columns regardless of the internal block size used.\n");
 
  fclose(OutFile);
  return(0);  /* on success */
} 
