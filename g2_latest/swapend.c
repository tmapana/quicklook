/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: swapend.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Swaps byte order for n-byte numbers in binary file.

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 1997-12-10
COMMENTS: 
Initial version based on rev4byte ver 0.1. Made more efficient
and add ability to cater for numbers other than 4-byte numbers.
Allow also operation in-place (overwrites the file so more dangerous).

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 2000-02-14
COMMENTS: 
Added fflush to fix the in-place operation. This required by ANSI C
when readgin and writing to same file (see fopen man page).

=========================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main (int argc, char *argv[])
{
  FILE *OutFile=NULL,*InFile=NULL,*msg=stdout;
  unsigned char *InArray=NULL,
                *OutArray=NULL;
  register long int line, val, byte, indx, indx2, BytesPerValLess1; 
  long int  Lines=0, LineByteSize=0, ValsPerLine=0, BytesPerVal=0, InPlace=0;
  time_t timeS,timeE; /* Start and end times */

  fprintf(msg,"\n-------\n"); 
  fprintf(msg,"SWAPEND (Ver. 0.1) - Code J.M. Horrell\n");

  if (argc < 6) { 
    fprintf(msg,"Swaps byte orders for n byte numbers in binary file.\n");
    fprintf(msg,"Can also perform in-place (specify OutFile as InFile - note\n");
    fprintf(msg,"results unpredicable, if interrupted before completion).\n");
    fprintf(msg,"\nUSAGE: swapend InFile OutFile Lines ValsPerLine BytesPerValue\n\n");
    fprintf(msg,"e.g. 1: swapend test.in test.out 10 10645 4\n");
    fprintf(msg,"e.g. 2 (in-place) : swapend test.in test.in 10 10645 4\n");
    fprintf(msg,"Note: more efficient with ValsPerLine > 1.\n\n");
    exit(1);
  }

  if (strcmp(argv[1],argv[2])==0)
     InPlace = 1;

  /* Open files and convert other parameters */
  if (!InPlace) {
    if ( (InFile = fopen (argv[1], "rb") ) == NULL ) {
      fprintf (msg,"ERROR: Input file %s not opened!\n",argv[1]);
      exit(-1);
    }
    if ( (OutFile = fopen (argv[2], "wb") ) == NULL ) {
      fprintf (msg,"ERROR: Output file %s not opened!\n",argv[2]);
      exit(-1);
    }
  }
  else {  /* for InPlace operation */
    if ( (InFile = fopen (argv[1], "r+b") ) == NULL ) {
      fprintf (msg,"ERROR: Input/output file %s not opened!\n",argv[1]);
      exit(-1);
    }
  }
    
  sscanf(argv[3], "%ld", &Lines);
  sscanf(argv[4], "%ld", &ValsPerLine);
  sscanf(argv[5], "%ld", &BytesPerVal);

  /* MESSAGES */
  fprintf(msg,"\nInput binary file         : %s\n",argv[1]);
  fprintf(msg,"Output binary file        : %s\n",argv[2]);
  fprintf(msg,"Number of lines to read   : %ld\n",Lines);
  fprintf(msg,"Number of values per line : %ld\n",ValsPerLine);
  fprintf(msg,"Bytes per value           : %ld\n",BytesPerVal);
  if (InPlace) {
    fprintf(msg,"In-place operation        : yes\n");
  }
  else {
    fprintf(msg,"In-place operation        : no\n");
  }


  /* Misc */
  LineByteSize = ValsPerLine*BytesPerVal;
  BytesPerValLess1 = BytesPerVal-1;
  timeS = time(NULL); /* get start time */
 
  /* Allocate memory */
  InArray = (unsigned char *)malloc(sizeof(unsigned char)*LineByteSize);
  OutArray = (unsigned char *)malloc(sizeof(unsigned char)*LineByteSize);
  if (InArray == NULL || OutArray == NULL) {
    fprintf(msg,"ERROR - in array mem allocation!\n");
    exit(-1);
  }

  /* Process */
  for (line=0; line<Lines; line++) {
     fread (InArray,LineByteSize,1,InFile);
     indx = 0;
     for (val=0; val<ValsPerLine; val++) {
       indx2 = indx+BytesPerValLess1;
       for (byte=0; byte<BytesPerVal; byte++) {
         OutArray[indx] = InArray[indx2-byte];
         indx++;
       }
     }
     if (!InPlace) {
       fwrite (OutArray,LineByteSize,1,OutFile);
     }
     else  {  /* if in-place operation */
       fseek(InFile,-LineByteSize,SEEK_CUR);
       fwrite(OutArray,LineByteSize,1,InFile);
       fflush(InFile);
     }

     if (line%500 == 0)
       fprintf (msg,"Line %ld done (%ld bytes)\n",line,ftell(InFile));
  }

  /* Get end time */
  timeE = time(NULL);
  fprintf(msg,"SWAP ENDIAN done - in %ld secs (%.2f min)\n",
	  timeE-timeS,(double)(timeE-timeS)/60.0);

  fclose(InFile);
  if (!InPlace) 
    fclose(OutFile);
  free(InArray);
  free(OutArray);

  return(0); /* for success */
 
 }  /* End main program */
