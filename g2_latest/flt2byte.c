/* Program to convert floats in binary file to unsigned chars in binary file.
   Logs overflows in the conversion (clips them), allows setting of a scale
   factor and keeps track of the max scaled float.

   J.Horrell 1998-10-02 : Creation
   Y.Tremeac 1998-11-05 : Add a DC offset to input values. The negative part
                          of the input values can then be taken into account.
                          Prints the min value.
   Y.Tremeac 1998-11-30 : The cast from byte to float is a truncation and
                          introduces a error of 0.5 in average compare to rounding.
                          This error has been fixed by addind a 0.5 offset
                          before casting/truncating.
   J.Horrell 1999-01-22 : (Ver 1.0) minor message cleanup
   J.Horrell 1999-08-03 : cosmetic message change                   
   J.Horrell 2000-02-16 : (Ver 1.1) Allow "math=" option before converting to byte.
                          This allows for raising to a power and/or log operations
                          prior to byte conversion (if both specified, raising to a
                          power performed first). Also change mandatory Scale and Offset
                          usage to "scale=" and "offset=" option.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <math.h>

#define Int4B long int

#define PROG_VERSION "1.1 - 2000-02-16"

int main(int argc, char *argv[])
  {
  FILE *InFile,*OutFile,*msg=stdout;
  float *InFloat;
  unsigned char *OutChar;
  float MaxScaled=0.0, MinScaled=255.0, Offset=0.0, Scale=1.0,tmp,Pow=1.0;
  Int4B Rows,
        Cols,
        Log = 0,
        OverFlows=0,
        i,j;

  printf ("\n-------------\n");
  printf ("Prog Flt2Byte (Ver. %s) - Code J.M.Horrell\n",PROG_VERSION);

  if (argc < 5)
  {
    printf ("Function : Converts floats to unsigned chars in binary files.\n");
    printf ("Tracks max scaled value and number of overflows (auto clipped).\n");
    printf("\nUSAGE: flt2byte [Infile] [OutFile] [Rows] [Colums]\n");
    printf ("\nwhere: \n");
    printf ("\tInFile  - the input file name\n");
    printf ("\tOutFile - the output file name\n");
    printf ("\tRows and Columns - no. of rows and cols of input floats\n");
    printf ("Optional params (use at end of line):\n");
    printf ("\t<math=t> - Math operation prior to scaling and byte conversion:\n");
    printf ("\t           powx - raise values to x power (e.g. pow0.5 for sqrt)\n");
    printf ("\t           log  - natural log of values\n");
    printf ("\t           logpowx - raise to x power, then log (e.g. logpow4.0)\n");
    printf ("\t<scale=n> - Scaling factor before byte conversion (default 1.0)\n");
    printf ("\t<offset=n> - Add a DC offset to the scaled data (e.g. 127.0)\n");
    printf ("e.g. flt2byte in.flt out.byt 2094 23 math=pow0.25 scale=1.8e-02\n"); 
    exit(1);
  }

  sscanf(argv[3], "%ld", &Rows);
  sscanf(argv[4], "%ld", &Cols);

  /* check for optional parameters */
  for (i=5; i<argc; i++) {
    if (strncmp(argv[i],"math=pow",8)==0) {
      sscanf(argv[i],"math=pow%f",&Pow);
    }
    else if (strncmp(argv[i],"math=log",8)==0) {
      Log = 1;
    }    
    if (strncmp(argv[i],"math=logpow",11)==0) {
      sscanf(argv[i],"math=logpow%f",&Pow);
      Log = 1;
    }
    else if (strncmp(argv[i],"scale=",6)==0) {
      sscanf(argv[i],"scale=%f",&Scale);
    }    
    else if (strncmp(argv[i],"offset=",7)==0) {
      sscanf(argv[i],"offset=%f",&Offset);
    }  
  }

  /* Allocate space for array */
  InFloat = (float *)malloc(sizeof(float)*Cols);
  OutChar = (unsigned char *)malloc(sizeof(unsigned char)*Cols);

  if (InFloat==NULL || OutChar==NULL) {
    fprintf(msg,"ERROR - in array memory allocation!\n"); exit(1);
  }

  /* messages */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"Input file   : %s\n",argv[1]);
  fprintf(msg,"Output file  : %s\n",argv[2]);
  fprintf(msg,"Rows         : %ld\n",Rows);
  fprintf(msg,"Columns      : %ld\n",Cols);
  if (Pow!=1.0)
    fprintf(msg,"Raise to pow : %f\n",Pow);
  if (Log)
    fprintf(msg,"Log          : yes\n");
  fprintf(msg,"Scale factor : %.4e\n",Scale);
  fprintf(msg,"DC offset    : %f\n",Offset);
  if (Pow!=1.0 && Log)
  fprintf(msg,"Note: pow performed before log operation\n");

  /* Open files */
  InFile  = fopen(argv[1],"rb");
  if (InFile==NULL) {
    fprintf(msg,"ERROR - Input file not found/opened!\n"); exit(1);
  }
  OutFile = fopen(argv[2],"wb");
  if (OutFile==NULL) {
    fprintf(msg,"ERROR opening output file\n"); exit(1);
  }

  /* Perform processing */
  for (i=0;i<Rows;i++) {

    /* Read in row of data */
    fread(InFloat, Cols, sizeof(float), InFile);

    for (j=0;j<Cols;j++) {
      tmp = InFloat[j];                   
      if (Pow!=1.0)
        tmp = (float)pow((double)tmp,(double)Pow);
      if (Log)
        tmp = (float)log((double)tmp);  
      tmp = tmp*Scale+Offset;  

      /* Check for max scaled value (note doesn't check for large negatives) */
      if (tmp > MaxScaled) {
        MaxScaled = tmp;
      }

      if (tmp < MinScaled) {
        MinScaled = tmp;
      }

      /* Clip output to within byte range */
      if (tmp < 0.0) {
        tmp = 0.0;
        OverFlows++;
      }
      else if (tmp > 255.0) {
        tmp = 255.0;
        OverFlows++;
      }

      /* Assign to output array */
      OutChar[j] = (unsigned char)(tmp+0.5);
      /* the 0.5 is here to compensate truncating effects */

    } /* end for j loop */

    /* Write to output file */
    fwrite(OutChar,Cols,sizeof(unsigned char),OutFile);
  }

  fprintf(msg,"\nOverflows    : %ld\n",OverFlows);
  fprintf(msg,"Max value    : %f (before clipping)\n",MaxScaled);
  fprintf(msg,"Min value    : %f (before clipping)\n",MinScaled);
  fprintf(msg,"Done!\n");

  /* Tidy up */
  free(InFloat);
  free(OutChar);
  fclose(OutFile);
  fclose(InFile);

  return(0);

  }  /* end main function */
