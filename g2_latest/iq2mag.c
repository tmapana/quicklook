/* Program to convert IQ unsigned char input to magnitude or power
   unsigned char output */
/* Revision History : 
   1998-09-30 Rewrote to take command line args (J. Horrell)
   1998-10-14 Introduce separate I- and Q DC offsets
   1998-11-30 Allow operation with float data 
   1998-12-10 Allow three params to be optional plus cosmetic changes.
   1999-08-03 cosmetic change to messages 
   1999-08-06 add StartRow option
*/
                    
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define Int4B long int
#define Int2B short int

int main (int argc, char *argv[])
{
  FILE	*InFile=NULL, *OutFile=NULL,*msg=stdout;
  unsigned char *InIQ=NULL, *OutDat=NULL;
  float *InIQFloat=NULL, *OutDatFloat=NULL;
  Int4B	i, line, Cols, Rows, Power=0,OverFlows=0,Float=0,StartRow=0;
  float Scale=1.0,DC_I=0.0,DC_Q=0.0,tmp;

  printf ("\n-----------\n");
  printf ("Prog IQ2MAG (Ver. 1998-12-10) - Code J.M.Horrell\n");

  if (argc < 7)
    { 
    printf ("Converts IQ input to magnitude or power.\n");
    printf ("Works for unsigned char and float binary data.\n");
    printf(
      "\nUSAGE: \n");
    printf(
     "  iq2mag InFile OutFile Rows Cols mag/pow float/char <optional params>\n");
    printf("\nwhere:\n"); 
    printf("    Rows and Cols are in points\n");
    printf("    mag/pow      - calculate magnitude or power\n");
    printf("    float/char   - floating point or unsigned char input and output\n");
    printf("    DC_I=[n]     - (opt) DC offset of I channel (default 0.0)\n");
    printf("    DC_Q=[m]     - (opt) DC offset of Q channel (default 0.0)\n");
    printf("    StartRow=[n] - (opt) The start row to use (default 0)\n");
    printf("    Scale=[s]    - (opt) mult output scale factor (default 1.0)\n");
    printf("\ne.g. 'iq2mag test.in test.out 3 45 pow char Scale=0.05'\n\n");
    exit(0); 
    }

  /* Read params */
  sscanf(argv[3], "%ld", &Rows);
  sscanf(argv[4], "%ld", &Cols);
  if ( strcmp(argv[5],"pow")==0 ) {Power = 1;}
  if ( strcmp(argv[6],"float")==0 ) {Float = 1;}

  /* Check for optional arguments */
  for (i=7;i<argc;i++) {
    if (strncmp(argv[i],"DC_I=",5)==0) { 
      sscanf(argv[i],"DC_I=%f",&DC_I);
    }
    else if (strncmp(argv[i],"DC_Q=",5)==0) { 
      sscanf(argv[i],"DC_Q=%f",&DC_Q);
    }
    else if (strncmp(argv[i],"Scale=",6)==0) { 
      sscanf(argv[i],"Scale=%f",&Scale); 
    }
    else if (strncmp(argv[i],"StartRow=",9)==0) { 
      sscanf(argv[i],"StartRow=%ld",&StartRow); 
    }
  }


  /* Open input file */
  if ( (InFile = fopen (argv[1], "rb") ) == NULL ) {
    printf ("ERROR: Input file %s not opened\n",argv[1]);
    exit(0);
  }

  /* Open output file */
  if ( (OutFile = fopen (argv[2], "wb") ) == NULL ) {
    printf ("ERROR: Output file %s not opened\n",argv[2]);
    exit(0);
  }


  /* Messages */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"Input / output files : %s / %s\n",argv[1],argv[2]);
  fprintf(msg,"Rows / Cols          : %ld / %ld\n",Rows,Cols);
  fprintf(msg,"Start row            : %ld\n",StartRow);
  if (Power) {
    fprintf(msg,"Magnitude or power   : power\n");
  }
  else {
    fprintf(msg,"Magnitude or power   : magnitude\n");
  }
  fprintf(msg,"DC offset (I/Q)      : %e / %e\n",DC_I,DC_Q);
  fprintf(msg,"Scale factor         : %e\n",Scale);
  if (Float) {
    fprintf(msg,"Data type            : float\n");
  }
  else {
    fprintf(msg,"Data type            : unsigned char\n");
  }


  if (Float) {

    /* PROCESS FOR FLOAT DATA */
    /*------------------------*/

    /* Allocate memory */
    InIQFloat = (float *)malloc(sizeof(float)*Cols*2);
    OutDatFloat = (float *)malloc(sizeof(float)*Cols);

    if (InIQFloat==NULL || OutDatFloat==NULL)
      fprintf(msg,"ERROR - in array memory allocation!\n");

    /* Move to start of data to be processed */
    fseek(InFile,StartRow*Cols*8,SEEK_SET);

     /* Process data */
    for (line = 0; line < Rows; line++) {

      fread (InIQFloat, sizeof(float), Cols*2, InFile);

      for (i = 0; i < Cols; i++) {
         /* Calc power */
         tmp =  (InIQFloat[i*2] - DC_I)*(InIQFloat[i*2] - DC_I)
	       +(InIQFloat[i*2+1] - DC_Q)*(InIQFloat[i*2+1] - DC_Q);

         /* If mag, take sqrt */
         if (!Power) tmp = (float)sqrt((double)tmp);       

         /* Multiply by scale */
         tmp *= Scale;       
         OutDatFloat[i] = tmp;
            
      }  /* end for i loop */

      /* Write to output file */
      fwrite (OutDatFloat,sizeof(float),Cols, OutFile);

    }  /* end for line loop */

    /* Tidy up */
    free(InIQFloat);
    free(OutDatFloat); 

  } /* end if Float */


  else { /* unsigned char data */

    /* PROCESS FOR UNSIGNED CHAR DATA */
    /*--------------------------------*/

    /* Allocate memory */
    InIQ = (unsigned char *)malloc(sizeof(unsigned char)*Cols*2);
    OutDat = (unsigned char *)malloc(sizeof(unsigned char)*Cols);
 
    if (InIQ==NULL || OutDat==NULL) {
      fprintf(msg,"ERROR - in array memory allocation!\n");
      exit(-1);
    }

    /* Move to start of data to be processed */
    fseek(InFile,StartRow*Cols*2,SEEK_SET);

    /* Process data */
    for (line = 0; line < Rows; line++) {

      fread (InIQ, sizeof(unsigned char), Cols*2, InFile);

      for (i = 0; i < Cols; i++) {
         /* Calc power */
         tmp = ((float) InIQ[i*2] - DC_I)*((float) InIQ[i*2] - DC_I)
	       +((float) InIQ[i*2+1] - DC_Q)*((float) InIQ[i*2+1] - DC_Q);

         /* If mag, take sqrt */
         if (!Power) tmp = (float)sqrt((double)tmp);       

         /* Multiply by scale */
         tmp *= Scale;       

         /* Check for overflows */
         if (tmp > 255.0) {
           OutDat[i] = 255;
           OverFlows++;
         }     
         else {
           OutDat[i] = (unsigned char)tmp;   
         }
            
      }  /* end for i loop */

      /* Write to output file */
      fwrite (OutDat,sizeof(unsigned char),Cols, OutFile);

    }  /* end for line loop */

    fprintf(msg,"\nNumber of overflows : %ld\n", OverFlows);

    /* Tidy up */
    free(InIQ);
    free(OutDat);

  }  /* end else unsigned char */
  
  /* Close files */
  fclose(InFile);
  fclose(OutFile);

  fprintf(msg,"\nIQ2MAG done!\n");

  return(0);  /* on success */
  
}
