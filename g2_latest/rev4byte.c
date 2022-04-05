/* Code: jmh     1/9/93   */
/* Ver (1997-12-10) - minor changes */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void main (int argc, char *argv[])
{
  FILE *outfile,*infile,*msg=stdout;
  unsigned char byte1,byte2,byte3,byte4;
  double i,points;

   fprintf(msg,"\nREV4BYTE Ver. 1997-12-10 (Code J.M. Horrell)\n");
   fprintf(msg,
    "Swaps byte orders for 4 byte numbers in binary file (not efficient)\n");

  if (argc < 4)
   { 
   fprintf(msg,"\nUSAGE: rev4byte InFile OutFile Points\n\n");
   fprintf(msg,"e.g. rev4byte test.in test.out 10645\n\n");
   exit(0); 
   }

  /* Open input file */
  fprintf (msg,"Opening file for input I and Q\n");
  if ( (infile = fopen (argv[1], "rb") ) == NULL )
    {
     fprintf (msg,"ERROR: Input file not opened\n");
     exit(0);
    }
  /* Open output file */
  fprintf (msg,"Opening file for output I and Q\n");
  if ( (outfile = fopen (argv[2], "wb") ) == NULL )
    {
     fprintf (msg,"ERROR: Output file not opened\n");
     exit(0);
    }

  sscanf(argv[3], "%lf", &points);

  for (i=0.0; i<points; i++)
    {
     fread (&byte1,sizeof(unsigned char),1,infile);
     fread (&byte2,sizeof(unsigned char),1,infile);
     fread (&byte3,sizeof(unsigned char),1,infile);
     fread (&byte4,sizeof(unsigned char),1,infile);
     fprintf (outfile, "%c%c%c%c",byte4,byte3,byte2,byte1);
     if ((long int)i%10000 == 0) fprintf (msg,"Point %10.0f done\r",i);
    }

  fprintf(msg,"Process complete!\n");

  fclose(infile);
  fclose(outfile);


 }  /* End main program */
