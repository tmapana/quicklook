#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main (int argc, char *argv[])
{
  FILE *File=NULL;
  unsigned char *Array=NULL;
  long int line, Lines=200, LineByteSize=10000;

  /* Open file for reading and writing */
  if ( (File = fopen ("test.dat", "r+b") ) == NULL ) {
      fprintf (stdout,"ERROR: Input/output file test.dat not opened!\n");
      exit(-1);
  }
    
  /* Allocate memory */
  Array = (unsigned char *)malloc(sizeof(unsigned char)*LineByteSize);
  if (Array == NULL) {
    fprintf(stdout,"ERROR - in array mem allocation!\n");
    exit(-1);
  }

  /* Process */
  for (line=0; line<Lines; line++) {
     fread (Array,1,LineByteSize,File);
     fseek(File,-LineByteSize,SEEK_CUR);
     fwrite(Array,1,LineByteSize,File);
     
     if (line%50 == 0){
       fprintf (stdout,"Line %ld done (%ld bytes)\n",line,ftell(File));
     }         
  }

  fprintf(stdout,"done!!\n");

  fclose(File);
  free(Array);

  return(0); /* for success */
 
 }  /* End main program */
