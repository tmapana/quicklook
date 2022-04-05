/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2b2tif.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Converts binary byte array to TIFF file.
Based on TIFF Revision 6.0 spec. Supports both little-endian (i386)
and big-endian (SUN) TIFF files. 

VERSION/AUTHOR/DATE : 1997-06-13 / Jasper Horrell / 1997-06-13
COMMENTS: 
Initial version.

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS:
Include orientation info in messages. Rename to g2b2tif.c and create 
g2b2tif.h. Allow run as function or standalone. Change width and 
height to rows and cols (order of params swapped). Clean up messages. 
Use Int2B and Int4B.  

VERSION/AUTHOR/DATE : 1999-01-27 / Jasper Horrell / 1999-01-27
COMMENTS:
Minor clean to usage info. Add this header info.
1999-08-03 - cosmetic message change

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/

#include "g2func.h"
#include "g2b2tif.h"

#define PROG_VERSION "1999-01-27"


#if FUNC_B2TIF
Int2B B2Tif (struct B2TifStruct Cmd)
#else
Int2B main(int argc, char *argv[])
#endif
{
  FILE   *InFile, *OutFile,*msg=stdout;
  char BigEndian = 'N',
       ProgVersion[STRING_SPACE]="",
       InputFileName[STRING_SPACE],
       OutputFileName[STRING_SPACE],
       Endian[STRING_SPACE]= {"little"}, /* 'little' (i386) or 'big' (SUN) */
       Header[] = { "B2TIF (Ver. 1999-01-26) - J.M. Horrell   "}; 
                    /* note length header fixed to 41 chars! */
  Int2B  Orient=1, /* init to default value */
         TmpS;  /* Temp 2 byte int */
  Int4B	 i = 0,
         StripsPerImage,  /* changed from short */
         NumStripBytes,
         NumRowsPerStrip,
         DataAddress,
         value,
         TmpL,  /* Temp long int */
         Cols,Rows;    /* changed from short */
  unsigned char *TIFFBuffer;

#if FUNC_B2TIF     

  msg = Cmd.msg;
  strcpy(ProgVersion,Cmd.Version);
  strcpy(InputFileName,Cmd.InputFileName);
  strcpy(OutputFileName,Cmd.OutputFileName);
  Rows = Cmd.Rows;
  Cols = Cmd.Cols;
  strcpy(Endian,Cmd.Endian);
  Orient = Cmd.Orient;

  fprintf(msg,"\nBYTE TO TIFF STAGE...\n");

  /* Check version IDs match */
  if (strcmp(ProgVersion,PROG_VERSION)!=0) {
    fprintf(msg,
      "WARNING - command version (%s) not same as program (%s)!\n\n",
      ProgVersion,PROG_VERSION);
  } 

#else  /* called from command line */

  fprintf(msg,"\n----------\n");         
  fprintf(msg,"Prog B2TIF (Ver. %s) - Code: J.M. Horrell\n",PROG_VERSION); 

  if (argc < 5) {
    fprintf(msg,"Converts byte file to TIFF (8-bit) grayscale image\n\n");
    fprintf(msg,"USAGE: b2tif [InFile] [OutFile] [Rows] [Cols]\n\n");
    fprintf(msg,"Optional params (use at end of line):\n");
    fprintf(msg,
      "  <endian=little/big> - little-endian (i386, the default) or\n");
    fprintf(msg,
      "                        big-endian (SUN) TIFF format\n");
    fprintf(msg,
      "  <orient=n> - TIFF orientation:\n");
    fprintf(msg,"               1 - row 0 top, col 0 left (default)\n");       
    fprintf(msg,"               2 - row 0 top, col 0 right\n");       
    fprintf(msg,"               3 - row 0 bottom, col 0 right\n");       
    fprintf(msg,"               4 - row 0 bottom, col 0 left\n");       
    fprintf(msg,"               5 - row 0 left, col 0 top\n");       
    fprintf(msg,"               6 - row 0 right, col 0 top\n");       
    fprintf(msg,"               7 - row 0 right, col 0 bottom\n");       
    fprintf(msg,"               8 - row 0 left, col 0 bottom\n");       
    fprintf(msg,"e.g. 'b2tif test.in test.tif 34 400 orient=4'\n");
    exit(1);
  }

  strcpy(InputFileName,argv[1]);
  strcpy(OutputFileName,argv[2]);
  sscanf(argv[3],"%ld",&Rows);
  sscanf(argv[4],"%ld",&Cols);

  /* check for optional parameters */
  for (i=5; i<argc; i++) {
    if (strncmp(argv[i],"endian=",7)==0) {
      sscanf(argv[i],"endian=%s",Endian);
    }
    else if (strncmp(argv[i],"orient=",7)==0) {
      sscanf(argv[i],"orient=%hd",&Orient);
    }    
  }

#endif

  /* set big/little endian flag */
  if ( strcmp(Endian,"big")==0 || strcmp(Endian,"BIG")==0 ) {
    BigEndian = 'Y';
  }

  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"Input file       : %s\n",InputFileName);
  fprintf(msg,"Output file      : %s\n",OutputFileName);
  fprintf(msg,"Rows             : %ld\n",Rows);
  fprintf(msg,"Columns          : %ld\n",Cols);
  if (BigEndian == 'Y') {
    fprintf(msg,"TIFF format      : big-endian (SUN)\n");
  }
  else {
    fprintf(msg,"TIFF format      : little-endian (i386)\n");
  }
  if (Orient == 1) {
    fprintf(msg,"Orientation      : 1 - row 0 top, col 0 left\n\n"); 
  }
  else if (Orient == 2) {
    fprintf(msg,"Orientation      : 2 - row 0 top, col 0 right\n\n");
  }
  else if (Orient == 3) {
    fprintf(msg,"Orientation      : 3 - row 0 bottom, col 0 right\n\n");
  }
  else if (Orient == 4) {
    fprintf(msg,"Orientation      : 4 - row 0 bottom, col 0 left\n\n");
  }
  else if (Orient == 5) {
    fprintf(msg,"Orientation      : 5 - row 0 left, col 0 top\n\n");
  }
  else if (Orient == 6) {
    fprintf(msg,"Orientation      : 6 - row 0 right, col 0 top\n\n");
  }
  else if (Orient == 7) {
    fprintf(msg,"Orientation      : 7 - row 0 right, col 0 bottom\n\n"); 
  }
  else if (Orient == 8) {
    fprintf(msg,"Orientation      : 8 - row 0 left, col 0 bottom\n\n"); 
  }          
  else {
    fprintf(msg,"ERROR - unknown orientation %hd!\n",Orient);
    exit(1);
  }

  if(( InFile = fopen(InputFileName,"rb")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open input file %s!\n",
            InputFileName); 
    exit(1); 
  }
  if(( OutFile = fopen(OutputFileName,"wb")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open output file %s!\n",
            OutputFileName); 
    exit(1); 
  }

  NumRowsPerStrip = Rows;
  StripsPerImage = 1;
  NumStripBytes = NumRowsPerStrip * Cols;
  DataAddress = 240L;
  
  /* Allocate mem */
  TIFFBuffer = (unsigned char *)malloc(sizeof(unsigned char)*Cols);
  if (TIFFBuffer == NULL) { 
    printf("Error in array memory allocation!\n"); exit(1); 
  }

  /*******************/
  /* WRITE TIFF FILE */
  /*******************/  

  /* TIFF HEADER */
  if (BigEndian != 'Y')
    fprintf( OutFile, "%c%c", 0x49,0x49);    /* little-endian indicator */
  else 
    fprintf( OutFile, "%c%c", 0x4d,0x4d);    /* big-endian indicator */
  TmpS=42;  fwrite( &TmpS,2,1,OutFile);  /* TIFF file indicator */
  TmpL= 8;  fwrite( &TmpL,4,1,OutFile);  /* first IFD byte offset */
   
  /* IMAGE FILE DIRECTORY (IFD) */
  /* Consists of 2 bytes for the no. of field entries, followed by a */ 
  /* sequence of 12 byte field entries, followed by a 4-byte offset of */
  /* the next IFD (zero if none). Each 12 byte sequence contains the field */
  /* tag (2 bytes), the field type (2 bytes - e.g. 3 = SHORT, 4=LONG), */
  /* the number of values of the */ 
  /* indicated type (4 bytes), and the value offset which points to the */ 
  /* position in the file of the value. If and only if the value fits into */
  /* 4 bytes, the value offset contains the value itself (used here) */  
  TmpS=14;  fwrite( &TmpS,2,1,OutFile);	/* no. of IFD field entries */
  
  TmpS=256; fwrite( &TmpS,2,1,OutFile);  /* ImageWidth Tag */
  TmpS=4;   fwrite( &TmpS,2,1,OutFile);     /* (field type "long") */
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);     /* (no. of values) */
            fwrite( &Cols,4,1,OutFile);    /* (value or value offset) */
  
  TmpS=257; fwrite( &TmpS,2,1,OutFile);  /* ImageLength Tag */
  TmpS=4;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
            fwrite( &Rows,4,1,OutFile); 

  TmpS=258; fwrite( &TmpS,2,1,OutFile);  /* BitsPerSample Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);     /* (field type "short") */
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=8;   fwrite( &TmpS,2,1,OutFile);
  TmpS=0;   fwrite( &TmpS,2,1,OutFile); 

  TmpS=259; fwrite( &TmpS,2,1,OutFile);  /* Compression Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=1;   fwrite( &TmpS,2,1,OutFile); 
  TmpS=0;   fwrite( &TmpS,2,1,OutFile); 

  TmpS=262; fwrite( &TmpS,2,1,OutFile);  /* PhotometricInterpretation Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=1;   fwrite( &TmpS,2,1,OutFile); 
  TmpS=0;   fwrite( &TmpS,2,1,OutFile); 

  TmpS=273; fwrite( &TmpS,2,1,OutFile);  /* StripOffsets Tag */
  TmpS=4;   fwrite( &TmpS,2,1,OutFile);
            fwrite( &StripsPerImage,4,1,OutFile);
            fwrite( &DataAddress,4,1,OutFile);  

  TmpS=274; fwrite( &TmpS,2,1,OutFile);  /* Orientation Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);     /* (field type "short") */
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=Orient; fwrite( &TmpS,2,1,OutFile);
  TmpS=0;   fwrite( &TmpS,2,1,OutFile); 
	    
  TmpS=277; fwrite( &TmpS,2,1,OutFile);  /* SamplesPerPixel Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=1;   fwrite( &TmpS,2,1,OutFile);
  TmpS=0;   fwrite( &TmpS,2,1,OutFile);     

  TmpS=278; fwrite( &TmpS,2,1,OutFile);  /* RowsPerStrip Tag */
  TmpS=4;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
            fwrite( &NumRowsPerStrip,4,1,OutFile);

  TmpS=279; fwrite( &TmpS,2,1,OutFile);  /* StripBytesCounts Tag */
  TmpS=4;   fwrite( &TmpS,2,1,OutFile);
            fwrite( &StripsPerImage,4,1,OutFile);
            fwrite( &NumStripBytes,4,1,OutFile);  
 
  TmpS=282; fwrite( &TmpS,2,1,OutFile);  /* XResolution Tag */
  TmpS=5;   fwrite( &TmpS,2,1,OutFile);     /* (field type "rational") */
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpL=224; fwrite( &TmpL,4,1,OutFile);    /* (address for Xres value) */

  TmpS=283; fwrite( &TmpS,2,1,OutFile);  /* YResolution Tag */
  TmpS=5;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpL=232; fwrite( &TmpL,4,1,OutFile);    /* (address for Yes value) */


  TmpS=296; fwrite( &TmpS,2,1,OutFile);  /* ResolutionUnit Tag */
  TmpS=3;   fwrite( &TmpS,2,1,OutFile);
  TmpL=1;   fwrite( &TmpL,4,1,OutFile);
  TmpS=2;   fwrite( &TmpS,2,1,OutFile);     /* (select inches) */
  TmpS=0;   fwrite( &TmpS,2,1,OutFile);   

  TmpS=305; fwrite( &TmpS,2,1,OutFile);  /* Software Tag */
  TmpS=2;   fwrite( &TmpS,2,1,OutFile);     /* (field type "ASCII") */
  TmpL=42;  fwrite( &TmpL,4,1,OutFile);     /* (length of ASCII header) */
  TmpL=182; fwrite( &TmpL,4,1,OutFile);     /* (address of ASCII header)*/

  TmpL=0;   fwrite( &TmpL,4,1,OutFile);  /* End of IFD Mark */
 
  /* WRITE TAG VALUES (> 4 bytes) & IMAGE DATA */
   
   fseek( OutFile, 182L, 0);             /* Write Software Label */
   fwrite( Header, 42, 1, OutFile);	 

   fseek( OutFile, 224L, 0);             /* Write XRes */
   TmpL=4;  fwrite( &TmpL,4,1,OutFile);
   if (Cols < Rows) { 
     value = Cols; 
   }
   else { 
     value = Rows; 
   }
   fwrite( &value, 4, 1, OutFile);			    

   fseek( OutFile, 232L, 0);             /* Write YRes */
   TmpL=4;  fwrite( &TmpL,4,1,OutFile);   
   fwrite( &value, 4, 1, OutFile);
 
   /* Write image data to TIFF file */
   fseek( OutFile, DataAddress, 0);
   for ( i = 0; i < Rows; i++){     /* Copy rows */
     fread( TIFFBuffer, Cols, 1, InFile);
     fwrite( TIFFBuffer, Cols, 1, OutFile);
     if (i%50 == 0) 
       printf("Row %ld written\r",i);
   }
   fprintf(msg,"\nBYTE TO TIFF DONE!\n");

   /* Tidy up */
   free (TIFFBuffer);
   fclose (InFile);
   fclose (OutFile);

   return(0);
}

