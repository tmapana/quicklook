/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2imgfu.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Functions to perform ground range projection 
of data and write data to AIRSAR and SUN Raster File formats.

VERSION/AUTHOR/DATE : g2imgf1b / Jasper Horrell / 1997-11-28
COMMENTS: 
Option of making Sun Raster File function a standalone program.

VERSION/AUTHOR/DATE : G2ImgFu / Jasper Horrell / 1997-12-08
COMMENTS:
Rename. Allow power in dB option for SunRas and AIRSAR functions, 
tidy up AIRSAR func code.

VERSION/AUTHOR/DATE : 1997-12-15 / Jasper Horrell  / 1997-12-15
COMMENTS:
Change ground range projection function so that it reads data rng line 
by rng line (i.e. does not need to perform corner turn).

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell  / 1999-01-26
COMMENTS:
Move FUNC_SUN_RAS to g2func.h file.

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/


#include"g2func.h"
#include"g2imgfu.h"   /* includes structure defns */

/* Special defines */
#define ROUND_TOL 1.0e-8

/************************************************************************/
/* FUNCTION: UnBlockAzcFile() */
/* Notes: This function is required when more than one azimuth block is
   processed. In this case, the output azimuth compressed file has been 
   appended with each azimuth block. This function rewrites the file so that 
   the image is correctly written as if processed as one large azimuth block */
Int2B UnBlockAzcFile(FILE *msg,
                     char InFileName[],
                     char OutFileName[],
                     Int4B DetectMethod,
                     Int4B DataType,
                     Int4B RngBins,
                     Int4B PRIsPerAzBlock,
                     Int4B NumAzBlocks )
{
  FILE *InFile,*OutFile;
  time_t TimeS, TimeE;
  Int2B Error=0;
  Int4B X2=1,AzBlock,RngBin,Vals,BytesPerAzBlock;
  unsigned char *Data0=NULL;
  unsigned Int2B *Data1=NULL;
  Int4B *Data2=NULL;
  float *Data3=NULL;
  double *Data4=NULL;

  fprintf(msg,"\nUNBLOCKING AZ COMPRESSED FILE...\n");
  fprintf(msg,"Input / output files           : %s / %s\n",
	  InFileName,OutFileName);  
  
  /* Get start time */
  TimeS  = time(NULL);

  /* Misc */
  if (DetectMethod==0) { X2 = 2; }  

  Vals = PRIsPerAzBlock*X2;

  /* Open files */
  if ( (InFile = fopen (InFileName,"rb")) == NULL )
    { fprintf(msg,"ERROR: Input file %s not found/opened\n",
                 InFileName); return(0); }
  if ( (OutFile = fopen (OutFileName,"wb")) == NULL )
    { fprintf(msg,"ERROR: Ouput file %s not opened\n",
                 OutFileName); return(0); }

  /* Allocate mem, etc */
  if (DataType==0)
    { 
    BytesPerAzBlock = PRIsPerAzBlock*sizeof(unsigned char)*X2;
    Data0 = (unsigned char *)malloc(sizeof(unsigned char)*Vals); 
    if (Data0==NULL) Error=1;
    }
  else if (DataType==1)
    {  
    BytesPerAzBlock = PRIsPerAzBlock*sizeof(unsigned Int2B)*X2;
    Data1 = (unsigned Int2B *)malloc(sizeof(unsigned Int2B)*Vals); 
    if (Data1==NULL) Error=1;
    }
  else if (DataType==2)
    {  
    BytesPerAzBlock = PRIsPerAzBlock*sizeof(Int4B)*X2;
    Data2 = (Int4B *)malloc(sizeof(Int4B)*Vals); 
    if (Data2==NULL) Error=1;
    }
  else if (DataType==3)
    {  
    BytesPerAzBlock = PRIsPerAzBlock*sizeof(float)*X2;
    Data3 = (float *)malloc(sizeof(float)*Vals); 
    if (Data3==NULL) Error=1;
    }
  else if (DataType==4)
    {  
    BytesPerAzBlock = PRIsPerAzBlock*sizeof(double)*X2;
    Data4 = (double *)malloc(sizeof(double)*Vals); 
    if (Data4==NULL) Error=1;
    }
  else 
    {
    fprintf(msg,"ERROR - unknown data type %ld!\n",DataType);
    return(0);
    }  

  if (Error)
    {
    fprintf(msg,"ERROR - in array mem allocation!\n");
    return(0);
    }

  /* Read and write data */
  fseek(InFile,0,0); /* move to start */
  for (RngBin=0;RngBin<RngBins;RngBin++)
    {
    fseek(InFile,RngBin*BytesPerAzBlock,0);  /* move to start posn for rbin */
    for (AzBlock=0; AzBlock<NumAzBlocks; AzBlock++)
      {
      if (DataType==0)
        { 
        fread(Data0,sizeof(unsigned char),Vals,InFile);
        fwrite(Data0,sizeof(unsigned char),Vals,OutFile);
        }
      else if (DataType==1)
        { 
        fread(Data1,sizeof(unsigned Int2B),Vals,InFile);
        fwrite(Data1,sizeof(unsigned Int2B),Vals,OutFile);
        }
      else if (DataType==2)
        { 
        fread(Data2,sizeof(Int4B),Vals,InFile);
        fwrite(Data2,sizeof(Int4B),Vals,OutFile);
        }
      else if (DataType==3)
        { 
        fread(Data3,sizeof(float),Vals,InFile);
        fwrite(Data3,sizeof(float),Vals,OutFile);
        }
      else if (DataType==4)
        { 
        fread(Data4,sizeof(double),Vals,InFile);
        fwrite(Data4,sizeof(double),Vals,OutFile);    
        }
      if (AzBlock != NumAzBlocks-1)
        fseek(InFile,(RngBins-1)*BytesPerAzBlock,1); /* move to next posn */
      }  /* end for AzBlock loop */

    } /* end for RngBin loop  */

  /* Get end time */
  TimeE = time(NULL);
  fprintf(msg,"UNBLOCKING done - in %ld secs (%.2f min)\n",
	  TimeE-TimeS,(double)(TimeE-TimeS)/60.0);  
  
  /* Tidy up */
  fclose(InFile);
  fclose(OutFile);
  if (DataType==0) { free(Data0); }
  else if (DataType==1) { free(Data1); }
  else if (DataType==2) { free(Data2); }
  else if (DataType==3) { free(Data3); }
  else if (DataType==4) { free(Data4); }


  return(1);
}


/************************************************************************/
/* FUNCTION: ProjImageToGrndRng() */
/* Notes: Requires power detected input file written rng line by rng line 
   (i.e. as if out of radar). Output file is also written rng line by rng 
   line. 
   For now, the number of output ground rng bins points is made
   the same as the number of input slant range bins. However, two separate 
   arrays such as OutFloat and InFloat are defined here so as to enable an 
   easier upgrade. The ground range spacing to use is calculated using the 
   ground range spacing at far swath. However, the output rng bins are 
   taken first from near swath out. This implies that some of the far 
   swath is discarded.
   Data types: 0 - unsigned char, 1 - unsigned Int2B, 2 - Int4B,
                     3 - float, 4 - double */

Int2B ProjImageToGrndRng(FILE *msg,
			 char InFileName[],
                         char GRngFileName[],
                         Int4B AzPixels,
                         Int4B RngBins,
                         Int4B DataType,
			 double FirstRngBinSlDist,    /* (m) */
			 double SlRngBinSpacing,      /* (m) */
                         double *RetGrndRngBinSpacing, /* (m) */
			 double Height,               /* (m) */
                         Int4B InterpPts )  /* even */
{
  FILE *InFile,*GRngFile;
  Int4B StartIndx=0,RngBin=0,RngLine,
        n,InterpPtsD2;
  double arg, 
         FirstRngBinGrndDist,
         InterpPosn=0.0,
         sinc=1.0,
         tmpDouble=0.0;
  unsigned char *InUnChar=NULL,
                *OutUnChar=NULL;
  unsigned Int2B *InUnInt2B=NULL,
                 *OutUnInt2B=NULL;
  time_t TimeS, TimeE;
  Int4B *InInt4B=NULL,
        *OutInt4B=NULL;
  float *InFloat=NULL,
        *OutFloat=NULL; 
  double *GrndRngLine=NULL,
         *SlRngLine=NULL,
         *InterpPt=NULL;

  fprintf(msg,"\nGROUND RANGE PROJECTION OF IMAGE...\n");
  fprintf(msg,"Input / output files           : %s / %s\n",
	  InFileName,GRngFileName);
  
  /* Get start time */
  TimeS  = time(NULL);

  /* Checks */
  if (InterpPts % 2 != 0 && InterpPts!=1)
    {
    fprintf(msg,"ERROR - Interpolation kernel must be 1 or even!\n");
    fprintf(msg,"No ground range projection performed!\n");
    return(0);
    }
  if (RngBins < InterpPts)
    {
    fprintf(msg,"WARNING - Fewer rng bins (%ld) than interp kernel (%ld)!\n",
	    RngBins, InterpPts);
    fprintf(msg,"No ground range projection performed.\n");
    return(0);
    }
  if (DataType != 0 && DataType != 1 && DataType !=2 && 
      DataType != 3 && DataType != 4)
    {
    fprintf(msg,"ERROR - Unknown data type %ld!\n",DataType); return(0);
    }
  
  /* Misc */
  InterpPtsD2 = (Int4B)(InterpPts/2);
  FirstRngBinGrndDist = sqrt( pow(FirstRngBinSlDist,2.0)-pow(Height,2.0) );
  *RetGrndRngBinSpacing =
       (sqrt(pow(FirstRngBinSlDist+(RngBins-1)*SlRngBinSpacing,2.0)
	  - pow(Height,2.0) ) )
     - (sqrt(pow(FirstRngBinSlDist+(RngBins-2)*SlRngBinSpacing,2.0)
	  - pow(Height,2.0) ) ); 
  fprintf(msg,"Slant / ground rng bin spacing : %.4f / %.4f m\n",
	  SlRngBinSpacing,*RetGrndRngBinSpacing);
  fprintf(msg,"Input az width / rng height    : %ld / %ld pixels\n",
	  AzPixels,RngBins); 
  fprintf(msg,"Output rng width / az height   : %ld / %ld pixels\n",
	  RngBins,AzPixels); 
  if (InterpPts == 1)
    {fprintf(msg,"Interpolation kernel           : nearest neighbour\n");}
  else
    {fprintf(msg,"Interpolation kernel           : %ld pts\n",InterpPts);}

  
  /* Open input file (slant range, power image) */
  if ( (InFile = fopen (InFileName,"rb")) == NULL )
    { fprintf(msg,"ERROR: Slant rng image file %s not found/opened\n",
                 InFileName); return(0); }

  /* Open output ground rng image file */
  if ( (GRngFile = fopen (GRngFileName,"wb")) == NULL )
    { fprintf(msg,"ERROR: Ground rng image file %s not opened\n",
                 GRngFileName); return(0); }

  /* Allocate array memory */
  SlRngLine = (double *)malloc(sizeof(double)*RngBins);
  GrndRngLine = (double *)malloc(sizeof(double)*RngBins);
  InterpPt = (double *)malloc(sizeof(double)*RngBins);
  if (SlRngLine==NULL || GrndRngLine==NULL || InterpPt==NULL)
    {
    fprintf (msg,"ERROR - in grnd rng projection array mem alloc!\n");
    return(0);
    }

  if (DataType == 0)
    { 
    InUnChar  = (unsigned char *)malloc(sizeof(unsigned char)*RngBins);
    OutUnChar = (unsigned char *)malloc(sizeof(unsigned char)*RngBins);
    }
  else if (DataType == 1)
    { 
    InUnInt2B = (unsigned Int2B *)malloc(sizeof(unsigned Int2B)*RngBins);
    OutUnInt2B = (unsigned Int2B *)malloc(sizeof(unsigned Int2B)*RngBins);
    }
  else if (DataType == 2)
    { 
    InInt4B = (Int4B *)malloc(sizeof(Int4B)*RngBins);
    OutInt4B = (Int4B *)malloc(sizeof(Int4B)*RngBins);
    }
  else if (DataType == 3) /* note - mem already allocated for data type 4 */
    { 
    InFloat = (float *)malloc(sizeof(float)*RngBins);
    OutFloat = (float *)malloc(sizeof(float)*RngBins);
    }

  /* Calc interp positions required */
  for (RngBin=0; RngBin<RngBins;RngBin++)
    {
    InterpPt[RngBin] =
      (sqrt(pow(FirstRngBinGrndDist+RngBin*(*RetGrndRngBinSpacing),2.0)
           + pow(Height,2.0) ) - FirstRngBinSlDist ) / SlRngBinSpacing;
    /* fprintf(msg,"InterpPt[%ld] = %f\n",RngBin,InterpPt[RngBin]); */
    }  
  fprintf(msg,"Sl rbins used in grnd rng proj : %ld - %ld\n",
	  (Int4B)InterpPt[0],(Int4B)(InterpPt[RngBins-1]+0.5));

  /*------------------------------------*/   
  /* REPEAT FOR EACH RANGE LINE IN FILE */  
  for (RngLine = 0; RngLine<AzPixels; RngLine++)
    { 
    /* Read in range line from image file and convert to double */
    if (DataType == 0)
      { 
      fread(InUnChar,sizeof(unsigned char),RngBins,InFile);
      for (RngBin = 0; RngBin<RngBins; RngBin++)       
        SlRngLine[RngBin] = (double)InUnChar[RngBin];
      }
    else if (DataType == 1)
      { 
      fread(InUnInt2B,sizeof(unsigned Int2B),RngBins,InFile);
      for (RngBin = 0; RngBin<RngBins; RngBin++)
        SlRngLine[RngBin] = (double)InUnInt2B[RngBin];
      }
    else if (DataType == 2)
      { 
      fread(InInt4B,sizeof(Int4B),RngBins,InFile); 
      for (RngBin = 0; RngBin<RngBins; RngBin++)
        SlRngLine[RngBin] = (double)InInt4B[RngBin];
      }
    else if (DataType == 3)
      { 
      fread(InFloat,sizeof(float),RngBins,InFile); 
      for (RngBin = 0; RngBin<RngBins; RngBin++)
        SlRngLine[RngBin] = (double)InFloat[RngBin]; 
      }
    else 
      { 
      for (RngBin = 0; RngBin<RngBins; RngBin++)
        fread(&SlRngLine[RngBin],sizeof(double),RngBins,InFile);
      }
         
    /* Project to ground range sampling */
    for (RngBin=0;RngBin<RngBins;RngBin++)  /* output rng bin */
      {
      InterpPosn = InterpPt[RngBin];
      if (InterpPts == 1)  /* nearest neighbour */
        { GrndRngLine[RngBin] = SlRngLine[(Int4B)(InterpPosn+0.5)];} 
      else if ( modf(InterpPosn+ROUND_TOL,&tmpDouble) < 2.0*ROUND_TOL )
        { GrndRngLine[RngBin]
	    = SlRngLine[(Int4B)(InterpPosn+ROUND_TOL)]; }
      else /* Shannon interpolation */
	{
        StartIndx = (Int4B)(InterpPosn+1.0) - InterpPtsD2;
	if (StartIndx < 0)
	  { StartIndx = 0; }
	if (StartIndx+InterpPts>RngBins)
	  { StartIndx = RngBins - InterpPts; }
        GrndRngLine[RngBin] = 0.0;  /* initialise */
 	for (n=StartIndx;n<StartIndx+InterpPts;n++)
          {
          arg = PI*(InterpPosn-(double)n);
          if (arg < ROUND_TOL && arg > -ROUND_TOL) 
            { sinc = 1.0; }
          else 
            { sinc = sin(arg)/arg; }
          GrndRngLine[RngBin] += SlRngLine[n]*sinc;
          }
	}  /* end else Shannon */
      } /* end for RngBin */


    /* Write output line (note - written rng line by rng line), convert back
	 to orig data type */  
    if (DataType == 0)
      {
      for (RngBin = 0; RngBin<RngBins; RngBin++)
	{ 
        if (GrndRngLine[RngBin] > 255.0)  /* check not out of range */
          { OutUnChar[RngBin] = 255; }
        else if (GrndRngLine[RngBin] < 0.0)
          { OutUnChar[RngBin] = 0; }  
        else
          { OutUnChar[RngBin] = (unsigned char)GrndRngLine[RngBin]; }
        }
      fwrite(OutUnChar,sizeof(unsigned char),RngBins,GRngFile);
      }
    else if (DataType == 1)
      {
      for (RngBin = 0; RngBin<RngBins; RngBin++)
	{ 
        if (GrndRngLine[RngBin] > 65535.0) /* check in range */
          { OutUnInt2B[RngBin] = 65535; }
        else if (GrndRngLine[RngBin] < 0.0)
          { OutUnInt2B[RngBin] = 0; }   
        else
          { OutUnInt2B[RngBin] = (unsigned Int2B)GrndRngLine[RngBin]; }
        }
      fwrite(OutUnInt2B,sizeof(unsigned Int2B),RngBins,GRngFile);
      }
    else if (DataType == 2)
      {
      for (RngBin = 0; RngBin<RngBins; RngBin++)
	{ OutInt4B[RngBin] = (Int4B)GrndRngLine[RngBin]; }
      fwrite(OutInt4B,sizeof(Int4B),RngBins,GRngFile);
      }
    else if (DataType == 3)
      {
      for (RngBin = 0; RngBin<RngBins; RngBin++)
	{ OutFloat[RngBin] = (float)GrndRngLine[RngBin]; }
      fwrite(OutFloat,sizeof(float),RngBins,GRngFile);
      }
    else 
      {
      fwrite(GrndRngLine,sizeof(double),RngBins,GRngFile);
      }

    }  /* end for RngLine */
  /*----------------------*/


  /* Tidy up */
  if (DataType == 0)
    { free(InUnChar); free(OutUnChar); }
  else if (DataType == 1)
    { free(InUnInt2B); free(OutUnInt2B); }
  else if (DataType == 2)
    { free(InInt4B); free(OutInt4B); } 
  else if (DataType == 3)
    { free(InFloat); free(OutFloat); }

  free (SlRngLine);
  free (GrndRngLine);
  free (InterpPt);
  fclose(InFile);
  fclose(GRngFile);

  /* Get end time */
  TimeE = time(NULL);
  fprintf(msg,"GRND RANGE PROJECTION done - in %ld secs (%.2f min)\n",
	  TimeE-TimeS,(double)(TimeE-TimeS)/60.0);  
  
  return(1); /* to indicate successful operation */
} /* end function ProjImageToGrndRng */


/**************************************************************/
/* Function/Prog WriteSUNRasterFile() */
/* Writes image data to SUN Raster File Format
    (format ver.  DLR GENESIS Sun Raster File Format, unless marked
    below as "RRSG addition"). Note that the SUN raster file is
    written as if on a Sun machine (big-endian), but the input file is
    assumed to have been written on the same machine as running this
    function (re: byte order).  The width parameter is the length of 
    a row in pixels (rows read from input file). The rasterfile standard 
    requires that the bytes per row be even (rows padded here by one byte,
    if necessary).
    Input data types: 0 - unsigned char, 1 - unsigned Int2B,
                      2 - Int4B, 3 - float, 4 - double
   DetectMethod: 0 - complex, 1 - magnitude, 2 - power, 3 - power in dB */
#if FUNC_SUN_RASTER
Int2B WriteSUNRasterFile(FILE *msg,
			 char InFileName[],
			 char RasFileName[],
			 Int4B DataType,
			 Int4B DetectMethod,
			 Int4B Width,  /* pixels */
			 Int4B Height) /* pixels */
#else
Int2B main (int argc,char *argv[])
#endif
{
#if !FUNC_SUN_RASTER
  FILE *msg=NULL;
  char InFileName[STRING_SPACE],RasFileName[STRING_SPACE];
  Int4B DataType=0,DetectMethod=2;
  Int4B Width=0,Height=0;  /* pixels */
#endif

  FILE *InFile=NULL,*RasFile=NULL;
  Int4B X2=1,PixelBitSize=0,SRFType=0,Row,ValWidth,
        Val,indx;
  Int2B BigEndian=0;
  time_t TimeE, TimeS;
  
  unsigned char *DataUnChar=NULL,*ptr=NULL;
  unsigned Int2B *DataUnInt2B=NULL;
  Int4B *DataInt4B=NULL;
  float *DataFloat=NULL;
  double *DataDouble=NULL;

#if FUNC_SUN_RASTER  
  fprintf(msg,"\nWRITING IMAGE TO SUN RASTERFILE FORMAT...\n");
#else
  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\nProg: SUNRAS (Ver. 1997-12-08)\n");
  fprintf(msg,"Code: J.M. Horrell (Copyright UCT)\n");

  if (argc<7)
  {
  fprintf(msg,"Writes data to Sun Raster File format.\n\n");
  fprintf(msg,
  "USAGE : sunras [InFile] [OutFile] [Type] [Detection] [Width] [Height]");
  fprintf(msg,"\n\nType (data type) : 0 - unsigned char\n");
  fprintf(msg,"                 : 1 - unsigned short integer (2 bytes)\n");
  fprintf(msg,"                 : 2 - integer (4 bytes)\n");
  fprintf(msg,"                 : 3 - float (4 bytes)\n");
  fprintf(msg,"                 : 4 - double (8 bytes)\n");
  fprintf(msg,"Detection : 0 - none (complex data)\n");
  fprintf(msg,"          : non zero - mag, power, etc (real data)\n");
  fprintf(msg,"Width and height : in pixels\n\n"); 
  exit(1);
  }

  fprintf(msg,"\n");

  strcpy(InFileName,argv[1]);
  strcpy(RasFileName,argv[2]);
  sscanf(argv[3],"%ld",&DataType);
  sscanf(argv[4],"%ld",&DetectMethod);  
  sscanf(argv[5],"%ld",&Width);
  sscanf(argv[6],"%ld",&Height);

#endif

  fprintf(msg,"Input / output files            : %s / %s\n",
	  InFileName,RasFileName);
  fprintf(msg,"Image width / height            : %ld / %ld pixels\n",
	  Width,Height);
  if (DetectMethod==0)
    { fprintf(msg,"Image type                      : complex\n");}
  else
    { fprintf(msg,"Image type                      : real\n");}    

  /* Get start time */
  TimeS = time(NULL);
  
  /* Misc */
  BigEndian = HostBigEndian();  /* returns 1 if host big-endian (Sun) */
  
  if (DetectMethod == 0)
    {
    X2 = 2;  /* scale arrays by factor of two for complex data */
    switch(DataType)
      {
      case 0:   /* unsigned char */
        fprintf(msg,"Data type                       : unsigned char\n");
	PixelBitSize = 16;
	SRFType = 88;     /* RRSG addition */
        break;
      case 1:   /* unsigned Int2B */
        fprintf(msg,
	   "Data type                       : unsigned 2 byte int\n");
        PixelBitSize = 32;
	SRFType = 97;     
        break;
      case 2:   /* Int4B */
        fprintf(msg,"Data type                       : 4 byte int\n");
        PixelBitSize = 64;
	SRFType = 87;     /* RRSG addition */
        break;
      case 3:   /* float */
        fprintf(msg,"Data type                       : float\n");
	PixelBitSize = 64;
	SRFType = 96;
        break;
      case 4:   /* double */
        fprintf(msg,"Data type                       : double\n");
	PixelBitSize = 128;
	SRFType = 95;
        break; 
      } /* end switch DataType */
    } /* end if detect method complex */
  else /* data not complex */
    {

    switch(DataType)
      {
      case 0:   /* unsigned char */
        PixelBitSize = 8;
	SRFType = 1;
        break;
      case 1:   /* unsigned Int2B */
        PixelBitSize = 16;
	SRFType = 1;     
        break;
      case 2:   /* Int4B */
        PixelBitSize = 32;
	SRFType = 89;     /* RRSG addition */
        break;
      case 3:   /* float */
        PixelBitSize = 32;
	SRFType = 99;
        break;
      case 4:   /* double */
        PixelBitSize = 64;
	SRFType = 98;
        break; 
      } /* end switch DataType */
    }
  ValWidth = Width*X2; /* number of values per row (I and Q
			    separate) */

  fprintf(msg,"Bits per pixel                  : %ld\n",PixelBitSize);
  fprintf(msg,"Sun rasterfile type             : %ld\n",SRFType);
  
  if (PixelBitSize != 8 || SRFType != 1)
    fprintf(msg,
      "WARNING -  Rasterfile not standard type (may be of limited use)!\n");
  
  /* Open files */
  if ( (InFile = fopen (InFileName,"rb")) == NULL )
    { fprintf(msg,"ERROR: input image file %s not found/opened!\n",
                 InFileName); return(0); }
  if ( (RasFile = fopen (RasFileName,"wb")) == NULL )
    { fprintf(msg,"ERROR: output image file %s not opened!\n",
                 RasFileName); return(0); }

  /* Allocate mem */

  switch(DataType)
    {
    case 0:
      DataUnChar = (unsigned char *)malloc(sizeof(unsigned char)*ValWidth);
      break;
    case 1:
      DataUnInt2B = (unsigned Int2B *)malloc(sizeof(unsigned Int2B)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(unsigned Int2B)*ValWidth);	
      break;
    case 2:
      DataInt4B = (Int4B *)malloc(sizeof(Int4B)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(Int4B)*ValWidth);	
      break;
    case 3:
      DataFloat = (float *)malloc(sizeof(float)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(float)*ValWidth);	
      break;
    case 4:
      DataDouble = (double *)malloc(sizeof(double)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(double)*ValWidth);	
      break;
    default:
      fprintf(msg,"ERROR - Unknown data type %ld!\n",DataType);
      return(0);
    } /* end switch DataType */

  if (DataUnChar==NULL && DataUnInt2B==NULL && DataInt4B==NULL
      && DataFloat==NULL && DataDouble==NULL)
    {
    fprintf (msg,"ERROR - in raster file conversion array mem alloc!\n");
    return(0);
    }

  /* Write SUN raster header bytes */
  BigEndWriteInt4B(RasFile,0x59a66a95);
  BigEndWriteInt4B(RasFile,Width);
  BigEndWriteInt4B(RasFile,Height);
  BigEndWriteInt4B(RasFile,PixelBitSize);
  BigEndWriteInt4B(RasFile,(Int4B)(Width*Height*PixelBitSize/8));
  BigEndWriteInt4B(RasFile,SRFType);  
  BigEndWriteInt4B(RasFile,0);  /* colour map type */
  BigEndWriteInt4B(RasFile,0);  /* coulour map size  */
  

  /* Repeat for each row */
  for (Row=0;Row<Height;Row++)
    {
    indx = 0;
    /* Read and write image row */
    switch(DataType)
      {
      case 0:
        fread(DataUnChar,sizeof(unsigned char),ValWidth,InFile);
	fwrite(DataUnChar,sizeof(unsigned char),ValWidth,RasFile);
        if (ValWidth % 2 != 0)
          { fprintf(RasFile,"%c",0); } /* pad out to nearest 16 bits */
        break;
      case 1:
        fread(DataUnInt2B,sizeof(unsigned Int2B),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataUnInt2B,sizeof(unsigned Int2B),ValWidth,RasFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataUnInt2B[Val];
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;
	    }  
	  fwrite(DataUnChar,1,sizeof(unsigned Int2B)*ValWidth,RasFile);
	  }  /* end else */
        break;
      case 2:
        fread(DataInt4B,sizeof(Int4B),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataInt4B,sizeof(Int4B),ValWidth,RasFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataInt4B[Val];
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(Int4B)*ValWidth,RasFile);
	  }  /* end else */  
        break;
      case 3:
        fread(DataFloat,sizeof(float),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataFloat,sizeof(float),ValWidth,RasFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataFloat[Val];
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(float)*ValWidth,RasFile);
	  }  /* end else */  
        break;
      case 4:
        fread(DataDouble,sizeof(double),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataDouble,sizeof(double),ValWidth,RasFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataDouble[Val];
	    DataUnChar[indx++] = *(ptr+7);
	    DataUnChar[indx++] = *(ptr+6);
	    DataUnChar[indx++] = *(ptr+5);	    
	    DataUnChar[indx++] = *(ptr+4);
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(double)*ValWidth,RasFile);
	  }  /* end else */  
        break;
      }  /* end switch data type */
    
    }  /* end for row */
    
  /* Tidy up */
  fclose (InFile);
  fclose (RasFile);

  switch(DataType)
    {
    case 0:
      free (DataUnChar);
      break;
    case 1:
      free (DataUnInt2B);
      if (!BigEndian) free(DataUnChar);	
      break;
    case 2:
      free (DataInt4B);
      if (!BigEndian) free(DataUnChar);
      break;
    case 3:
      free(DataFloat);
      if (!BigEndian) free(DataUnChar);	
      break;
    case 4:
      free(DataDouble);
      if (!BigEndian) free(DataUnChar);	
      break;
    } /* end switch DataType */

  /* Get end time */
  TimeE = time(NULL);
  fprintf(msg,"SUN RASTER FILE written - in %ld secs (%.2f min)\n",
	  TimeE-TimeS,(double)(TimeE-TimeS)/60.0);  
  
  return(1);
}
  

/**************************************************************/
/* Function WriteAIRSARFormat() */
/* Writes image data to AIRSAR Format
   (format ver. 0.01 - July 17 1996). 
   For now, only single polarization without using compressed Stokes
   matrix format. Big-endian (Sun) byte order. All records must be 
   of equal length (longest header 5000 bytes - parameter header).
   At present only new header and parameter header used.
   Input data types: 0 - unsigned char, 1 - unsigned Int2B,
                     2 - Int4B, 3 - float, 4 - double
   DetectMethod: 0 - complex, 1 - magnitude, 2 - power, 3 - power in dB */
Int2B WriteAIRSARFormat(FILE *msg,
			char InFileName[],
			char AIRFileName[],
                        Int4B DataType,
			Int4B DetectMethod,
			Int4B WidthPixels,
			Int4B HeightPixels,
		       	struct AIRSARHead1Struct *Head1,
		        struct AIRSARParamHeadStruct *Param)
  
{

  FILE *InFile, *AIRFile;
  time_t TimeE,TimeS;
  char tmps[STRING_SPACE];
  Int2B BigEndian;
  Int4B X2=1,Val,ValWidth,Row,PixelBitSize=0,indx,RecordByteLen=0,
        DataByteWidth,i,HEADER_RECORDS=2;

  unsigned char *DataUnChar=NULL,*ptr=NULL;
  unsigned Int2B *DataUnInt2B=NULL;
  Int4B *DataInt4B=NULL;
  float *DataFloat=NULL;
  double *DataDouble=NULL;

 
  fprintf(msg,"\nWRITING IMAGE TO AIRSAR FORMAT...\n");
  fprintf(msg,"Input / output files            : %s / %s\n",
	  InFileName,AIRFileName);
  fprintf(msg,"Image width / height            : %ld / %ld pixels\n",
	  WidthPixels,HeightPixels);
  if (DetectMethod==0)
    { fprintf(msg,"Image type                      : complex\n");}
  else if (DetectMethod==1)
    { fprintf(msg,"Image type                      : magnitude\n");}    
  else if (DetectMethod==2)
    { fprintf(msg,"Image type                      : power\n");}
  else if (DetectMethod==3)
    { fprintf(msg,"Image type                      : power in dB\n");}
  else
    { fprintf(msg,"Image type                      : unknown!!\n");}


  /* Get start time */
  TimeS = time(NULL);
  
  /* Misc */
  BigEndian = HostBigEndian();  /* returns 1 if host big-endian (Sun) */
 
  if (DetectMethod == 0)
    { 
    strcpy(Head1->DataType,"COMPLEX_");
    strcpy(Param->DataType,"COMPLEX_");
    }
  else if (DetectMethod == 1)
    {
    strcpy(Head1->DataType,"MAG_");
    strcpy(Param->DataType,"MAG_");
    }
  else if (DetectMethod == 2)
    {
    strcpy(Head1->DataType,"POWER_");
    strcpy(Param->DataType,"POWER_");
    }
  else if (DetectMethod == 3)
    {
    strcpy(Head1->DataType,"POWER_DB_");
    strcpy(Param->DataType,"POWER_DB_");
    }
  else
    {
    strcpy(Head1->DataType,"UNKNOWN_");
    strcpy(Param->DataType,"UNKNOWN_");
    }

  switch(DataType)
    {
    case 0:   /* unsigned char */
      PixelBitSize = 8;
      strcat(Head1->DataType,"UNCHAR");
      strcat(Param->DataType,"UNCHAR");
      fprintf(msg,
        "Data type                       : unsigned char\n");
      break;
    case 1:   /* unsigned Int2B */
      PixelBitSize = 16; 
      strcat(Head1->DataType,"UNINTEGER*2");
      strcat(Param->DataType,"UNINTEGER*2");
      fprintf(msg,
	"Data type                       : unsigned 2 byte int\n");
      break;
    case 2:   /* Int4B */
      PixelBitSize = 32;
      strcat(Head1->DataType,"INTEGER*4");
      strcat(Param->DataType,"INTEGER*4");
      fprintf(msg,
        "Data type                       : 4 byte int\n");
      break;
    case 3:   /* float */
      PixelBitSize = 32;
      strcat(Head1->DataType,"FLOAT");
      strcat(Param->DataType,"FLOAT");
       fprintf(msg,
        "Data type                       : float\n");
      break;
    case 4:   /* double */
      PixelBitSize = 64;
      strcat(Head1->DataType,"DOUBLE");
      strcat(Param->DataType,"DOUBLE"); 
      fprintf(msg,
        "Data type                       : double\n");
      break; 
    } /* end switch DataType */

  
  if (DetectMethod == 0)  /* complex data */
    {
    X2 = 2;  /* scale arrays by factor of two for complex data */
    switch(DataType)
      {
      case 0:   /* unsigned char */
        fprintf(msg,"DC offset                       : 127\n");
        break;
      case 1:   /* unsigned Int2B */
        fprintf(msg,"DC offset                       : 32767\n");
        break;
      case 2:   /* Int4B */
      case 3:   /* float */
      case 4:   /* double */
        fprintf(msg,"DC offset                       : none\n");
        break; 
      } /* end switch DataType */
    } /* end if detect method complex */

  PixelBitSize = PixelBitSize*X2;
  ValWidth = WidthPixels*X2; /* number of values per row (I and Q
		  	        separate) */
  DataByteWidth = (Int4B)(PixelBitSize/8)*WidthPixels;
  RecordByteLen = DataByteWidth;
  if (RecordByteLen < 5000) 
    RecordByteLen = 5000; /* ensure all records equal length */
  Head1->RecordLen = RecordByteLen;
  Head1->NumBytesPerSample = (Int4B)(PixelBitSize/8);
  Head1->ByteOffset1stDataRecord = HEADER_RECORDS*RecordByteLen;
  Head1->ByteOffsetParamHeader = RecordByteLen;
 
  fprintf(msg,"Bits per pixel                  : %ld\n",PixelBitSize);
  fprintf(msg,"Record length                   : %ld bytes\n",
          RecordByteLen);  

  /* Open files */
  if ( (InFile = fopen (InFileName,"rb")) == NULL )
    { fprintf(msg,"ERROR: input image file %s not found/opened\n",
                 InFileName); return(0); }
  if ( (AIRFile = fopen (AIRFileName,"wb")) == NULL )
    { fprintf(msg,"ERROR: output image file %s not opened\n",
                 AIRFileName); return(0); }


  /* Allocate mem */
  switch(DataType)
    {
    case 0:
      DataUnChar = (unsigned char *)malloc(sizeof(unsigned char)*ValWidth);
      break;
    case 1:
      DataUnInt2B = (unsigned Int2B *)malloc(sizeof(unsigned Int2B)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(unsigned Int2B)*ValWidth);	
      break;
    case 2:
      DataInt4B = (Int4B *)malloc(sizeof(Int4B)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(Int4B)*ValWidth);	
      break;
    case 3:
      DataFloat = (float *)malloc(sizeof(float)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(float)*ValWidth);	
      break;
    case 4:
      DataDouble = (double *)malloc(sizeof(double)*ValWidth);
      if (!BigEndian)
        DataUnChar =
	  (unsigned char *)malloc(sizeof(double)*ValWidth);	
      break;
    default:
      fprintf(msg,"ERROR - Unknown data type %ld!\n",DataType);
      return(0);
    } /* end switch DataType */

  if (DataUnChar==NULL && DataUnInt2B==NULL && DataInt4B==NULL
      && DataFloat==NULL && DataDouble==NULL)
    {
    fprintf (msg,"ERROR - in AIRSAR file conversion array mem alloc!\n");
    return(0);
    }


  /* WRITE AIRSAR IMAGE ASCII HEADERS (probably best to make all
     fields check for the default values - update later!)*/

  /* Write first header */
  sprintf(tmps,"RECORD LENGTH IN BYTES =");
  fprintf(AIRFile,"%-.24s",tmps);
  sprintf(tmps,"%ld",Head1->RecordLen);
  fprintf(AIRFile,"%26s",tmps);
  sprintf(tmps,"NUMBER OF HEADER RECORDS =");
  fprintf(AIRFile,"%-.26s",tmps);
  sprintf(tmps,"%ld",Head1->NumHeaderRecords);
  fprintf(AIRFile,"%24s",tmps);  
  sprintf(tmps,"NUMBER OF SAMPLES PER RECORD =");
  fprintf(AIRFile,"%-.30s",tmps);
  sprintf(tmps,"%ld",Head1->SamplesPerRecord);
  fprintf(AIRFile,"%20s",tmps);  
  sprintf(tmps,"NUMBER OF LINES IN IMAGE =");
  fprintf(AIRFile,"%-.26s",tmps);
  sprintf(tmps,"%ld",Head1->LinesInImage);
  fprintf(AIRFile,"%24s",tmps);  
  sprintf(tmps,"NUMBER OF BYTES PER SAMPLE ="); /* field 5 */
  fprintf(AIRFile,"%-.28s",tmps);
  sprintf(tmps,"%ld",Head1->NumBytesPerSample);
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"SASAR PROCESSOR VERSION =");
  fprintf(AIRFile,"%-.25s",tmps);
  sprintf(tmps,"%s",Head1->SASARProcVersion);
  fprintf(AIRFile,"%25s",tmps);  
  sprintf(tmps,"DATA TYPE =");
  fprintf(AIRFile,"%-.11s",tmps);
  sprintf(tmps,"%s",Head1->DataType);
  fprintf(AIRFile,"%39s",tmps);
  sprintf(tmps,"RANGE PROJECTION =");
  fprintf(AIRFile,"%-.18s",tmps);
  sprintf(tmps,"%s",Head1->RangeProjection);
  fprintf(AIRFile,"%32s",tmps);
  sprintf(tmps,"RANGE PIXEL SPACING (METRES) =");
  fprintf(AIRFile,"%-.30s",tmps);
  sprintf(tmps,"%.4f",Head1->RangePixelSpacing);
  fprintf(AIRFile,"%20s",tmps);
  sprintf(tmps,"AZIMUTH PIXEL SPACING (METRES) ="); /* field 10 */
  fprintf(AIRFile,"%-.32s",tmps);
  sprintf(tmps,"%.4f",Head1->AzPixelSpacing);
  fprintf(AIRFile,"%18s",tmps);  
  sprintf(tmps,"BYTE OFFSET OF OLD HEADER =");
  fprintf(AIRFile,"%-.27s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffsetOldHeader);
  fprintf(AIRFile,"%23s",tmps);    
  sprintf(tmps,"BYTE OFFSET OF USER HEADER =");
  fprintf(AIRFile,"%-.28s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffsetUserHeader);
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"BYTE OFFSET OF FIRST DATA RECORD =");
  fprintf(AIRFile,"%-.34s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffset1stDataRecord);
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"BYTE OFFSET OF PARAMETER HEADER =");
  fprintf(AIRFile,"%-.33s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffsetParamHeader);
  fprintf(AIRFile,"%17s",tmps);
  sprintf(tmps,"LINE FORMAT OF DATA =");  /* field 15 */
  fprintf(AIRFile,"%-.21s",tmps);
  sprintf(tmps,"%s",Head1->DataLineFormat);
  fprintf(AIRFile,"%29s",tmps);
  sprintf(tmps,"BYTE OFFSET OF CALIBRATION HEADER =");
  fprintf(AIRFile,"%-.35s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffsetCalHeader);
  fprintf(AIRFile,"%15s",tmps);
  sprintf(tmps,"BYTE OFFSET OF DEM HEADER =");
  fprintf(AIRFile,"%-.27s",tmps);
  sprintf(tmps,"%ld",Head1->ByteOffsetDEMHeader);
  fprintf(AIRFile,"%23s",tmps);  
  fprintf(AIRFile,"%150s"," ");  /* blank fill fields 18-20 of header */

  for (i=0;i<RecordByteLen-1000;i++)
    fprintf(AIRFile,"%s"," "); /*blank pad to record length*/ 

    
  /* Write parameter header */
  sprintf(tmps,"NAME OF HEADER");
  fprintf(AIRFile,"%-.14s",tmps);
  sprintf(tmps,"%s",Param->HeaderName);
  fprintf(AIRFile,"%36s",tmps);  
  sprintf(tmps,"SITE NAME");
  fprintf(AIRFile,"%-.9s",tmps);
  sprintf(tmps,"%30s",Param->SiteName);
  fprintf(AIRFile,"%41s",tmps);
  sprintf(tmps,"LATITUDE OF SITE (DEGREES)");
  fprintf(AIRFile,"%-.26s",tmps);
  sprintf(tmps,"%.4f",Param->SiteLatitude);
  fprintf(AIRFile,"%24s",tmps);    
  sprintf(tmps,"LONGITUDE OF SITE (DEGREES)");
  fprintf(AIRFile,"%-.27s",tmps);
  sprintf(tmps,"%.4f",Param->SiteLongitude);
  fprintf(AIRFile,"%23s",tmps);    
  sprintf(tmps,"IMAGE TITLE");  /* field 5 */
  fprintf(AIRFile,"%-.11s",tmps);
  sprintf(tmps,"%s",Param->ImageTitle);
  fprintf(AIRFile,"%39s",tmps);  
  sprintf(tmps,"CDROM ID");
  fprintf(AIRFile,"%-.8s",tmps);
  if (Param->CDROMId != DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->CDROMId); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%42s",tmps);      
  sprintf(tmps,"FREQUENCY");
  fprintf(AIRFile,"%-.9s",tmps);
  sprintf(tmps,"%s",Param->Freq);
  fprintf(AIRFile,"%41s",tmps);
  sprintf(tmps,"POLARIZATION");
  fprintf(AIRFile,"%-.12s",tmps);
  sprintf(tmps,"%s",Param->Polarization);
  fprintf(AIRFile,"%38s",tmps);
  sprintf(tmps,"DATA TYPE");
  fprintf(AIRFile,"%-.9s",tmps);
  sprintf(tmps,"%s",Param->DataType);
  fprintf(AIRFile,"%41s",tmps);
  sprintf(tmps,"DATA ID");   /* field 10 */
  fprintf(AIRFile,"%-.7s",tmps);
  if (Param->DataId != DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->DataId); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%43s",tmps);
  sprintf(tmps,"ARCHIVAL FLAG");
  fprintf(AIRFile,"%-.13s",tmps);
  sprintf(tmps,"%ld",Param->ArchivalFlg);
  fprintf(AIRFile,"%37s",tmps);
  sprintf(tmps,"CDROM START PRI ID");
  fprintf(AIRFile,"%-.18s",tmps);
  sprintf(tmps,"%ld",Param->CDROMStartPRIId);
  fprintf(AIRFile,"%32s",tmps);
  sprintf(tmps,"PROC SCENE START PRI ID");
  fprintf(AIRFile,"%-.23s",tmps);
  sprintf(tmps,"%ld",Param->ProcStartPRIId);
  fprintf(AIRFile,"%27s",tmps);
  sprintf(tmps,"LATITUDE AT START OF SCENE (DEGREES)");
  fprintf(AIRFile,"%-.36s",tmps);
  sprintf(tmps,"%.4f",Param->LatitudeStartScene);
  fprintf(AIRFile,"%14s",tmps);
  sprintf(tmps,"LONGITUDE AT START OF SCENE (DEGREES)"); /* field 15 */
  fprintf(AIRFile,"%-.37s",tmps);
  sprintf(tmps,"%.4f",Param->LongitudeStartScene);
  fprintf(AIRFile,"%13s",tmps);
  sprintf(tmps,"LATITUDE AT END OF SCENE (DEGREES)");
  fprintf(AIRFile,"%-.34s",tmps);
  sprintf(tmps,"%.4f",Param->LatitudeEndScene);
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"LONGITUDE AT END OF SCENE (DEGREES)");
  fprintf(AIRFile,"%-.35s",tmps);
  sprintf(tmps,"%.4f",Param->LongitudeEndScene);
  fprintf(AIRFile,"%15s",tmps);
  sprintf(tmps,"APPROXIMATE STARTING HDDT FOOTAGE");
  fprintf(AIRFile,"%-.33s",tmps);
  if (Param->ApproxStartHDDTFootage!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->ApproxStartHDDTFootage); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%17s",tmps);
  sprintf(tmps,"TIME OF ACQUISITION: UTC YEAR");
  fprintf(AIRFile,"%-.29s",tmps);
  sprintf(tmps,"%s",Param->AcquisitionDate);
  fprintf(AIRFile,"%21s",tmps);
  sprintf(tmps,"TIME OF ACQUISITION: UTC DAY"); /* field 20 */
  fprintf(AIRFile,"%-.28s",tmps);
  sprintf(tmps,"%ld",Param->AcquisitionDay);
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"TIME OF ACQUISITION: SECONDS OF DAY");
  fprintf(AIRFile,"%-.35s",tmps);
  sprintf(tmps,"%.1f",Param->AcquisitionSec);
  fprintf(AIRFile,"%15s",tmps);
  sprintf(tmps,"RECORD WINDOW DURATION (MICROSECONDS)");
  fprintf(AIRFile,"%-.37s",tmps);
  sprintf(tmps,"%ld",Param->RecordWindowDuration);
  fprintf(AIRFile,"%13s",tmps);
  sprintf(tmps,"FREQUENCIES COLLECTED");
  fprintf(AIRFile,"%-.21s",tmps);
  sprintf(tmps,"%s",Param->FreqsCollected);
  fprintf(AIRFile,"%29s",tmps); 
  sprintf(tmps,"DIGITAL DELAY (MICROSECONDS)");
  fprintf(AIRFile,"%-.28s",tmps);
  sprintf(tmps,"%.1f",Param->DigitalDelay);
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"CHIRP DELAY (MICROSECONDS)"); /* field 25 */
  fprintf(AIRFile,"%-.26s",tmps);
  if (Param->ChirpDelay!=DEFAULT_VAL_F) 
    {  sprintf(tmps,"%.1f",Param->ChirpDelay); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%24s",tmps);
  sprintf(tmps,"PROCESSOR DELAY (MICROSECONDS)");
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->ProcDelay!=DEFAULT_VAL_I) 
    {  sprintf(tmps,"%ld",Param->ProcDelay); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps); 
  sprintf(tmps,"PRF AT START OF TRANSFER (HZ)");
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->PRFStartTransfer!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->PRFStartTransfer); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps);
  sprintf(tmps,"SAMPLING RATE (MHZ)");
  fprintf(AIRFile,"%-.19s",tmps);
  if (Param->SamplingRate!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->SamplingRate); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%31s",tmps);
  sprintf(tmps,"CENTRE FREQUENCY AT VIDEO (MHZ)");
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->CentreFreq!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CentreFreq); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps); 
  sprintf(tmps,"CHIRP BANDWIDTH (MHZ)"); /* field 30 */
  fprintf(AIRFile,"%-.21s",tmps);
  if (Param->ChirpBandwidth!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->ChirpBandwidth); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%29s",tmps);
  sprintf(tmps,"TYPE OF CHIRP USED (ANALOG OR DIGITAL)");
  fprintf(AIRFile,"%-.38s",tmps);
  sprintf(tmps,"%s",Param->TypeOfChirp);
  fprintf(AIRFile,"%12s",tmps);
  sprintf(tmps,"PULSE LENGTH (MICROSECONDS)");
  fprintf(AIRFile,"%-.27s",tmps);
  if (Param->PulseLen!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->PulseLen); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%23s",tmps); 
  sprintf(tmps,"PROCESSOR WAVELENGTH (METRES)");
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->ProcWavelen!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.5f",Param->ProcWavelen); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps);
  sprintf(tmps,"BAROMETRIC ALTITUDE (METRES)");
  fprintf(AIRFile,"%-.28s",tmps);
  if (Param->BaroAltitude!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->BaroAltitude); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"RADAR ALTIMETER ALTITUDE (METRES)");  /* field 35 */
  fprintf(AIRFile,"%-.33s",tmps);
  if (Param->RadarAltimeterAlt!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->RadarAltimeterAlt); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%17s",tmps); 
  sprintf(tmps,"HEIGHT AGL USED IN PROCESSOR (METRES)");
  fprintf(AIRFile,"%-.37s",tmps);
  if (Param->ProcAltitude!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->ProcAltitude); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%13s",tmps);
  sprintf(tmps,"ELEVATION OF INVESTIGATOR SITE (METRES)");
  fprintf(AIRFile,"%-.39s",tmps);
  if (Param->ElevInvestigatorSite!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%f",Param->ElevInvestigatorSite); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%11s",tmps);
  sprintf(tmps,"AIRCRAFT TRACK ANGLE (DEGREES)");
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->AircraftTrackAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->AircraftTrackAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps); 
  sprintf(tmps,"AIRCRAFT YAW ANGLE (DEGREES)");
  fprintf(AIRFile,"%-.28s",tmps);
  if (Param->AircraftYawAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->AircraftYawAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"AIRCRAFT PITCH ANGLE (DEGREES)"); /* field 40 */
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->AircraftPitchAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->AircraftPitchAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps);
  sprintf(tmps,"AIRCRAFT ROLL ANGLE (DEGREES)");
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->AircraftRollAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->AircraftRollAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps); 
  sprintf(tmps,"PROCESSOR YAW ANGLE USED (DEGREES)");
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->ProcYawAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->ProcYawAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"PROCESSOR PITCH ANGLE USED (DEGREES)");
  fprintf(AIRFile,"%-.36s",tmps);
  if (Param->ProcPitchAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->ProcPitchAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%14s",tmps);
  sprintf(tmps,"PROCESSOR ROLL ANGLE USED (DEGREES)");
  fprintf(AIRFile,"%-.35s",tmps);
  if (Param->ProcRollAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->ProcRollAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%15s",tmps); 
  sprintf(tmps,"NOMINAL PRF RATIO (HZ/KNOT)"); /* field 45 */
  fprintf(AIRFile,"%-.27s",tmps);
  if (Param->NomPRFRatioHzPerKnot!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->NomPRFRatioHzPerKnot); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%23s",tmps);
  sprintf(tmps,"NOMINAL PRF RATIO (1/METRES)");
  fprintf(AIRFile,"%-.28s",tmps);
  if (Param->NomPRFRatioPerMetre!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->NomPRFRatioPerMetre); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%22s",tmps);
  sprintf(tmps,"PRF RATIO CORRECTION FACTOR USED");
  fprintf(AIRFile,"%-.32s",tmps);
  if (Param->PRFCorrectionFactor!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->PRFCorrectionFactor); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%18s",tmps); 
  sprintf(tmps,"RANGE FFT SIZE");
  fprintf(AIRFile,"%-.14s",tmps);
  if (Param->RngFFTSize!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->RngFFTSize); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%36s",tmps);
  sprintf(tmps,"AZIMUTH FFT SIZE");
  fprintf(AIRFile,"%-.16s",tmps);
  if (Param->AzFFTSize!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->AzFFTSize); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%34s",tmps);
  sprintf(tmps,"FRAME SIZE (RANGE LINES)"); /* field 50 */
  fprintf(AIRFile,"%-.24s",tmps);
  if (Param->FrameSize!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->FrameSize); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%26s",tmps); 
  sprintf(tmps,"NUMBER OF FRAMES PROCESSED");
  fprintf(AIRFile,"%-.26s",tmps);
  if (Param->NumFramesProcessed!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->NumFramesProcessed); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%24s",tmps);
  sprintf(tmps,"RANGE ALIGNMENT DELAY USED, HH (MICROSEC)");
  fprintf(AIRFile,"%-.41s",tmps);
  if (Param->RngAlignDelayHH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->RngAlignDelayHH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%9s",tmps);
  sprintf(tmps,"RANGE ALIGNMENT DELAY USED, HV (MICROSEC)");
  fprintf(AIRFile,"%-.41s",tmps);
  if (Param->RngAlignDelayHV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->RngAlignDelayHV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%9s",tmps); 
  sprintf(tmps,"RANGE ALIGNMENT DELAY USED, VH (MICROSEC)");
  fprintf(AIRFile,"%-.41s",tmps);
  if (Param->RngAlignDelayVH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->RngAlignDelayVH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%9s",tmps);
  sprintf(tmps,"RANGE ALIGNMENT DELAY USED, VV (MICROSEC)"); /* field 55 */
  fprintf(AIRFile,"%-.41s",tmps);
  if (Param->RngAlignDelayVV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->RngAlignDelayVV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%9s",tmps);
  sprintf(tmps,"NEAR SLANT RANGE (METERS)");
  fprintf(AIRFile,"%-.25s",tmps);
  if (Param->NearSlantRng!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->NearSlantRng); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%25s",tmps); 
  sprintf(tmps,"FAR SLANT RANGE (METERS)");
  fprintf(AIRFile,"%-.24s",tmps);
  if (Param->FarSlantRng!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->FarSlantRng); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%26s",tmps);
  sprintf(tmps,"NEAR LOOK ANGLE (DEGREES)");
  fprintf(AIRFile,"%-.25s",tmps);
  if (Param->NearLookAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->NearLookAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%25s",tmps); 
  sprintf(tmps,"FAR LOOK ANGLE (DEGREES)");
  fprintf(AIRFile,"%-.24s",tmps);
  if (Param->FarLookAngle!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->FarLookAngle); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%26s",tmps);
  sprintf(tmps,"NUMBER OF LOOKS PROCESSED IN AZIMUTH"); /* field 60 */
  fprintf(AIRFile,"%-.36s",tmps);
  if (Param->NumProcAzLooks!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->NumProcAzLooks); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%14s",tmps); 
  sprintf(tmps,"NUMBER OF LOOKS PROCESSED IN RANGE");
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->NumProcRngLooks!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%f",Param->NumProcRngLooks); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"RANGE WEIGHTING USED");
  fprintf(AIRFile,"%-.20s",tmps);
  sprintf(tmps,"%s",Param->RngWeighting);
  fprintf(AIRFile,"%30s",tmps); 
  sprintf(tmps,"RANGE WEIGHTING COEFFICIENT");
  fprintf(AIRFile,"%-.27s",tmps);
  if (Param->RngWeightingCoeff!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->RngWeightingCoeff); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%23s",tmps);
  sprintf(tmps,"AZIMUTH WEIGHTING USED");
  fprintf(AIRFile,"%-.22s",tmps);
  sprintf(tmps,"%s",Param->AzWeighting);
  fprintf(AIRFile,"%28s",tmps); 
  sprintf(tmps,"AZIMUTH WEIGHTING COEFFICIENT"); /* field 65 */ 
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->AzWeightingCoeff!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->AzWeightingCoeff); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps);
  sprintf(tmps,"PERCENT OF PRF BANDWIDTH USED");
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->PercentPRFBWProc!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->PercentPRFBWProc); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps); 
  sprintf(tmps,"DESKEW FLAG (1=DESKEWED, 2=NOT DESKEWED)");
  fprintf(AIRFile,"%-.40s",tmps);
  if (Param->DeskewFlg!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->DeskewFlg); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%10s",tmps);
  sprintf(tmps,"SLANT RANGE SAMPLE SPACING (METERS)");
  fprintf(AIRFile,"%-.35s",tmps);
  if (Param->SlantRngSampleSpacing!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->SlantRngSampleSpacing); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%15s",tmps); 
  sprintf(tmps,"NOMINAL SLANT RANGE RESOLUTION (METERS)");
  fprintf(AIRFile,"%-.39s",tmps);
  if (Param->NomSlantRngRes!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->NomSlantRngRes); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%11s",tmps);
  sprintf(tmps,"AZIMUTH SAMPLE SPACING (METERS)"); /* field 70 */
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->AzSampleSpacing!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->AzSampleSpacing); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps); 
  sprintf(tmps,"NOMINAL AZIMUTH RESOLUTION (METERS)");
  fprintf(AIRFile,"%-.35s",tmps);
  if (Param->NomAzRes!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->NomAzRes); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%15s",tmps);
  sprintf(tmps,"NUMBER OF INTERPOLATION POINTS USED IN RMC");
  fprintf(AIRFile,"%-.42s",tmps);
  if (Param->NumInterpPtsRMC!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->NumInterpPtsRMC); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%8s",tmps);
  sprintf(tmps,"AZIMUTH REFERENCE SIZE/LOOK, NEAR RANGE");
  fprintf(AIRFile,"%-.39s",tmps);
  if (Param->AzRefSizePerLookNear!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->AzRefSizePerLookNear); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%11s",tmps); 
  sprintf(tmps,"AZIMUTH REFERENCE SIZE/LOOK, FAR RANGE");
  fprintf(AIRFile,"%-.38s",tmps);
  if (Param->AzRefSizePerLookFar!=DEFAULT_VAL_I) 
    { sprintf(tmps,"%ld",Param->AzRefSizePerLookFar); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%12s",tmps);
  sprintf(tmps,"IMAGE CENTRE LATITUDE (DEGREES)"); /* field 75 */
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->ImageCentreLat!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->ImageCentreLat); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps); 
  sprintf(tmps,"IMAGE CENTRE LONGITUDE (DEGREES)");
  fprintf(AIRFile,"%-.32s",tmps);
  if (Param->ImageCentreLong!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->ImageCentreLong); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%18s",tmps);
  sprintf(tmps,"CALTONE VIDEO FREQUENCY (MHZ)");
  fprintf(AIRFile,"%-.29s",tmps);
  if (Param->CalToneVideoFreq!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.4f",Param->CalToneVideoFreq); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%21s",tmps); 
  sprintf(tmps,"CALTONE POWER MEASURED, DB, HH");
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->CalTonePowerHH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalTonePowerHH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps);
  sprintf(tmps,"CALTONE POWER MEASURED, DB, HV");
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->CalTonePowerHV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalTonePowerHV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps); 
  sprintf(tmps,"CALTONE POWER MEASURED, DB, VH");/* field 80 */
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->CalTonePowerVH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalTonePowerVH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps);
  sprintf(tmps,"CALTONE POWER MEASURED, DB, VV");
  fprintf(AIRFile,"%-.30s",tmps);
  if (Param->CalTonePowerVV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalTonePowerVV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%20s",tmps); 
  sprintf(tmps,"CALIBRATION FACTOR APPLIED, DB, HH");
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->CalibFactorHH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalibFactorHH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"CALIBRATION FACTOR APPLIED, DB, HV");
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->CalibFactorHV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalibFactorHV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps); 
  sprintf(tmps,"CALIBRATION FACTOR APPLIED, DB, VH");
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->CalibFactorVH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalibFactorVH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps);
  sprintf(tmps,"CALIBRATION FACTOR APPLIED, DB, VV");/* field 85 */
  fprintf(AIRFile,"%-.34s",tmps);
  if (Param->CalibFactorVV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.2f",Param->CalibFactorVV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%16s",tmps); 
  sprintf(tmps,"MEASURED AND CORRECTED HV/VH POWER RATIO");
  fprintf(AIRFile,"%-.40s",tmps);
  if (Param->MeasCorrHVVHPowerRatio!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.3f",Param->MeasCorrHVVHPowerRatio); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%10s",tmps);
  sprintf(tmps,"MEASURED AND CORRECTED HV/VH PHASE (DEG)");
  fprintf(AIRFile,"%-.40s",tmps);
  if (Param->MeasCorrHVVHPhase!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->MeasCorrHVVHPhase); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%10s",tmps); 
  sprintf(tmps,"CALTONE PHASE MEASURED, DEG, HH");
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->CalTonePhaseMeasHH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->CalTonePhaseMeasHH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps);
  sprintf(tmps,"CALTONE PHASE MEASURED, DEG, HV");
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->CalTonePhaseMeasHV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->CalTonePhaseMeasHV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps); 
  sprintf(tmps,"CALTONE PHASE MEASURED, DEG, VH");/* field 90 */
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->CalTonePhaseMeasVH!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->CalTonePhaseMeasVH); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps);
  sprintf(tmps,"CALTONE PHASE MEASURED, DEG, VV");
  fprintf(AIRFile,"%-.31s",tmps);
  if (Param->CalTonePhaseMeasVV!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->CalTonePhaseMeasVV); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%19s",tmps); 
  sprintf(tmps,"GENERAL SCALE FACTOR");
  fprintf(AIRFile,"%-.20s",tmps);
  if (Param->GenScaleFactor!=DEFAULT_VAL_F) 
    { sprintf(tmps,"%.1f",Param->GenScaleFactor); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%30s",tmps);
  sprintf(tmps,"GPS ALTITUDE, M");
  fprintf(AIRFile,"%-.15s",tmps);
  if (Param->GPSAlt!=DEFAULT_VAL_F)
    {  sprintf(tmps,"%.2f",Param->GPSAlt); }
  else 
    { sprintf(tmps,"%s"," "); }
  fprintf(AIRFile,"%35s",tmps); 
  fprintf(AIRFile,"%350s"," "); /* blank fields 94-100 */

  for (i=0;i<RecordByteLen-5000;i++)
    fprintf(AIRFile,"%s"," "); /*blank pad to record length*/ 

  /* WRITE IMAGE DATA */
  /* Repeat for each row */
  for (Row=0;Row<HeightPixels;Row++)
    {
    indx = 0;
    /* Read and write image row (record) */
    switch(DataType)
      {
      case 0:
        fread(DataUnChar,sizeof(unsigned char),ValWidth,InFile);
	fwrite(DataUnChar,sizeof(unsigned char),ValWidth,AIRFile);
        break;
      case 1:
        fread(DataUnInt2B,sizeof(unsigned Int2B),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataUnInt2B,sizeof(unsigned Int2B),ValWidth,AIRFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataUnInt2B[Val];
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;
	    }  
	  fwrite(DataUnChar,1,sizeof(unsigned Int2B)*ValWidth,AIRFile);
	  }  /* end else */
        break;
      case 2:
        fread(DataInt4B,sizeof(Int4B),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataInt4B,sizeof(Int4B),ValWidth,AIRFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataInt4B[Val];
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(Int4B)*ValWidth,AIRFile);
	  }  /* end else */  
        break;
      case 3:
        fread(DataFloat,sizeof(float),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataFloat,sizeof(float),ValWidth,AIRFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataFloat[Val];
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(float)*ValWidth,AIRFile);
	  }  /* end else */  
        break;
      case 4:
        fread(DataDouble,sizeof(double),ValWidth,InFile);
        if (BigEndian)
	  {fwrite(DataDouble,sizeof(double),ValWidth,AIRFile);}
        else /* for small endian host */
	  {
          for (Val=0;Val<ValWidth;Val++)
	    {
            ptr = (unsigned char *)&DataDouble[Val];
	    DataUnChar[indx++] = *(ptr+7);
	    DataUnChar[indx++] = *(ptr+6);
	    DataUnChar[indx++] = *(ptr+5);	    
	    DataUnChar[indx++] = *(ptr+4);
	    DataUnChar[indx++] = *(ptr+3);
	    DataUnChar[indx++] = *(ptr+2);
	    DataUnChar[indx++] = *(ptr+1);
	    DataUnChar[indx++] = *ptr;	    
	    }  
	  fwrite(DataUnChar,1,sizeof(double)*ValWidth,AIRFile);
	  }  /* end else */  
        break;
      }  /* end switch data type */

    /* Blank fill to record length */          
    for (i=0;i<RecordByteLen-DataByteWidth;i++)
      fprintf(AIRFile,"%s"," "); /*blank pad to record length*/ 

    
    }  /* end for row (record) */



  /* Tidy up */
  fclose (InFile);
  fclose (AIRFile);

  switch(DataType)
    {
    case 0:
      free (DataUnChar);
      break;
    case 1:
      free (DataUnInt2B);
      if (!BigEndian) free(DataUnChar);	
      break;
    case 2:
      free (DataInt4B);
      if (!BigEndian) free(DataUnChar);
      break;
    case 3:
      free(DataFloat);
      if (!BigEndian) free(DataUnChar);	
      break;
    case 4:
      free(DataDouble);
      if (!BigEndian) free(DataUnChar);	
      break;
    } /* end switch DataType */

  /* Get end time */
  TimeE = time(NULL);
  fprintf(msg,"AIRSAR FORMAT IMAGE written - in %ld secs (%.2f min)\n",
	  TimeE-TimeS,(double)(TimeE-TimeS)/60.0);  
  
  return(1);
}

