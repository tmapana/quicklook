/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2proc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Unauthorized use prohibited!!
Part of G2 SAR processor. Main program for integrated G2 SAR Processor.

VERSION/AUTHOR/DATE : pre 1998-01-25 / Jasper Horrell / pre 1999-01-25
COMMENTS: 
1a - 19970821
1b - 19970926 major additions.
1c - 19971009 major additions incl ground rng projection function.
1d - (1997-12-05) - print template command file. Accomodate changes in 
     the modules such as program version checking etc.
G2Proc (1997-12-05) - rename.
(1997-12-08) - rename rng focus updates to RngFocSegments, remove G2 tag to 
      variables in command file and add dollar symbol to SeekP searches.
(1997-12-15) - add in call to corner turn function before ground range 
      projection. Also change groung range projection function so that no 
      corner turn is performed by it.
(1998-01-21) - fixed bug with ground range projection input file name.

VERSION/AUTHOR/DATE : 1999-01-25 / Jasper Horrell / 1999-01-25
COMMENTS:
Started changing to parse processing parameters directly from SASAR lbr 
file instead of from SASAR CEOS format files.

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#include"g2func.h"
#include"g2parse.h"
#include"g2sniffdc.h"
#include"g2rnc.h"
#include"g2cor.h"
#include"g2azc.h"
#include"g2imgfu.h"

/* Misc defns limited to this file */
#define PROG_VERSION "G2 - 1999-01-25"
#define SNIFFDC_VERSION "1999-01-25"
#define RNC_VERSION "1999-01-22"
#define COR_VERSION "1999-01-20"
#define AZC_VERSION "1999-01-20"
#define RNG_COMPRESSED_FILE "tmpg2.rnc"
#define CORNER_TURNED_FILE "tmpg2.cor"
#define AZ_COMPRESSED_FILE "tmpg2.azc"
#define UNBLOCK_AZC_FILE "tmpg2unb.azc"
#define GRNG_FILE "tmpg2.grg"
#define MOCOMP_FILE "tmpg2.moc"
#define STC_FILE "tmpg2.stc"

#define NOM_HEIGHT 10000.0
    
/* Function prototypes */
Int2B WriteG2TmplCmdFile(char CmdFileName[]);

/************************************************************************/
/* MAIN PROGRAM */
void main(Int2B argc, char *argv[])
{
  FILE *cmdfp,*logfp=NULL,*msg=stdout,*LBRFile,*RawFile;
  char tmps[STRING_SPACE],RunID[STRING_SPACE],Messages[STRING_SPACE]="displ",
       LogName[STRING_SPACE],
       LBRName[STRING_SPACE],RawName[STRING_SPACE],
       SALSiteName[STRING_SPACE],
       MocompFlg,SUNRasterFlg,AIRSARFileFlg,
       RncCmdName[]= RNC_CMD_FILE,
       CorCmdName[]= COR_CMD_FILE,
       AzcCmdName[]= AZC_CMD_FILE,
       RasterFileName[STRING_SPACE],AIRSARFileName[STRING_SPACE],
       IMOPolarization[3],AzcFileName[STRING_SPACE],
       ProgVersion[STRING_SPACE];
  Int2B DataChannel = 0, /* only one data channel possible now, but keep for 
                            future upgrade possibility. Note independent of
                            HBR number. */
        HBR = DEFAULT_VAL_I;
  Int4B AzBlock,
        AzFFTSize,
        AzLooks,
        AzProcStartRngBin = 0,
        AzProcRngBins = 0,
        AzRefFuncSign,
        DetectMethod,
        Errors = 0,  
        GrndRngInterpPts,
        GrndRngProjected=0,
        HeightPixels=0,
        i,
        IMOAcqDay,
        IMOAcqYear,
        IMOFirstPRIID,
        IMODataRecords,
        IMODCOffset,
        IMORecordLen,
        IMORngBins,
        MocompPhaseSign,
        NumAzBlocks=1, 
        OutputDataType,
        PRIsPerAzBlock,
        RngCurvBatchSize,
        RngCurvInterpSize,
        RngFocSegments,
        SALPresumRatio,
        StartPRIOffset, 
        tmpl,
        WidthPixels=0;
  double AzScale,
         AzWinConstFreq,
         AzWinConstTime,
         AzWinTimeBroadenFactor=1.0, 
         GrndRngBinSpacing=0.0,
         IMOAcqSec,
         IMOCarrierFreq,
         IMODelay0,
         IMOEndAlt=0.0,    
         IMOEndLat=0.0,
         IMOEndLong=0.0,
         IMOMidAlt=0.0,   
         IMOMidLat=0.0,
         IMOMidLong=0.0,
         IMOStartAlt=0.0, 
         IMOStartLat=0.0,
         IMOStartLong=0.0,
         IMOPulseLen,
         LookOverlapFrac,
         NomAzRes,
         SALADFreq,
         SALAveTerrainAlt,
         SALAzBeamwidth,
         SALNomGroundSpeed,
         SALNomPRF;
  struct RncFileStruct RncCmd;
  struct CorFileStruct CorCmd;
  struct AzcFileStruct AzcCmd;
  struct AIRSARHead1Struct AIRHead1;
  struct AIRSARParamHeadStruct AIRParam;
  struct CntrlStruct Gen; /* for data channel info */
  struct DataChannelStruct Chan[MAX_DATA_CHANNELS];

  if (argc != 2)
    {
    fprintf(msg,"\nProg: G2 SAR Processor (Ver. %s)\n",PROG_VERSION);  
    fprintf(msg,"Code: J.M. Horrell (Copyright UCT 1997)\n");
    fprintf(msg,"\nUSAGE : g2 [Cmd File]\n\n");
    fprintf(msg,"To see command file structure, type `g2 -tmpl'\n"); 
    fprintf(msg,"(generates a template command file `g2.tmp')\n\n");  
    exit(1);
    } 

  /* Write template command file, if requested */ 
  if (argc==2 && strcmp(argv[1],"-tmpl")==0)
    {
     if (WriteG2TmplCmdFile("g2.tmp"))
       fprintf(msg,"Template command file `g2.tmp' written!\n");
     exit(1);
    }
  
  /* Open command file (as binary) */  
  if ( (cmdfp = fopen (argv[1], "rb") ) == NULL )
    {
    fprintf(msg,"ERROR: G2 cmd file %s not found/opened\n",argv[1]);
    exit(1);
    }

  /* Check version IDs match */
  if (!SeekP(cmdfp,"$ProgramVersion","=>",-1,0)) 
    { fprintf(msg,
      "WARNING - in parsing command file (prog version not found)!\n\n"); } 
  else
    {
    ReadStr(cmdfp,ProgVersion,STRING_SPACE); 
    if (strcmp(ProgVersion,PROG_VERSION)!=0)
      {
      fprintf(msg,
        "\nWARNING - command file ver. (%s) not as for prog. (%s)!\n",
        ProgVersion,PROG_VERSION);
      }
    } /* end else */

  /* Read RunID from command file  */
  if (!SeekP(cmdfp,"$RunID","=>",-1,0))
    {
    fprintf(msg,"ERROR : No RunID found in G2 cmd file %s!!\n",argv[1]);
    exit(1);
    } 
  ReadStr(cmdfp,RunID,STRING_SPACE);

  /* Read message output destination */
  if (!SeekP(cmdfp,"$Messages","=>",-1,0))
    { fprintf(msg,"WARNING - messages not specfied in command file!\n"); }
  else
    { ReadStr(cmdfp,Messages,STRING_SPACE); }
  if (strcmp(Messages,"file")==0)
    {
    /* Open log file */
    strcpy(LogName,RunID);
    strcat(LogName,".log");
    if ( (logfp = fopen (LogName, "wt") ) == NULL )
      {
      fprintf(msg,
         "\nERROR : Could not open log file for RunID %s!!\n",RunID );
      exit(1);
      }
    msg = logfp;  /* assign messages to the log file */
    }
  else
    { msg = stdout; }
    
  /* print header info now that message destination fixed. */
  fprintf(msg,"\nProg: G2 SAR Processor (Ver. %s)\n",PROG_VERSION);  
  fprintf(msg,"Code: J.M. Horrell (Copyright UCT 1997-1999)\n");
  fprintf(msg,"\nMessage Log File:\n");
  fprintf(msg,"=================\n\n");  


  /* Read LBR file name and open file */
  if (!SeekP(cmdfp,"$LBRFile","=>",-1,0)) { 
    fprintf(msg,"ERROR - LBR file name not found in command file!\n"); 
    exit(1); 
  }
  else { 
    ReadStr(cmdfp,LBRName,STRING_SPACE); 
  }
  if ( (LBRFile = fopen (LBRName, "rb") ) == NULL ) { 
    fprintf(msg,"ERROR: LBRFile %s not found/opened!\n",LBRName); exit(1); 
  }
  else { 
    fprintf(msg,"LBR file                       : %s\n",LBRName); 
  }  

  /* Read raw data file name and open file */
  if (!SeekP(cmdfp,"$RawDataFile","=>",-1,0)) { 
    fprintf(msg,"ERROR - Raw data file name not found in cmd file!\n");
    exit(1);
  }
  else { 
    ReadStr(cmdfp,RawName,STRING_SPACE); 
  }
  if ( (RawFile = fopen (RawName, "rb") ) == NULL ) { 
    fprintf(msg,"ERROR: RawFile %s not found/opened!\n",RawName); exit(1); 
  }
  else { 
    fprintf(msg,"Raw data file                  : %s\n",RawName); 
  }    

  /* Read HBR to use for the processing */
  if (!SeekP(cmdfp,"$HBR","=>",-1,0)) { 
    fprintf(msg,"ERROR - HBR to process not found in cmd file!\n"); 
    exit(1);
  }
  else { 
    fscanf(cmdfp,"%hd",&HBR);
    fprintf(msg,"HBR to process                 : %d\n",HBR);
  }

  /* Init errors for remaining parsing and initialization */
  Errors=0;

  /* Read StartPRIOffset */
  if (!SeekP(cmdfp,"$StartPRIOffset","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&StartPRIOffset);  
  fprintf(msg,"Spec start PRI offset          : %ld\n",StartPRIOffset);  


  /* PARSE LBR FILE FOR MAIN PARAMETERS */
  ParseLBRFile(msg,LBRFile,&Gen,Chan);


  /* READ PARAMS FROM CEOS FILES COMMON TO ALL PROCESSING */
  fseek(RawFile,8,0); 
  IMORecordLen = BigEndReadInt4B(RawFile);

  fseek(RawFile,IMORecordLen+12,0); 
  IMOFirstPRIID = BigEndReadInt4B(RawFile);
  fprintf(msg,"First PRI ID (IMOP file)       : %ld\n",IMOFirstPRIID);

  fseek(RawFile,180,0); 
  fscanf(RawFile,"%6ld",&IMODataRecords);
  fprintf(msg,"PRIs (data records IMOP file)  : %ld\n",IMODataRecords);  

  fseek(RawFile,216,0); 
  fscanf(RawFile,"%4ld",&tmpl);
  if (tmpl != 8)
    {
    fprintf(msg,
       "ERROR - only 8 bits per I/Q input sample supported (not %ld)!\n"
          ,tmpl);
    exit(1); 
    }

  fseek(RawFile,IMORecordLen+64,0); 
  if (BigEndReadInt2B(RawFile) != 1)
    {
    fprintf(msg,
       "ERROR - G2 software requires range compressed input data!\n");
    exit(1);
    }

  fseek(RawFile,248,0); 
  fscanf(RawFile,"%8ld",&IMORngBins);
  fprintf(msg,"Range bins (IMOP file)         : %ld\n",IMORngBins);  

  fseek(RawFile,IMORecordLen+120,0); 
  IMODelay0 = ((double)BigEndReadInt4B(RawFile))*1.0e-9; /*(sec)*/
  fprintf(msg,"Delay 1st rng bin (IMOP file)  : %.3f usec\n",
	  IMODelay0*1.0e+6);
  
  fseek(RawFile,IMORecordLen+84,0);
  IMOCarrierFreq = ((double)BigEndReadInt4B(RawFile))*1000.0; /*(Hz)*/
  fprintf(msg,"Carrier freq (IMOP file)       : %.3f MHz\n",
    IMOCarrierFreq*1.0e-6);

  fseek(LBRFile,720+934,0);
  fscanf(LBRFile,"%16lf",&SALNomPRF);
  fprintf(msg,"Nom PRF (SARL file)            : %.3f Hz\n",SALNomPRF);

  fseek(LBRFile,720+1782,0);
  fscanf(LBRFile,"%8ld",&SALPresumRatio);
  fprintf(msg,"Presum ratio (SARL file)       : %ld\n",SALPresumRatio);

  fseek(LBRFile,720+1766,0);
  fscanf(LBRFile,"%16lf",&SALNomGroundSpeed);
  fprintf(msg,"Nom grnd speed (SARL file)     : %.3f m/s\n",
          SALNomGroundSpeed);

  fseek(RawFile,640,0);
  fscanf(RawFile,"%8ld",&IMODCOffset);
  fprintf(msg,"DC offset of data (IMOP file)  : %ld\n",IMODCOffset);
  
  fseek(LBRFile,720+710,0);
  fscanf(LBRFile,"%16lf",&SALADFreq);
  SALADFreq *= 1.0e+6;  /* convert to Hz */
  fprintf(msg,"AD freq (SARL file)            : %.3f MHz\n",
          SALADFreq*1.0e-6);

  fseek(RawFile,IMORecordLen+68,0);
  IMOPulseLen =  ((double)BigEndReadInt4B(RawFile))*1.0e-09; /* convert
							       to sec */
  fprintf(msg,"Pulse length (IMO file)        : %.3f nsec\n",
	  IMOPulseLen*1.0e+9);
  
  fseek(LBRFile,720+966,0);
  fscanf(LBRFile,"%16lf",&SALAzBeamwidth);
  fprintf(msg,"Az beamwidth (SARL file)       : %.3f deg\n",SALAzBeamwidth);

  fseek(LBRFile,720+36,0);
  ReadStr2(LBRFile,SALSiteName,32);
  fprintf(msg,"Site name (SARL file)          : %s\n",SALSiteName);  

  fseek(LBRFile,720+308,0);
  fscanf(LBRFile,"%16lf",&SALAveTerrainAlt);  
  fprintf(msg,"Ave terrain altitude (SAL file): %.3f m\n",
         SALAveTerrainAlt);

    

  /* READ PARAMS SPECIFIED IN COMMAND FILE FOR PRODUCTS  */
  if (!SeekP(cmdfp,"$PRIsPerAzBlock","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&PRIsPerAzBlock);
  fprintf(msg,"Spec PRIs per az block         : %ld\n",
            PRIsPerAzBlock);
  if (PRIsPerAzBlock > IMODataRecords)
    {
    PRIsPerAzBlock = IMODataRecords;
    fprintf(msg,
	"WARNING - PRIs per az block reset to num data records!\n");
    } /* end if SynPRIs... */
  
  if (!SeekP(cmdfp,"$NumAzBlocks","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&NumAzBlocks);
  if (NumAzBlocks!=0)
    { fprintf(msg,"Spec az blocks                 : %ld\n",NumAzBlocks); }
  else
    {
    fprintf(msg,"Spec az processing             : to end IMOP file\n");
    NumAzBlocks = (long int)(IMODataRecords/PRIsPerAzBlock);
    fprintf(msg,"Calculated az blocks required  : %ld\n",NumAzBlocks); 
    }

  if (!SeekP(cmdfp,"$MocompFlg","=>",-1,0)) Errors++;
  ReadChar(cmdfp,&MocompFlg);
  if (MocompFlg == 'Y')
    { fprintf(msg,"Spec motion compensation       : yes\n"); }
  else
    {
    MocompFlg = 'N';
    fprintf(msg,"Spec motion compensation       : none\n");
    }

  if (!SeekP(cmdfp,"$MocompPhaseSign","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&MocompPhaseSign);
  if (MocompFlg == 'Y')
    fprintf(msg,"Spec mocomp phase sign         : %ld\n",
	    MocompPhaseSign);  

  if (!SeekP(cmdfp,"$AzProcStartRngBin","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&AzProcStartRngBin);
  if ( (AzProcStartRngBin > IMORngBins-1) || (AzProcStartRngBin < 0) )
    { 
    AzProcStartRngBin = 0; 
    fprintf(msg,
       "WARNING - AzProcStartRngBin out of range - reset to zero!\n");
    }

  if (!SeekP(cmdfp,"$AzProcRngBins","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&AzProcRngBins);
  if ( AzProcRngBins <= 0 )
    { 
    AzProcRngBins = IMORngBins - AzProcStartRngBin; 
    }
  else if ( AzProcStartRngBin+AzProcRngBins > IMORngBins )
    {
    AzProcRngBins = IMORngBins - AzProcStartRngBin; 
    fprintf(msg,"WARNING - AzProcRngBins too large - reset!\n");
    }
  fprintf(msg,"Spec az processing of rng bins : %ld - %ld\n",
      AzProcStartRngBin,AzProcStartRngBin+AzProcRngBins-1);

  if (!SeekP(cmdfp,"$AzFFTSize","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&AzFFTSize);
  fprintf(msg,"Spec az FFT size               : %ld\n",AzFFTSize);

  if (!SeekP(cmdfp,"$RngFocSegments","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&RngFocSegments);
  if (RngFocSegments <= 0) RngFocSegments = AzProcRngBins;
  fprintf(msg,"Spec rng focus updates         : %ld\n",RngFocSegments);

  if (!SeekP(cmdfp,"$AzRefFuncSign","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&AzRefFuncSign);
  fprintf(msg,"Spec az ref func sign          : %ld\n",AzRefFuncSign);  
  
  if (!SeekP(cmdfp,"$NomAzRes","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%lf",&NomAzRes);
  fprintf(msg,"Spec nominal az resolution     : %.3f m\n",NomAzRes);

  if (!SeekP(cmdfp,"$AzWinConstTime","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%lf",&AzWinConstTime);
  fprintf(msg,"Spec az win const (time)       : %.3f\n",
       AzWinConstTime);  

  if (!SeekP(cmdfp,"$AzWinConstFreq","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%lf",&AzWinConstFreq);
  fprintf(msg,"Spec az win const (freq)       : %.3f\n",
       AzWinConstFreq);  

  if (!SeekP(cmdfp,"$AzLooks","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&AzLooks);
  fprintf(msg,"Spec azimuth looks             : %ld\n",AzLooks);

  if (!SeekP(cmdfp,"$LookOverlapFrac","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%lf",&LookOverlapFrac);
  fprintf(msg,"Spec look overlap fraction     : %.3f\n",LookOverlapFrac);
  
  if (!SeekP(cmdfp,"$RngCurvInterpSize","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&RngCurvInterpSize);
  fprintf(msg,"Spec rng curv interp size      : %ld\n",RngCurvInterpSize);

  if (!SeekP(cmdfp,"$RngCurvBatchSize","=>",-1,0)) Errors++; 
  fscanf(cmdfp,"%ld",&RngCurvBatchSize);
  fprintf(msg,"Spec rng curv batch size       : %ld\n",RngCurvBatchSize);

  if (!SeekP(cmdfp,"$DetectMethod","=>",-1,0)) Errors++;  
  fscanf(cmdfp,"%ld",&DetectMethod);
  if (DetectMethod==0)
    { fprintf(msg,"Spec detection method          : none (complex)\n"); }
  else if (DetectMethod==1)   
    { fprintf(msg,"Spec detection method          : magnitude\n"); }
  else if (DetectMethod==2)
    { fprintf(msg,"Spec detection method          : power\n"); }
  else if (DetectMethod==3)
    { fprintf(msg,"Spec detection method          : power in dB\n"); }
  else
    { fprintf(msg,"ERROR - unsupported detection method %ld\n",
	      DetectMethod); exit(1); }
  
  if (!SeekP(cmdfp,"$OutputDataType","=>",-1,0)) Errors++;    
  fscanf(cmdfp,"%ld",&OutputDataType);
  switch(OutputDataType)
    {
    case 0:
      fprintf(msg,
         "Spec output data type          : unsigned char (0 - 255)\n");
      break;
    case 1:
      fprintf(msg,
         "Spec output data type          : unsigned short int (0 - 65535)\n");
      break;
    case 2:
      fprintf(msg,
	 "Spec output data type          : long int (+- 2 147 483 647)\n");
      break;
    case 3:
      fprintf(msg,
         "Spec output data type          : float (+-3.4x10^(+-28))\n");
      break;
    case 4:
      fprintf(msg,
	 "Spec output data type          : double (+- 1.7*10^(+-308))\n");
    }  /* end switch OutputDataType */


  if (!SeekP(cmdfp,"$AzScale","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%lf",&AzScale);
  fprintf(msg,"Spec az scale factor           : %6.4e\n",AzScale);

  if (!SeekP(cmdfp,"$GrndRngInterpPts","=>",-1,0)) Errors++;
  fscanf(cmdfp,"%ld",&GrndRngInterpPts);
  fprintf(msg,"Spec grnd rng interp kernel    : %ld pts\n",
	    GrndRngInterpPts);
   
  if (!SeekP(cmdfp,"$SUNRasterFlg","=>",-1,0)) Errors++;
  ReadChar(cmdfp,&SUNRasterFlg);
  if (SUNRasterFlg == 'Y')
    { fprintf(msg,"Spec SUN raster file output    : yes\n"); }
  else
    {
    SUNRasterFlg = 'N';
    fprintf(msg,"Spec SUN raster file output    : no\n");    
    }

  if (!SeekP(cmdfp,"$AIRSARFileFlg","=>",-1,0)) Errors++;
  ReadChar(cmdfp,&AIRSARFileFlg);   
  if (AIRSARFileFlg == 'Y')
    { fprintf(msg,"Spec AIRSAR format file output : yes\n"); }
  else
    {
    AIRSARFileFlg = 'N';
    fprintf(msg,"Spec AIRSAR format file output : no\n");    
    }  

  /* Close spec file */
  fclose(cmdfp);

  if (Errors!=0)
    {fprintf(msg,"WARNING - %ld errors in parsing command file!\n\n",Errors);}


  /* MISC DATA for processed image annotation (AIRSAR format image) */    

  fseek(RawFile,0,0);
  for (i=0;i<1+StartPRIOffset;i++)
    { fseek(RawFile,IMORecordLen,1);} /* find start of image */
  fseek(RawFile,36,1);
  IMOAcqYear = BigEndReadInt4B(RawFile); /* IMU INT values */
  IMOAcqDay = BigEndReadInt4B(RawFile);
  IMOAcqSec = (double)BigEndReadInt4B(RawFile)*1.0e-03;
  fseek(RawFile,84,1);
  IMOStartLat = (double)BigEndReadInt4B(RawFile)*1.0e-06;     
  IMOStartLong = (double)BigEndReadInt4B(RawFile)*1.0e-06;
  IMOStartAlt = (double)BigEndReadInt4B(RawFile)*1.0e-03;           
  fprintf(msg,"Start proc lat/long (IMO file) : %.3f / %.3f deg\n",
          IMOStartLat,IMOStartLong);
  fprintf(msg,"Start proc aircraft alt (AMSL) : %.3f m\n",
          IMOStartAlt); 

  fseek(RawFile,0,0);
  for (i=0;i<1+StartPRIOffset+(Int4B)(NumAzBlocks*PRIsPerAzBlock/2);i++)
    { fseek(RawFile,IMORecordLen,1);} /* find centre of image */
  fseek(RawFile,132,1);
  IMOMidLat = (double)BigEndReadInt4B(RawFile)*1.0e-06;     
  IMOMidLong = (double)BigEndReadInt4B(RawFile)*1.0e-06;
  IMOMidAlt = (double)BigEndReadInt4B(RawFile)*1.0e-03;      
  fprintf(msg,"Mid proc lat / long (IMO file) : %.3f / %.3f deg\n",
          IMOMidLat,IMOMidLong);
  fprintf(msg,"Mid proc aircraft alt (AMSL)   : %.3f m\n",
          IMOMidAlt); 

  fseek(RawFile,0,0);
  for (i=0;i<1+StartPRIOffset+NumAzBlocks*PRIsPerAzBlock-1;i++)
    { fseek(RawFile,IMORecordLen,1);} /* find end of image */
  fseek(RawFile,132,1);
  IMOEndLat = (double)BigEndReadInt4B(RawFile)*1.0e-06;     
  IMOEndLong = (double)BigEndReadInt4B(RawFile)*1.0e-06;
  IMOEndAlt = (double)BigEndReadInt4B(RawFile)*1.0e-03;
  fprintf(msg,"End proc lat / long (IMO file) : %.3f / %.3f deg\n",
          IMOEndLat,IMOEndLong);
  fprintf(msg,"End proc aircraft alt (AMSL)   : %.3f m\n",
          IMOEndAlt); 

  fseek(RawFile,IMORecordLen+52,0);  /* find polarization flags */
  if (BigEndReadInt2B(RawFile)==1)
    { strcpy(IMOPolarization,"V"); }
  else 
    { strcpy(IMOPolarization,"H"); }
  if (BigEndReadInt2B(RawFile)==1)
    { strcat(IMOPolarization,"V"); }
  else 
    { strcat(IMOPolarization,"H"); }

  if (AzWinConstTime==0.08) AzWinTimeBroadenFactor = 1.4; /* for now */


  /* SET UP DATA FOR COMMAND FILES  */

  /* Set up default rnc command file params */  
  RncCmd.msg = msg;
  RncCmd.ScreenUpdateRate=10;
  RncCmd.StartProcessPRI=1+StartPRIOffset; /* skip first record of IMO file*/
  RncCmd.InputPRF=SALNomPRF/(double)SALPresumRatio;
  RncCmd.DopCentroid=0.0;
  RncCmd.PreSumRatio=1;
  RncCmd.PreSummedPulsesToUse=PRIsPerAzBlock;
  RncCmd.HeaderBytes=0;
  RncCmd.FooterBytes=0;
  RncCmd.InputFileRngBins=IMORngBins;
  RncCmd.StartRngBin=0;
  RncCmd.RngBinsToProcess=IMORngBins; 
  strcpy(RncCmd.InputFileName,RawName);
  strcpy(RncCmd.OutputFileName,RNG_COMPRESSED_FILE );
  RncCmd.STCFlg = 'N';
  strcpy(RncCmd.STCFileName,STC_FILE);
  RncCmd.MocompFlg=MocompFlg;
  RncCmd.MocompRngShiftFlg='N';
  RncCmd.MocompRngShiftIndex=0;
  RncCmd.RngShiftInterpSize=8;
  RncCmd.MocompPhaseSign=MocompPhaseSign;
  strcpy(RncCmd.MocompFileName,MOCOMP_FILE);
  RncCmd.MocompRngUpdates=1;
  RncCmd.MocompStartPRI=0;
  RncCmd.RefFuncPhaseSign=-1;
  RncCmd.ChirpBandwidth=-1.0;   /* no pulse compression, for now */
  RncCmd.PulseLen=IMOPulseLen;  /* not actually used for now */
  RncCmd.SRCFlg='N';
  RncCmd.NomGroundSpeed=SALNomGroundSpeed;  
  RncCmd.SquintAngle=0.0; 
  RncCmd.SRCFocusRng=0.0; 
  RncCmd.RngWalkRngShiftFlg='N';
  RncCmd.RngWalkPhaseShiftFlg='N';
  RncCmd.AzBeamwidth=SALAzBeamwidth;
  RncCmd.CarrierFreq=IMOCarrierFreq;
  RncCmd.A2DFreq=SALADFreq;
  RncCmd.InputDCOffsetI=-999.0;
  RncCmd.InputDCOffsetQ=-999.0
  RncCmd.RngFFTSize=8192;
  RncCmd.WinConstTime=0.08;
  RncCmd.InputDataType=0;
  RncCmd.OutputDataType=3;
  RncCmd.Scale=1.0;
      

  /* Set up default corner turn command file params (assumes Rngcom used) */
  strcpy(CorCmd.InputFileName,RNG_COMPRESSED_FILE); 
  strcpy(CorCmd.OutputFileName,CORNER_TURNED_FILE);
  CorCmd.BytesPerValue = 2;
  CorCmd.RowHeaderBytesToSkip = 0;
  CorCmd.RowFooterBytesToSkip = 0;
  CorCmd.InputStartRow = 0;
  CorCmd.InputRowsToUse = PRIsPerAzBlock;
  CorCmd.MaxInputRows = PRIsPerAzBlock;
  CorCmd.MaxInputColumns = IMORngBins;
  CorCmd.InputStartColumn = 0;
  CorCmd.ColumnsPerProcSegment = IMORngBins; /*process all range bins*/
  CorCmd.NumProcSegments = 1;

  if (MocompFlg != 'Y') /* no mocomp, hence no rng compression function */
    {
    /* Alter corner turn command file params so as to skip rng compression
       stage if no mocomp */
    CorCmd.InputStartRow = 1+StartPRIOffset; /* skip 1st record of IMOP file*/ 
    CorCmd.RowHeaderBytesToSkip = 412;
    CorCmd.RowFooterBytesToSkip = IMORecordLen-412-2*IMORngBins;
    CorCmd.MaxInputRows = IMODataRecords+1; /* +1 for the file decrip rec*/
    CorCmd.MaxInputColumns = IMORngBins;
    strcpy(CorCmd.InputFileName,RawName);
    }
  
  /* Set up default azimuth compression command file params */
  AzcCmd.ScreenUpdateRate = 10;
  AzcCmd.InputStartSampleDelay = IMODelay0; /* Note for future:  change if rng
					       compression module used */
  AzcCmd.CarrierFreq = IMOCarrierFreq;
  AzcCmd.InputPRF = SALNomPRF/(double)SALPresumRatio;
  AzcCmd.NomGroundSpeed = SALNomGroundSpeed;
  AzcCmd.InputFileAzPts = PRIsPerAzBlock;
  AzcCmd.StartProcessAzPt = 0;
  AzcCmd.AzPtsToProcess = PRIsPerAzBlock;
  AzcCmd.InputFileRngBins = IMORngBins;
  AzcCmd.StartProcessRngBin = AzProcStartRngBin;
  AzcCmd.RngBinsToProcess = AzProcRngBins;
  AzcCmd.InputDCOffset = IMODCOffset;
  AzcCmd.FFTSize = AzFFTSize;
  AzcCmd.InvFFTSizeReduc = 1; /* no IFFT size reduc */
  strcpy(AzcCmd.InputFileName,CORNER_TURNED_FILE);
  strcpy(AzcCmd.OutputFileName,AZ_COMPRESSED_FILE);
  AzcCmd.AppendExistOutFileFlg = 'N';  
  AzcCmd.RngFocSegments = RngFocSegments;
  AzcCmd.RefFuncSign = AzRefFuncSign;
  AzcCmd.A2DFreq = SALADFreq;
  AzcCmd.NomAzRes = NomAzRes;
  AzcCmd.WinConstTime = AzWinConstTime;
  AzcCmd.NumLooks = AzLooks;
  AzcCmd.LookOverlapFrac = LookOverlapFrac;
  AzcCmd.WinConstFreq = AzWinConstFreq;
  AzcCmd.RngCurvInterpSize = RngCurvInterpSize;
  AzcCmd.RngCurvBatchSize = RngCurvBatchSize;
  AzcCmd.PostSumRatio = 1; /* no post summing */
  AzcCmd.DetectMethod = DetectMethod;
  AzcCmd.OutputDataType = OutputDataType;
  AzcCmd.Scale = AzScale;
  AzcCmd.ReportMax = 1; /* report max values */


  
  /***************************/
  /* PERFORM CORE PROCESSING */

  /*------------------------------------------------------------*/
  /* REPEAT PROCESSING FOR EACH AZ BLOCK (typically SYN product)
     else NumAzBlocks=1 (other products) */
  for (AzBlock=0; AzBlock<NumAzBlocks;AzBlock++)
    {
    if (NumAzBlocks != 1)
      fprintf(msg,"\nProcessing az block %ld [%ld remain(s)]...\n",
            AzBlock,NumAzBlocks-(AzBlock+1));

    /* ONLY PERFORM RNG COMPRESSION STAGE IF MOCOMP SELECTED */
    if (MocompFlg == 'Y')
      {

      /* Update rnc start posn, from second block */
      if (AzBlock != 0) RncCmd.StartProcessPRI += PRIsPerAzBlock;
	  
      /* Write rnc command file */
      if (!WriteRncCmdFile(&RncCmd,RncCmdName))
        {
        fprintf(msg,"ERROR - Rng compression command file not written!\n");
        exit(1);
        }

      /* PERFORM RNG COMPRESSION (or mocomp correction) */
      G2Rnc(msg,RncCmdName);

      }   /* end if mocomp */
    else  /* IF NO MOCOMP (i.e. RNGCOM NOT USED) */
      {
      /* Update corner start posn, from second block */
      if (AzBlock != 0) CorCmd.InputStartRow += PRIsPerAzBlock;
      }
	
    /* Write corner turn command file */
    if (!WriteCorCmdFile(&CorCmd,CorCmdName))
      {
      fprintf(msg,"ERROR - Corner turn command file not written!\n");
      exit(1);
      }

    /* PERFORM CORNER TURN FOR AZ BLOCK */
    G2Cor(msg,CorCmdName);

    /* Write azimuth compression command file */
    if (AzBlock!=0) AzcCmd.AppendExistOutFileFlg = 'Y';
    if (!WriteAzcCmdFile(&AzcCmd,AzcCmdName))
      {
      fprintf(msg,"ERROR - Az compression command file not written!\n");
      exit(1);
      }     

    /* PERFORM AZIMUTH COMPRESSION FOR AZ BLOCK*/
    if (!G2Azc(msg,AzcCmdName) )
      {
      fprintf(msg,"ERROR - Az compression abnormal exit!!\n");
      exit(1);
      };

    }  /* END FOR AZBLOCK LOOP */
  /*---------------------------*/


  /* UNSCRAMBLE AZC OUTPUT FILE */  
  if (NumAzBlocks > 1)
    {
    if (!UnBlockAzcFile(msg,
                        AZ_COMPRESSED_FILE,
                        UNBLOCK_AZC_FILE,
                        DetectMethod,
                        OutputDataType,
                        AzProcRngBins,
                        PRIsPerAzBlock,
                        NumAzBlocks ) )
      {
      fprintf(msg,"ERROR - in unscrambling the az compressed file!\n");
      exit(1);
      } 
    strcpy(AzcFileName,UNBLOCK_AZC_FILE);
    }
  else
    {
    strcpy(AzcFileName,AZ_COMPRESSED_FILE);
    }

  
  /* PROJECT IMAGE TO GROUND RANGE, if appropriate */
  if ( (DetectMethod==2) && (GrndRngInterpPts!=0) )
    {   
    /* Set up corner turn command file */
    strcpy(CorCmd.InputFileName,AzcFileName); 
    strcpy(CorCmd.OutputFileName,CORNER_TURNED_FILE);
    if (OutputDataType==0) 
      { CorCmd.BytesPerValue = 1; }
    else if (OutputDataType==1) 
      { CorCmd.BytesPerValue = 2; } 
    else if (OutputDataType==2) 
      { CorCmd.BytesPerValue = 4; } 
    else if (OutputDataType==3) 
      { CorCmd.BytesPerValue = 4; }
    else if (OutputDataType==4) 
      { CorCmd.BytesPerValue = 8; }  
    CorCmd.RowHeaderBytesToSkip = 0;
    CorCmd.RowFooterBytesToSkip = 0;
    CorCmd.InputStartRow = 0;
    CorCmd.InputRowsToUse = AzProcRngBins;
    CorCmd.MaxInputRows = AzProcRngBins;
    CorCmd.MaxInputColumns = PRIsPerAzBlock*NumAzBlocks;
    CorCmd.InputStartColumn = 0;
    CorCmd.ColumnsPerProcSegment = PRIsPerAzBlock;
    CorCmd.NumProcSegments = NumAzBlocks;
    /* Write corner turn command file */
    if (!WriteCorCmdFile(&CorCmd,CorCmdName))
      {fprintf(msg,"ERROR - Corner turn cmd file not written!\n"); exit(1);}
    /* Perform corner turn to arrange data for ground range projection */
    G2Cor(msg,CorCmdName);

    /* Perform ground range projection */
    if (ProjImageToGrndRng(msg,
		       CORNER_TURNED_FILE, /* fixed 1998-01-21 */
                       GRNG_FILE,
                       PRIsPerAzBlock*NumAzBlocks,
                       AzProcRngBins,
                       OutputDataType,
		       (double)(C*IMODelay0/2.0),
		       (double)(C/(2.0*SALADFreq)),
                       &GrndRngBinSpacing,
  /* TEMP, for ATP!!! */    NOM_HEIGHT,
			   /*     IMOMidAlt-SALAveTerrainAlt, */
		       GrndRngInterpPts))
      GrndRngProjected = 1;
    } /* end if grnd rng projection */
  else if ( (DetectMethod!=2) && (GrndRngInterpPts!=0) )
    { 
    fprintf(msg,"WARNING - correct grnd rng interp requires power image!\n");
    fprintf(msg,"No ground rng projection performed!\n");
    }


  /* WRITE SUN RASTER FILE, if required  */
  if (SUNRasterFlg == 'Y')
    { 
    strcpy(RasterFileName,RunID);
    strcat(RasterFileName,".ras");
    if (GrndRngProjected) /* then image transposed */
      {
      strcpy(tmps,GRNG_FILE);
      WidthPixels = AzProcRngBins;
      HeightPixels = PRIsPerAzBlock*NumAzBlocks;
      }
    else  /* image orientation as out of azc */
      {
      strcpy(tmps,AzcFileName);
      WidthPixels = PRIsPerAzBlock*NumAzBlocks;
      HeightPixels = AzProcRngBins;
      }
    if (!WriteSUNRasterFile(msg,tmps,RasterFileName,OutputDataType,
			  DetectMethod, WidthPixels,HeightPixels))
      fprintf(msg,"WARNING - SUN raster file not correctly written!\n");
    } /* end if SUNRasterflg =  Y */
    

  /* WRITE AIRSAR FORMAT IMAGE FILE, if required  */
  if (AIRSARFileFlg == 'Y')
    {
    strcpy(AIRSARFileName,RunID);
    strcat(AIRSARFileName,".air");
    if (GrndRngProjected) /* then image transposed */
      {
      strcpy(tmps,GRNG_FILE);
      WidthPixels = AzProcRngBins;
      HeightPixels = PRIsPerAzBlock*NumAzBlocks;
      }
    else  /* not grnd rng projected so image orientation as out of azc */
      {
      strcpy(tmps,AzcFileName);
      WidthPixels = PRIsPerAzBlock*NumAzBlocks;
      HeightPixels = AzProcRngBins;
      }

    /* Set up AIRSAR format header info */
    /* First header */
    AIRHead1.RecordLen = 0;    /* calc in function */
    AIRHead1.NumHeaderRecords = 2; /* new header and parameter header */
    AIRHead1.SamplesPerRecord = WidthPixels; 
    AIRHead1.LinesInImage = HeightPixels;
    AIRHead1.NumBytesPerSample = 0; /* calc in function */
    strcpy(AIRHead1.SASARProcVersion,PROG_VERSION);
    strcpy(AIRHead1.DataType," "); /* assign in function */
    if (GrndRngProjected)
      { strcpy(AIRHead1.RangeProjection,"GROUND"); }
    else 
      { strcpy(AIRHead1.RangeProjection,"SLANT"); }      
    if (GrndRngProjected)
      { AIRHead1.RangePixelSpacing = GrndRngBinSpacing; }
    else 
      { AIRHead1.RangePixelSpacing = (double)(C/(2.0*SALADFreq)); }
    AIRHead1.AzPixelSpacing = 
        SALNomGroundSpeed*(double)SALPresumRatio/SALNomPRF;
    AIRHead1.ByteOffsetOldHeader = 0;
    AIRHead1.ByteOffsetUserHeader = 0;
    AIRHead1.ByteOffset1stDataRecord = 0; /* calc in function */
    AIRHead1.ByteOffsetParamHeader = 0; /* calc in function */
    if (GrndRngProjected)
      { strcpy(AIRHead1.DataLineFormat,"RANGE"); }
    else 
      { strcpy(AIRHead1.DataLineFormat,"AZIMUTH"); }
    AIRHead1.ByteOffsetCalHeader = 0;    
    AIRHead1.ByteOffsetDEMHeader = 0;

    /* Parameter header */
    strcpy(AIRParam.HeaderName,"PARAMETER");
    strcpy(AIRParam.SiteName,SALSiteName);
    AIRParam.SiteLatitude = IMOMidLat;
    AIRParam.SiteLongitude = IMOMidLong;
    strcpy(AIRParam.ImageTitle," "); /* field 5 */
    AIRParam.CDROMId = DEFAULT_VAL_I;
    if (IMOCarrierFreq>100.0e+06 && IMOCarrierFreq<200.0e+06)
      { strcpy(AIRParam.Freq,"VHF"); }
    else 
      { strcpy(AIRParam.Freq," "); }
    strcpy(AIRParam.Polarization,IMOPolarization);
    strcpy(AIRParam.DataType," "); /* assign in function */
    AIRParam.DataId = DEFAULT_VAL_I; /* field 10 */
    AIRParam.ArchivalFlg = 0;
    AIRParam.CDROMStartPRIId = IMOFirstPRIID;
    AIRParam.ProcStartPRIId = IMOFirstPRIID+StartPRIOffset;
    AIRParam.LatitudeStartScene = IMOStartLat;
    AIRParam.LongitudeStartScene = IMOStartLong; /* field 15 */
    AIRParam.LatitudeEndScene = IMOEndLat;
    AIRParam.LongitudeEndScene = IMOEndLong;
    AIRParam.ApproxStartHDDTFootage = DEFAULT_VAL_I;
    sprintf(AIRParam.AcquisitionDate,"%ld",IMOAcqYear);
    AIRParam.AcquisitionDay = IMOAcqDay;  /* field 20 */
    AIRParam.AcquisitionSec = IMOAcqSec;  
    AIRParam.RecordWindowDuration = 1.0e+06*IMORngBins/SALADFreq;
    strcpy(AIRParam.FreqsCollected," ");
    AIRParam.DigitalDelay = IMODelay0*1.0e+06;
    AIRParam.ChirpDelay = DEFAULT_VAL_F;  /* field 25 */
    AIRParam.ProcDelay = AzProcStartRngBin;
    AIRParam.PRFStartTransfer = AzcCmd.InputPRF;;
    AIRParam.SamplingRate = SALADFreq*1.0e-06;
    AIRParam.CentreFreq = AzcCmd.CarrierFreq*1.0e-06;
    AIRParam.ChirpBandwidth = DEFAULT_VAL_F; /* field 30 */
    strcpy(AIRParam.TypeOfChirp,"NONE");
    AIRParam.PulseLen = IMOPulseLen*1.0e+06;
    AIRParam.ProcWavelen = C/AzcCmd.CarrierFreq;
    AIRParam.BaroAltitude = DEFAULT_VAL_F;
    AIRParam.RadarAltimeterAlt = DEFAULT_VAL_F; /* field 35 */
    AIRParam.ProcAltitude = IMOMidAlt-SALAveTerrainAlt;
    AIRParam.ElevInvestigatorSite = DEFAULT_VAL_F;
    AIRParam.AircraftTrackAngle = DEFAULT_VAL_F;
    AIRParam.AircraftYawAngle = DEFAULT_VAL_F;
    AIRParam.AircraftPitchAngle = DEFAULT_VAL_F;  /* field 40 */
    AIRParam.AircraftRollAngle = DEFAULT_VAL_F;
    AIRParam.ProcYawAngle = DEFAULT_VAL_F;
    AIRParam.ProcPitchAngle = DEFAULT_VAL_F;
    AIRParam.ProcRollAngle = DEFAULT_VAL_F;
    AIRParam.NomPRFRatioHzPerKnot = DEFAULT_VAL_F; /* field 45 */
    AIRParam.NomPRFRatioPerMetre = AzcCmd.InputPRF/AzcCmd.NomGroundSpeed;
    AIRParam.PRFCorrectionFactor = DEFAULT_VAL_F;
    AIRParam.RngFFTSize = DEFAULT_VAL_I;
    AIRParam.AzFFTSize = AzcCmd.FFTSize;
    AIRParam.FrameSize = DEFAULT_VAL_I;  /* field 50 */
    AIRParam.NumFramesProcessed = DEFAULT_VAL_F;
    AIRParam.RngAlignDelayHH = DEFAULT_VAL_F;
    AIRParam.RngAlignDelayHV = DEFAULT_VAL_F;
    AIRParam.RngAlignDelayVH = DEFAULT_VAL_F;
    AIRParam.RngAlignDelayVV = DEFAULT_VAL_F; /* field 55 */
    AIRParam.NearSlantRng = (double)
         (C*(AzcCmd.InputStartSampleDelay+
             AzcCmd.StartProcessRngBin/AzcCmd.A2DFreq)/2.0);
    AIRParam.FarSlantRng =  AIRParam.NearSlantRng+
           (double)(C*(AzcCmd.RngBinsToProcess-1)/(AzcCmd.A2DFreq*2.0));
    AIRParam.NearLookAngle = (180.0/PI)*acos(AIRParam.ProcAltitude/
                             AIRParam.NearSlantRng);
    AIRParam.FarLookAngle = (180.0/PI)*acos(AIRParam.ProcAltitude/
                             AIRParam.FarSlantRng);
    AIRParam.NumProcAzLooks =                /* field 60 */
      (double)AzLooks - ((double)AzLooks-1.0)*LookOverlapFrac; 
    AIRParam.NumProcRngLooks = DEFAULT_VAL_F;
    strcpy(AIRParam.RngWeighting," ");
    AIRParam.RngWeightingCoeff = DEFAULT_VAL_F;
    strcpy(AIRParam.AzWeighting," ");
    AIRParam.AzWeightingCoeff = DEFAULT_VAL_F; /* field 65 */
    AIRParam.PercentPRFBWProc = DEFAULT_VAL_F;
    AIRParam.DeskewFlg = DEFAULT_VAL_I;
    AIRParam.SlantRngSampleSpacing = C/(2.0*AzcCmd.A2DFreq);
    AIRParam.NomSlantRngRes = C*IMOPulseLen/2.0;;
    AIRParam.AzSampleSpacing = AIRHead1.AzPixelSpacing;  /* field 70 */
    AIRParam.NomAzRes = NomAzRes*AzWinTimeBroadenFactor;
    AIRParam.NumInterpPtsRMC = AzcCmd.RngCurvInterpSize;
    AIRParam.AzRefSizePerLookNear = DEFAULT_VAL_I;
    AIRParam.AzRefSizePerLookFar = DEFAULT_VAL_I;
    AIRParam.ImageCentreLat = DEFAULT_VAL_F;
    AIRParam.ImageCentreLong = DEFAULT_VAL_F;
    AIRParam.CalToneVideoFreq = DEFAULT_VAL_F;
    AIRParam.CalTonePowerHH = DEFAULT_VAL_F;
    AIRParam.CalTonePowerHV = DEFAULT_VAL_F;
    AIRParam.CalTonePowerVH = DEFAULT_VAL_F;
    AIRParam.CalTonePowerVV = DEFAULT_VAL_F;
    AIRParam.CalibFactorHH = DEFAULT_VAL_F;
    AIRParam.CalibFactorHV = DEFAULT_VAL_F;
    AIRParam.CalibFactorVH = DEFAULT_VAL_F;
    AIRParam.CalibFactorVV = DEFAULT_VAL_F;
    AIRParam.MeasCorrHVVHPowerRatio = DEFAULT_VAL_F;
    AIRParam.MeasCorrHVVHPhase = DEFAULT_VAL_F;
    AIRParam.CalTonePhaseMeasHH = DEFAULT_VAL_F;
    AIRParam.CalTonePhaseMeasHV = DEFAULT_VAL_F;
    AIRParam.CalTonePhaseMeasVH = DEFAULT_VAL_F;
    AIRParam.CalTonePhaseMeasVV = DEFAULT_VAL_F;
    AIRParam.GenScaleFactor = DEFAULT_VAL_F;
    AIRParam.GPSAlt = DEFAULT_VAL_F;     

    
    if (!WriteAIRSARFormat(msg,tmps,AIRSARFileName,OutputDataType,
			 DetectMethod,WidthPixels,HeightPixels,
			 &AIRHead1,&AIRParam))
      fprintf(msg,
      "WARNING - AIRSAR image format file not correctly written!\n");
    } /* end if AIRSARFileFlg = Y */
 
 
  /* Tidy up */
  fclose(LBRFile);
  fclose(RawFile);
  if (strcmp(Messages,"file")==0) fclose(logfp);
  
} /* end main prog */


/************************************************************************/
/* FUNCTION : WriteG2TmplCmdFile() */
Int2B WriteG2TmplCmdFile(char CmdFileName[])
{
  FILE *Out; 
  if ( (Out = fopen (CmdFileName, "wt") ) == NULL ) return(0);

  fprintf(Out,"G2 SAR Processor command file (SASAR System)\n");
  fprintf(Out,"$ProgramVersion (jmh) => %s\n\n\n",PROG_VERSION);

  fprintf(Out,
   "Processing run parameters (fields change from run to run):\n");
  fprintf(Out,
   "----------------------------------------------------------\n");
  fprintf(Out,"$RunID                       => sas01\n");
  fprintf(Out,"$Messages (file/displ)       => displ\n");
  fprintf(Out,"$LBRFile                     => lbr1b.atp\n");
  fprintf(Out,"$RawDataFile                 => ivh01-H4.000\n");
  fprintf(Out,"$HBR                         => 4\n");
  fprintf(Out,"$StartPRIOffset              => 0\n\n");

  fprintf(Out,"Fields specified for product (usually fixed):\n");
  fprintf(Out,"---------------------------------------------\n");
  fprintf(Out,"$NumAzBlocks [0-full data set]          => 1\n");
  fprintf(Out,"$PRIsPerAzBlock []                      => 1001\n");

  fprintf(Out,"

  fprintf(Out,"$MocompFlg [(Y/N)]                      => N\n");
  fprintf(Out,"$MocompPhaseSign [(+-1)]                => -1\n");
  
  fprintf(Out,"$AzProcStartRngBin []                   => 0\n");
  fprintf(Out,"$AzProcRngBins [0 - all]                => 0\n");   
  fprintf(Out,"$AzFFTSize [power of 2]                 => 2048\n");
  fprintf(Out,"$RngFocSegments [-1 for max]            => 128\n");
  fprintf(Out,"$AzRefFuncSign [(+-1)]                  => -1\n"); 
  fprintf(Out,"$NomAzRes [(m)- see note]               => 15.0\n");
  fprintf(Out,"$AzWinConstTime []                      => 0.08\n");
  fprintf(Out,"$AzWinConstFreq []                      => 1.0\n");
  fprintf(Out,"$AzLooks []                             => 1\n");
  fprintf(Out,"$LookOverlapFrac []                     => 0.5\n");
  fprintf(Out,"$RngCurvInterpSize []                   => 4\n");
  fprintf(Out,"$RngCurvBatchSize []                    => 256\n");
  fprintf(Out,"$DetectMethod [see note]                => 2\n");
  fprintf(Out,"$OutputDataType [see note]              => 0\n");  
  fprintf(Out,"$AzScale [multiplicative]               => 3.0e-7\n");
  fprintf(Out,"$GrndRngInterpPts [see note]            => 0\n"); 
  fprintf(Out,"$SUNRasterFlg [(Y/N)]                   => Y\n");
  fprintf(Out,"$AIRSARFileFlg [(Y/N) - see note]       => N\n\n");

  fprintf(Out,"Notes:\n");
  fprintf(Out,"------\n\n");
  fprintf(Out,
    "NomAzRes : The nominal azimuth resolution in metres assuming a\n");
  fprintf(Out,
    "           window broadening factor of unity.\n\n");
  fprintf(Out,"RngCurvInterpSize: 0 - no rng curv correction\n");
  fprintf(Out,"                   1 - nearest neighbour\n");
  fprintf(Out,"                   else, even positive\n\n");

  fprintf(Out,"DetectMethod: 0 - complex output\n");
  fprintf(Out,"              1 - magnitude\n");
  fprintf(Out,
   "              2 - power (obligatory for ground range projection)\n");
  fprintf(Out,
   "              3 - power in dB\n\n");
  fprintf(Out,
   "OutputDataType: 0 - unsigned char (1 byte) (DCOff 127 - complex)\n");
  fprintf(Out,
   "                1 - unsigned short int (2 bytes) (DCOff 32767 - complex)\n");
  fprintf(Out,"                2 - long int (4 bytes)\n");
  fprintf(Out,"                3 - float (4 bytes)\n");
  fprintf(Out,"                4 - double (8 bytes)\n\n");

  fprintf(Out,"GrndRngInterpPts: 0 - no grnd rng projection\n");
  fprintf(Out,"                  1 - nearest neighbour\n");
  fprintf(Out,
   "                  else, even positive (the more pts, the better)\n");
  fprintf(Out,"                  Note: requires power detection!\n\n");

  fclose(Out);
  return(1);
} 



