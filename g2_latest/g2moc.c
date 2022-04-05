/*==========================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2moc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Motion compensation module for SASAR data. 
Assumed Marconi FIN3110 IMU. Two major functions:

1) Print out summary log file of the motion data recorded in the input
LBR file for the data channel. Variables associated with this function
usually start with "LBR". Note that for some unknown reason, the 
"Sync IDA" of the last motion record of the  SASAR LBR file is not found. 
Seems to be a glitch in the way the record is written as, if the LBR file
is hacked to add another record at the end, it still misses the old last 
record, but picks up the new last record.

2) Calculate the motion comp shifts required for SASAR data. The 
motion corrections are calculated over a specified subset of the PRI
range of the data channel and the variables associated with this usually
start with "Proc". The motion corrections make use of the IMU Inert lat, 
Inert long and INT height output. The lat, long and height values are 
transformed to cartesian coords with some vector maths needed to obtain 
the ground reference track in the swath etc. The Maple package was used 
to perform the symbolic manipulations programmed here (i.e. don't expect 
to easily understand the maths straight from the program). Corrections 
set to zero if no data exists (at start say) and equal to final 
correction value at end, if processing past end of motion data. Else 
linear interpolation between the range shifts is used.

The range shift output is written to binary file in the form of an Int4B  
containing the G2 PRI ID (increments by one) and a float corresponding to 
the range shift required. A positive range shift implies that the data 
should be shifted away. In addition, there is an option to write the 
G2 PRI ID and range shift to a text file (for plotting purposes).

General philosophy is not to trust all the LBR file parameters. If 
a parameter has been already specified for this prog, then rather 
use that.

Note that the LBR PRI IDs for a data channel are actually disk sector
IDs (where each sector is 512 bytes) and thus do not increment by one
(typically 8 for 2048 rng bins, 16 for 4096 rng bins). In this file,
the PRI ID's follow the SASAR scheme. This is different to the rest
of the G2 SAR processor where PRI IDs increment by one for each stored 
range line. Thus, for the output of this prog to be compatible with
the rest of the processor, the PRI IDs written to the binary output
file increment by one and are referenced such that the LBR file start PRI
ID is zero.

VERSION/AUTHOR/DATE : 1999-02-04 / Jasper Horrell / 1999-02-04
COMMENTS: 
First working version. Change from use of INT Lat and Long to Inert Lat 
and Long (keep INT hgt). 

VERSION/AUTHOR/DATE : 1999-07-20 / Jasper Horrell / 1999-07-20
COMMENTS:
Cosmetic changes to output. Last version prior to incorporation of DGPS
data.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-07-23
COMMENTS:
Read in DGPS data (converted to readable format with correct LBR PRI IDs)
from file "DGPS.dat". Uses this lat,long, alt data instead of the IMU data.

VERSION/AUTHOR/DATE : 0.3 / Jasper Horrell / 1999-07-29
COMMENTS:
Read in G2-unpacked DGPS file name. Read in delay0 rather than rely on LBR file.

VERSION/AUTHOR/DATE : 0.4 / Jasper Horrell / 1999-08-02
COMMENTS:
Write log file. Remove LBR_MOT_SHIFT write compile option. Add in option to 
specify A/D frequency and change messages (removed some - not relevant here).
(1999-08-14) Compile using height info.

VERSION/AUTHOR/DATE : 0.5 / Jasper Horrell / 1999-08-26
COMMENTS:
Clean up and remove requirement to read LBR file. Passes parameters via command
instead. Intended for compatibility with G2 processor glued with python. Motion
records read from ASCII file (LBR PRI, UTCsecs, LatDeg, LongDeg, AltMetre). This
input motion record file probably created by combining IMU and DGPS data. Removed
possibility of compiling as a function.

VERSION/AUTHOR/DATE : 0.6 / Jasper Horrell / 1999-09-15
COMMENTS:
Add in ability to write out the image corner positions in geodetic coordinates.
Add checks in function to calculate ground ref point. Handle special case of
calculating a zero look angle and also calculating a "ground" reference point 
when there is a zero look angle.

2000-03-07 Change to track time from G2 PRI in output text file (for EUSAR 2000
paper plot).

=========================================*/

#include "g2func.h"

#define PROG_VERSION "0.6"
#define DEFAULT_OUT_EXT "_moc"
#define DEFAULT_LOG_EXT "_mlog"
#define WGS84_a 6378137.0
#define WGS84_b 6356752.314
#define WGS84_e_sq 6.694380066e-03
#define WGS84_ep_sq 6.73949681994e-03
#define TO_RAD 0.01745329252
#define TO_DEG 57.29577951


/* Structure definitions for this file */
struct cartes {
  double x;
  double y;
  double z;
 };

struct geodetic {
  double latrad;   /* latitude in radians */
  double longrad;  /* longitude in radians */
  double alt; /* altitude above ellipsoid in metres */
};

struct MotStruct {
  Int4B FilePosn;     /* byte position of start of record */
  Int4B LP;          /* LBR PRI ID (i.e. increments by 16 for 4096 rng bins) */
  double UTCsecs;
  double Lat;  /* in degrees */
  double Long; /* in degrees */
  double Alt;  /* in metres */
};


/* Function prototypes for functions in this file */
Int2B ConvToCartesECEF(double LatRad, double LongRad, double h,
                   struct cartes *Pt); 
Int2B ConvToGeodetic(struct cartes Pt,struct geodetic *GPosn);
double CalcLookAngleRad(double heightAGL, double slantrng);               
double CalcCartesDist(struct cartes p1,struct cartes p2);
Int2B CalcGrndRefPt(char Dirn, 
                       double NomSlantRng,
                       double LookAngleRad, 
                       struct cartes FlightRef, 
                       struct cartes NomFlight,
                       struct cartes *GrndRefPt);
Int2B QuadRoots(double a, double b, double c, 
                double *root1, double *root2);
Int2B CalcR(double h, double k, double m, 
            double q, double r, double s,
            double p, double u, 
            struct cartes *R);   
Int2B CrossDot(struct cartes FlightRef,
               struct cartes R,
               struct cartes NomFlight,
               double *CrossCalc);
Int2B ReadMotRec(FILE *MotFile,
                       struct MotStruct *MotData);                       
Int2B MovingAveSmooth(double In[], double Out[], 
                      Int4B ArrSize, Int4B KernSize);

/*===========MAIN PROG/FUNCTION===========*/

Int2B main (int argc,char *argv[])
{
  FILE *msg=stdout,*OutFile=NULL,*OutTextFile=NULL,*MotFile=NULL,*LogFile;
  char OutFileName[STRING_SPACE]="",OutTextFileName[STRING_SPACE]="",
       MotFileName[STRING_SPACE]="",LogFileName[STRING_SPACE]="",
       OutTextFileFlg = 'N', /* for range shifts */
       AntennaDirn = 'R'; /* default value */
  Int4B DataStartLP=DEFAULT_VAL_I,
        EndG2PRI, /* referenced from start of raw data file */
        errors,
        FarProcRngBin=4095,
        i,
        KernSize=1,  /* moving average kernel size */
        LPIncr=DEFAULT_VAL_I,
        MaxRngShiftIndx=DEFAULT_VAL_I,
        MotRec,
        NearProcRngBin=0,
        OutG2PRI,
        ProcEndFound=0,
        ProcEndLP=DEFAULT_VAL_I,
        *ProcMotShiftLP, /* LBR PRI IDs for elements of ProcMotShift array */
        ProcNumMotRecords=0,
        ProcStartFound=0,       
        ProcStartLP=DEFAULT_VAL_I,
        PRIChange,
        RefRngBin=DEFAULT_VAL_I,
        StartG2PRI; /* referenced from start of raw data file */
  float MaxRngShift=0.0,
        RngShift;      
  double A2DFreq=DEFAULT_VAL_F,
         Delay0=DEFAULT_VAL_F,
         FarLookAngleRad,
         FarProcSlantRng,
         NearLookAngleRad,
         NearProcSlantRng,
         NomLookAngleRad,
         NomSlantRng,
         PRF=DEFAULT_VAL_F, /* of stored data (Hz) */
         ProcAveAlt,
         ProcMotDist,
         *ProcMotShift, /* holds the rng shifts for the motion records */
         *ProcSmoothShift, /* smoothed range shift for the mot records */  
         ProcMotTime,
         ProcAveGrndSpeed,
         RngShiftCur,
         RngShiftIncr,
         RngShiftNext,    
         TerrainAlt=DEFAULT_VAL_F;
  struct MotStruct	ProcMotDataStart,
  					ProcMotDataEnd,
  					MotDataCur,
  					MotDataNext;                 
  struct cartes p1,p2,
                AircraftPosn,
                GrndRefPt,
                GrndRefStartPt,
                ImFarEndPt,
                ImFarStartPt,
                ImNearEndPt,
                ImNearStartPt,
                ImRefEndPt,
                ImRefStartPt,
                NomFlight,
                RefEndPt,
                RefStartPt,
                RefM;

  struct geodetic GrndRefEndGPosn,
                  GrndRefStartGPosn,
                  ImFarEndGPosn,
                  ImFarStartGPosn,
                  ImNearEndGPosn,
                  ImNearStartGPosn,
                  ImRefEndGPosn,    /* extrapolated platform posn for image */
                  ImRefStartGPosn;  /* extrapolated platform posn for image */
           
  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n------------\n");
  fprintf(msg,"Prog: MOCOMP (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1999\n");

  if (argc < 8) {
    fprintf(msg,
      "Calculates motion comp range shifts for SASAR data.\n\n");
    fprintf(msg,"USAGE : mocomp [required params] <opt params>\n");
    fprintf(msg,
      "e.g. 'mocomp MotF=mot.txt Delay0=2.0e-6 etc...'\n\n");
    fprintf(msg,"Required:\n");
    fprintf(msg,
      "<MotF=m>     - motion file of LBR PRI, UTCSec, LatDeg, LongDeg, AltM\n");
    fprintf(msg,
      "<Delay0=m>   - delay to first range sample (secs)\n");
    fprintf(msg,
      "<RefBin=n>   - the range bin to which mocomp is calculated\n");
    fprintf(msg,
      "<A2DFreq=m>  - A/D frequency (Hz)\n");  
    fprintf(msg,
      "<PRF=m>      - Pulse repetition frequency of stored raw data (Hz)\n");  
    fprintf(msg,
      "<DataStartLP=n> - LBR PRI corresponding to start of raw data\n");
    fprintf(msg,
      "<ProcStartLP=n> - start LBR PRI to process (disk sector)\n");
    fprintf(msg,
      "<ProcEndLP=n>   - end LBR PRI to process (disk sector)\n");      
    fprintf(msg,
      "<LPIncr=n>      - LBR PRI increment\n\n");      
   
    fprintf(msg,"Optional:\n");      
    fprintf(msg,
      "<LogF=t>       - log file name (else default name)\n");
    fprintf(msg,
      "<OutF=t>       - the output rng shift file (else default name)\n");
    fprintf(msg,
      "<OutTxtF=t>    - the output rng shift text file (else not written)\n");
    fprintf(msg,
      "<KernS=n>      - rng shift moving ave smooth kernel size (default 1)\n");
    fprintf(msg,
      "<AntDirn=L/R>  - antenna direction (default R)\n");
    fprintf(msg,         
      "<TerAlt=m>     - average terrain altitude ASL (m) (default zero)\n");
    fprintf(msg,         
      "<NearRngBin=m> - near rng bin processed - for geocoding (default 0)\n");
    fprintf(msg,         
      "<FarRngBin=m>  - far rng bin processed - for geocoding (default 4095)\n");

    exit(1);
  } 

  /* read command line parameters */
  for (i=1; i<argc; i++) {
    if (strncmp(argv[i],"MotF=",5)==0) {
      sscanf(argv[i],"MotF=%s",MotFileName);
    }    
    else if (strncmp(argv[i],"Delay0=",7)==0) {
      sscanf(argv[i],"Delay0=%lf",&Delay0);
    }
    else if (strncmp(argv[i],"RefBin=",7)==0) {
      sscanf(argv[i],"RefBin=%ld",&RefRngBin);
    }
    else if (strncmp(argv[i],"A2DFreq=",8)==0) {
      sscanf(argv[i],"A2DFreq=%lf",&A2DFreq);
    }
    else if (strncmp(argv[i],"PRF=",4)==0) {
      sscanf(argv[i],"PRF=%lf",&PRF);
    }
    else if (strncmp(argv[i],"DataStartLP=",12)==0) {
      sscanf(argv[i],"DataStartLP=%ld",&DataStartLP);
    }
    else if (strncmp(argv[i],"ProcStartLP=",12)==0) {
      sscanf(argv[i],"ProcStartLP=%ld",&ProcStartLP);
    }
    else if (strncmp(argv[i],"ProcEndLP=",10)==0) {
      sscanf(argv[i],"ProcEndLP=%ld",&ProcEndLP);
    }    
    else if (strncmp(argv[i],"LPIncr=",7)==0) {
      sscanf(argv[i],"LPIncr=%ld",&LPIncr);
    }    

    else if (strncmp(argv[i],"LogF=",5)==0) {
      sscanf(argv[i],"LogF=%s",LogFileName);
    }    
    else if (strncmp(argv[i],"OutF=",5)==0) {
      sscanf(argv[i],"OutF=%s",OutFileName);
    }    
    else if (strncmp(argv[i],"OutTxtF=",8)==0) {
      sscanf(argv[i],"OutTxtF=%s",OutTextFileName);
    } 
    else if (strncmp(argv[i],"KernS=",6)==0) {
      sscanf(argv[i],"KernS=%ld",&KernSize);
    }  
    else if (strncmp(argv[i],"AntDirn=",8)==0) {
      sscanf(argv[i],"AntDirn=%c",&AntennaDirn);
    } 
    else if (strncmp(argv[i],"TerAlt=",7)==0) {
      sscanf(argv[i],"TerAlt=%lf",&TerrainAlt);
    } 
    else if (strncmp(argv[i],"NearRngBin=",11)==0) {
      sscanf(argv[i],"NearRngBin=%ld",&NearProcRngBin);
    }
    else if (strncmp(argv[i],"FarRngBin=",10)==0) {
      sscanf(argv[i],"FarRngBin=%ld",&FarProcRngBin);
    }
    else {
      fprintf(msg,"ERROR - parameter %s unknown!\n",argv[i]);
      exit(1);
    }      
  }  /* end for i loop */

   
  /* Check that required params have been specified */
  errors = 0;
  if (strcmp(MotFileName,"")==0) {
    fprintf(msg,"ERROR - MotF not specified!\n"); errors++;
  }
  if (Delay0 == DEFAULT_VAL_F) {
    fprintf(msg,"ERROR - Delay0 not specified!\n"); errors++;
  }
  if (RefRngBin == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - RefBin not specified!\n"); errors++;
  }
  if (ProcStartLP == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - ProcStartLP not specified!\n"); errors++;
  }
  if (ProcEndLP == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - ProcEndLP not specified!\n"); errors++;
  } 
  if (A2DFreq == DEFAULT_VAL_F) {
    fprintf(msg,"ERROR - A2DFreq not specified!\n"); errors++;
  } 
  if (A2DFreq == DEFAULT_VAL_F) {
    fprintf(msg,"ERROR - PRF not specified!\n"); errors++;
  } 
  if (DataStartLP == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - DataStartLP not specified!\n"); errors++;
  } 
  if (ProcStartLP == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - ProcStartLP not specified!\n"); errors++;
  } 
  if (ProcEndLP == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - ProcEndLP not specified!\n"); errors++;
  } 
  if (LPIncr == DEFAULT_VAL_I) {
    fprintf(msg,"ERROR - LPIncr not specified!\n"); errors++;
  } 
  if (errors != 0) exit(1);

  /* Create out file name and log file name and open files */
  if (strcmp(OutFileName,"")==0) {
    strcpy(OutFileName,MotFileName);
    strcat(OutFileName,DEFAULT_OUT_EXT);
  }  
  if (strcmp(LogFileName,"")==0) {  /* if log file name not defined */
    strcpy(LogFileName,MotFileName);
    strcat(LogFileName,DEFAULT_LOG_EXT);
  }
  if(( LogFile = fopen(LogFileName,"wt")) == NULL) { 
    fprintf(msg,
      "ERROR - Unable to open output log file %s!\n",LogFileName); 
    exit(1); 
  }  
  if(( MotFile = fopen(MotFileName,"r")) == NULL) { 
    fprintf(msg,
      "ERROR - Unable to input motion file %s!\n",MotFileName); 
    exit(1); 
  }  
  if(( OutFile = fopen(OutFileName,"wb")) == NULL) { 
    fprintf(msg,
      "ERROR - Unable to open output file %s!\n",OutFileName); 
    exit(1); 
  }
  if (strcmp(OutTextFileName,"")!=0) {
    if(( OutTextFile = fopen(OutTextFileName,"w")) == NULL) { 
      fprintf(msg,
        "ERROR - Unable to open output text file %s!\n",OutTextFileName); 
      exit(1); 
    }
    OutTextFileFlg = 'Y';
  } /* end if strcmp */
  
  /* Reset any values which have been set to type defaults. */
  if (KernSize == DEFAULT_VAL_I) {
    KernSize = 1;
  }
  else if (KernSize % 2 == 0) {  /* ensure odd */
    KernSize++;
  }
  if (TerrainAlt == DEFAULT_VAL_F) {
    TerrainAlt = 0.0;
  }  
  
  /* Assign messages to the log file */
  msg = LogFile;
  fprintf(msg,"Prog: MOCOMP (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1999\n");

  /* MESSAGES (more later) */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"Input motion file name        : %s\n",MotFileName);            
  fprintf(msg,"Output rng shift binary file  : %s\n",OutFileName);
  if (OutTextFileFlg == 'Y') {
    fprintf(msg,"Output rng shift text file    : %s\n",OutTextFileName);
  }
  else {
    fprintf(msg,"Output rng shift text file    : [not written]\n");
  }
  fprintf(msg,"Smoothing kernel size         : %ld\n",KernSize);
  fprintf(msg,"Antenna direction             : %c\n",AntennaDirn);
  fprintf(msg,"Delay0                        : %.5e sec\n",Delay0);
  fprintf(msg,"Reference range bin           : %ld\n",RefRngBin);
  fprintf(msg,"A/D freq                      : %.5e Hz\n",A2DFreq);
  fprintf(msg,"PRF                           : %.5e Hz\n",PRF);
  fprintf(msg,"Terrain altitude              : %f m\n",TerrainAlt);
  fprintf(msg,"LBR PRI ID increment          : %ld\n",LPIncr);
  fprintf(msg,"LBR PRI for start of raw data : %ld\n",DataStartLP);   

  /********* START PROCESSING ********/

  /* First need to find the start and end record to use for the processing.
   * This will be used to set up the reference track and also to calc the
   * average ground speed. The start record will be first record >=
   * ProcStartLP and the end record will be first record >= ProcEndLP if 
   * enough records exist, else the last available record. */

  ProcStartFound=0;
  ProcEndFound=0;
  ProcNumMotRecords=0;
  ProcAveAlt = 0.0;

  while ( !ReadMotRec(MotFile,&MotDataNext) && !ProcEndFound ) {
    /* Check for start record to process */
    if ( !ProcStartFound && MotDataNext.LP>=ProcStartLP ) {
      ProcMotDataStart = MotDataNext;
      ProcStartFound = 1;
    }
    /* Keep some stats of the processed section */
    if (ProcStartFound && !ProcEndFound) {
      ProcNumMotRecords++;
      ProcAveAlt += MotDataNext.Alt;
    }
    /* Check for end record to process.  */
    if ( !ProcEndFound && MotDataNext.LP>=ProcEndLP ) {
      ProcEndFound = 1;
    }
    /* Reassigns end processed record until while exits */
    ProcMotDataEnd = MotDataNext;
  }  /* end while */

  if (!ProcStartFound) {
    fprintf(msg,"ERROR - no motion records found within proc range!\n");
    exit(1);
  }
  if (ProcNumMotRecords==1) {
    fprintf(msg,"ERROR - only one valid motion record found!\n");
    exit(1);
  }

  /* Calc processed average altitude and average ground speed */
  ProcAveAlt = ProcAveAlt / ProcNumMotRecords;
  ConvToCartesECEF(ProcMotDataEnd.Lat * TO_RAD,
               ProcMotDataEnd.Long * TO_RAD,
               ProcMotDataStart.Alt, /* use same height as start for calc */
               &p2);
  ConvToCartesECEF(ProcMotDataStart.Lat * TO_RAD,
               ProcMotDataStart.Long * TO_RAD,
               ProcMotDataStart.Alt,
               &p1);                            
  ProcMotDist = CalcCartesDist(p1,p2);
  ProcMotTime = ((ProcMotDataEnd.LP - ProcMotDataStart.LP)/LPIncr)/PRF;
 
  if (ProcMotTime != 0.0) {
    ProcAveGrndSpeed = ProcMotDist / ProcMotTime;
  }
  else {
    fprintf(msg,"ERROR - processed motion duration zero!\n");
    exit(1);
  }  

  StartG2PRI = (Int4B)((ProcStartLP-DataStartLP)/LPIncr);
  EndG2PRI = (Int4B)((ProcEndLP-DataStartLP)/LPIncr);

  /* more messages */
  fprintf(msg,"\nPARAMETERS FOR MOTION RECORDS PROCESSED:\n"); 
  fprintf(msg,"Proc G2 / LBR start PRI                : %ld / %ld\n",
          StartG2PRI,ProcStartLP);
  fprintf(msg,"Proc G2 / LBR end PRI                  : %ld / %ld\n",
          EndG2PRI,ProcEndLP);
  fprintf(msg,"Proc start mot record LBR PRI          : %ld\n",
          ProcMotDataStart.LP);
  fprintf(msg,"Proc end mot record LBR PRI            : %ld\n",
          ProcMotDataEnd.LP);
  fprintf(msg,"Proc platform start DGPS lat / long    : %.10f / %.10f deg\n",
          ProcMotDataStart.Lat,ProcMotDataStart.Long);
  fprintf(msg,"Proc platform end DGPS lat / long      : %.10f / %.10f deg\n",
          ProcMotDataEnd.Lat,ProcMotDataEnd.Long);
  fprintf(msg,"Proc start DGPS UTC                    : %.1f sec\n",
          ProcMotDataStart.UTCsecs);
  fprintf(msg,"Proc end DGPS UTC                      : %.1f sec\n",
          ProcMotDataEnd.UTCsecs);
  fprintf(msg,"Proc start DGPS altitude               : %.2f m\n",
          ProcMotDataStart.Alt);                  
  fprintf(msg,"Proc end DGPS altitude                 : %.2f m\n",
          ProcMotDataEnd.Alt);
  fprintf(msg,"Proc num motion records       : %ld\n",ProcNumMotRecords);
  fprintf(msg,"Proc average ground speed     : %f m/s\n",ProcAveGrndSpeed);
  fprintf(msg,"Proc horizontal motion dist.  : %f m\n",ProcMotDist);  
  fprintf(msg,"Proc motion duration          : %f secs\n",ProcMotTime);   
  fprintf(msg,"Proc average altitude (AMSL)  : %f m\n", ProcAveAlt);


  /***************************
    Perform motion comp calcs 
  ****************************/

  /* Calculate gradients for equation of reference flight path using the
   * constant start height for reference track. */
  ConvToCartesECEF(ProcMotDataEnd.Lat * TO_RAD,
               ProcMotDataEnd.Long * TO_RAD,
               ProcMotDataStart.Alt, /* use same height as start for ref calc */
               &RefEndPt);
  ConvToCartesECEF(ProcMotDataStart.Lat * TO_RAD,
               ProcMotDataStart.Long * TO_RAD,
               ProcMotDataStart.Alt,
               &RefStartPt);         
  RefM.x = (RefEndPt.x - RefStartPt.x) / ProcMotTime;   /* calc gradients */
  RefM.y = (RefEndPt.y - RefStartPt.y) / ProcMotTime;
  RefM.z = (RefEndPt.z - RefStartPt.z) / ProcMotTime;    

  /* Calculate start position of ground reference track */
  NomFlight.x = RefEndPt.x - RefStartPt.x;  /* vector along ref flight path */
  NomFlight.y = RefEndPt.y - RefStartPt.y;
  NomFlight.z = RefEndPt.z - RefStartPt.z;
  NomSlantRng = (C/2.0)*(Delay0 + (double)RefRngBin/A2DFreq);
  NomLookAngleRad = CalcLookAngleRad(ProcMotDataStart.Alt-TerrainAlt,NomSlantRng);
  CalcGrndRefPt(AntennaDirn, 
                   NomSlantRng,
                   NomLookAngleRad, 
                   RefStartPt,
                   NomFlight,
                   &GrndRefStartPt);

  ConvToGeodetic(GrndRefStartPt,&GrndRefStartGPosn); /* used for log only */
  

  /****** Calc the geodetic coords of the image corners *****/ 
  /* This only used for logs. Note assumes zero-Doppler coord processing. */

  /* First we need to find the position on the extrapolated reference track
   * corresponding to the start and end of the image. Extrapolation may be
   * required as the first and last motion records processed may not 
   * necessarily correspond exactly to the start and end of the image to
   * be processed.*/
  ImRefStartPt.x = RefStartPt.x - 
                 RefM.x*(ProcMotDataStart.LP-ProcStartLP)/(LPIncr*PRF);
  ImRefStartPt.y = RefStartPt.y - 
                 RefM.y*(ProcMotDataStart.LP-ProcStartLP)/(LPIncr*PRF);   
  ImRefStartPt.z = RefStartPt.z - 
                 RefM.z*(ProcMotDataStart.LP-ProcStartLP)/(LPIncr*PRF);
  ImRefEndPt.x = RefEndPt.x - 
               RefM.x*(ProcMotDataEnd.LP-ProcEndLP)/(LPIncr*PRF);
  ImRefEndPt.y = RefEndPt.y - 
               RefM.y*(ProcMotDataEnd.LP-ProcEndLP)/(LPIncr*PRF);
  ImRefEndPt.z = RefEndPt.z - 
               RefM.z*(ProcMotDataEnd.LP-ProcEndLP)/(LPIncr*PRF);

  /* Now to find the four corner pts in cartesian ECEF coords */ 
  NearProcSlantRng = (C/2.0)*(Delay0 + (double)NearProcRngBin/A2DFreq);
  FarProcSlantRng = (C/2.0)*(Delay0 + (double)FarProcRngBin/A2DFreq); 
  NearLookAngleRad = CalcLookAngleRad(ProcMotDataStart.Alt-TerrainAlt,
                                      NearProcSlantRng);
  FarLookAngleRad = CalcLookAngleRad(ProcMotDataStart.Alt-TerrainAlt,
                                     FarProcSlantRng);

  CalcGrndRefPt(AntennaDirn,NearProcSlantRng,NearLookAngleRad, 
                ImRefStartPt,NomFlight,&ImNearStartPt);
  CalcGrndRefPt(AntennaDirn,NearProcSlantRng,NearLookAngleRad, 
                ImRefEndPt,NomFlight,&ImNearEndPt);
  CalcGrndRefPt(AntennaDirn,FarProcSlantRng,FarLookAngleRad, 
                ImRefStartPt,NomFlight,&ImFarStartPt);
  CalcGrndRefPt(AntennaDirn,FarProcSlantRng,FarLookAngleRad, 
               ImRefEndPt,NomFlight,&ImFarEndPt);

  /* Finally, convert the image corner positions and the image reference track
   * to geodetic coords */
  ConvToGeodetic(ImNearStartPt,&ImNearStartGPosn);   
  ConvToGeodetic(ImNearEndPt,&ImNearEndGPosn);  
  ConvToGeodetic(ImFarStartPt,&ImFarStartGPosn);  
  ConvToGeodetic(ImFarEndPt,&ImFarEndGPosn);                     
  ConvToGeodetic(ImRefStartPt,&ImRefStartGPosn);
  ConvToGeodetic(ImRefEndPt,&ImRefEndGPosn);   

  /*******CALCULATE THE RANGE SHIFT FOR PROCESSED MOTION RECORDS********/
  /* In the first step, the range shift is calculated only at the LBR PRIs
   * corresponding to the motion records to be processed. After this step,
   * the range shifts have been stored in the ProcMotShift array and the 
   * corresponding LBR PRIs in the ProcMotShiftLP array. */

  /* Allocate mem for the rng shifts */
  ProcMotShift = (double *)malloc(sizeof(double)*ProcNumMotRecords);
  ProcSmoothShift = (double *)malloc(sizeof(double)*ProcNumMotRecords);
  ProcMotShiftLP = (Int4B *)malloc(sizeof(Int4B)*ProcNumMotRecords);
  if (ProcMotShift == NULL || ProcSmoothShift==NULL 
      || ProcMotShiftLP == NULL) {
    fprintf(msg,"ERROR - in shift array allocation!\n");
    exit(1);
  } 
 
  /* Reset file pointer to start of motion records and search for start proc
   * motion record. After this step, MotDataNext contains the motion record
   * which immediately follows the first mot record processed. */
  fseek(MotFile,0,SEEK_SET);
  MotDataNext.LP = -999;
  while (MotDataNext.LP <= ProcMotDataStart.LP) {
    if (ReadMotRec(MotFile,&MotDataNext)) {
      fprintf(msg,"ERROR - motion record not found!\n");
      exit(1);      
    }
  }  /* end while */

  /* Assign first element of range shift arrays. Shift for first processed 
   * record set to zero. */
  ProcMotShiftLP[0] = ProcMotDataStart.LP;
  ProcMotShift[0] = 0.0;

  /* Perform motion calculation for each processed motion record */ 
  for (MotRec=0; MotRec< ProcNumMotRecords-1; MotRec++) {

    /* Calculate the GrndRefPt at PRINext */
    PRIChange = (Int4B)( (MotDataNext.LP - ProcMotDataStart.LP)/LPIncr ); 
    GrndRefPt.x = GrndRefStartPt.x + RefM.x*PRIChange/PRF;
    GrndRefPt.y = GrndRefStartPt.y + RefM.y*PRIChange/PRF;
    GrndRefPt.z = GrndRefStartPt.z + RefM.z*PRIChange/PRF;

    /* Calc distance from aircraft to GrndRefPt at PRINext*/
    ConvToCartesECEF(MotDataNext.Lat * TO_RAD,
                 MotDataNext.Long * TO_RAD,
                 MotDataNext.Alt,		/* if using actual reported DGPS alt*/
//     			 ProcMotDataStart.Alt,  /* if using constant alt */
                 &AircraftPosn);  
                       
    RngShiftNext = NomSlantRng - CalcCartesDist(AircraftPosn,GrndRefPt);
    ProcMotShiftLP[MotRec+1] = MotDataNext.LP;
    ProcMotShift[MotRec+1] = RngShiftNext;
    
    /* Read next PRI */
    if (ReadMotRec(MotFile,&MotDataNext)) break;   
     
  } /* end for MotRec loop */

  ConvToGeodetic(GrndRefPt,&GrndRefEndGPosn); /* used for log only */

  /* Smooth the range shift data */
  if (MovingAveSmooth(ProcMotShift,
                       ProcSmoothShift,
                       ProcNumMotRecords,
                       KernSize) ) {
    fprintf(msg,"ERROR - in moving ave function call!\n");
    exit(1);
  }

  /******Write out range shifts to file*******/
  /* Linearly interpolate between available rng shifts. */

  /* Initialise output (G2) PRI ID */
  OutG2PRI = (Int4B)( (ProcStartLP - DataStartLP)/LPIncr );

  /* Write zero range shift for LBR PRIs <= ProcMotDataStart.LP */
  RngShift = 0.0;
  for (i=0; i<=(Int4B)( (ProcMotDataStart.LP-ProcStartLP)/
                        LPIncr ); i++) {
    if (OutTextFileFlg == 'Y')
      fprintf(OutTextFile,"%f\t%f\n",(double)OutG2PRI/PRF,RngShift); 
    fwrite(&OutG2PRI,sizeof(Int4B),1,OutFile);
    OutG2PRI++;
    fwrite(&RngShift,sizeof(float),1,OutFile);
  }

  /* Write out rng shifts over region of valid mot records. Note that
   * the zero shift corresponding to ProcMotDataStart has already been
   * written. */
  for (MotRec=0; MotRec< ProcNumMotRecords-1; MotRec++) {

    MotDataCur.LP = ProcMotShiftLP[MotRec];
    RngShiftCur = ProcSmoothShift[MotRec];
    MotDataNext.LP = ProcMotShiftLP[MotRec+1];
    RngShiftNext = ProcSmoothShift[MotRec+1];

    PRIChange = (Int4B)((MotDataNext.LP - MotDataCur.LP)/LPIncr); 
    RngShiftIncr = (RngShiftNext - RngShiftCur)/PRIChange;

    /* Write out for each G2 PRI ID */
    for (i=0; i < PRIChange; i++) {

      /* Calc range shift */
      RngShift = (float)(RngShiftCur + (i+1)*RngShiftIncr);
      if (fabs(RngShift) > fabs(MaxRngShift)) {
         MaxRngShift = RngShift;
         MaxRngShiftIndx = OutG2PRI;          
      }
  
      /* Write range shift to output file */
     if (OutTextFileFlg == 'Y')
        fprintf(OutTextFile,"%f\t%f\n",(double)OutG2PRI/PRF,RngShift); 
      fwrite(&OutG2PRI,sizeof(Int4B),1,OutFile);
      OutG2PRI++;
      fwrite(&RngShift,sizeof(float),1,OutFile);

    } /* end for i loop */

 } /* end for MotRec */

  /* Write extra range shifts at end, if necessary, repeat last value */
  for (i=0; i<(Int4B)((ProcEndLP-ProcMotDataEnd.LP)/LPIncr); i++) {
    if (OutTextFileFlg == 'Y')
      fprintf(OutTextFile,"%f\t%f\n",(double)OutG2PRI/PRF,RngShift);  
    fwrite(&OutG2PRI,sizeof(Int4B),1,OutFile);
    OutG2PRI++;
    fwrite(&RngShift,sizeof(float),1,OutFile);

  } 

  fprintf(msg,"Reference range bin           : %ld\n",RefRngBin);
  fprintf(msg,"NomSlantRng                   : %f m\n",NomSlantRng);
  fprintf(msg,"NomLookAngle                  : %f rad (%f deg)\n",
          NomLookAngleRad, NomLookAngleRad*TO_DEG);
  fprintf(msg,"Max range shift               : %f m\n",MaxRngShift);
  fprintf(msg,"Max range shift G2 PRI        : %ld\n",MaxRngShiftIndx);
  fprintf(msg,"Ground ref track start lat    : %.10f deg\n",
          GrndRefStartGPosn.latrad*TO_DEG);
  fprintf(msg,"Ground ref track start long   : %.10f deg\n",
          GrndRefStartGPosn.longrad*TO_DEG);          
  fprintf(msg,"Ground ref start alt          : %.2f m (AMSL)\n",
          GrndRefStartGPosn.alt);
  fprintf(msg,"Ground ref track end lat      : %.10f deg\n",
          GrndRefEndGPosn.latrad*TO_DEG);
  fprintf(msg,"Ground ref track end long     : %.10f deg\n",
          GrndRefEndGPosn.longrad*TO_DEG);          
  fprintf(msg,"Ground ref end alt            : %.2f m (AMSL)\n",
          GrndRefEndGPosn.alt);

  fprintf(msg,"\nIMAGE GEOCODING:\n");
  fprintf(msg,"Image range bins processed          : %ld - %ld\n",
          NearProcRngBin,FarProcRngBin); 
  fprintf(msg,"Image near swath start latitude     : %.10f deg\n",
          ImNearStartGPosn.latrad*TO_DEG);
  fprintf(msg,"Image near swath start longitude    : %.10f deg\n",
          ImNearStartGPosn.longrad*TO_DEG);
  fprintf(msg,"Image near swath start altitude     : %.2f m (AMSL)\n",
          ImNearStartGPosn.alt);          
  fprintf(msg,"Image far swath start latitude      : %.10f deg\n",
          ImFarStartGPosn.latrad*TO_DEG);
  fprintf(msg,"Image far swath start longitude     : %.10f deg\n",
          ImFarStartGPosn.longrad*TO_DEG);
  fprintf(msg,"Image far swath start altitude      : %.2f m (AMSL)\n",
          ImFarStartGPosn.alt);     
  fprintf(msg,"Image near swath end latitude       : %.10f deg\n",
          ImNearEndGPosn.latrad*TO_DEG);
  fprintf(msg,"Image near swath end longitude      : %.10f deg\n",
          ImNearEndGPosn.longrad*TO_DEG);
  fprintf(msg,"Image near swath end altitude       : %.2f m (AMSL)\n",
          ImNearEndGPosn.alt);          
  fprintf(msg,"Image far swath end latitude        : %.10f deg\n",
          ImFarEndGPosn.latrad*TO_DEG);
  fprintf(msg,"Image far swath end longitude       : %.10f deg\n",
          ImFarEndGPosn.longrad*TO_DEG);
  fprintf(msg,"Image far swath end altitude        : %.2f m (AMSL)\n",
          ImFarEndGPosn.alt);   
  fprintf(msg,"Image platform ref track start lat  : %.10f deg\n",
          ImRefStartGPosn.latrad*TO_DEG);
  fprintf(msg,"Image platform ref track start long : %.10f deg\n",
          ImRefStartGPosn.longrad*TO_DEG);
  fprintf(msg,"Image platform ref track start alt  : %.2f m (AMSL)\n",
          ImRefStartGPosn.alt);
  fprintf(msg,"Image platform ref track end lat    : %.10f deg\n",
          ImRefEndGPosn.latrad*TO_DEG);
  fprintf(msg,"Image platform ref track end long   : %.10f deg\n",
          ImRefEndGPosn.longrad*TO_DEG);
  fprintf(msg,"Image platform ref track end alt    : %.2f m (AMSL)\n",
          ImRefEndGPosn.alt);

  fprintf(stdout,"DONE!\n");
  
  /* Clean up */
  free(ProcMotShift);
  free(ProcMotShiftLP);

  if (OutTextFileFlg == 'Y')
    fclose(OutTextFile);

  fclose(LogFile);
  fclose(MotFile);
  fclose(OutFile);

  return(0);  /* on success */

}  /***** end main prog or function ******/



/*===============FUNCTIONS================*/

/*============ReadMotRec============*/
/* Function to parse DGPS file from current posn to read 
   next motion record (incl LBR PRI ID). To use only the unpacked IMU
   data, uncomment the two lines marked and remove the first text
   line from the unpacked IMU file (also rename to *_IMU+DGPS) */

Int2B ReadMotRec(FILE *MotFile,struct MotStruct *MotData)
{
  Int4B FilePosn; 
  // if only using IMU data, uncomment next line 
  // double t1,t2,t3,t4,t5,t6;
  FilePosn = ftell(MotFile);
  if (fscanf(MotFile,"%ld %lf %lf %lf %lf",&MotData->LP,&MotData->UTCsecs,
                   &MotData->Lat,&MotData->Long,&MotData->Alt) 
     != 5) { return(1); }
  // if only using IMU data, uncomment next line   
  // fscanf(MotFile,"%lf %lf %lf %lf %lf %lf\n",&t1,&t2,&t3,&t4,&t5,&t6);
  MotData->FilePosn = FilePosn; /* assign after fscanf in case exits early */
  return(0);
}                       


/*============ConvToCartesECEF================*/
/* Function to convert from lat,long and ellipsoidal height to
   3-D Cartesian coords (Earth centred, Earth fixed) */
Int2B ConvToCartesECEF(double LatRad, double LongRad, double h,
                   struct cartes *Pt)
{
  double N;  

  N = WGS84_a / sqrt(1 - WGS84_e_sq * sin(LatRad)*sin(LatRad));
  Pt->x = (N+h)*cos(LatRad)*cos(LongRad);
  Pt->y = (N+h)*cos(LatRad)*sin(LongRad);
  Pt->z = (N*(1-WGS84_e_sq) + h)*sin(LatRad);

  return(0);
}                      

/*============ConvToGeodetic================*/
/* Function to convert from cartesian ECEF (earth centred, 
   earth fixed coords) to geodetic lat,long and ellipsoidal height */
Int2B ConvToGeodetic(struct cartes Pt,struct geodetic *GPosn)   
{
  double N,p,x,y,z,theta,phi;  

  x = Pt.x;
  y = Pt.y;
  z = Pt.z; 
  p = sqrt(x*x + y*y);
  theta = atan(z*WGS84_a/(p*WGS84_b));
  phi = atan((z+WGS84_ep_sq*WGS84_b*sin(theta)*sin(theta)*sin(theta))/
             (p-WGS84_e_sq*WGS84_a*cos(theta)*cos(theta)*cos(theta)));
  N = WGS84_a / sqrt(1 - WGS84_e_sq * sin(phi)*sin(phi));
  GPosn->latrad = phi;
  GPosn->longrad = atan2(y,x);
  GPosn->alt = p/cos(phi) - N;  

  return(0);
}     

/*===========CalcLookAngleRad===============*/
/* Function to calc look angle in radians to data point */
double CalcLookAngleRad(double heightAGL, double slantrng)
{
  if (heightAGL > slantrng)
    return(0.0);
  else
    return(acos(heightAGL/slantrng));
}

/*===========CalcCartesDist===============*/
/* Function to calculate the straight line distance between two
   points specified in cartesian coords */

double CalcCartesDist(struct cartes p1,struct cartes p2)
{
  double dist;
  dist = (p2.x - p1.x)*(p2.x - p1.x) +
         (p2.y - p1.y)*(p2.y - p1.y) +
         (p2.z - p1.z)*(p2.z - p1.z);
  dist = sqrt(dist);
  return (dist);

}
   
/*============CalcGrndRefPt============*/
/* Function to calc a position in the ground swath in ECEF cartesian coords.
 * The pt calculated is at a position on the ground orthogonal to the
 * nominal flight direction at the FlightRef pt and at specified slant
 * range and look angle. This function is used to find the start of the
 * ground reference track and also the position of the image corners for
 * geocoding purposes. The Maple-calculated solution is used here. */
Int2B CalcGrndRefPt(char Dirn, /* antenna dirn (L/R) */
                       double NomSlantRng,
                       double LookAngleRad, 
                       struct cartes FlightRef, /* coords of ref start pt */
                       struct cartes NomFlight, /* vector along ref flight path */
                       struct cartes *GrndRefPt) /* coords of grnd ref start */
{

  double a,b,c,    /* coeffs of quadratic to solve */
         h,k,m,    /* FlightRef (x,y,z) */
         q,r,s,    /* NomFlight (x,y,z) */
         t,        /* NomSlantRng */
         u,        /* R dot FlightRef */
         p1,p2,    /* possible z components for R */
         CrossCalc1,
         CrossCalc2;
  struct cartes  GPt,
                 R,R1,R2; /* vector from FlightRef posn to GrndRefPt 
                            and two possibles */
  struct geodetic FlightRefGPosn,GrndRefGPosn;

  /* handle special case and then return where LookAngle zero
   * (i.e. want pt directly below aircraft flight posn. - sample 
   * before nadir return) */
  if (LookAngleRad == 0.0) {
    ConvToGeodetic(FlightRef,&FlightRefGPosn);
    GrndRefGPosn.latrad = FlightRefGPosn.latrad;
    GrndRefGPosn.longrad = FlightRefGPosn.longrad;
    GrndRefGPosn.alt = FlightRefGPosn.alt - NomSlantRng;
    ConvToCartesECEF(GrndRefGPosn.latrad,GrndRefGPosn.longrad,
                     GrndRefGPosn.alt,&GPt);
    GrndRefPt->x = GPt.x;
    GrndRefPt->y = GPt.y;
    GrndRefPt->z = GPt.z;
    return(0);                    
  }

                 
  h = FlightRef.x;
  k = FlightRef.y;
  m = FlightRef.z;
  q = NomFlight.x;
  r = NomFlight.y;
  s = NomFlight.z;
  t = NomSlantRng;

  /* Calc u (R.FlightRef) */
  u  = sqrt(FlightRef.x*FlightRef.x + FlightRef.y*FlightRef.y +
            FlightRef.z*FlightRef.z);
  u *= NomSlantRng*cos(PI - LookAngleRad);
  
  /* calc coeffs of quadratic used to find z component of R */
  a = s*s*k*k + m*m*q*q + q*q*k*k - 2.0*r*s*m*k + r*r*m*m +
      h*h*s*s - 2.0*h*s*m*q + h*h*r*r - 2.0*h*r*q*k;
  b = 2.0*r*s*u*k + 2.0*h*s*u*q - 2.0*u*m*r*r - 2.0*u*q*q*m;
  c = r*r*u*u + u*u*q*q - t*t*h*h*r*r + 2.0*t*t*q*h*r*k -
      t*t*q*q*k*k;    

  /* Calc roots of quadratic and the two possible R's */  
  QuadRoots(a,b,c,&p1,&p2);
  CalcR(h,k,m,q,r,s,p1,u,&R1);
  CalcR(h,k,m,q,r,s,p2,u,&R2);

  /* Check which R gives correct direction (L or R) */
  CrossDot(FlightRef,R1,NomFlight,&CrossCalc1);
  CrossDot(FlightRef,R2,NomFlight,&CrossCalc2);
  
  if ( (Dirn == 'R' && CrossCalc1 >= 0.0) || 
       (Dirn == 'L' && CrossCalc1  < 0.0) )  {
    R.x = R1.x;
    R.y = R1.y;
    R.z = R1.z;
  }
  else if ( (Dirn == 'R' && CrossCalc2 >= 0.0) ||
       (Dirn == 'L' && CrossCalc2 < 0.0) ) {
    R.x = R2.x;
    R.y = R2.y;
    R.z = R2.z;
  }

  /* Calc GrndRefPt */
  GrndRefPt->x = FlightRef.x + R.x;
  GrndRefPt->y = FlightRef.y + R.y;
  GrndRefPt->z = FlightRef.z + R.z;

  return(0);

} /* end function CalcGrndRefPt */            


/*===========QuadRoots===============*/

/* Calcs the roots of a quadratic */
Int2B QuadRoots(double a, double b, double c, 
                double *root1, double *root2) 
{
  double tmp = b*b-4*a*c;
  if (a != 0.0) {
    if (tmp >= 0.0) {
      *root1 = (-b+sqrt(tmp))/(2*a);
      *root2 = (-b-sqrt(tmp))/(2*a);
    }
    else {
      printf("ERROR - sqrt of neg in QuadRoots routine!\n");
      exit(1);
    }
  }
  else {
    printf("ERROR - denom zero in QuadRoots routine!\n");
    exit(1);
  }

  return(0);
}  /* end function QuadRoots */                


/*==========CalcR  ==================*/

/* Part of Maple-calculated solution for the SASAR motion 
   comp */
Int2B CalcR(double h, double k, double m, 
            double q, double r, double s,
            double p, double u, 
            struct cartes *R)            
{
  double tmp;
  tmp = h*r-q*k;
  if (tmp != 0.0) { 
    R->x = - (r*p*m - r*u - p*s*k)/tmp;
    R->y = - (h*p*s - p*m*q + u*q)/tmp;
    R->z = p;
  }
  else {
    printf("ERROR - denom zero in CalcR routine!\n");
    exit(1);
  }

  return(0);
}  /* end function CalcR */

/*==========CrossDot=================*/

/* Calcs (FlightRef x R) . NomFlight */
Int2B CrossDot(struct cartes FlightRef,
               struct cartes R,
               struct cartes NomFlight,
               double *CrossCalc)
{

  double a,b,c,
         d,e,f;

  a = FlightRef.x;
  b = FlightRef.y;
  c = FlightRef.z;
  d = R.x;
  e = R.y;
  f = R.z;

  *CrossCalc  = (b*f-c*e)*NomFlight.x;
  *CrossCalc += (c*d-a*f)*NomFlight.y;
  *CrossCalc += (a*e-b*d)*NomFlight.z;

  return(0);

}  /* end func CrossDot */
               

/*==========MovingAveSmooth===========*/

/* Performs a moving average smoothing operation on real data.
The kernel size should be an odd integer. It is assumed that the
start and end values of the array are repeated beyond the array limits. */ 
Int2B MovingAveSmooth(double In[],double Out[],
                      Int4B ArrSize, Int4B KernSize)
{
  Int4B i,indx,KernSizeD2;
  double *InEx;

  if ( KernSize % 2 == 0 ) KernSize++;
  KernSizeD2 = (Int4B)(KernSize/2);
  
  InEx = (double *)malloc(sizeof(double)*(ArrSize+KernSize-1));
  if (InEx==NULL) return(1);

  /* set up extended input array */
  indx = 0;
  for (i=0; i< KernSizeD2; i++) {
    InEx[indx++] = In[0];
  }
  for (i=0; i< ArrSize; i++) {
    InEx[indx++] = In[i];
  }
  for (i=0; i< KernSizeD2; i++) {
    InEx[indx++] = In[ArrSize-1];
  }

  /* perform moving average calc */
  for (indx=0; indx<ArrSize;indx++) {
    Out[indx] = 0.0;
    for (i=0; i<KernSize; i++) {
      Out[indx] += InEx[indx+i];   
    }
    Out[indx] /= KernSize;
  }

  free(InEx);
  return(0);  /* on success */

} /* end function MovingAveSmooth */

  