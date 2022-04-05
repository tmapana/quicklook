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

VERSION/AUTHOR/DATE :
COMMENTS:

=========================================*/

#include "g2func.h"
#include "g2moc.h"
#include "g2parse.h"

#define PROG_VERSION "1999-02-04"

#define DEFAULT_OUT_EXT "_moc"
#define DEFAULT_LOG_EXT "_mlog"
#define WGS84_a 6378137.0
#define WGS84_b 6356752.314
#define WGS84_e_sq 6.694380066e-03
#define TO_RAD 0.01745329252
#define TO_DEG 57.29577951

#define WRITE_LBR_MOT_SHIFTS 0

/* Structure definitions for this file */
struct cartes {
  double x;
  double y;
  double z;
 };

/* Function prototypes for functions in this file */
Int2B ConvToCartes(double LatRad, double LongRad, double h,
                   struct cartes *Pt); 
double CalcCartesDist(struct cartes p1,struct cartes p2);  
double CalcTimeSecs(struct IMUStruct Data);
Int2B CalcStartGrndRef(char Dirn, 
                       double NomSlantRng,
                       double LookAngleRad, 
                       struct cartes StartRef, 
                       struct cartes NomFlight,
                       struct cartes *StartGrndRef);
Int2B QuadRoots(double a, double b, double c, 
                double *root1, double *root2);
Int2B CalcR(double h, double k, double m, 
            double q, double r, double s,
            double p, double u, 
            struct cartes *R);   
Int2B CrossDot(struct cartes StartRef,
               struct cartes R,
               struct cartes NomFlight,
               double *CrossCalc);
Int2B ReadLBRMotionRec(FILE *LBRFile,
                       char Radar,
                       Int4B *PRI,
                       struct IMUStruct *IMUData);
Int2B MovingAveSmooth(double In[], double Out[], 
                      Int4B ArrSize, Int4B KernSize);

/*===========MAIN PROG/FUNCTION===========*/

#if FUNC_MOC
Int2B G2Moc (struct MocCmdStruct Cmd) 
#else
Int2B main (int argc,char *argv[])
#endif
{
  FILE *msg=stdout, *LBRFile=NULL, *OutFile=NULL,*LogFile=NULL,*OutTextFile=NULL;
  char ProgVersion[STRING_SPACE]="",LBRFileName[STRING_SPACE],
       OutFileName[STRING_SPACE]="",OutTextFileName[STRING_SPACE]="",
       LogFileName[STRING_SPACE]="",
       OutTextFileFlg = 'N', /* for range shifts */
       AntennaDirn = 'R', /* default value */
       HBR='3',   /* default value */
       Radar='B'; /* default value */
  Int4B i,
        EndG2PRI=DEFAULT_VAL_I, /* zero corresponds to LBR start PRI, incr 1 */ 
        KernSize=1,  /* moving average kernel size */
        LBRMotEndPRI,
        LBRMotStartPRI,
        LBRNumMotRecords=0,
        MaxRngShiftIndx=DEFAULT_VAL_I,
        MotRec,
        MotRecLength,
        OutPRI,
        ProcEndPRI=DEFAULT_VAL_I,
        ProcMotEndPRI=DEFAULT_VAL_I,
        ProcMotEndPRIFound=0,
        *ProcMotShiftPRI, /* LBR PRI IDs for elements of ProcMotShift array */
        ProcMotStartPRI=DEFAULT_VAL_I,
        ProcMotStartPRIFound=0,
        ProcNumMotRecords=0,
        ProcStartPRI=DEFAULT_VAL_I,
        PRIChange,
        PRICur,
        PRINext,
        RefRngBin=DEFAULT_VAL_I,
        StartG2PRI=DEFAULT_VAL_I; /* zero corresponds to LBR start PRI, incr 1 */ 
  float MaxRngShift=0.0,
        RngShift;      
  double LBRAveGPSHgt,
         LBRDist,
         LBRTime,
         LBRAveGrndSpeed,
         NomLookAngleRad,
         NomSlantRng,
         ProcAveGPSHgt,
         ProcMotDist,
         *ProcMotShift, /* holds the rng shifts for the motion records */
         *ProcSmoothShift, /* smoothed range shift for the mot records */  
         ProcMotTime,
         ProcAveGrndSpeed,
         NomAlt=DEFAULT_VAL_F, /* Nominal flight altitude (include for now 
                          as height info not available at present) */
         RngShiftCur,
         RngShiftIncr,
         RngShiftNext,    
         TerrainAlt=DEFAULT_VAL_F, 
         TimeIncr;    
  struct CntrlStruct ParseCntrl;
  struct DataChannelStruct Chan[1];
  struct IMUStruct LBRMotDataEnd,
                   LBRMotDataStart,
                   ProcMotDataEnd,
                   ProcMotDataStart,
                   IMUDataCur,
                   IMUDataNext;
  struct cartes p1,p2,
                AircraftPosn,
                GrndRefPt,
                NomFlight,
                RefEndPt,
                RefStartPt,
                RefM;
           
#if FUNC_MOC

  /* Assign variables from the structure passed */
  msg = Cmd.msg;
  strcpy(ProgVersion,Cmd.Version);
  strcpy(LBRFileName,Cmd.LBRFileName);
  strcpy(OutFileName,Cmd.OutFileName);
  strcpy(OutTextFileName,Cmd.OutTextFileName);
  strcpy(LogFileName,Cmd.LogFileName);
  RefRngBin    = Cmd.RefRngBin;
  ProcStartPRI = Cmd.ProcStartPRI;
  ProcEndPRI   = Cmd.ProcEndPRI;
  StartG2PRI   = Cmd.StartG2PRI;
  EndG2PRI     = Cmd.EndG2PRI;
  KernSize     = Cmd.KernSize;
  AntennaDirn  = Cmd.AntennaDirn;
  Radar        = Cmd.Radar;
  HBR          = Cmd.HBR;
  TerrainAlt   = Cmd.TerrainAlt;
  NomAlt       = Cmd.NomAlt;
  

  fprintf(msg,"\nMOTION COMP CALCULATION STAGE...\n");
  
  /* Check version IDs match */
  if (strcmp(ProgVersion,PROG_VERSION)!=0) {
    fprintf(msg,
      "WARNING - command version (%s) not same as program (%s)!\n\n",
      ProgVersion,PROG_VERSION);
  }

#else /* called from command line */

  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\nProg: MOCOMP (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1999\n");

  if (argc < 2 || argc > 9) {
    fprintf(msg,
      "Calcs motion comp range shifts for SASAR data.\n");
    fprintf(msg,
      "Also writes out ASCII summary motion log file.\n\n");
    fprintf(msg,"USAGE : mocomp [LBR_file]\n");
    fprintf(msg,"Optional params (use at end of line):\n");
    fprintf(msg,
      "<OutF=t>      - the output rng shift file (else default name)\n");
    fprintf(msg,
      "<OutTxtF=t>   - the output rng shift text file (else not written)\n");
    fprintf(msg,
      "<LogF=t>      - the LBR motion log file (else default name)\n");
    fprintf(msg,
      "<RefBin=n>    - sets the range bin to which mocomp is calculated\n");
    fprintf(msg,
      "<StartP=n>    - start LBR PRI (disk sector) (default 1st PRI)\n");
    fprintf(msg,
      "<EndP=n>      - end LBR PRI (disk sector) (default last PRI)\n");
    fprintf(msg,
      "<StartG2P=n>  - start G2 PRI (0 is start data, incr is 1)\n");
    fprintf(msg,
      "<EndG2P=n>    - end G2 PRI (0 is start data, incr is 1)\n");     
    fprintf(msg,
      "<KernS=n>     - rng shift moving ave smooth kernel size (default 1)\n");
    fprintf(msg,
      "<AntDirn=L/R> - antenna direction, R (right) is the default\n");
    fprintf(msg,
      "<Radar=A/B>   - B is the default\n");
    fprintf(msg,
      "<HBR=n>       - in range from 1 - 6, 3 is default\n");
    fprintf(msg,         
      "<TerAlt=m>    - in metres ASL (default zero if unspec in LBR file\n");
    fprintf(msg,
      "<NomAlt=m>    - in metres ASL (include since no heights for now)\n");  
    fprintf(msg,
      "e.g. 'mocomp' lbrfile.001 OutF=mocfile.0 Radar=A NomAlt=3000'\n");
    exit(0);
  } 

  strcpy(LBRFileName,argv[1]);

  /* check for optional parameters */
  for (i=2; i<argc; i++) {
    if (strncmp(argv[i],"OutF=",5)==0) {
      sscanf(argv[i],"OutF=%s",OutFileName);
    }    
    else if (strncmp(argv[i],"OutTxtF=",8)==0) {
      sscanf(argv[i],"OutTxtF=%s",OutTextFileName);
    } 
    else if (strncmp(argv[i],"LogF=",5)==0) {
      sscanf(argv[i],"LogF=%s",LogFileName);
    }    
    else if (strncmp(argv[i],"RefBin=",7)==0) {
      sscanf(argv[i],"RefBin=%ld",&RefRngBin);
    }
    else if (strncmp(argv[i],"StartP=",7)==0) {
      sscanf(argv[i],"StartP=%ld",&ProcStartPRI);
    }
    else if (strncmp(argv[i],"EndP=",5)==0) {
      sscanf(argv[i],"EndP=%ld",&ProcEndPRI);
    }    
    else if (strncmp(argv[i],"StartG2P=",9)==0) {
      sscanf(argv[i],"StartG2P=%ld",&StartG2PRI);
    }
    else if (strncmp(argv[i],"EndG2P=",7)==0) {
      sscanf(argv[i],"EndG2P=%ld",&EndG2PRI);
    }  
    else if (strncmp(argv[i],"KernS=",6)==0) {
      sscanf(argv[i],"KernS=%ld",&KernSize);
    }  
    else if (strncmp(argv[i],"AntDirn=",8)==0) {
      sscanf(argv[i],"AntDirn=%c",&AntennaDirn);
    } 
    else if (strncmp(argv[i],"Radar=",6)==0) {
      sscanf(argv[i],"Radar=%c",&Radar);
    } 
    else if (strncmp(argv[i],"HBR=",4)==0) {
      sscanf(argv[i],"HBR=%c",&HBR);
    } 
    else if (strncmp(argv[i],"TerAlt=",7)==0) {
      sscanf(argv[i],"TerAlt=%lf",&TerrainAlt);
    } 
    else if (strncmp(argv[i],"NomAlt=",7)==0) {
      sscanf(argv[i],"NomAlt=%lf",&NomAlt);
    }
    else {
      fprintf(msg,"ERROR - parameter %s unknown!\n",argv[i]);
      exit(0);
    }      
  }  /* end for i loop */


#endif

  /* Create out file name and log file name and open files */

  if (strcmp(OutFileName,"")==0) {
    strcpy(OutFileName,LBRFileName);
    strcat(OutFileName,DEFAULT_OUT_EXT);
  }  
  if (strcmp(LogFileName,"")==0) {
    strcpy(LogFileName,LBRFileName);
    strcat(LogFileName,DEFAULT_LOG_EXT);
  }  
  
  if (strcmp(OutTextFileName,"")!=0) {
    if(( OutTextFile = fopen(OutTextFileName,"w")) == NULL) { 
      fprintf(msg,"ERROR - Unable to open output text file %s!\n",
              OutTextFileName); 
      exit(0); 
    }
    OutTextFileFlg = 'Y';
  } /* end if strcmp */
  
  if(( LBRFile = fopen(LBRFileName,"rb")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open LBR file %s!\n",
            LBRFileName); 
    exit(0); 
  }
  if(( OutFile = fopen(OutFileName,"wb")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open output file %s!\n",
            OutFileName); 
    exit(0); 
  }
  if(( LogFile = fopen(LogFileName,"w")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open log file %s!\n",
            LogFileName); 
    exit(0); 
  }

  /* Scan through LBR file for relevant params using g2parse.c 
     functions */
  ParseCntrl.NumDataChannels=1;
  ParseCntrl.CurrentChannel=0;
  Chan[0].Radar = Radar;
  Chan[0].HBR = HBR;
  ParseDataChannels(msg,LBRFile,&ParseCntrl,Chan);

  /* Reset any values which have been set to type defaults. Note that G2
     PRI's take priority over the LBR PRI's, if both specified. */
  if (ProcStartPRI == DEFAULT_VAL_I) {
    ProcStartPRI = Chan[0].StartPRI;
  }
  if (ProcEndPRI == DEFAULT_VAL_I) {
    ProcEndPRI = Chan[0].EndPRI;
  } 

  if (StartG2PRI != DEFAULT_VAL_I) {
    ProcStartPRI = Chan[0].StartPRI + StartG2PRI*Chan[0].PRIIDIncr;
  }
  else {  /* set for message purposes */
    StartG2PRI = (Int4B)((ProcStartPRI - Chan[0].StartPRI)/Chan[0].PRIIDIncr);
  }

  if (EndG2PRI != DEFAULT_VAL_I) {
    ProcEndPRI = Chan[0].StartPRI + EndG2PRI*Chan[0].PRIIDIncr;
  }
  else { /* set for message display purposes */ 
    EndG2PRI = (Int4B)((ProcEndPRI-Chan[0].StartPRI)/Chan[0].PRIIDIncr);
  }

  if (RefRngBin == DEFAULT_VAL_I) {
    RefRngBin = (Int4B)(0.5*Chan[0].RngBins);
  }
  if (KernSize == DEFAULT_VAL_I) {
    KernSize = 1;
  }
  else if (KernSize % 2 == 0) {  /* ensure odd */
    KernSize++;
  }
  
  if (TerrainAlt == DEFAULT_VAL_F && 
      Chan[0].AveTerrainHeight!=DEFAULT_VAL_F) {
    TerrainAlt = Chan[0].AveTerrainHeight;
  }
  else if (TerrainAlt == DEFAULT_VAL_F) {
    TerrainAlt = 0.0;
  }  

  /* Misc */
  TimeIncr = (double)Chan[0].PresumRatio / Chan[0].MasterPRF;

  /* MESSAGES (more later) */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"LBR file name                 : %s\n",LBRFileName);
  fprintf(msg,"LBR motion log file           : %s\n",LogFileName);
  fprintf(msg,"Output rng shift binary file  : %s\n",OutFileName);
  if (OutTextFileFlg == 'Y') {
    fprintf(msg,"Output rng shift text file    : %s\n",OutTextFileName);
  }
  else {
    fprintf(msg,"Output rng shift text file    : [not written]\n");
  }
  fprintf(msg,"Smoothing kernel size         : %ld\n",KernSize);
  fprintf(msg,"Antenna direction             : %c\n",AntennaDirn);
  fprintf(msg,"Radar                         : %c\n",Radar);
  fprintf(msg,"HBR                           : %c\n",HBR);
  fprintf(msg,"Pulse Length                  : %.5e sec\n",Chan[0].PulseLen);
  fprintf(msg,"Delay0 (meas)                 : %.5e sec\n",Chan[0].Delay0Meas);
  fprintf(msg,"Range bins                    : %ld\n",Chan[0].RngBins);
  fprintf(msg,"Terrain altitude              : %f m\n",TerrainAlt); 

  /********* START PROCESSING ********/

  /* Read first two motion records */
  
  SeekP(LBRFile,":EndTimingCFG","",SEEKP_NOLIM,SEEKP_SET);
  if (!ReadLBRMotionRec(LBRFile,Radar,&PRICur,&IMUDataCur)) {
    fprintf(msg,"ERROR - no motion records in LBR file!\n");
    exit(0);
  }
  
  MotRecLength = ftell(LBRFile);

  if (!ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
    fprintf(msg,"ERROR - only one motion record found in LBR file!\n");
    exit(0);
  }  
  MotRecLength = ftell(LBRFile) - MotRecLength;

  /* Carry on reading until record found closest to start LBR PRI (note may
  be repeated PRI's at start) */
  while (PRINext <= Chan[0].StartPRI) {
    PRICur = PRINext;
    IMUDataCur = IMUDataNext;
    if (!ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      fprintf(msg,"ERROR - no motion record found for LBR start PRI!\n");
      exit(0);
    }    
  }/* end while */

  if (PRICur < Chan[0].StartPRI) {  /* ensure PRICur >= StartPRI */
    PRICur = PRINext;
    IMUDataCur = IMUDataNext;    
  }

  /* At this point the next PRI is > Chan[0].StartPRI and the 
  current PRI is >= Chan[0].StartPRI. Need to move file pointer back
  a record so that next 'current PRI' read follows the prev current PRI. */
  fseek(LBRFile,-(MotRecLength+50),SEEK_CUR); /* add 50 to be sure */
  LBRMotStartPRI = PRICur;
  LBRMotDataStart = IMUDataCur;


  /* READ and WRITE LOG of LBR motion records while current PRI <= 
  LBR endPRI.*/
  LBRNumMotRecords = 0;
  LBRAveGPSHgt = 0.0;
  ProcNumMotRecords = 0;
  ProcAveGPSHgt = 0.0;
  fprintf(LogFile,
    "Lat (deg)\tLong (deg)\tHeight (m)\tPRI ID\tGPS UTC (hh:min:sec)\n");

  while (PRICur <= Chan[0].EndPRI) {

    /* Mark motion record to be used as the first processed */
    if (!ProcMotStartPRIFound && PRICur >= ProcStartPRI && PRICur<=ProcEndPRI) {
      ProcMotStartPRI = PRICur;
      ProcMotDataStart = IMUDataCur;
      ProcMotStartPRIFound = 1;
    }

    /* Increment counters */
    if (ProcMotStartPRIFound && !ProcMotEndPRIFound) {
      ProcNumMotRecords++;
      ProcAveGPSHgt += IMUDataCur.GPSHgt;
    }
    LBRNumMotRecords++;
    LBRAveGPSHgt += IMUDataCur.GPSHgt;
 
    /* Write out IMU (INT) data to log file */
    fprintf(LogFile,"%f\t%f\t%f\t",IMUDataCur.INTLat,
            IMUDataCur.INTLong, IMUDataCur.INTHgt);        
    fprintf(LogFile,"%ld\t", PRICur);
    fprintf(LogFile,"%02ld:%02ld:%02ld\n",IMUDataCur.GPS_UTC_hrs,
            IMUDataCur.GPS_UTC_min,IMUDataCur.GPS_UTC_sec);

    /* Read next motion record */
    if (!ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      break;
    }    
    /* Mark motion record to be used as last processed */
    if (!ProcMotEndPRIFound && PRINext > ProcEndPRI) {
      ProcMotEndPRI = PRICur;
      ProcMotDataEnd = IMUDataCur;
      ProcMotEndPRIFound = 1;
    }    
    if (PRINext > Chan[0].EndPRI) {
      break;
    }
    else {
      PRICur = PRINext;
      IMUDataCur = IMUDataNext;
    }

  } /* end while PRICur */


  /* Record ProcEndPRI if not found */
  if (!ProcMotEndPRIFound) {
    ProcMotEndPRI = PRICur;
    ProcMotDataEnd = IMUDataCur;
    ProcMotEndPRIFound = 1;
  }

  if (!ProcMotStartPRIFound) {
    fprintf(msg,"ERROR - no motion records found within PRI range!\n");
    exit(0);
  }

  /* Record last current motion record */
  LBRMotEndPRI = PRICur;      
  LBRMotDataEnd = IMUDataCur;


  
  /* Calc LBR average GPS height and average ground speed. Assumes here that
  the INT value are more accurate over a long period (although not smooth). */
  LBRAveGPSHgt = LBRAveGPSHgt / LBRNumMotRecords;
  ConvToCartes(LBRMotDataEnd.INTLat * TO_RAD,
               LBRMotDataEnd.INTLong * TO_RAD,
               LBRMotDataStart.INTHgt,  /* use same height as start for calc */
               &p2);
  ConvToCartes(LBRMotDataStart.INTLat * TO_RAD,
               LBRMotDataStart.INTLong * TO_RAD,
               LBRMotDataStart.INTHgt,
               &p1);             
  LBRDist = CalcCartesDist(p1,p2);
  LBRTime = ((LBRMotEndPRI - LBRMotStartPRI)/Chan[0].PRIIDIncr)*TimeIncr;
  if (LBRTime != 0.0) {
    LBRAveGrndSpeed = LBRDist / LBRTime;
  }
  else {
    LBRAveGrndSpeed = DEFAULT_VAL_F;
  }  

  /* Calc Proc average GPS height and average ground speed */
  ProcAveGPSHgt = ProcAveGPSHgt / ProcNumMotRecords;
  ConvToCartes(ProcMotDataEnd.INTLat * TO_RAD,
               ProcMotDataEnd.INTLong * TO_RAD,
               ProcMotDataStart.INTHgt, /* use same height as start for calc */
               &p2);
  ConvToCartes(ProcMotDataStart.INTLat * TO_RAD,
               ProcMotDataStart.INTLong * TO_RAD,
               ProcMotDataStart.INTHgt,
               &p1);    
                         
  ProcMotDist = CalcCartesDist(p1,p2);
  ProcMotTime = ((ProcMotEndPRI - ProcMotStartPRI)/Chan[0].PRIIDIncr)*TimeIncr;  
  if (ProcMotTime != 0.0) {
    ProcAveGrndSpeed = ProcMotDist / ProcMotTime;
  }
  else {
    ProcAveGrndSpeed = DEFAULT_VAL_F;
  }  

  /* more messages */
  fprintf(msg,"LBR num motion records        : %ld\n",LBRNumMotRecords);  
  fprintf(msg,"LBR PRI ID increment          : %ld\n",Chan[0].PRIIDIncr);
  fprintf(msg,"LBR start PRI (all data)      : %ld\n",Chan[0].StartPRI);
  fprintf(msg,"LBR end PRI (all data)        : %ld\n\n",Chan[0].EndPRI);

  fprintf(msg,"G2 / LBR start PRI to process : %ld / %ld\n",
          StartG2PRI,ProcStartPRI);
  fprintf(msg,"G2 / LBR end PRI to process   : %ld / %ld\n",
          EndG2PRI,ProcEndPRI);
  fprintf(msg,"Start IMU record processed    : %ld\n",ProcMotStartPRI);
  fprintf(msg,"End IMU record processed      : %ld\n",ProcMotEndPRI);
  fprintf(msg,"Num motion records processed  : %ld\n",ProcNumMotRecords);
  fprintf(msg,"Average ground speed          : %f m/s (INT data)\n",
          ProcAveGrndSpeed);
  fprintf(msg,"Horizontal motion distance    : %f m (INT data)\n",
          ProcMotDist);  
  fprintf(msg,"Motion duration               : %f secs\n",ProcMotTime);   
  fprintf(msg,"Average GPS height            : %f m\n", ProcAveGPSHgt);

  /* Write out extra info to end of log file */
  fprintf(LogFile,"\nRadar            : %c\n",Radar);
  fprintf(LogFile,"HBR              : %c\n",HBR);
  fprintf(LogFile,"CarrierFreq      : %.5e Hz\n",Chan[0].CarrierFreq);
  fprintf(LogFile,"Polarization     : %s\n",Chan[0].Polarization);
  fprintf(LogFile,"Pulse Length     : %.5e sec\n",Chan[0].PulseLen);
  fprintf(LogFile,"Delay0 (meas)    : %.5e sec\n",Chan[0].Delay0Meas);
  fprintf(LogFile,"Range bins       : %ld\n",Chan[0].RngBins);
  fprintf(LogFile,"PRI ID increment : %ld\n\n",Chan[0].PRIIDIncr);

  fprintf(LogFile,
    "LBR average ground speed    : %f m/s (INT data)\n", LBRAveGrndSpeed);
  fprintf(LogFile,
    "LBR horizontal motion dist. : %f m (INT data)\n",LBRDist);  
  fprintf(LogFile,
    "LBR motion duration         : %f secs\n",LBRTime);   
  fprintf(LogFile,
    "LBR average GPS height      : %f m\n", LBRAveGPSHgt);
  fprintf(LogFile,
    "LBR num motion records      : %ld\n",LBRNumMotRecords);  
  fprintf(LogFile,
    "LBR start PRI               : %ld\n",Chan[0].StartPRI);
  fprintf(LogFile,
    "LBR end PRI                 : %ld\n\n",Chan[0].EndPRI);


  /********************************************************************
    Perform motion comp calcs using inert data (smoother than INT data)
  ********************************************************************/

  /* Calculate gradients for equation of reference flight path */
  ConvToCartes(ProcMotDataEnd.INTLat * TO_RAD,  /* jmh from inert */
               ProcMotDataEnd.INTLong * TO_RAD, /* jmh from inert */
               ProcMotDataStart.INTHgt, /* use same height as start for calc */
               &RefEndPt);
  ConvToCartes(ProcMotDataStart.INTLat * TO_RAD, /* jmh from inert */
               ProcMotDataStart.INTLong * TO_RAD, /* jmh from inert */
               ProcMotDataStart.INTHgt,
               &RefStartPt);         
  RefM.x = (RefEndPt.x - RefStartPt.x) / ProcMotTime;   /* calc gradients */
  RefM.y = (RefEndPt.y - RefStartPt.y) / ProcMotTime;
  RefM.z = (RefEndPt.z - RefStartPt.z) / ProcMotTime;    

  /* Calculate start position of ground reference track */
  NomFlight.x = RefEndPt.x - RefStartPt.x;  /* vector along ref flight path */
  NomFlight.y = RefEndPt.y - RefStartPt.y;
  NomFlight.z = RefEndPt.z - RefStartPt.z;
  NomSlantRng = (C/2.0)*(Chan[0].Delay0Meas + 
                (double)RefRngBin/Chan[0].ADFreq);
  NomLookAngleRad = acos((ProcMotDataStart.INTHgt-TerrainAlt)/NomSlantRng);
  CalcStartGrndRef(AntennaDirn, 
                   NomSlantRng,
                   NomLookAngleRad, 
                   RefStartPt,
                   NomFlight,
                   &GrndRefPt);

/*
 printf("RefStart LatDeg / LongDeg / Hgt = %.9f/ %.9f/ %.9f\n",
        ProcMotDataStart.INTLat,ProcMotDataStart.INTLong,ProcMotDataStart.INTHgt); 
 printf("RefEnd   LatDeg / LongDeg / Hgt = %.9f/ %.9f/ %.9f\n",
        ProcMotDataEnd.INTLat,ProcMotDataEnd.INTLong,ProcMotDataStart.INTHgt); 
 printf("RefStartPt x / y / z = %f / %f / %f\n",
        RefStartPt.x,RefStartPt.y,RefStartPt.z);
 printf("RefEndPt   x / y / z = %f / %f / %f\n",
        RefEndPt.x,RefEndPt.y,RefEndPt.z);
 printf("RefM  x / y / z = %f / %f / %f\n",
        RefM.x,RefM.y,RefM.z);
 printf("NomFlight  x / y / z = %f / %f / %f\n",
        NomFlight.x,NomFlight.y,NomFlight.z);        
 printf("Proc mot distance = %f\n",ProcMotDist);
 printf("C = %f\n",C);
 printf("Chan[0].Delay0Meas = %.8f secs\n",Chan[0].Delay0Meas);
 printf("Chan[0].ADFreq = %f Hz\n",Chan[0].ADFreq); 
 printf("GrndRefPt  x / y / z = %f / %f / %f\n",
        GrndRefPt.x,GrndRefPt.y,GrndRefPt.z);
 */


  /*******CALCULATE THE RANGE SHIFTS*********/

  /* Allocate mem for the rng shifts */
  ProcMotShift = (double *)malloc(sizeof(double)*ProcNumMotRecords);
  ProcSmoothShift = (double *)malloc(sizeof(double)*ProcNumMotRecords);
  ProcMotShiftPRI = (Int4B *)malloc(sizeof(Int4B)*ProcNumMotRecords);
  if (ProcMotShift == NULL || ProcSmoothShift==NULL 
     || ProcMotShiftPRI == NULL) {
    fprintf(msg,"ERROR - in shift array allocation!\n");
    exit(0);
  } 
 
  /* Reset LBR file pointer to start of motion records and search for 
  first motion record (the range shift to the ProMotStartPRI should
  come out to be zero). */
  SeekP(LBRFile,":EndTimingCFG","",SEEKP_NOLIM,SEEKP_SET);

  PRICur = ProcMotStartPRI;  /* needed for first GrndRefPt calc */
  PRINext = -999;
  while (PRINext <= ProcMotStartPRI) {
    if (!ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      fprintf(msg,"ERROR - mocomp calc record not found!\n");
      exit(0);      
    }
  }  /* end while */

  ProcMotShiftPRI[0] = ProcMotStartPRI;
  ProcMotShift[0] = 0.0;

  RngShiftCur = 0.0;     
  /* Perform motion calculations (interpolate between rng shifts) */ 
  for (MotRec=0; MotRec< ProcNumMotRecords-1; MotRec++) {

    /* Calculate the GrndRefPt at PRINext */
    PRIChange = (Int4B)( (PRINext - PRICur)/Chan[0].PRIIDIncr ); 
    GrndRefPt.x = GrndRefPt.x + RefM.x*TimeIncr*PRIChange;
    GrndRefPt.y = GrndRefPt.y + RefM.y*TimeIncr*PRIChange;
    GrndRefPt.z = GrndRefPt.z + RefM.z*TimeIncr*PRIChange;

    /* Calc distance from aircraft to GrndRefPt at PRINext*/
    ConvToCartes(IMUDataNext.INTLat * TO_RAD, /* jmh from inert */
                 IMUDataNext.INTLong * TO_RAD, /* jmh from inert */
                 IMUDataNext.INTHgt,
                 &AircraftPosn);  
                       
    RngShiftNext = NomSlantRng - CalcCartesDist(AircraftPosn,GrndRefPt);
    ProcMotShiftPRI[MotRec+1] = PRINext;
    ProcMotShift[MotRec+1] = RngShiftNext;
    
    /* Read next PRI */
    PRICur = PRINext;
    RngShiftCur = RngShiftNext; 
    if (!ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      break;      
    }   
     
  } /* end for MotRec loop */

  /* Smooth the range shift data */
  if (!MovingAveSmooth(ProcMotShift,
                       ProcSmoothShift,
                       ProcNumMotRecords,
                       KernSize) ) {
    fprintf(msg,"ERROR - in moving ave function call!\n");
    exit(0);
  }


  /******Write out range shifts to file (first block is for output of actual
  range shifts for each motion record) *********/

#if WRITE_LBR_MOT_SHIFTS
  for (MotRec=0; MotRec < ProcNumMotRecords; MotRec++) {
    if (OutTextFileFlg == 'Y')
      fprintf(OutTextFile,"%ld\t%f\n",ProcMotShiftPRI[MotRec],
              ProcSmoothShift[MotRec]);  
  }
  fprintf(msg,"Early exit after writing LBR motion shifts - recompile\n");
  fprintf(msg,"for full operation with WRITE_LBR_MOT_SHIFTS disabled\n");
  exit(0);
#endif


  /* Initialise output PRI ID */
  OutPRI = (Int4B)( (ProcStartPRI - Chan[0].StartPRI)/Chan[0].PRIIDIncr );

  /* Write zero range shift for PRIs < ProcMotStartPRI */
  RngShift = 0.0;
  for (i=0; i<=(Int4B)( (ProcMotStartPRI-ProcStartPRI)/Chan[0].PRIIDIncr ); i++) {
    if (OutTextFileFlg == 'Y')
      fprintf(OutTextFile,"%f\n",RngShift); /* jmh from inert */
      /* fprintf(OutTextFile,"%ld\t%f\n",OutPRI,RngShift);*/ /* jmh from inert */
    fwrite(&OutPRI,sizeof(Int4B),1,OutFile);
    OutPRI++;
    fwrite(&RngShift,sizeof(float),1,OutFile);
  }

  /* Write out rng shifts over region of valid mot records */
  for (MotRec=0; MotRec< ProcNumMotRecords-1; MotRec++) {

    PRICur = ProcMotShiftPRI[MotRec];
    RngShiftCur = ProcSmoothShift[MotRec];
    PRINext = ProcMotShiftPRI[MotRec+1];
    RngShiftNext = ProcSmoothShift[MotRec+1];

    PRIChange = (Int4B)( (PRINext - PRICur)/Chan[0].PRIIDIncr ); 
    RngShiftIncr = (RngShiftNext - RngShiftCur)/PRIChange;

    /* Write out for each G2 PRI ID */
    for (i=0; i < PRIChange; i++) {

      /* Calc range shift */
      RngShift = (float)(RngShiftCur + (i+1)*RngShiftIncr);
      if (fabs(RngShift) > fabs(MaxRngShift)) {
         MaxRngShift = RngShift;
         MaxRngShiftIndx = OutPRI;          
      }
  
      /* Write range shift to output file */
      if (OutTextFileFlg == 'Y')
        fprintf(OutTextFile,"%f\n",RngShift); /* jmh from inert */
        /* fprintf(OutTextFile,"%ld\t%f\n",OutPRI,RngShift);*/ /* jmh from inert */
      fwrite(&OutPRI,sizeof(Int4B),1,OutFile);
      OutPRI++;
      fwrite(&RngShift,sizeof(float),1,OutFile);

    } /* end for i loop */

 } /* end for MotRec */

  /* Write extra range shifts at end, if necessary, repeat last value */
  for (i=0; i<(Int4B)((ProcEndPRI-ProcMotEndPRI)/Chan[0].PRIIDIncr); i++) {
    if (OutTextFileFlg == 'Y')
      fprintf(OutTextFile,"%f\n",RngShift); /* jmh from inert */
      /* fprintf(OutTextFile,"%ld\t%f\n",OutPRI,RngShift);*/ /* jmh from inert */
    fwrite(&OutPRI,sizeof(Int4B),1,OutFile);
    OutPRI++;
    fwrite(&RngShift,sizeof(float),1,OutFile);
  } 

  fprintf(msg,"Reference range bin           : %ld\n",RefRngBin);
  fprintf(msg,"NomSlantRng                   : %f m\n",NomSlantRng);
  fprintf(msg,"NomLookAngle                  : %f rad (%f deg)\n",
          NomLookAngleRad, NomLookAngleRad*TO_DEG);
  fprintf(msg,"Max range shift               : %f m\n",MaxRngShift);
  fprintf(msg,"Max range shift G2 PRI        : %ld\n",MaxRngShiftIndx);
  
  /* Clean up */
  free(ProcMotShift);
  free(ProcMotShiftPRI);

  if (OutTextFileFlg == 'Y')
    fclose(OutTextFile);
  
  fclose(LBRFile);
  fclose(OutFile);
  fclose(LogFile);

  return(1);

}  /***** end main prog or function ******/



/*===============FUNCTIONS================*/

/*============ReadLBRMotionRec============*/
/* Function to parse LBR file from current posn to read 
   next motion record and PRI ID. Removed \n from SeekP call as not
   picking up first PRI. Note it returns without altering params, if no
   "Sync IDA" found */

Int2B ReadLBRMotionRec(FILE *LBRFile,
                       char Radar,
                       Int4B *PRI,
                       struct IMUStruct *IMUData)
{

  fseek(LBRFile,-1,SEEK_CUR);  /* move back one byte, otherwise misses 
                                  alternate records for some reason */   
  if (!SeekP(LBRFile,"Sync IDA","",SEEKP_NOLIM,SEEKP_CUR)) /* no IMU data */ 
    return(0);  
  if (Radar == 'B') 
    SeekP(LBRFile,"IDB","",50,-1); /* find 'IDB' for Radar B */
  fscanf(LBRFile,"%ld",PRI);  /* Read in PRI (is ptr)*/  
  SeekP(LBRFile,"\n","",-1,-1);
  ReadIMUMsg5(LBRFile,IMUData); /* Read and convert IMU record (is ptr)*/

  return(1);
}                       


/*============ConvToCartes================*/
/* Function to convert from lat,long and ellipsoidal height to
   3-D Cartesian coords */
Int2B ConvToCartes(double LatRad, double LongRad, double h,
                   struct cartes *Pt)
{
  double N;  

  N = WGS84_a / sqrt(1 - WGS84_e_sq * sin(LatRad)*sin(LatRad));
  Pt->x = (N+h)*cos(LatRad)*cos(LongRad);
  Pt->y = (N+h)*cos(LatRad)*sin(LongRad);
  Pt->z = (N*(1-WGS84_e_sq) + h)*sin(LatRad);

  return(1);
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

/*===========CalcTimeSecs=================*/
/* Converts GPS hours, min, secs in IMUStruct to secs*/

double CalcTimeSecs(struct IMUStruct Data)
{
 double secs;
 secs = 3600.0*Data.GPS_UTC_hrs + 60.0*Data.GPS_UTC_min + Data.GPS_UTC_sec;
 return (secs);
}
   
/*============CalcStartGrndRef============*/
/* Function to calc the start reference position in the ground swath.
   Uses the Maple-calculated solution. */
Int2B CalcStartGrndRef(char Dirn, /* antenna dirn (L/R) */
                       double NomSlantRng,
                       double LookAngleRad, 
                       struct cartes StartRef, /* coords of ref start pt */
                       struct cartes NomFlight, /* vector along ref flight path */
                       struct cartes *StartGrndRef) /* coords of grnd ref start */
{

  double a,b,c,    /* coeffs of quadratic to solve */
         h,k,m,    /* StartRef (x,y,z) */
         q,r,s,    /* NomFlight (x,y,z) */
         t,        /* NomSlantRng */
         u,        /* R dot StartRef */
         p1,p2,    /* possible z components for R */
         CrossCalc1,
         CrossCalc2;
  struct cartes R,R1,R2; /* vector from StartRef posn to StartGrndRef 
                            and two possibles */
                 
  h = StartRef.x;
  k = StartRef.y;
  m = StartRef.z;
  q = NomFlight.x;
  r = NomFlight.y;
  s = NomFlight.z;
  t = NomSlantRng;

  /* Calc u (R.StartRef) */
  u  = sqrt(StartRef.x*StartRef.x + StartRef.y*StartRef.y +
            StartRef.z*StartRef.z);
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
  CrossDot(StartRef,R1,NomFlight,&CrossCalc1);
  CrossDot(StartRef,R2,NomFlight,&CrossCalc2);
  
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

  /* Calc StartGrndRef */
  StartGrndRef->x = StartRef.x + R.x;
  StartGrndRef->y = StartRef.y + R.y;
  StartGrndRef->z = StartRef.z + R.z;

  return(1);

} /* end function CalcStartGrndRef */            


/*===========QuadRoots===============*/

/* Calcs the roots of a quadratic */
Int2B QuadRoots(double a, double b, double c, 
                double *root1, double *root2) 
{
  *root1 = (-b+sqrt(b*b-4*a*c))/(2*a);
  *root2 = (-b-sqrt(b*b-4*a*c))/(2*a);

  return(1);
}  /* end function QuadRoots */                


/*==========CalcR  ==================*/

/* Part of Maple-calculated solution for the SASAR motion 
   comp */
Int2B CalcR(double h, double k, double m, 
            double q, double r, double s,
            double p, double u, 
            struct cartes *R)            
{

  R->x = - (r*p*m - r*u - p*s*k)/(h*r - q*k);
  R->y = - (h*p*s - p*m*q + u*q)/(h*r - q*k);
  R->z = p;

  return(1);
}  /* end function CalcR */

/*==========CrossDot=================*/

/* Calcs (StartRef x R) . NomFlight */
Int2B CrossDot(struct cartes StartRef,
               struct cartes R,
               struct cartes NomFlight,
               double *CrossCalc)
{

  double a,b,c,
         d,e,f;

  a = StartRef.x;
  b = StartRef.y;
  c = StartRef.z;
  d = R.x;
  e = R.y;
  f = R.z;

  *CrossCalc  = (b*f-c*e)*NomFlight.x;
  *CrossCalc += (c*d-a*f)*NomFlight.y;
  *CrossCalc += (a*e-b*d)*NomFlight.z;

  return(1);

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
  if (InEx==NULL) return(0);

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
  return(1);

} /* end function MovingAveSmooth */

  