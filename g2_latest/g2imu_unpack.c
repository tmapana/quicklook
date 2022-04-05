/*==========================================
COPYRIGHT: UCT Radar Remote Sensing Group 1999
FILE NAME: g2moc.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Extract motion records from Marconi FIN3110
IMU and write part of the data to ASCII file. 

Note that for some unknown reason, the "Sync IDA" of the last motion
record of the  SASAR LBR file is not found. Seems to be a glitch in the
way the record is written as, if the LBR fileis hacked to add another
record at the end, it still misses the old last record, but picks up
the new last record.

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 1999-07-26
COMMENTS:
Extracted/hacked from the mocomp program.

=========================================*/

#include "g2func.h"
#include "g2parse.h"

#define PROG_VERSION "1.0"

#define DEFAULT_OUT_EXT "_unpack"
#define WGS84_a 6378137.0
#define WGS84_b 6356752.314
#define WGS84_e_sq 6.694380066e-03
#define TO_RAD 0.01745329252
#define TO_DEG 57.29577951


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
Int2B ReadLBRMotionRec(FILE *LBRFile,
                       char Radar,
                       Int4B *PRI,
                       struct IMUStruct *IMUData);


/*===========MAIN PROG/FUNCTION===========*/

Int2B main (int argc,char *argv[])
{
  FILE *msg=stdout, *LBRFile=NULL, *OutFile=NULL;
  char ProgVersion[STRING_SPACE]="",LBRFileName[STRING_SPACE],
       OutFileName[STRING_SPACE]="",
       HBR='3',   /* default value */
       Radar='B'; /* default value */
  Int4B i,
        LBRMotEndPRI,
        LBRMotStartPRI,
        LBRNumMotRecords=0,
        MotRecLength,
        PRICur,
        PRINext;
  double LBRAveGPSHgt,
         LBRDist,
         LBRTime,
         LBRAveGrndSpeed,
         TimeIncr;    
  struct CntrlStruct ParseCntrl;
  struct DataChannelStruct Chan[1];
  struct IMUStruct LBRMotDataEnd,
                   LBRMotDataStart,
                   IMUDataCur,
                   IMUDataNext;
  struct cartes p1,p2;
                              
        
  /* Assign message output */
  msg = stdout;

  fprintf(msg,"\n----------------\n");
  fprintf(msg,"Prog: imu_unpack (Ver. %s)\n",PROG_VERSION);
  fprintf(msg,"Code: J.M. Horrell (C) UCT Radar Remote Sensing Group 1999\n");

  if (argc < 2) {
    fprintf(msg,
      "Writes out IMU data in ASCII format after reading from SASAR LBR file.\n\n");
    fprintf(msg,"USAGE : imu_unpack [LBR_file]\n");
    fprintf(msg,"Optional params (use at end of line):\n");
    fprintf(msg,"<OutF=t>      - t is the output file name (else default name)\n");
    fprintf(msg,"<Radar=A/B>   - B is the default\n");
    fprintf(msg,"<HBR=n>       - in range from 1 - 6, 3 is default\n");
    fprintf(msg,"e.g. 'imu_unpack' lbrfile.001 HBR=4'\n");
    exit(1);
  } 


  /* check for optional parameters */
  for (i=2; i<argc; i++) {
    if (strncmp(argv[i],"OutF=",5)==0) {
      sscanf(argv[i],"OutF=%s",OutFileName);
    } 
    else if (strncmp(argv[i],"Radar=",6)==0) {
      sscanf(argv[i],"Radar=%c",&Radar);
    } 
    else if (strncmp(argv[i],"HBR=",4)==0) {
      sscanf(argv[i],"HBR=%c",&HBR);
    } 
    else {
      fprintf(msg,"ERROR - parameter %s unknown!\n",argv[i]);
      exit(1);
    }      
  }  /* end for i loop */


  /* Create file names and log file name and open files */
  strcpy(LBRFileName,argv[1]);
  if (strcmp(OutFileName,"")==0) {   /* if not exist, create outfile name */ 
    strcpy(OutFileName,LBRFileName);
    strcat(OutFileName,DEFAULT_OUT_EXT);
  }  
  
  if(( LBRFile = fopen(LBRFileName,"rb")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open LBR file %s!\n",
            LBRFileName); 
    exit(1); 
  }
  if(( OutFile = fopen(OutFileName,"w")) == NULL) { 
    fprintf(msg,"ERROR - Unable to open output file %s!\n",
            OutFileName); 
    exit(1); 
  }

  /* Scan through LBR file for relevant params using g2parse.c 
     functions */
  ParseCntrl.NumDataChannels=1;
  ParseCntrl.CurrentChannel=0;
  Chan[0].Radar = Radar;
  Chan[0].HBR = HBR;
  ParseDataChannels(msg,LBRFile,&ParseCntrl,Chan);

  /* Misc */
  TimeIncr = (double)Chan[0].PresumRatio / Chan[0].MasterPRF;

  /* MESSAGES (more later) */
  fprintf(msg,"\nMESSAGES:\n");
  fprintf(msg,"LBR file name                 : %s\n",LBRFileName);
  fprintf(msg,"Output file name              : %s\n",OutFileName);
  fprintf(msg,"Radar                         : %c\n",Radar);
  fprintf(msg,"HBR                           : %c\n",HBR);
  fprintf(msg,"Pulse Length                  : %.5e sec\n",Chan[0].PulseLen);
  fprintf(msg,"Delay0 (meas)                 : %.5e sec\n",Chan[0].Delay0Meas);
  fprintf(msg,"Range bins                    : %ld\n",Chan[0].RngBins);

  /********* START PROCESSING ********/

  /* Read first two motion records */
  
  SeekP(LBRFile,":EndTimingCFG","",SEEKP_NOLIM,SEEKP_SET);
  if (ReadLBRMotionRec(LBRFile,Radar,&PRICur,&IMUDataCur)) {
    fprintf(msg,"ERROR - no motion records in LBR file!\n");
    exit(1);
  }

  MotRecLength = ftell(LBRFile);

  if (ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
    fprintf(msg,"ERROR - only one motion record found in LBR file!\n");
    exit(1);
  }  
  MotRecLength = ftell(LBRFile) - MotRecLength;

  /* Carry on reading until record found closest to start LBR PRI (note may
  be repeated PRI's at start) */
  while (PRINext <= Chan[0].StartPRI) {
    PRICur = PRINext;
    IMUDataCur = IMUDataNext;
    if (ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      fprintf(msg,"ERROR - no motion record found for LBR start PRI!\n");
      exit(1);
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

  /* Write header for output file */
  WriteMotionHeader(OutFile);

  /* READ and WRITE LBR motion records while current PRI <= LBR endPRI.*/
  LBRNumMotRecords = 0;
  LBRAveGPSHgt = 0.0;

 while (PRICur <= Chan[0].EndPRI) {

    LBRNumMotRecords++;
    LBRAveGPSHgt += IMUDataCur.GPSHgt;

    /* Write part of motion record to ASCII file */  
    WriteMotion(OutFile,IMUDataCur,PRICur);

    /* Read next motion record */
    if (ReadLBRMotionRec(LBRFile,Radar,&PRINext,&IMUDataNext)) {
      break;
    }    
    /* Mark motion record to be used as last processed */
    if (PRINext > Chan[0].EndPRI) {
      break;
    }
    else {
      PRICur = PRINext;
      IMUDataCur = IMUDataNext;
    }

  } /* end while PRICur */

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

  /* more messages */
  fprintf(msg,"LBR num motion records        : %ld\n",LBRNumMotRecords);  
  fprintf(msg,"LBR PRI ID increment          : %ld\n",Chan[0].PRIIDIncr);
  fprintf(msg,"LBR start PRI (all data)      : %ld\n",Chan[0].StartPRI);
  fprintf(msg,"LBR end PRI (all data)        : %ld\n\n",Chan[0].EndPRI);
  fprintf(msg,"LBR average ground speed      : %f m/s (INT data)\n",
          LBRAveGrndSpeed);
  fprintf(msg,"LBR horizontal motion dist.   : %f m (INT data)\n",LBRDist);  
  fprintf(msg,"LBR motion duration           : %f secs\n",LBRTime);   
  fprintf(msg,"LBR average GPS height        : %f m\n", LBRAveGPSHgt);
  fprintf(msg,"LBR num motion records        : %ld\n",LBRNumMotRecords);  
  fprintf(stdout,"DONE!\n");

  /* clean up */
  fclose(LBRFile);
  fclose(OutFile);

  return(0);

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
  if (SeekP(LBRFile,"Sync IDA","",SEEKP_NOLIM,SEEKP_CUR)) /* no IMU data */ 
    return(1);  
  if (Radar == 'B') 
    SeekP(LBRFile,"IDB","",50,-1); /* find 'IDB' for Radar B */
  fscanf(LBRFile,"%ld",PRI);  /* Read in PRI (is ptr)*/  
  SeekP(LBRFile,"\n","",-1,-1);
  ReadIMUMsg5(LBRFile,IMUData); /* Read and convert IMU record (is ptr)*/

  return(0);
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

  return(0);
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
   

  