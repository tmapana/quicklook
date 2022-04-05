/* Code to write SASAR VHF data to SASAR CEOS Format */
/* Author: J.M. Horrell */

/* Ver. 1a - 19970805 */
/* Ver. 1b - 19970819 - tidy up for DJGPP compiler. Major
   additions everywhere. Parse LBR file to find data */
/* Ver 1c - major additions and changes */
/* Ver 1d - 19970925 - finish off major additions */
/* Ver G2Ceos (19971205) - rename to as part of the new g2 code
   (post ATP) - different included header file name, PROG_VERSION added */
/* Ver 1998-11-24 - option to extract inertial data to ASCII file and to skip
   writing the bulk of the IMO file */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#include"g2func.h"

/* Misc defns limited to this file  */
#define PROG_VERSION "1997-12-05"
#define MAX_DATA_CHANNELS 6
#define COMMENT_SPACE     200
#define FEET2METRES       0.3048
#define SOFTWARE_ID       "G2CEOS"
#define CREATING_COUNTRY  "RSA"
#define CREATING_AGENCY   "UCT"
#define CREATING_FACILITY "RRSG"

/* Fixed fields, for now */
#define PHYS_VOL_ID       " "
#define LOGICAL_VOL_ID    " "
#define VOL_SET_ID        " "
#define LOGICAL_VOL       1
#define DEBUG             0
#define READSIM           0
#define SIMFILE           "atpsim1a.bin"
#define FILL_IMO          0
#define EXTRACT_MOTION    1

/* STRUCTURES */

/* Stores misc control info for transcription */
struct CntrlStruct {
  FILE *LBRFile;
  FILE *msg;
  Int2B NumDataChannels;
  Int2B CurrentChannel;
  Int2B LogicalVol;
  Int2B FileNumOfType;
  Int2B CEOSFileNum;
};

/* Stores info for data channel */
struct DataChannelStruct {
  /* General config info */
  char Radar;
  char HBR;
  char Sensor[STRING_SPACE];
  char Location[STRING_SPACE];
  char Mission[STRING_SPACE];
  char Track[STRING_SPACE];
  char Date[11];
  double NomGroundSpeed;    /*(m/s)*/
  double AveTerrainHeight;  /*(m)*/
  Int4B GenCommentLen;
  char GenComment[COMMENT_SPACE];
  /* LBR config info */
  char PRICorrectionFlg;
  /* TMC config info */
  double ADFreq;            /*(MHz)*/   
  Int4B RngBins;
  /* HBR config info */
  Int2B Bits;
  double MasterPRF;         /*(Hz)*/
  double CarrierFreq;       /*(MHz)*/
  char Polarization[3];
  char RngCompressedFlg;
  char Modulation[STRING_SPACE];
  double PulseLen;          /*(usec)*/
  double PulseBandwidth;    /*(MHz)*/
  double PulsePower;        /*(kW)*/
  double AzBeamwidth;       /*(deg)*/
  double ElevBeamwidth;     /*(deg)*/    
  double DepressAngle;      /*(deg)*/ 
  double SquintAngle;       /*(deg)*/
  double CalToneLevel;
  double CalFreq;           /*(MHz)*/
  double SysLosses;         /*(dB)*/
  double Delay0Meas;        /*(usec)*/
  double Delay0Calc;        /*(usec)*/
  char STCFlg;
  Int2B DGCChannel;
  double RxNoiseFig;        /*(dB)*/
  Int2B DCOffset;
  Int2B PresumRatio;
  Int4B DigGain;
  char HBRComment[COMMENT_SPACE];
  /* Misc */
  Int4B StartPRI;
  Int4B EndPRI;
  Int4B PRIs;
  Int4B SceneCenUTC_hrs;    /*(hrs)*/
  Int4B SceneCenUTC_min;    /*(min)*/
  Int4B SceneCenUTC_sec;    /*(sec)*/
  double SceneCenLat;       /*(deg)*/
  double SceneCenLong;      /*(deg)*/
  double SceneCenTrueHead;  /*(deg)*/
};

/* Structure for IMU/GPS raw binary data */
struct Msg5Struct
  {
  Int4B GPS_UTC;
  Int4B InertLat;
  Int4B InertLong;
  Int4B BaroInertHgt;
  Int4B InertVelNorth;
  Int4B InertVelEast;
  Int4B InertVelVert;
  Int4B INTLat;
  Int4B INTLong;
  Int4B INTHgt;
  Int4B INTVelNorth;
  Int4B INTVelEast;
  Int4B INTVelVert;
  Int4B GPSLat;
  Int4B GPSLong;
  Int4B GPSHgt;
  Int4B GPSVelNorth;
  Int4B GPSVelEast;
  Int4B GPSVelVert;
  unsigned Int2B  GPSRxStatus;
  Int4B GPSEstHorizErr;
  Int4B GPSEstVertErr;
  Int4B InertRollAngle;
  Int4B InertPitchAngle;
  Int4B InertHeadingTrue;
  Int4B INTRollAngle;
  Int4B INTPitchAngle;
  Int4B INTHeadingTrue;
  };


/* Structure for IMU/GPS data, converted to readable format  */
struct IMUStruct
  {
  Int4B GPS_UTC_hrs;     /*(hrs)*/
  Int4B GPS_UTC_min;     /*(min)*/
  Int4B GPS_UTC_sec;     /*(sec)*/ 
  double InertLat;       /*(deg)*/
  double InertLong;      /*(deg)*/
  double BaroInertHgt;   /*(m)*/
  double InertVelNorth;  /*(m/s)*/
  double InertVelEast;   /*(m/s)*/
  double InertVelVert;   /*(m/s)*/
  double INTLat;         /*(deg)*/
  double INTLong;        /*(deg)*/
  double INTHgt;         /*(m)*/
  double INTVelNorth;    /*(m/s)*/
  double INTVelEast;     /*(m/s)*/ 
  double INTVelVert;     /*(m/s)*/
  double GPSLat;         /*(deg)*/
  double GPSLong;        /*(deg)*/
  double GPSHgt;         /*(m)*/  
  double GPSVelNorth;    /*(m/s)*/
  double GPSVelEast;     /*(m/s)*/
  double GPSVelVert;     /*(m/s)*/
  unsigned Int2B  GPSRxStatus;
  double GPSEstHorizErr;    /*(m)*/
  double GPSEstVertErr;     /*(m)*/
  double InertRollAngle;    /*(deg)*/
  double InertPitchAngle;   /*(deg)*/
  double InertHeadingTrue;  /*(deg)*/
  double INTRollAngle;      /*(deg)*/
  double INTPitchAngle;     /*(deg)*/
  double INTHeadingTrue;    /*(deg)*/
  };

/**********************************************************************/
void ConstructFileName(char FileType[],Int2B LogicalVol,
		      Int2B FileNumOfType,char FileName[])
{

  char tmps[80];

  /* Construct  file name */
  sprintf(FileName,"%-3s","sas"); /* add SASAR designator */
  sprintf(tmps,"%02d",LogicalVol);
  strcat(FileName,tmps);  /* add logical vol no. */
  if (strcmp(FileType,"VDF")==0)
    strcat(FileName,"vdf.");  /* add VDF designator */
  else if (strcmp(FileType,"SAL")==0)
    strcat(FileName,"sal.");  /* add SARL designator */
  else if (strcmp(FileType,"IMO")==0)
    strcat(FileName,"imo.");  /* add IMOP designator */
  else if (strcmp(FileType,"NVD")==0)
    strcat(FileName,"nvd.");  /* add SARL designator */
  sprintf(tmps,"%03d",FileNumOfType);
  strcat(FileName,tmps);  /* add file number of this type */

  return;
}  /* end ConstructFileName function */


/*****************************************/
/* Converts from raw Msg 5 format to understandable */
double ConvMsg5Val(Int4B RawVal, Int4B MSB)
{
  return ( ((double)ConvertLong(RawVal)) / (pow(2,31)/(double)MSB) );
}  

/***************************************************************/
/* Function to read in binary data message from IMU/GPS
   FIN 3110 system of Marconi and convert to readable format */
void ReadIMUMsg5(struct CntrlStruct *Gen,
	         struct IMUStruct *IMU)
{
  struct Msg5Struct Msg5;
  Int4B tmp4b;
  
  /* read in raw binary IMU values */
  fread( &Msg5.GPS_UTC, sizeof(Msg5), 1, Gen->LBRFile);

  /* Convert to readable format and in units of hrs,min,sec,deg,metres,secs */
  tmp4b = (Int4B) ConvMsg5Val(Msg5.GPS_UTC,131072 );
  IMU->GPS_UTC_hrs  = tmp4b / 3600.0;
  tmp4b -=  IMU->GPS_UTC_hrs * 3600.0;
  IMU->GPS_UTC_min = tmp4b / 60.0;
  tmp4b -= IMU->GPS_UTC_min * 60.0;
  IMU->GPS_UTC_sec = tmp4b;

  IMU->InertLat         = ConvMsg5Val(Msg5.InertLat,180);
  IMU->InertLong        = ConvMsg5Val(Msg5.InertLong,180);
  IMU->BaroInertHgt     = FEET2METRES*ConvMsg5Val(Msg5.BaroInertHgt,81920); 
  IMU->InertVelNorth    = ConvMsg5Val(Msg5.InertVelNorth,1000); 
  IMU->InertVelEast     = ConvMsg5Val(Msg5.InertVelEast,1000); 
  IMU->InertVelVert     = ConvMsg5Val(Msg5.InertVelVert,1000); 
  IMU->INTLat           = ConvMsg5Val(Msg5.INTLat,180); 
  IMU->INTLong          = ConvMsg5Val(Msg5.INTLong,180); 
  IMU->INTHgt           = FEET2METRES*ConvMsg5Val(Msg5.INTHgt,81920); 
  IMU->INTVelNorth      = ConvMsg5Val(Msg5.INTVelNorth,1000); 
  IMU->INTVelEast       = ConvMsg5Val(Msg5.INTVelEast,1000); 
  IMU->INTVelVert       = ConvMsg5Val(Msg5.INTVelVert,1000); 
  IMU->GPSLat           = ConvMsg5Val(Msg5.GPSLat,180); 
  IMU->GPSLong          = ConvMsg5Val(Msg5.GPSLong,180); 
  IMU->GPSHgt           = FEET2METRES*ConvMsg5Val(Msg5.GPSHgt,81920); 
  IMU->GPSVelNorth      = ConvMsg5Val(Msg5.GPSVelNorth,1000); 
  IMU->GPSVelEast       = ConvMsg5Val(Msg5.GPSVelEast,1000); 
  IMU->GPSVelVert       = ConvMsg5Val(Msg5.GPSVelVert,1000); 
  IMU->GPSEstHorizErr   = ConvMsg5Val(Msg5.GPSEstHorizErr,32768); 
  IMU->GPSEstVertErr    = ConvMsg5Val(Msg5.GPSEstVertErr,32768); 
  IMU->InertRollAngle   = ConvMsg5Val(Msg5.InertRollAngle,180); 
  IMU->InertPitchAngle  = ConvMsg5Val(Msg5.InertPitchAngle,180); 
  IMU->InertHeadingTrue = ConvMsg5Val(Msg5.InertHeadingTrue,180); 
  IMU->INTRollAngle     = ConvMsg5Val(Msg5.INTRollAngle,180); 
  IMU->INTPitchAngle    = ConvMsg5Val(Msg5.INTPitchAngle,180); 
  IMU->INTHeadingTrue   = ConvMsg5Val(Msg5.INTHeadingTrue,180); 

  return ;
}  


#if EXTRACT_MOTION

/***************************************************************/
/* Function to write out IMU/GPS data to ASCII file (appends)*/
void WriteMotion(struct CntrlStruct *Gen,
	         struct IMUStruct *IMU,
                 Int4B PRI)
{
  FILE *OutFile=NULL;
  
  OutFile = fopen("IMU.txt","a");
  fprintf(OutFile,
    "%ld\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
           PRI,IMU->GPS_UTC_hrs,IMU->GPS_UTC_min,IMU->GPS_UTC_sec,
           IMU->BaroInertHgt,
           IMU->InertVelNorth,IMU->InertVelEast,IMU->InertVelVert,
           IMU->INTLat,IMU->INTLong,IMU->INTHgt,
           IMU->INTVelNorth,IMU->INTVelEast,IMU->INTVelVert,
           IMU->InertRollAngle,IMU->InertPitchAngle,IMU->InertHeadingTrue,
           IMU->INTRollAngle,IMU->INTPitchAngle,IMU->INTHeadingTrue);

  fclose (OutFile);

} /* end function WriteMotion() */


#endif

/*********************************************/
void ParseLBRFile(struct CntrlStruct *Gen,
		  struct DataChannelStruct Chan[])
{
  Int2B ii,i;
  Int4B Fatals=0,
        MaxO=0, /* MaxOffset extent of search (bytes) for SeekP func */
        PrevWarns=0,
        PRICen,
        PRINext,
        St,   /* start (Whence)  posn for SeekP function */
        Warns=0;
  struct IMUStruct IMUDat;

  /* Repeat parsing for each valid data channel */
  for (ii=0;ii<Gen->NumDataChannels;ii++)
    {
    Warns = 0; 

    fprintf(Gen->msg,"Parsing data channel %d config info...\n",ii);
    
    /* Parse general configuration info */
    PrevWarns = Warns;
    MaxO = 1000; 
    if (!SeekP(Gen->LBRFile,"// General Configuration","",-1,0))
      { fprintf(Gen->msg,
	  "ERROR - general config header not found!\n"); exit(1);}
    else
      { St = ftell(Gen->LBRFile); }    
    if (!SeekP(Gen->LBRFile,"Sensor",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Sensor,STRING_SPACE)) Warns++; 
    if (!SeekP(Gen->LBRFile,"Location",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Location,STRING_SPACE)) Warns++;  
    if (!SeekP(Gen->LBRFile,"Mission Number",":",MaxO,St)) Warns++; 
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Mission,STRING_SPACE)) Warns++; 
    if (!SeekP(Gen->LBRFile,"Track Number",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Track,STRING_SPACE)) Warns++;  
    if (!SeekP(Gen->LBRFile,"Date",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Date,STRING_SPACE)) Warns++;
    if (!SeekP(Gen->LBRFile,"Nominal Ground Speed",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].NomGroundSpeed)) Warns++;      
    if (!SeekP(Gen->LBRFile,"Average Terrain Height",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].AveTerrainHeight))Warns++;

    Chan[ii].GenCommentLen = 0;
    sprintf(Chan[ii].GenComment,"%s","");
    if (!SeekP(Gen->LBRFile,"Length of comments",":",MaxO,St))
      { Warns++; }
    else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].GenCommentLen))
      { Warns++; }  
    else if (!SeekP(Gen->LBRFile,"Comments",":",MaxO,St))
      { Warns++; }
    else
      {
      if (Chan[ii].GenCommentLen > COMMENT_SPACE)
	Chan[ii].GenCommentLen = COMMENT_SPACE;
      for (i=0;i<Chan[ii].GenCommentLen;i++)
	Chan[ii].GenComment[i]=fgetc(Gen->LBRFile);  /* read gen comment */
      Chan[ii].GenComment[Chan[ii].GenCommentLen] = '\0';
      }
    
    
    if (Warns > PrevWarns)
      fprintf(Gen->msg,
	 "WARNINGS - during gen config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Parse LBR configuration info */
    PrevWarns = Warns;
    MaxO = 500;
    if (!SeekP(Gen->LBRFile,"// LBR Configuration","",-1,0))
      { fprintf(Gen->msg,
	  "ERROR - LBR configuration header not found!\n"); exit(1);}
    else
      { St = ftell(Gen->LBRFile); }
    if (!SeekP(Gen->LBRFile,"Enable PRI",":",MaxO,St)) Warns++;
    else if (!ReadChar(Gen->LBRFile, &Chan[ii].PRICorrectionFlg)) Warns++;

    if (Warns > PrevWarns)
      fprintf(Gen->msg,
	 "WARNINGS - during LBR config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Parse TMC configuration info */
    PrevWarns = Warns;
    MaxO = 1000;
    if (!SeekP(Gen->LBRFile,"// TMC Configuration","",-1,0))
      { fprintf(Gen->msg,
	  "ERROR - TMC configuration header not found!\n"); exit(1);}
    else
      { St = ftell(Gen->LBRFile);}

    if (Chan[ii].Radar == 'A')
      {
      if (!SeekP(Gen->LBRFile,"RADAR A Sampling Rate",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].ADFreq))
	Warns++; 
      if (!SeekP(Gen->LBRFile,"RADAR A Sampling Bins",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].RngBins))
        { fprintf(Gen->msg,
	  "ERROR - Num rng bins not read!\n"); exit(1); } 
      if (!SeekP(Gen->LBRFile,"RADAR A Start PRIID",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].StartPRI))
        { fprintf(Gen->msg,
	  "ERROR - Start PRI ID not read!\n"); exit(1); } 
      }
    else if (Chan[ii].Radar == 'B')
      {
      if (!SeekP(Gen->LBRFile,"RADAR B Sampling Rate",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].ADFreq))
	Warns++; 
      if (!SeekP(Gen->LBRFile,"RADAR B Sampling Bins",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].RngBins))
        { fprintf(Gen->msg,
	  "ERROR - Num rng bins not read!\n"); exit(1);} 
      if (!SeekP(Gen->LBRFile,"RADAR B Start PRIID",":",MaxO,St))
	Warns++;
      else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].StartPRI))
        { fprintf(Gen->msg,
	  "ERROR - Start PRI ID not read!\n"); exit(1);} 
      }
    else
      { fprintf(Gen->msg,
	   "ERROR - Radar %c unknown!\n",Chan[ii].Radar); exit(1);}	

    if (Chan[ii].RngBins < 115)
      {
      fprintf(Gen->msg,
	 "ERROR - current data format requires at least 115 range bins\n");
      fprintf(Gen->msg,
	 "Could be fixed in later software versions with zero padding\n");
      exit(1);
      }
    
    if (Warns > PrevWarns)
      fprintf(Gen->msg,
	 "WARNINGS - during TMC config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Parse HBR configuration info */
    PrevWarns = Warns;
    Fatals = 0;
    MaxO = 1100;
    if (Chan[ii].HBR == '1')   
      { if (!SeekP(Gen->LBRFile,"// HBR 1 Configuration","",-1,0))
	  Fatals++; }
    else if (Chan[ii].HBR == '2')   
      { if (!SeekP(Gen->LBRFile,"// HBR 2 Configuration","",-1,0))
	  Fatals++; }
    else if (Chan[ii].HBR == '3')   
      { if (!SeekP(Gen->LBRFile,"// HBR 3 Configuration","",-1,0))
	  Fatals++; }
    else if (Chan[ii].HBR == '4')   
      { if (!SeekP(Gen->LBRFile,"// HBR 4 Configuration","",-1,0))
	  Fatals++; }
    else if (Chan[ii].HBR == '5')   
      { if (!SeekP(Gen->LBRFile,"// HBR 5 Configuration","",-1,0))
	  Fatals++; }
    else if (Chan[ii].HBR == '6')   
      { if (!SeekP(Gen->LBRFile,"// HBR 6 Configuration","",-1,0))
	  Fatals++; }
    else
      { Fatals++; }

    if (Fatals != 0)
      { fprintf(Gen->msg,"ERROR - HBR config header not found!\n");
        exit(1); }
    else
      { St = ftell(Gen->LBRFile); }
   
    if (!SeekP(Gen->LBRFile,"Bits per I and Q",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%hd",&Chan[ii].Bits)) Warns++; 
    if (!SeekP(Gen->LBRFile,"Master PRF",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].MasterPRF)) Warns++;
    if (!SeekP(Gen->LBRFile,"Carrier Frequency",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].CarrierFreq)) Warns++;
    if (!SeekP(Gen->LBRFile,"Polarization",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Polarization,3)) Warns++;
    if (!SeekP(Gen->LBRFile,"Range Compression",":",MaxO,St)) Warns++;
    else if (!ReadChar(Gen->LBRFile, &Chan[ii].RngCompressedFlg)) Warns++;
    if (!SeekP(Gen->LBRFile,"Pulse Modulation",":",MaxO,St)) Warns++;
    else if (!ReadStr(Gen->LBRFile,Chan[ii].Modulation,STRING_SPACE)) Warns++;
    if (!SeekP(Gen->LBRFile,"Pulse Length",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].PulseLen)) Warns++;
    if (!SeekP(Gen->LBRFile,"Pulse Bandwidth",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].PulseBandwidth)) Warns++;
    if (!SeekP(Gen->LBRFile,"Pulse Power",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].PulsePower)) Warns++;
    if (!SeekP(Gen->LBRFile,"Antenna azimuth b",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].AzBeamwidth)) Warns++;
    if (!SeekP(Gen->LBRFile,"Antenna elevation b",":",MaxO,St)) Warns++ ;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].ElevBeamwidth)) Warns++;
    if (!SeekP(Gen->LBRFile,"Antenna depression a",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].DepressAngle)) Warns++;
    if (!SeekP(Gen->LBRFile,"Antenna squint",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].SquintAngle)) Warns++;
    if (!SeekP(Gen->LBRFile,"Calibration tone",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].CalToneLevel)) Warns++;
    if (!SeekP(Gen->LBRFile,"Calibration frequency",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].CalFreq)) Warns++;
    if (!SeekP(Gen->LBRFile,"System Losses",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].SysLosses)) Warns++;
    if (!SeekP(Gen->LBRFile,"Measured Delay0",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].Delay0Meas)) Warns++;
    if (!SeekP(Gen->LBRFile,"Calculated delay0",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].Delay0Calc)) Warns++;
    if (!SeekP(Gen->LBRFile,"Calculated STC en",":",MaxO,St)) Warns++;
    else if (!ReadChar(Gen->LBRFile,&Chan[ii].STCFlg)) Warns++;
    if (!SeekP(Gen->LBRFile,"STC/DGC Assigned ch",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%hd",&Chan[ii].DGCChannel)) Warns++;
    if (!SeekP(Gen->LBRFile,"Receiver noise figure",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%lf",&Chan[ii].RxNoiseFig)) Warns++;
    if (!SeekP(Gen->LBRFile,"DC offset",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%hd",&Chan[ii].DCOffset)) Warns++;
    if (!SeekP(Gen->LBRFile,"Presum Ratio",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%hd",&Chan[ii].PresumRatio)) Warns++;
    if (!SeekP(Gen->LBRFile,"Digital Gain",":",MaxO,St)) Warns++;
    else if (!fscanf(Gen->LBRFile,"%ld",&Chan[ii].DigGain)) Warns++;
    
    if (Warns > PrevWarns)
      fprintf(Gen->msg,
	 "WARNINGS - during HBR config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Find end PRIID */
    fseek(Gen->LBRFile,-500,2); /* Move close to end of file */
    St = ftell(Gen->LBRFile);
    MaxO = 510;
    Fatals = 0;
    if (Chan[ii].Radar == 'A')
      { if (!SeekP(Gen->LBRFile,"RADAR A Stop PRIID Number",":",MaxO,St))
	Fatals++; }
    else if (Chan[ii].Radar == 'B')
      { if (!SeekP(Gen->LBRFile,"RADAR B Stop PRIID Number",":",MaxO,St))
	Fatals++; }   
    if (!fscanf( Gen->LBRFile,"%ld",&Chan[ii].EndPRI)) Fatals++;  

    if (Fatals > 0)
      { fprintf(Gen->msg,
	   "ERROR - Radar stop PRI ID not read!\n"); exit(1); }

#if DEBUG
    Chan[ii].PRIs = 50;   /* Arb, for testing only */
    Chan[ii].EndPRI = Chan[ii].StartPRI+Chan[ii].PRIs-1;
#else
    Chan[ii].PRIs = (Chan[ii].EndPRI-Chan[ii].StartPRI) + 1; /* calc PRIs */
#endif
    
    /* Find IMU data for PRIID close to scene centre, if exists */
    if (!SeekP(Gen->LBRFile,"\nSync IDA","",-1,0))
      {
      fprintf(Gen->msg,"WARNING - no IMU/GPS data found!\n");
      Chan[ii].SceneCenUTC_hrs = DEFAULT_VAL_I;
      Chan[ii].SceneCenUTC_min = DEFAULT_VAL_I;
      Chan[ii].SceneCenUTC_sec = DEFAULT_VAL_I;;
      Chan[ii].SceneCenLat = DEFAULT_VAL_F;
      Chan[ii].SceneCenLong = DEFAULT_VAL_F;
      Chan[ii].SceneCenTrueHead = DEFAULT_VAL_F;
      }
    else
      {
      PRICen = Chan[ii].StartPRI +
	        (Int4B)(( (double)Chan[ii].EndPRI-
		          (double)Chan[ii].StartPRI )/2.0);

      if (Chan[ii].Radar == 'B')
	SeekP(Gen->LBRFile,"IDB","",50,-1); /* move for Radar B */
      fscanf(Gen->LBRFile,"%ld",&PRINext);  /* Read in PRI */          

      while ( PRINext < PRICen )
	{
	if (!SeekP(Gen->LBRFile,"\nSync IDA","",-1,-1))
	  { break; } 
	else     
	  {
          if (Chan[ii].Radar == 'B')
            SeekP(Gen->LBRFile,"IDB","",50,-1); /* move for Radar B */
          fscanf(Gen->LBRFile,"%ld",&PRINext);  /* Read in PRI */     
	  }

	} /* end while */

      /* Read in IMU data line (note fp not updated in while loop
	 if Sync not found). If Sync found, the line read here is the
	 one equal to or immediately after the centre PRI  */	
      SeekP(Gen->LBRFile,"\n","",-1,-1);
      ReadIMUMsg5(Gen,&IMUDat); /* Read and convert IMU record */

      /* Assign values */
      Chan[ii].SceneCenUTC_hrs = IMUDat.GPS_UTC_hrs;
      Chan[ii].SceneCenUTC_min = IMUDat.GPS_UTC_min;
      Chan[ii].SceneCenUTC_sec = IMUDat.GPS_UTC_sec;
      Chan[ii].SceneCenLat = IMUDat.GPSLat;
      Chan[ii].SceneCenLong = IMUDat.GPSLong; 
      Chan[ii].SceneCenTrueHead = IMUDat.INTHeadingTrue;
      
      fprintf(Gen->msg,
	 "Scene centre GPS UTC (hh:min:sec): %02ld:%02ld:%02ld\n",
	 Chan[ii].SceneCenUTC_hrs,Chan[ii].SceneCenUTC_min,
	 Chan[ii].SceneCenUTC_sec);
	      
      }  /* end else */
 
       
    /* Display total warnings */
    fprintf(Gen->msg,
       "Errors in parsing config info (channel %d) : %ld\n",ii,Warns);
 
    
    } /* End for ii (DataChannel) loop */

} /* end ParseLBRFile function */



/*********************************************/
void WriteVDF(struct CntrlStruct *Gen,
	      struct DataChannelStruct Chan[])
{
  FILE *VDFFile;
  char VDFFileName[MAX_FILE_NAME]="";
  char tmps[STRING_SPACE],TheDate[STRING_SPACE],TheTime[STRING_SPACE];
  unsigned char unchar;
  Int2B channel;
  Int4B tmp4b,Bytes=0;


  /* Construct VDF output file name */
  ConstructFileName("VDF",Gen->LogicalVol,1,VDFFileName);

  /* Open output file */
  if ( (VDFFile = fopen (VDFFileName, "wb") ) == NULL )
    { 
    fprintf (Gen->msg,"ERROR: Output file %s not opened\n",VDFFileName); 
    exit(0); 
    }

  /* Write VDF file entries */

  /* Write VDF Volume Descriptor Record */
  tmp4b = 1;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  unchar = 192;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 192;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(VDFFile,"%c",unchar);
  tmp4b = 360;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(VDFFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%-2s",tmps);
  sprintf(tmps,"CCB-CCT-0002");
  Bytes += fprintf(VDFFile,"%12s",tmps);
  sprintf(tmps,"E");                       /* field 10 */
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"A");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"%s",SOFTWARE_ID);
  Bytes += fprintf(VDFFile,"%-12s",tmps);
  sprintf(tmps,"%s",PHYS_VOL_ID);
  Bytes += fprintf(VDFFile,"%-16s",tmps);
  sprintf(tmps,"%s",LOGICAL_VOL_ID);
  Bytes += fprintf(VDFFile,"%-16s",tmps);
  sprintf(tmps,"%s",VOL_SET_ID);             /* field 15 */
  Bytes += fprintf(VDFFile,"%-16s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");                      /* field 20 */
  Bytes += fprintf(VDFFile,"%4s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%4s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%4s",tmps);

  FindDateTime(TheDate,TheTime);
  fprintf(Gen->msg,"CEOS transcription date: %s, time: %s\n",TheDate,TheTime);

  Bytes += fprintf(VDFFile,"%8s",TheDate);     /* field 23 */
  Bytes += fprintf(VDFFile,"%8s",TheTime);     /* field 24 */
  sprintf(tmps,"%s",CREATING_COUNTRY);         /* field 25 */
  Bytes += fprintf(VDFFile,"%-12s",tmps);
  sprintf(tmps,"%s",CREATING_AGENCY);
  Bytes += fprintf(VDFFile,"%-8s",tmps);
  sprintf(tmps,"%s",CREATING_FACILITY);
  Bytes += fprintf(VDFFile,"%-12s",tmps);
  sprintf(tmps,"%d",2*Gen->NumDataChannels);
  Bytes += fprintf(VDFFile,"%4s",tmps);
  sprintf(tmps,"%d",Gen->NumDataChannels+2);
  Bytes += fprintf(VDFFile,"%4s",tmps);
  sprintf(tmps,"1");                     /* field 30 */
  Bytes += fprintf(VDFFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%88s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%100s",tmps);

  /* Repeat for each data channel in logical volume */
  for (channel=0; channel<Gen->NumDataChannels; channel++)
  {

  /* Write VDF file pointer record to SAR leader file */
  tmp4b = 2 + 3*channel;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  unchar = 219;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 192;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(VDFFile,"%c",unchar);
  tmp4b = 360;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(VDFFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"%d",1+2*channel);
  Bytes += fprintf(VDFFile,"%4s",tmps);
  ConstructFileName("SAL",Gen->LogicalVol,channel+1,tmps);  /* field 10 */
  Bytes += fprintf(VDFFile,"%-16s",tmps);
  sprintf(tmps,"SARLEADER FILE");
  Bytes += fprintf(VDFFile,"%-28s",tmps);
  sprintf(tmps,"SARL");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  sprintf(tmps,"MIXED BINARY AND ASCII");
  Bytes += fprintf(VDFFile,"%-28s",tmps);
  sprintf(tmps,"MBAA");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  sprintf(tmps,"2");                           /* field 15 */
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"720");
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"1886");
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"VARIABLE LEN");
  Bytes += fprintf(VDFFile,"%-12s",tmps);
  sprintf(tmps,"VARE");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  sprintf(tmps,"1");                          /* field 20 */
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"2");
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%100s",tmps);
  sprintf(tmps," ");                         /* field 25 */
  Bytes += fprintf(VDFFile,"%100s",tmps);

  /*Write VDF file pointer record to imagery options file */
  tmp4b = 3 + 3*channel;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  unchar = 219;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 192;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(VDFFile,"%c",unchar);
  tmp4b = 360;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(VDFFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"%d",2+2*channel);
  Bytes += fprintf(VDFFile,"%4s",tmps);
  ConstructFileName("IMO",Gen->LogicalVol,channel+1,tmps);  /* field 10 */
  Bytes += fprintf(VDFFile,"%-16s",tmps);
  sprintf(tmps,"IMAGERY OPTIONS FILE");
  Bytes += fprintf(VDFFile,"%-28s",tmps);
  sprintf(tmps,"IMOP");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  sprintf(tmps,"MIXED BINARY AND ASCII");
  Bytes += fprintf(VDFFile,"%-28s",tmps);
  sprintf(tmps,"MBAA");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  tmp4b = Chan[channel].PRIs + 1; /*records*/
  sprintf(tmps,"%ld",tmp4b);      /* field 15 */
  Bytes += fprintf(VDFFile,"%8s",tmps);
  tmp4b = (Chan[channel].RngBins*2)+412; /* 1st record length */
  sprintf(tmps,"%ld",tmp4b);
  Bytes += fprintf(VDFFile,"%8s",tmps);  
  tmp4b = (Chan[channel].RngBins*2)+412; /* max record length */
  sprintf(tmps,"%ld",tmp4b);
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"FIXED LENGTH");
  Bytes += fprintf(VDFFile,"%-12s",tmps);
  sprintf(tmps,"FIXD");
  Bytes += fprintf(VDFFile,"%-4s",tmps);
  sprintf(tmps,"1");                          /* field 20 */
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps,"%ld",Chan[channel].PRIs);
  Bytes += fprintf(VDFFile,"%8s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%100s",tmps);
  sprintf(tmps," ");                         /* field 25 */
  Bytes += fprintf(VDFFile,"%100s",tmps);

  /* Write VDF text record */
  tmp4b = 4 + 3*channel;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  unchar = 18;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 63;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(VDFFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(VDFFile,"%c",unchar);
  tmp4b = 360;
  Bytes += BigEndWriteInt4B(VDFFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(VDFFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(VDFFile,"%2s",tmps);
  sprintf(tmps,"SASAR RAW SAR DATA");
  Bytes += fprintf(VDFFile,"%-40s",tmps);
  sprintf(tmps,"PRODUCED: RSA/UCT/RRSG. DATE: %s. TIME: %s.",
		TheDate,TheTime);         /* field 10 */
  Bytes += fprintf(VDFFile,"%-60s",tmps);
  sprintf(tmps,"CD-ROM ID: %s",PHYS_VOL_ID);
  Bytes += fprintf(VDFFile,"%-40s",tmps);

  sprintf(tmps,"Mission: ");           /* field 12 */
  strcat(tmps, Chan[channel].Mission);
  strcat(tmps, " Track: ");
  strcat(tmps, Chan[channel].Track); 
  Bytes += fprintf(VDFFile,"%-40s",tmps);

  sprintf(tmps,"%s",Chan[channel].Location);
  Bytes += fprintf(VDFFile,"%-40s",tmps);
  sprintf(tmps,"%s VHF %s",Chan[channel].Sensor,Chan[channel].Polarization);
  Bytes += fprintf(VDFFile,"%-20s",tmps);
  sprintf(tmps," ");                        /* field 15 */
  Bytes += fprintf(VDFFile,"%-104s",Chan[channel].GenComment);

  }  /* end channel for loop */

  fprintf(Gen->msg,"Bytes written of VDF File = %ld\n",Bytes);

  /* Tidy up */
  fclose(VDFFile);
  return;
}

/**********************************************************/
void WriteSARL(struct CntrlStruct *Gen,
	       struct DataChannelStruct Chan[])
{
  FILE *SALFile;
  unsigned char unchar;
  char SALFileName[MAX_FILE_NAME]="";
  char tmps[STRING_SPACE];
  Int4B tmp4b,Bytes=0,i;
  Int2B ii;

  /* Misc */
  ii = Gen->CurrentChannel;
  
  /* Construct SARL output file name */
  ConstructFileName("SAL",Gen->LogicalVol,Gen->FileNumOfType,SALFileName);

  /* Open output SAR Leader File */
  if ( (SALFile = fopen (SALFileName, "wb") ) == NULL )
    { 
    fprintf (Gen->msg,"ERROR: Output file %s not opened\n",SALFileName); 
    exit(0); 
    }

  /* Write SARL file entries */

  /* Write SARL file descriptor record */
  tmp4b = 1;
  Bytes += BigEndWriteInt4B(SALFile,tmp4b);
  unchar = 63;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 192;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(SALFile,"%c",unchar);
  tmp4b = 720;
  Bytes += BigEndWriteInt4B(SALFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(SALFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%-2s",tmps);
  sprintf(tmps,"CEOS-SAR-CCT");
  Bytes += fprintf(SALFile,"%12s",tmps);
  sprintf(tmps,"B");                        /* field 10 */
  Bytes += fprintf(SALFile,"%2s",tmps);
  sprintf(tmps,"B");
  Bytes += fprintf(SALFile,"%2s",tmps);   
  sprintf(tmps,SOFTWARE_ID);
  Bytes += fprintf(SALFile,"%-12s",tmps);
  sprintf(tmps,"%d",Gen->CEOSFileNum);
  Bytes += fprintf(SALFile,"%4s",tmps);     
  ConstructFileName("SAL",Gen->LogicalVol,Gen->FileNumOfType,tmps);
  Bytes += fprintf(SALFile,"%-16s",tmps);
  sprintf(tmps,"FSEQ");                      /* field 15 */
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(SALFile,"%8s",tmps);   
  sprintf(tmps,"4");
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps,"FTYP");
  Bytes += fprintf(SALFile,"%-4s",tmps);
  sprintf(tmps,"5");
  Bytes += fprintf(SALFile,"%8s",tmps);   
  sprintf(tmps,"4");                      /* field 20 */
  Bytes += fprintf(SALFile,"%4s",tmps);   
  sprintf(tmps,"FLGT");
  Bytes += fprintf(SALFile,"%-4s",tmps);
  sprintf(tmps,"9");
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"4");
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps," ");                      /* fields 24-27 */ 
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%64s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(SALFile,"%6s",tmps);
  sprintf(tmps,"1886");                   /* field 30 */
  Bytes += fprintf(SALFile,"%6s",tmps);  
  for (i=0;i<28;i++) {
     sprintf(tmps,"0");                    /* fields 31-58 */
     Bytes += fprintf(SALFile,"%6s",tmps);  
     }
  sprintf(tmps," ");                    /* fields 59-68 */
  Bytes += fprintf(SALFile,"%60s",tmps);
  sprintf(tmps,"0");
  Bytes += fprintf(SALFile,"%6s",tmps);
  sprintf(tmps,"0");                      /* field 70 */
  Bytes += fprintf(SALFile,"%6s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%288s",tmps);

  /* Write data set summary record */
  tmp4b = 2;
  Bytes += BigEndWriteInt4B(SALFile,tmp4b);
  unchar = 18;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 10;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 18;
  Bytes += fprintf(SALFile,"%c",unchar);
  unchar = 20;                             /* field 5 */
  Bytes += fprintf(SALFile,"%c",unchar);
  tmp4b = 1886;
  Bytes += BigEndWriteInt4B(SALFile,tmp4b);
  sprintf(tmps,"1");
  Bytes += fprintf(SALFile,"%4s",tmps);

  if (Chan[ii].Polarization == "HH") { sprintf(tmps,"00"); }
  else if (Chan[ii].Polarization == "HV") { sprintf(tmps,"01"); }
  else if (Chan[ii].Polarization == "VH") { sprintf(tmps,"10"); }
  else if (Chan[ii].Polarization == "VV") { sprintf(tmps,"11"); }
  else { sprintf(tmps," "); }
  Bytes += fprintf(SALFile,"%4s",tmps);

  sprintf(tmps,Chan[ii].Mission);
  Bytes += fprintf(SALFile,"%-16s",tmps);
  sprintf(tmps,"%s",Chan[ii].Location);            /* field 10 */
  Bytes += fprintf(SALFile,"%-32s",tmps);
  sprintf(tmps,"%10s/%2ld:%2ld:%2ld.000",Chan[ii].Date,
	 Chan[ii].SceneCenUTC_hrs,Chan[ii].SceneCenUTC_min,
	 Chan[ii].SceneCenUTC_sec);  
  Bytes += fprintf(SALFile,"%-32s",tmps);
  sprintf(tmps," ");  /* Spare */
  Bytes += fprintf(SALFile,"%-16s",tmps);
  sprintf(tmps,"%16s"," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16s"," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16s"," ");                 /* field 15 */
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");                                 /* fields 16-24 */
  Bytes += fprintf(SALFile,"%144s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].AveTerrainHeight);  /* field 25 */
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%8ld",
	  Chan[ii].StartPRI + (Int4B)(Chan[ii].PRIs/2));
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"%8ld",(Int4B)(Chan[ii].RngBins/2)); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  if (Chan[ii].MasterPRF > 0.0)
    sprintf(tmps,"%16.7f",
	 0.001*(Chan[ii].NomGroundSpeed/Chan[ii].MasterPRF)*Chan[ii].PRIs );
  else sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  if (Chan[ii].ADFreq > 0.0)
    sprintf(tmps,"%16.7f",
	   0.001*Chan[ii].RngBins*1.5e+08/Chan[ii].ADFreq);
  else sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");                                /* field 30 */
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"1"); 
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%4s",tmps);
  sprintf(tmps,"%16s",Chan[ii].Mission); 
  Bytes += fprintf(SALFile,"%-16s",tmps);
  sprintf(tmps,"%32s",Chan[ii].Sensor); 
  Bytes += fprintf(SALFile,"%-32s",tmps);
  sprintf(tmps,"%8s",Chan[ii].Track);     /* field 35 */
  Bytes += fprintf(SALFile,"%-8s",tmps);
  sprintf(tmps,"%8.3f",Chan[ii].SceneCenLat); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"%8.3f",Chan[ii].SceneCenLong); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"%8.3f",Chan[ii].SceneCenTrueHead);
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"-90"); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps," ");   /* field 40 */
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"%8.3f",0.001*Chan[ii].CarrierFreq); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"%16.7f",(C/(Chan[ii].CarrierFreq*1.0e+6)) ); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"00");
  Bytes += fprintf(SALFile,"%2s",tmps);
  sprintf(tmps,"UNMODULATED");
  Bytes += fprintf(SALFile,"%-16s",tmps);
  sprintf(tmps," ");                                /* fields 45-54 */
  Bytes += fprintf(SALFile,"%160s",tmps);
  sprintf(tmps," ");                               /* field 55 */
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%-8s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].ADFreq); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  if (Chan[ii].Delay0Meas > 0.0)
    sprintf(tmps,"%16.7f",Chan[ii].Delay0Meas);
  else
    sprintf(tmps,"%16.7f",Chan[ii].Delay0Calc);
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].PulseLen); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"YES");                              /* field 60 */ 
  Bytes += fprintf(SALFile,"%-4s",tmps);
  sprintf(tmps,"YES");
  Bytes += fprintf(SALFile,"%-4s",tmps);
  sprintf(tmps," ");  
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");  
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"8");
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps,"UNIFORM");                     /* field 65 */
  Bytes += fprintf(SALFile,"%-12s",tmps);
  sprintf(tmps," ");  
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");                                  /* field 70 */
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16.7f",90.0-Chan[ii].DepressAngle);
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%-4s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].MasterPRF); 
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].ElevBeamwidth); /* field 75 */
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%16.7f",Chan[ii].AzBeamwidth);  
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps," ");   
  Bytes += fprintf(SALFile,"%-32s",tmps);
  sprintf(tmps," "); 
  Bytes += fprintf(SALFile,"%8s",tmps);
  sprintf(tmps," ");                        /* field 80 */
  Bytes += fprintf(SALFile,"%-8s",tmps);
  sprintf(tmps," ");                        /* fields 81-125 */
  Bytes += fprintf(SALFile,"%-720s",tmps);

  sprintf(tmps,"%16.7f",Chan[ii].NomGroundSpeed);
  Bytes += fprintf(SALFile,"%16s",tmps);
  sprintf(tmps,"%8d",Chan[ii].PresumRatio);
  Bytes += fprintf(SALFile,"%8s",tmps);

  if (Chan[ii].PRICorrectionFlg == 'Y')     /* field 128 */
    { sprintf(tmps,"YES"); }
  else
    { sprintf(tmps,"NOT"); }
  Bytes += fprintf(SALFile,"%-4s",tmps);

  sprintf(tmps," ");
  Bytes += fprintf(SALFile,"%-92s",tmps);

  fprintf(Gen->msg,"Bytes written of SARL File = %ld\n",Bytes);

  /* Tidy up */
  fclose(SALFile);

  return;
}


/********************************************/
void WriteIMOP(struct CntrlStruct *Gen,
	       struct DataChannelStruct Chan[])
{
  FILE *IMOFile;
#if READSIM
  FILE *SimFile;
  unsigned char *SimData;
#endif
  unsigned char unchar;
  char IMOFileName[MAX_FILE_NAME]="";
  char tmps[STRING_SPACE];
  Int2B tmp2b,ii=0,TxPolar,RxPolar,RngCompress,PulseModulation=3;
  Int4B tmp4b,i,Record,CurrentIMUPRI,NextIMUPRI,PRI,DataRecordLen,
	SensorParamsUpdateFlg,PlatformPosnUpdateFlg,LookAngle,
        StartSlantRange,IMUValidFlg=0,Year,Month,Day,Delay0;
  double InstPRF,Bytes=0.0;	
  struct IMUStruct IMU;
  
  /* Misc */
#if READSIM
  SimData = (unsigned char *)malloc(sizeof(unsigned char)*2*Chan[ii].RngBins);
  if (SimData==NULL)
    { fprintf(Gen->msg,"ERROR - in SimData array allocation!\n"); exit(1); }
  if ((SimFile = fopen(SIMFILE,"rb")) == NULL)
    { fprintf(Gen->msg,"ERROR - sim file %s not opened!\n",SIMFILE); exit(1); }
  fprintf(Gen->msg,"Using raw data from sim file : %s\n",SIMFILE);
#endif 

  ii = Gen->CurrentChannel;

  DataRecordLen = 412 + 2*Chan[ii].RngBins;
  if (DataRecordLen<648) DataRecordLen = 648; /* pad out records if fewer 
                                                 than 118 rng bins */

  tmps[0] = Chan[ii].Date[0]; tmps[1] = Chan[ii].Date[1];
  tmps[2] = Chan[ii].Date[2]; tmps[3] = Chan[ii].Date[3];
  tmps[4] = '\0';
  sscanf(tmps,"%4ld",&Year);
  tmps[0] = Chan[ii].Date[6]; tmps[1] = Chan[ii].Date[7];
  tmps[2] = '\0';
  sscanf(tmps,"%2ld",&Month);
  tmps[0] = Chan[ii].Date[9]; tmps[1] = Chan[ii].Date[10];
  tmps[2] = '\0';
  sscanf(tmps,"%2ld",&Day);  

  if (Month == 2)
    { Day += 31; }
  else if (Month == 3)
    { Day += 59; }
  else if (Month == 4)
    { Day += 90; }  
  else if (Month == 5)
    { Day += 120; }    
  else if (Month == 6)
    { Day += 151; }
  else if (Month == 7)
    { Day += 181; }
  else if (Month == 8)
    { Day += 212; }
  else if (Month == 9)
    { Day += 243; }
  else if (Month == 10)
    { Day += 273; }
  else if (Month == 11)
    { Day += 304; }
  else if (Month == 12)
    { Day += 334; }
  if ( (Year % 4 == 0) && (Month > 2) )
    Day++;
  
  if (Chan[ii].Polarization[0] == 'H')
    { TxPolar = 0; }
  else
    { TxPolar = 1; }
  if (Chan[ii].Polarization[1] == 'H')
    { RxPolar = 0; }
  else
    { RxPolar = 1; }

  if (Chan[ii].RngCompressedFlg == 'Y')
    { RngCompress = 1; }
  else
    { RngCompress = 0; }

  if ( (strcmp(Chan[ii].Modulation,"MONO") == 0) ||
       (strcmp(Chan[ii].Modulation,"Mono") == 0) )		 
    PulseModulation = 2;
  else
    PulseModulation = 3;
    
  LookAngle = (Int4B)((90.0-Chan[ii].DepressAngle)*1.0e+06);

  if (Chan[ii].Delay0Meas > 0.0) /* val in struc is usec */
    {
    StartSlantRange =   /* need value in cm */
      (Int4B)((C/2.0)*Chan[ii].Delay0Meas*1.0e-4);
    Delay0 =            /* need value in nsec */
      (Int4B)(Chan[ii].Delay0Meas*1000.0);
    }
  else
    {
    StartSlantRange = 
      (Int4B)((C/2.0)*Chan[ii].Delay0Calc*1.0e-3);
    Delay0 = (Int4B)(Chan[ii].Delay0Calc*1000.0);
    }  
   
  /* Construct IMOP output file name */
  ConstructFileName("IMO",Gen->LogicalVol,Gen->FileNumOfType,IMOFileName);
  
  /* Open output file */  
  if ( (IMOFile = fopen (IMOFileName, "wb") ) == NULL )
    { 
    fprintf (Gen->msg,"ERROR: Output file %s not opened\n",IMOFileName); 
    exit(0); 
    }

  /* Read first IMU Record, if exists */
  if (!SeekP(Gen->LBRFile,"\nSync IDA","",-1,0))
    {
    IMUValidFlg = 0; 
    /* set to default value, if no valid data */
    IMU.GPS_UTC_hrs = DEFAULT_VAL_F;
    IMU.GPS_UTC_min = DEFAULT_VAL_F;
    IMU.GPS_UTC_sec = DEFAULT_VAL_F;    
    IMU.InertLat = DEFAULT_VAL_F;
    IMU.InertLong = DEFAULT_VAL_F;
    IMU.BaroInertHgt = DEFAULT_VAL_F;
    IMU.InertVelNorth = DEFAULT_VAL_F;
    IMU.InertVelEast = DEFAULT_VAL_F;
    IMU.InertVelVert = DEFAULT_VAL_F; 
    IMU.INTLat = DEFAULT_VAL_F;
    IMU.INTLong = DEFAULT_VAL_F;
    IMU.INTHgt = DEFAULT_VAL_F;
    IMU.INTVelNorth = DEFAULT_VAL_F;
    IMU.INTVelEast = DEFAULT_VAL_F;
    IMU.INTVelVert = DEFAULT_VAL_F;
    IMU.GPSLat = DEFAULT_VAL_F;
    IMU.GPSLong = DEFAULT_VAL_F;
    IMU.GPSHgt = DEFAULT_VAL_F;
    IMU.GPSVelNorth = DEFAULT_VAL_F;
    IMU.GPSVelEast = DEFAULT_VAL_F;
    IMU.GPSVelVert = DEFAULT_VAL_F;
    IMU.GPSRxStatus= DEFAULT_VAL_F;
    IMU.GPSEstHorizErr = DEFAULT_VAL_F;
    IMU.GPSEstVertErr = DEFAULT_VAL_F;
    IMU.InertRollAngle = DEFAULT_VAL_F;
    IMU.InertPitchAngle = DEFAULT_VAL_F;
    IMU.InertHeadingTrue = DEFAULT_VAL_F;
    IMU.INTRollAngle = DEFAULT_VAL_F;
    IMU.INTPitchAngle = DEFAULT_VAL_F;
    IMU.INTHeadingTrue = DEFAULT_VAL_F;;
    }  
  else /* valid IMU data exists */
    {
    IMUValidFlg = 1;

    /* Read in first IMU PRI ID, Inst PRI and data */
    if (Chan[ii].Radar == 'B')
      { SeekP(Gen->LBRFile,"IDB","",50,-1); } /* move slightly for Radar B */
    fscanf(Gen->LBRFile,"%ld",&CurrentIMUPRI);  /* Read in PRI ID*/        
    SeekP(Gen->LBRFile,"PRI","",50,-1);
    fscanf(Gen->LBRFile,"%lf",&InstPRF); /* Read in Inst PRI */
    if (InstPRF != 0.0) InstPRF = 1.0/InstPRF;  /* Convert to Inst PRF */
    SeekP(Gen->LBRFile,"\n","",-1,-1);
    ReadIMUMsg5(Gen,&IMU); /* Read and convert first IMU record */
#if EXTRACT_MOTION
    WriteMotion(Gen,&IMU,CurrentIMUPRI);
#endif    

    /* Read in second IMU PRI */
    if (!SeekP(Gen->LBRFile,"\nSync IDA","",500,-1))  
      { NextIMUPRI = Chan[ii].EndPRI + 2; }
    else
      {
      if (Chan[ii].Radar == 'B')
        SeekP(Gen->LBRFile,"IDB","",50,-1); /* move slightly for Radar B */
      fscanf(Gen->LBRFile,"%ld",&NextIMUPRI);  /* Read in next IMU PRI */ 
      }

    } /* end else valid IMU data exists */
  

  /* Write IMOP file descriptor record */
  tmp4b = 1;
  Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);    	 
  unchar = 63;
  Bytes += (double)fprintf(IMOFile,"%c",unchar);
  unchar = 192;
  Bytes += (double)fprintf(IMOFile,"%c",unchar); 
  unchar = 18; 
  Bytes += (double)fprintf(IMOFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += (double)fprintf(IMOFile,"%c",unchar);
  tmp4b = DataRecordLen; 
  Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += (double)fprintf(IMOFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-2s",tmps);
  sprintf(tmps,"CEOS-SAR-CCT");
  Bytes += (double)fprintf(IMOFile,"%12s",tmps);
  sprintf(tmps,"B");                        /* field 10 */
  Bytes += (double)fprintf(IMOFile,"%2s",tmps);
  sprintf(tmps,"B");
  Bytes += (double)fprintf(IMOFile,"%2s",tmps);
  sprintf(tmps,SOFTWARE_ID);
  Bytes += (double)fprintf(IMOFile,"%-12s",tmps);
  sprintf(tmps,"%d",Gen->CEOSFileNum);
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);     
  ConstructFileName("IMO",Gen->LogicalVol,Gen->FileNumOfType,tmps);
  Bytes += (double)fprintf(IMOFile,"%-16s",tmps);
  sprintf(tmps,"FSEQ");                      /* field 15 */
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"1");
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);   
  sprintf(tmps,"4");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"FTYP");
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"5");
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps,"4");                      /* field 20 */
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"FLGT");
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"9");
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps,"4");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps," ");                      /* fields 24-27 */ 
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);  
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-64s",tmps);
  sprintf(tmps,"%6ld",Chan[ii].PRIs);
  Bytes += (double)fprintf(IMOFile,"%6s",tmps);
  sprintf(tmps,"%6ld",DataRecordLen); /* field 30 */
  Bytes += (double)fprintf(IMOFile,"%6s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-24s",tmps);
  sprintf(tmps,"8");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"2");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"2");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"IQ");                             /* field 35 */
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"1");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"%8ld",Chan[ii].PRIs); 
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"%8ld",Chan[ii].RngBins); 
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);  
  sprintf(tmps," ");                               /* field 40 */
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"BSQ");
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"1");
  Bytes += (double)fprintf(IMOFile,"%2s",tmps);
  sprintf(tmps," ");                               /* field 45 */
  Bytes += (double)fprintf(IMOFile,"%2s",tmps);
  sprintf(tmps,"412");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps,"%8ld",2*Chan[ii].RngBins);
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps,"0");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps," ");                               /* field 50 */
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                    
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                   
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                    
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                   
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                              /* field 55 */ 
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-28s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                            
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                            
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps," ");                                /* field 60 */
  Bytes += (double)fprintf(IMOFile,"%-8s",tmps);
  sprintf(tmps,"COMPLEX UNSIGNED INTEGER");                            
  Bytes += (double)fprintf(IMOFile,"%-28s",tmps);
  sprintf(tmps,"CIU2");                            
  Bytes += (double)fprintf(IMOFile,"%-4s",tmps);
  sprintf(tmps,"0");
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);    
  sprintf(tmps,"0");                            
  Bytes += (double)fprintf(IMOFile,"%4s",tmps);  
  sprintf(tmps,"255");                           /* field 65 */
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps," ");
  Bytes += (double)fprintf(IMOFile,"%-192s",tmps);
  sprintf(tmps,"%8ld",(Int4B)Chan[ii].DCOffset);
  Bytes += (double)fprintf(IMOFile,"%8s",tmps);
  sprintf(tmps," ");
  for (i=648;i<DataRecordLen;i++)
    Bytes += (double)fprintf(IMOFile,"%1s",tmps);

  /* Write IMOP signal data records */
  PRI = Chan[ii].StartPRI;
  for (Record=2;Record<Chan[ii].PRIs+2;Record++)
    {
     
    tmp4b = Record;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    unchar = 50;
    Bytes += (double)fprintf(IMOFile,"%c",unchar);
    unchar = 10;
    Bytes += (double)fprintf(IMOFile,"%c",unchar);
    unchar = 18; 
    Bytes += (double)fprintf(IMOFile,"%c",unchar);
    unchar = 20;                             /* field 5 */
    Bytes += (double)fprintf(IMOFile,"%c",unchar);
    tmp4b = DataRecordLen;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = Record-2+Chan[ii].StartPRI;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = Chan[ii].RngBins;               /* field 10 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);    
    tmp4b = (Int4B)((DataRecordLen/2)-(206+Chan[ii].RngBins));
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);

    /* Note re IMU data (if IMU/GPS data exists):
       New IMU data is written if this data exists for a particular radar
       PRI. If new IMU data does not exist, the old data is carried over.
       At the start, the closest IMU data to the radar data is used
       (i.e. carried back) if no IMU data exists for the first few PRI. */
    while ((PRI >= NextIMUPRI) && (IMUValidFlg)) 
      {
      CurrentIMUPRI = NextIMUPRI;

      SeekP(Gen->LBRFile,"PRI","",50,-1);
      fscanf(Gen->LBRFile,"%lf",&InstPRF); /* Read in Inst PRI */
      if (InstPRF != 0.0) InstPRF = 1.0/InstPRF;  /* Convert to Inst PRF */

      SeekP(Gen->LBRFile,"\n","",-1,-1);
      ReadIMUMsg5(Gen,&IMU); /* Read and convert  IMU record */
#if EXTRACT_MOTION
      WriteMotion(Gen,&IMU,CurrentIMUPRI);
#endif        
      /* Read in next IMU PRI ID */
      if (!SeekP(Gen->LBRFile,"\nSync IDA","",500,-1))  
        { NextIMUPRI = Chan[ii].EndPRI + 2; }
      else
        {
        if (Chan[ii].Radar == 'B')
          SeekP(Gen->LBRFile,"IDB","",50,-1); /* move slightly for Radar B */
        fscanf(Gen->LBRFile,"%ld",&NextIMUPRI);  /* Read in next IMU PRI */ 
        }

      } /* end while loop */

    if (CurrentIMUPRI != PRI)
      {
      SensorParamsUpdateFlg = 0;
      PlatformPosnUpdateFlg = 0;
      }
    else
      {
      SensorParamsUpdateFlg = 1;
      PlatformPosnUpdateFlg = 1;
      }

#if FILL_IMO    
    Bytes += (double)BigEndWriteInt4B(IMOFile,SensorParamsUpdateFlg);             
    Bytes += (double)BigEndWriteInt4B(IMOFile,Year);  
    Bytes += (double)BigEndWriteInt4B(IMOFile,Day);
    tmp4b = IMU.GPS_UTC_hrs*3600000+IMU.GPS_UTC_min*60000+
            IMU.GPS_UTC_sec*1000 ;          /* field 15 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp2b = 0;            /* tmp, for now */
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);  
    tmp2b = 0;            /* tmp, for now */
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
    tmp2b = TxPolar;
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
    tmp2b = RxPolar;   
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
    tmp4b = (Int4B)(InstPRF*1000);             /* field 20 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = CurrentIMUPRI;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp2b = RngCompress;
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
    tmp2b = PulseModulation;
    Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
    tmp4b = (Int4B)(Chan[ii].PulseLen*1000.0);  
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = 0;                                 /* field 25 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0; 
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = 0;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(Chan[ii].CarrierFreq*1000.0);            
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = 0;            
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b); 
    tmp4b = 0;                               /* field 30 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = 0;   
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;          
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = LookAngle;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = 0;                             /* field 35 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = StartSlantRange;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);  
    tmp4b = Delay0;       
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0; 
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    Bytes += (double)BigEndWriteInt4B(IMOFile,PlatformPosnUpdateFlg);    
    tmp4b = (Int4B)(IMU.INTLat*1.0e+06);   /* field 40 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTLong*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTHgt*1.0e+03);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b); 
    tmp4b = (Int4B)(IMU.INTVelNorth*1.0e+03); /* field 44 (1) */       
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTVelEast*1.0e+03);  /* field 44 (2) */       
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTVelVert*1.0e+03);  /* field 44 (3) */       
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    for (i=0;i<3;i++) {
      tmp4b = 0;                            /* field 45 */
      Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
      }  
    tmp4b = 0;            
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTHeadingTrue*1.0e+06);;          
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTPitchAngle*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.INTRollAngle*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = 0;                                /* field 50 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);

    tmp4b = (Int4B)(IMU.InertLat*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);    
    tmp4b = (Int4B)(IMU.InertLong*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.InertVelNorth*1.0e+03); /* field 53 (1) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.InertVelEast*1.0e+03); /* field 53 (2) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.InertVelVert*1.0e+03); /* field 53 (3) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   

    tmp4b = (Int4B)(IMU.InertRollAngle*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.InertPitchAngle*1.0e+06); /* field 55 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.InertHeadingTrue*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);

    tmp4b = (Int4B)(IMU.BaroInertHgt*1.0e+03);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
 
    tmp4b = (Int4B)(IMU.GPSLat*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);    
    tmp4b = (Int4B)(IMU.GPSLong*1.0e+06);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.GPSHgt*1.0e+03);    /* field 60 */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b); 

    tmp4b = (Int4B)(IMU.GPSVelNorth*1.0e+03); /* field 61 (1) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.GPSVelEast*1.0e+03); /* field 61 (2) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   
    tmp4b = (Int4B)(IMU.GPSVelVert*1.0e+03); /* field 61 (3) */
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);   

    tmp4b = (Int4B)(IMU.GPSRxStatus);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.GPSEstHorizErr);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
    tmp4b = (Int4B)(IMU.GPSEstVertErr);
    Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);

    for (i=0;i<37;i++) 
      {
      tmp4b = 0;                             /* field 65 */ 
      Bytes += (double)BigEndWriteInt4B(IMOFile,tmp4b);
      }
#endif

    /* Write radar data */
#if READSIM
    fread(SimData,sizeof(unsigned char),2*Chan[ii].RngBins,SimFile);
    Bytes += fwrite(SimData,sizeof(unsigned char),2*Chan[ii].RngBins,IMOFile);
#elseif FILL_IMO
    for (i=0;i<Chan[ii].RngBins/2;i++)
      {
      tmp2b = i%255;            /* tmp, for now */
      Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
      tmp2b = -(i%255);            /* tmp, for now */
      Bytes += (double)BigEndWriteInt2B(IMOFile,tmp2b);
      }
#endif

#if FILL_IMO   
    /* Pad out record if fewer than 118 rng bins */
    sprintf(tmps," ");
    for (i=(412+2*Chan[ii].RngBins);i<DataRecordLen;i++)
      Bytes += (double)fprintf(IMOFile,"%1s",tmps);
#endif
    PRI++;  /* increment PRI ID counter */
    
    }  /* end for Record loop */

  fprintf(Gen->msg,"Bytes written of IMOP File = %.0f (%.1f MB)\n",
        Bytes,Bytes/1048576.0);

  /* Tidy up */
  fclose(IMOFile);
#if READSIM
  fclose(SimFile);
  free(SimData);
#endif
  return;
}

/*******************************************/
void WriteNVDF(struct CntrlStruct *Gen)
{
  FILE *NVDFile;
  unsigned char unchar;
  char NVDFileName[MAX_FILE_NAME]="";
  char tmps[STRING_SPACE];
  Int4B Bytes=0,tmp4b;

  /* Construct NVDF output file name */
  ConstructFileName("NVD",Gen->LogicalVol,1,NVDFileName);

  /* Open output file */
  if ( (NVDFile = fopen (NVDFileName, "wb") ) == NULL )
    { 
    fprintf (Gen->msg,"ERROR: Output file %s not opened\n",NVDFileName); 
    exit(1); 
    }

  /* Write VDF Volume Descriptor Record */
  tmp4b = 1;
  Bytes += BigEndWriteInt4B(NVDFile,tmp4b);
  unchar = 192;
  Bytes += fprintf(NVDFile,"%c",unchar);
  unchar = 192;
  Bytes += fprintf(NVDFile,"%c",unchar);
  unchar = 63;
  Bytes += fprintf(NVDFile,"%c",unchar);
  unchar = 18;                             /* field 5 */
  Bytes += fprintf(NVDFile,"%c",unchar);
  tmp4b = 360;
  Bytes += BigEndWriteInt4B(NVDFile,tmp4b);
  sprintf(tmps,"A");
  Bytes += fprintf(NVDFile,"%-2s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(NVDFile,"%-2s",tmps);
  sprintf(tmps,"CCB-CCT-0002");
  Bytes += fprintf(NVDFile,"%12s",tmps);
  sprintf(tmps,"E");             /* field 10 */
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps,"A");
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps,SOFTWARE_ID);
  Bytes += fprintf(NVDFile,"%-12s",tmps);   
  sprintf(tmps,PHYS_VOL_ID); 
  Bytes += fprintf(NVDFile,"%-16s",tmps);     
  sprintf(tmps,LOGICAL_VOL_ID); 
  Bytes += fprintf(NVDFile,"%-16s",tmps);
  sprintf(tmps,VOL_SET_ID);     /* field 15 */
  Bytes += fprintf(NVDFile,"%-16s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps,"1");
  Bytes += fprintf(NVDFile,"%2s",tmps);
  sprintf(tmps," ");           /* field 20 */
  Bytes += fprintf(NVDFile,"%4s",tmps);
  sprintf(tmps,"2");
  Bytes += fprintf(NVDFile,"%4s",tmps);
  sprintf(tmps,"2");
  Bytes += fprintf(NVDFile,"%4s",tmps);
  sprintf(tmps," ");          /* fields 23-30 */
  Bytes += fprintf(NVDFile,"%60s",tmps);
  sprintf(tmps," ");
  Bytes += fprintf(NVDFile,"%88s",tmps);
  sprintf(tmps," ");         /* field 32 */
  Bytes += fprintf(NVDFile,"%100s",tmps);

  fprintf(Gen->msg,"Bytes written of NVD File = %ld\n",Bytes);

  /* Tidy up */
  fclose(NVDFile);
  return;
}


/*****************************************************/
void main(int argc, char *argv[])
{
  Int2B DataChannel=0;
  Int4B TimeStart,TimeEnd;
  double DiskSpace=0;
  struct CntrlStruct Gen; /* current channel and general info */
  struct DataChannelStruct Chan[MAX_DATA_CHANNELS];
 

  /* Assign message output */
  Gen.msg = stdout;

  fprintf(Gen.msg,"\nProg: G2 CEOS Formatter (Ver. %s)\n",PROG_VERSION);
  fprintf(Gen.msg,"Code: J.M. Horrell (Copyright UCT)\n");
  fprintf(Gen.msg,"[Writes native SASAR output to SASAR CEOS Format]\n\n");
  
  if (argc < 2)
    {fprintf(Gen.msg,"Usage: ceos [lbr_file_name]\n"); exit(1);}

  /* Messages about compile-time options */
#if DEC_CC
  fprintf(Gen.msg,"[Code compiled with DEC CC compiler options]\n");
#elif GNUC
  fprintf(Gen.msg,"[Code compiled with GNU C compiler options]\n");
#elif BC31
  fprintf(Gen.msg,"[Code compiled with Borland C 3.1 compiler options]\n");
#else
  fprintf(Gen.msg,"\nERROR - no compiler specified at compile time!!\n");
#endif

  /* Open LBR file */
  if ( (Gen.LBRFile = fopen (argv[1], "rb") ) == NULL )
    {
    fprintf (Gen.msg,"ERROR: Input LBR file %s not opened\n",argv[1]);
    exit(1);
    }

  /* Initialise */
  Gen.FileNumOfType = 1;
  Gen.CEOSFileNum = 1;
  Gen.LogicalVol = LOGICAL_VOL;
  Gen.NumDataChannels = 0;

  /* Checks */
  if (Gen.LogicalVol > 99 || Gen.LogicalVol < 1) {
    fprintf(Gen.msg,"ERROR - Logical volume %d out of range (1-99)\n",
            Gen.LogicalVol);
    exit(1); }


  /* Read Radars and HBRs used */
  if (!SeekP(Gen.LBRFile,"// General Configuration","",-1,0))
    {
    fprintf(Gen.msg,"ERROR at start - General config header not found!\n");
    exit(1);
    }
 
  if (!SeekP(Gen.LBRFile,"Systems enabled",":",200,-1))
    {
    fprintf(Gen.msg,"ERROR - 'Systems enabled' parameter not found!\n");
    exit(1);
    }
    
  DataChannel = 0;
  while (DataChannel < MAX_DATA_CHANNELS)
    {
    if (SeekP(Gen.LBRFile,"RADAR_","",100,-1))
      {
      Chan[DataChannel].Radar = getc(Gen.LBRFile);      
      if (SeekP(Gen.LBRFile,"_HBR","",100,-1))
	{
        Chan[DataChannel].HBR = getc(Gen.LBRFile);
        Gen.NumDataChannels++;
        }  /* end if FindStr _HBR */
      else
	DataChannel = MAX_DATA_CHANNELS;
      }  /* end if FindStr RADAR_ */
    else
      DataChannel = MAX_DATA_CHANNELS;
    DataChannel++;
    }  /* end while DataChannel */
    
  if (Gen.NumDataChannels==0)
    { fprintf(Gen.msg,"ERROR - No valid data channels!\n"); exit(1);}
  else
    fprintf(Gen.msg,"Number of data channels = %d\n",Gen.NumDataChannels);
 

  /* Parse LBR File for relevant radar and HBR data */
  ParseLBRFile(&Gen,Chan);

  /* Calc disk space needed for transcription (incl. LBR file) */
  DiskSpace = 0.0;
  fseek(Gen.LBRFile,0,2);
  DiskSpace += ftell(Gen.LBRFile);                /* LBR output file */
  fseek(Gen.LBRFile,0,0); /* reset pointer, just in case */
  DiskSpace += 360.0 + (double)Gen.NumDataChannels*(3.0*360.0); /* VDF file */ 
  DiskSpace += (double)Gen.NumDataChannels*(720.0+1886.0);    /* SARL file */
  for (DataChannel=0;DataChannel<Gen.NumDataChannels;DataChannel++)
    {
    DiskSpace += (((double)Chan[DataChannel].EndPRI-
		   (double)Chan[DataChannel].StartPRI)+2.0)*
                   (412.0+(double)Chan[DataChannel].RngBins*2.0);
    fprintf(Gen.msg,"Data channel %d : PRIS %ld - %ld, RngBins %ld\n",
	DataChannel, Chan[DataChannel].StartPRI,
	    Chan[DataChannel].EndPRI, Chan[DataChannel].RngBins );     

    }
  DiskSpace += 360.0;                               /* NVDF file */
  fprintf(Gen.msg,
     "\nTotal space required for output : %.0f bytes (%.1f MB)\n",
     DiskSpace,(DiskSpace/1048576.0) ); 

  /* Get start time */  
  TimeStart = (Int4B)time(NULL);

  /* Call main functions to write files */
  WriteVDF(&Gen,Chan);      /* VDF File */
  for (DataChannel=0;DataChannel<Gen.NumDataChannels;DataChannel++)
    {
    Gen.CurrentChannel = DataChannel;
    fprintf(Gen.msg,"Writing data channel %d: radar = %c, HBR = %c\n",
          DataChannel,Chan[DataChannel].Radar,Chan[DataChannel].HBR);
    WriteSARL(&Gen,Chan);   /* SARL File */
    Gen.CEOSFileNum++;   /* incr */
    WriteIMOP(&Gen,Chan);   /* IMOP File */
    Gen.CEOSFileNum++;   /* incr */
    Gen.FileNumOfType++; /* incr */
    } /* end for DataChannel loop */
  WriteNVDF(&Gen);    /* NVDF File */

  /* Display time taken for transcription */
  TimeEnd = (Int4B)time(NULL);
  fprintf(Gen.msg,"Time for transcription : %ld sec (%.2f min)\n",
  TimeEnd-TimeStart,(double)((TimeEnd-TimeStart)/60.0));

  /* Tidy up */
  fclose(Gen.LBRFile);

}  /* end main program */
