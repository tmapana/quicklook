/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2parse.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2parse.c. 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Add this header info. Convert to Hz, sec, etc.

VERSION/AUTHOR/DATE : 0.2 /jmh / 1999-07-22
COMMENTS:
(IMUStruct) GPS_UTC_sec changed to double from Int4B. Changed and 
added WriteMotion prototype.

=========================================*/


#define COMMENT_SPACE 200
#define MAX_DATA_CHANNELS 6
#define BYTES_PER_SECTOR 512

/* Stores misc control info for transcription */
struct CntrlStruct {
  Int2B NumDataChannels;
  Int2B CurrentChannel;
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
  double ADFreq;            /*(Hz)*/   
  Int4B RngBins;
  /* HBR config info */
  Int2B Bits;
  double MasterPRF;         /*(Hz)*/
  double CarrierFreq;       /*(Hz)*/
  char Polarization[3];
  char RngCompressedFlg;
  char Modulation[STRING_SPACE];
  double PulseLen;          /*(sec)*/
  double PulseBandwidth;    /*(Hz)*/
  double PulsePower;        /*(W)*/
  double AzBeamwidth;       /*(deg)*/
  double ElevBeamwidth;     /*(deg)*/    
  double DepressAngle;      /*(deg)*/ 
  double SquintAngle;       /*(deg)*/
  double CalToneLevel;
  double CalFreq;           /*(MHz)*/
  double SysLosses;         /*(dB)*/
  double Delay0Meas;        /*(sec)*/
  double Delay0Calc;        /*(sec)*/
  char STCFlg;
  Int2B DGCChannel;
  double RxNoiseFig;        /*(dB)*/
  Int2B DCOffset;
  Int2B PresumRatio;
  Int4B DigGain;
  char HBRComment[COMMENT_SPACE];
  /* Misc */
  Int4B PRIIDIncr; /* increment of the PRI ID per pulse */
  Int4B StartPRI;
  Int4B EndPRI;
  Int4B PRIs;
  Int4B SceneCenUTC_hrs;    /*(hrs)*/
  Int4B SceneCenUTC_min;    /*(min)*/
  double SceneCenUTC_sec;    /*(sec)*/
  double SceneCenLat;       /*(deg)*/
  double SceneCenLong;      /*(deg)*/
  double SceneCenHgt;       /*(m)*/
  double SceneCenTrueHead;  /*(deg)*/
};

/* Structure for IMU/GPS raw binary data */
struct Msg5Struct {
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
struct IMUStruct {
  Int4B GPS_UTC_hrs;     /*(hrs)*/
  Int4B GPS_UTC_min;     /*(min)*/
  double GPS_UTC_sec;     /*(sec)*/ 
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


/* Prototypes for functions in g2parse.c file */
void ParseLBRFile(FILE *msg,
                  FILE *LBRFile,
                  struct CntrlStruct *Gen,
		          struct DataChannelStruct Chan[]);

void ParseDataChannels(FILE *msg,
                       FILE *LBRFile,
                       struct CntrlStruct *Gen,
		               struct DataChannelStruct Chan[]);		  

void ReadIMUMsg5(FILE *LBRFile,
	             struct IMUStruct *IMU);

void WriteMotionHeader(FILE *OutF);

void WriteMotion(FILE *OutF, struct IMUStruct IMU,Int4B PRI);
