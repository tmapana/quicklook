/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2rnc.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for range compression.

VERSION/AUTHOR/DATE : 1999-01-28 / Jasper Horrell / 1999-01-28
COMMENTS: 
Incorporates intereference suppression code.

VERSION/AUTHOR/DATE : 0.2 / Jasper Horrell / 1999-09-14
COMMENTS:
Changed update rates for interference suppression to Int4B. 

VERSION/AUTHOR/DATE : 1.2 / Jasper Horrell / 1999-09-16
COMMENTS:
Add LogFileName to struct.

=========================================*/

struct RncFileStruct
  {
  /* general */
  FILE   *msg;
  char   Version[STRING_SPACE];
  Int4B  ScreenUpdateRate;
  char   LogFileName[STRING_SPACE];
  char   InputFileName[STRING_SPACE];
  char   OutputFileName[STRING_SPACE];
  Int4B  InputDataType;
  Int4B  OutputDataType;
  Int4B  StartProcessPRI; 
  Int4B  PreSumRatio;
  Int4B  PreSummedPulsesToUse;
  Int4B  InputFileRngBins;
  Int4B  StartRngBin;
  Int4B  RngBinsToProcess; 
  Int4B  HeaderBytes;
  Int4B  FooterBytes;
  double InputDCOffsetI;
  double InputDCOffsetQ;
  double InputIQRatio;
  /* misc (used throughout) */
  double CarrierFreq;        /* SRC, range walk, mocomp calcs */
  double A2DFreq;
  Int4B  RngFFTSize;
  Int4B  RngShiftInterpSize; /* mocomp rng shift, rng walk rng shift */
  double Scale;
  /* range compression specific */
  char   RngComFlg;
  Int4B  RngComRefFuncPhaseSign;  
  double RngComChirpBandwidth;
  double RngComPulseLen;
  double RngComWinConstTime;
  /* motion comp specific */
  char   MoCompFlg;
  char   MoCompFileName[STRING_SPACE];
  char   MoCompRngShiftFlg;
  Int4B  MoCompRngShiftSign;
  Int4B  MoCompRngShiftIndex;
  Int4B  MoCompPhaseSign;
  Int4B  MoCompRngUpdates;
  Int4B  MoCompStartPRI;
  /* LMS interference suppression */ /* Richard */
  char   LmsFlg;
  Int4B  LmsUpdateRate;
  Int4B  LmsNumWeights;
  Int4B  LmsSidelobeOrder;
  /* Notch interference suppression */ /* Richard */
  char   NotchFlg;
  Int4B  NotchUpdateRate;
  Int4B  NotchNumFFTLines;
  double NotchCutoff;
  Int4B  NotchMedianKernLen;
  /* STC specific */
  char   STCFlg;
  char   STCFileName[STRING_SPACE];
  /* SRC, Doppler centroid, and range walk specific */
  char   SRCFlg;            /* requires range compress */
  double DopCentroid;
  char   RngWalkRngShiftFlg;
  char   RngWalkPhaseShiftFlg;
  double RngWalkAzBeamwidth; /* range walk calcs */
  double SRCFocusRng;        /* SRC */
  double NomGroundSpeed;     /* SRC and range walk calcs */
  double SquintAngle;        /* SRC and range walk calcs */
  double InputPRF;           /* rng walk and DopCentroid calcs */
  };

/* Header file with rng compression function prototype contained in 'g2rnc.c'*/
Int2B G2Rnc (struct RncFileStruct Cmd);
