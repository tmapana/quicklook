/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2azc.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2azc.c 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Updated with this header info.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 1999-09-16
COMMENTS:
Add LogFileName to struct.

=========================================*/


struct AzcFileStruct
  {
  FILE *msg;
  char Version[STRING_SPACE];
  Int4B ScreenUpdateRate;
  char LogFileName[STRING_SPACE];
  double InputStartSampleDelay;
  double CarrierFreq;
  double InputPRF;
  double NomGroundSpeed;
  Int4B InputFileAzPts;
  Int4B StartProcessAzPt;
  Int4B AzPtsToProcess;
  Int4B InputFileRngBins;
  Int4B StartProcessRngBin;
  Int4B RngBinsToProcess;
  double InputDCOffsetI; /* new */
  double InputDCOffsetQ; /* new */
  Int4B FFTSize;
  Int4B InvFFTSizeReduc;
  char InputFileName[STRING_SPACE];
  char OutputFileName[STRING_SPACE];
  char AppendExistOutFileFlg;  
  Int4B RngFocSegments;
  Int4B RefFuncSign;
  double A2DFreq;
  double NomAzRes; 
  double WinConstTime;
  Int4B NumLooks;
  double LookOverlapFrac;
  double WinConstFreq;
  Int4B RngCurvInterpSize;
  Int4B RngCurvBatchSize;
  Int4B PostSumRatio;
  Int4B DetectMethod;
  Int4B InputDataType; /* new */
  Int4B OutputDataType; 
  double Scale;
  Int4B ReportMax;
  };


/* Prototypes for functions in g2azc.c file */

Int2B RngCurvStraighten(double   A2DFreq,
		      float    **CorrectedBatchFreq,
		      Int4B FFTSize,
		      Int4B MasterRefFreqPtsD2,
		      Int4B *NextRngBinToStraighten,
		      double   RngBinSize,
		      Int4B RngCurvBatchSize,
		      double   RngCurvCalc,
		      float    **RngCurvedBatchFreq,
		      Int4B RngCurvInterpSize,
		      double   *RngCurvRefRng,
		      Int4B *StraightenedRngBins,
		      Int4B *StraightenedRngBinsUsed
		      );
		      
Int2B G2Azc (struct AzcFileStruct Cmd);
