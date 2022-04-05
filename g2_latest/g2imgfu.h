/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2imgfu.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2imgfu.c. 

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1999-01-26
COMMENTS: 
Add this header info.

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/


/* Structs for AIRSAR format data */
struct AIRSARHead1Struct
  {
  Int4B RecordLen;
  Int4B NumHeaderRecords;
  Int4B SamplesPerRecord;
  Int4B LinesInImage;
  Int4B NumBytesPerSample;
  char SASARProcVersion[STRING_SPACE];
  char DataType[STRING_SPACE];
  char RangeProjection[STRING_SPACE];
  double RangePixelSpacing;
  double AzPixelSpacing;
  Int4B ByteOffsetOldHeader;
  Int4B ByteOffsetUserHeader;
  Int4B ByteOffset1stDataRecord;
  Int4B ByteOffsetParamHeader;
  char DataLineFormat[STRING_SPACE];
  Int4B ByteOffsetCalHeader;    
  Int4B ByteOffsetDEMHeader;
  };

struct AIRSARParamHeadStruct
  {
  char HeaderName[STRING_SPACE];
  char SiteName[STRING_SPACE];
  double SiteLatitude;
  double SiteLongitude;
  char ImageTitle[STRING_SPACE];
  Int4B CDROMId;
  char Freq[STRING_SPACE];
  char Polarization[STRING_SPACE];
  char DataType[STRING_SPACE];
  Int4B DataId;
  Int4B ArchivalFlg;
  Int4B CDROMStartPRIId;
  Int4B ProcStartPRIId;
  double LatitudeStartScene;
  double LongitudeStartScene;
  double LatitudeEndScene;
  double LongitudeEndScene;
  Int4B ApproxStartHDDTFootage;
  char AcquisitionDate[STRING_SPACE];
  Int4B AcquisitionDay;
  double AcquisitionSec;  
  Int4B RecordWindowDuration;
  char FreqsCollected[STRING_SPACE];
  double DigitalDelay;
  double ChirpDelay;
  Int4B ProcDelay;
  double PRFStartTransfer;
  double SamplingRate;
  double CentreFreq;
  double ChirpBandwidth;
  char TypeOfChirp[STRING_SPACE];
  double PulseLen;
  double ProcWavelen;
  double BaroAltitude;
  double RadarAltimeterAlt;
  double ProcAltitude;
  double ElevInvestigatorSite;
  double AircraftTrackAngle;
  double AircraftYawAngle;
  double AircraftPitchAngle;
  double AircraftRollAngle;
  double ProcYawAngle;
  double ProcPitchAngle;
  double ProcRollAngle;
  double NomPRFRatioHzPerKnot;
  double NomPRFRatioPerMetre;
  double PRFCorrectionFactor;
  Int4B RngFFTSize;
  Int4B AzFFTSize;
  Int4B FrameSize;
  Int4B NumFramesProcessed;
  double RngAlignDelayHH;
  double RngAlignDelayHV;
  double RngAlignDelayVH;
  double RngAlignDelayVV;
  double NearSlantRng;
  double FarSlantRng;
  double NearLookAngle;
  double FarLookAngle;
  double NumProcAzLooks;
  double NumProcRngLooks;
  char RngWeighting[STRING_SPACE];
  double RngWeightingCoeff;
  char AzWeighting[STRING_SPACE];
  double AzWeightingCoeff;
  double PercentPRFBWProc;
  Int4B DeskewFlg;
  double SlantRngSampleSpacing;
  double NomSlantRngRes;
  double AzSampleSpacing;
  double NomAzRes;
  Int4B NumInterpPtsRMC;
  Int4B AzRefSizePerLookNear;
  Int4B AzRefSizePerLookFar;
  double ImageCentreLat;
  double ImageCentreLong;
  double CalToneVideoFreq;
  double CalTonePowerHH;
  double CalTonePowerHV;
  double CalTonePowerVH;
  double CalTonePowerVV;
  double CalibFactorHH;
  double CalibFactorHV;
  double CalibFactorVH;
  double CalibFactorVV;
  double MeasCorrHVVHPowerRatio;
  double MeasCorrHVVHPhase;
  double CalTonePhaseMeasHH;
  double CalTonePhaseMeasHV;
  double CalTonePhaseMeasVH;
  double CalTonePhaseMeasVV;
  double GenScaleFactor;
  double GPSAlt;  
  };  


/* Prototypes for functions in g2imgfu.c file */
Int2B UnBlockAzcFile(FILE *msg,
                     char InFileName[],
                     char OutFileName[],
                     Int4B DetectMethod,
                     Int4B DataType,
                     Int4B RngBins,
                     Int4B PRIsPerAzBlock,
                     Int4B NumAzBlocks );


Int2B ProjImageToGrndRng(FILE *msg,
			 char AzcFileName[],
                         char GRngFileName[],
                         Int4B AzPixels,
                         Int4B RngBins,
                         Int4B DataType,
			 double FirstRngBinSlDist,     /* (m) */
			 double SlRngBinSpacing,       /* (m) */
                         double *RetGrndRngBinSpacing, /* (m) */
			 double Height,                /* (m) */
                         Int4B InterpPts ); 

Int2B WriteSUNRasterFile(FILE *msg,
			 char InFileName[],
			 char RasFileName[],
			 Int4B DataType,
			 Int4B DetectMethod,
			 Int4B Width,
			 Int4B Height);

Int2B WriteAIRSARFormat(FILE *msg,
			char InFileName[],
			char AIRFileName[],
                        Int4B DataType,
			Int4B DetectMethod,
			Int4B WidthPixels,
			Int4B HeightPixels,
		       	struct AIRSARHead1Struct *Head1,
			struct AIRSARParamHeadStruct *Param);
