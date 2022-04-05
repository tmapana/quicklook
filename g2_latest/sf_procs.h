/*===================================================================
COPYRIGHT: UCT Radar Remote Sensing Group 2000
FILE NAME: sf_procs.h
CODE CONTROLLER: Richard Lord

DESCRIPTION: 
Header file for sf_procs.c.

Function prototypes:
G2AddSubSpectrum
G2FFTRangeLine
G2GetRangeLine
G2GetPointSim
G2PhaseComp
G2PutRangeLine
G2RangeCompress
G2GetRangeProfile


VERSION/AUTHOR/DATE : 0.1 / Richard Lord / 2000-03-31
COMMENTS: 

=====================================================================*/


int G2AddSubSpectrum(
              FILE *msg,
              Int4B RngFFTSize,
              Int4B NewRngBinsToProcess,
              Int4B step,
              double fshift,
              double rline[],
              double newrline[]
              );

int G2FFTRangeLine(
              FILE *msg,
              Int4B RngFFTSize,
              double InputA2DFreq,
              double InputStartSampleDelay,
              double rline[]
              );

int G2GetRangeLine( 
              FILE *msg,
              FILE *infile,
              Int4B RngBinsToProcess,
              Int4B RngFFTSize,
              Int4B InputDataType,
              double InputDCOffsetI,
              double InputDCOffsetQ,
              double rline[]
              );

int G2GetPointSim(
              FILE *msg,
              Int4B RngFFTSize,
              Int4B upsample,
              double fcenter,
              double NarrowChirpBandwidth,
              double InputA2DFreq,
              double NarrowPulseLen,
              double InputStartSampleDelay,
              double RangeTarget,
              double rline[]
              );

int G2PhaseComp(
              FILE *msg,
              Int4B RngFFTSize,
              Int4B NewRngBinsToProcess,
              double InputA2DFreq,
              double tshift,
              double newrline[]
              );

int G2PutRangeLine( 
              FILE *msg,
              FILE *outfile,
              Int4B RngBinsToProcess,
              double rline[]
              );

int G2RangeCompress(
              FILE *msg,
              Int4B RngFFTSize,
              Int4B RngComRefFuncPhaseSign,
              double NarrowChirpBandwidth,
              double InputA2DFreq,
              double NarrowPulseLen,
              double rline[]
              );

int G2GetRangeProfile(
              FILE *msg,
              Int4B NewRngBinsToProcess,
              double newrline[],
              double Hsys[]
              );




