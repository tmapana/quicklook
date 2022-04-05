/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2func.h
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Header file for g2func.c. Also included widely
throughout G2 processor modules.

VERSION/AUTHOR/DATE : 1999-01-16 / Jasper Horrell / 1999-01-16
COMMENTS: 
Changed default values from 999 to -999.

VERSION/AUTHOR/DATE : 1999-01-26 / Jasper Horrell / 1997-01-26
COMMENTS:
Include FUNC_ defines here rather than in each file to make compiling 
main prog easier to set up.

VERSION/AUTHOR/DATE : 1999-01-27 / Richard Lord, Jasper Horrell / 1999-01-27
COMMENTS:
Add Richard's median function prototype. 

VERSION/AUTHOR/DATE : 1999=01=31 / Jasper Horrell / 1999-01-31
COMMENTS:
Add FUNC_MOC define. Added SEEKP_SET, SEEKP_CUR, SEEKP_NOLIM

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 1999-08-04
COMMENTS:
Add nextPowerOfTwo function prototype.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 2000-04-04
COMMENTS:
Added CFFT_float function prototype.

=========================================*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>


/* SET COMPILER TYPE (choose one only) */
#define DEC_CC 0
#define GNUC 1
#define BC31 0

/* Set up data types for compiler */
#if DEC_CC
#define Int2B short int
#define Int4B int
#elif GNUC
#define Int2B short int
#define Int4B long int
#elif BC31
#define Int2B short int
#define Int4B long int
#endif

/* Compile programs as functions (set to 1) or standalones (set to 0).
   Set all to 1 if compling main G2 program. */
#define FUNC_SNIFFDC 0
#define FUNC_MOC 0
#define FUNC_RNGCOM 0
#define FUNC_CORNER 0
#define FUNC_AZCOM 0
#define FUNC_B2TIF 0
#define FUNC_SUN_RASTER 0

/* Misc defns */
#define C 2.99792e+08
#ifndef PI
#define PI 3.141592654
#endif
#define STRING_SPACE 150
#define MAX_FILE_NAME 150
#define DEFAULT_VAL_I -999
#define DEFAULT_VAL_F -999.0
#define SEEKP_SET 0
#define SEEKP_CUR -1
#define SEEKP_NOLIM -1

/* Function prototypes for g2func.c file */
Int2B HostBigEndian();
Int2B ReadStr(FILE *fp, char Str[], Int2B MaxLen);
Int2B ReadStr2(FILE *fp, char Str[], Int2B MaxLen);
Int2B ReadChar(FILE *fp, char *Ch);
Int2B FindStr(FILE *fp, char SearchStr[]);
Int2B SeekP(FILE *fp,char Str1[],char Str2[],Int4B MaxOffset,Int4B Whence);
Int2B FindParam(FILE *fp, char ParamStr[]);
Int2B seekarrow(FILE *fp);
Int2B skipchar(FILE *fp,char ch);
Int2B GetString(FILE *fp,char *string);
void CFFT(double data[],Int4B nn,Int4B isign);
void CFFT_float(float data[],Int4B nn,Int4B sign);
Int2B Cmult(double z1[2],double z2[2], double *z);
Int2B IQ2power(double pow_array[],double iq_array[], Int4B size);
Int4B CorrelMaxIndex(double array1[],double array2[],Int4B size);
Int4B ConvertLong(Int4B data);
void FindDateTime(char DateNow[], char TimeNow[]);
Int4B BigEndWriteInt4B(FILE *fp,Int4B value);
Int2B BigEndWriteInt2B(FILE *fp,Int2B value);
Int4B BigEndReadInt4B(FILE *fp);
Int2B BigEndReadInt2B(FILE *fp);
double median(double array[], Int4B n);
Int4B nextPowerOfTwo(Int4B MinSize);
