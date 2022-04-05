/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 2000
FILE NAME: pga.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION:
Possible part of G2 SAR processor. Phase gradient autofocus
for fine tuning SAR mocomp. 

VERSION/AUTHOR/DATE : 0.1 / Jasper Horrell / 2000-04-03
COMMENTS:
Initial version based on description of PGA in book by Jacowatz -
"Spotlight Mode SAR - Signal Processing Approach" (Kluwer Acedemic
Publishers).
2000-04-05 - Fixed window FFT so that image not shifted in azimuth.
             Fixed bug in ConjMultSum indexing. Also, minor cleanup.

============================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "g2func.h"

#define Int4B long int
#define DEF_VAL_I -999999

/***********
MAIN PROGRAM
************/

int main(int argc, char *argv[])
{
  FILE *InFile=NULL,
       *msg=stdout,
       *OutFile=NULL,
       *TmpFile=NULL;

  char tmpstr[STRING_SPACE]="",
       CycShiftFileName[STRING_SPACE]="",
       PhErrorFileName[STRING_SPACE]="",
       RowSumDBFileName[STRING_SPACE]="";

  Int4B BytesPerVal,
        BytesToSkip,
        Cols=-DEF_VAL_I,
        ColsToProc=DEF_VAL_I,
        ColsToProcD2,
        ColsToProcX2,
        FFTSize,
        FFTSizeX2,
        Iteration=0,
        MaxCol=DEF_VAL_I, // column with max value for row
        MaxIterations=DEF_VAL_I,
        ProgressiveWin=0, // default is off
        ShiftX2=0,
        StartCol=0, // to process
        StartRow=0, // to process
        Rows=DEF_VAL_I,
        RowsToProc=DEF_VAL_I,
        WinColSize=DEF_VAL_I,
        WinColSizeX2=DEF_VAL_I,
        WinLeftCol,  // column marking left side of window
        WinRightCol; // column marking right side of window
        
  register Int4B i,j,
                 indx,indxi,indxq,
                 row,col;

  float *ConjMultSum=NULL,
        **InData=NULL,
        *InDataBlock=NULL,
        *FFTData=NULL,
        *PhError=NULL,
        *RowData=NULL,
        *RowPow=NULL,
        **ShiftData=NULL,
        *ShiftDataBlock=NULL,
        MaxVal,  // for row (power)
        PhaseErrorRMS,  // RMS phase error
        PhaseErrorTgt=20.0,
        tmpi,tmpq,tmpi2,tmpq2,
        WinThresh=-10.0,  // window threshold value in dB for size determination
        WinThreshPow;  // WinThresh (power) times peak power
    
  time_t TimeEnd,
         TimeStart;
         

  printf("\nProgram: PGA - Ver. 0.1 (Phase Gradient Autofocus for SAR)\n");
  printf("Code: J.M. Horrell - (C) UCT RRSG (2000)\n\n");

  if (argc < 5 || argc > 13) {
    printf("1) Read azimuth line format complex float data into 2-D array.\n");
    printf("2) Find for each az line position of max power value.\n");
    printf("3) Cyclic shift data, range sum power, and determine window width.\n");
    printf("4) Apply window to each shifted line and FFT\n");
    printf("5) Estimate phase error for each az position.\n");
    printf("6) Correct each unshifted line for the phase error in az freq domain.\n");
    printf("7) Iterate (2-6) until RMS phase error below threshold.\n");
    printf("8) Write final focused image to output file.\n");
    printf("\nUSAGE:\n");
    printf(
     "  pga [InFile] [OutFile] [Rows] [Cols] <optional params>\n");
    printf("\nwhere:\n");
    printf("  InFile    - (req) input file name\n");
    printf("  OutFile   - (req) output file name\n");
    printf("  Rows      - (req) input file max rows\n");
    printf("  Cols      - (req) input file max columns (complex points per row)\n");
    printf("  WinThresh=f     - (opt) window threshold in dB,'prog' for progressive)\n");
    printf("                    (default=-10.0).\n");
    printf("  StartRow=n      - (opt) start row to use for processing (default=0)\n");
    printf("  RowsToProc=n    - (opt) num of rows to process (default=Rows-StartRow)\n"); 
    printf("  StartCol=n      - (opt) start column to process (default=0)\n");
    printf("  ColsToProc=n    - (opt) num of columns to process (default=Cols-StartCol)\n");
    printf("  PhaseErrorTgt=f - (opt) max RMS phase error in deg (default=20.0)\n"); 
    printf("  MaxIterations=n - (opt) max iterations (default use PhaseErrorTgt)\n");
    printf("  RowSumDBFile=s  - (opt) file name for row power dB sum (default none)\n");
    printf("  CycShiftFile=s  - (opt) file output of cyclic shifted data (default none)\n");
    printf("  PhErrorFile=s   - (opt) file output of deg phase error (default none)\n");
    printf("\ne.g. pga in.dat out.dat 1024 512 StartRow=42 PhaseErrorTgt=5.7\n\n");
    exit(1);
    }

  sscanf(argv[3],"%ld",&Rows);
  sscanf(argv[4],"%ld",&Cols);

  for (i=5; i<argc;i++) {
    if (strncmp(argv[i],"WinThresh=",10)==0) {
      sscanf(argv[i],"WinThresh=%s",tmpstr);
      if (strcmp(tmpstr,"prog")==0) {
        ProgressiveWin=1;
      }
      else { 
        sscanf(tmpstr,"%f",&WinThresh);
      }
    }
    if (strncmp(argv[i],"StartRow=",9)==0) { 
      sscanf(argv[i],"StartRow=%ld",&StartRow);
    }  
    if (strncmp(argv[i],"RowsToProc=",11)==0) { 
      sscanf(argv[i],"RowsToProc=%ld",&RowsToProc);
    }
    if (strncmp(argv[i],"StartCol=",9)==0) { 
      sscanf(argv[i],"StartCol=%ld",&StartCol);
    }  
    if (strncmp(argv[i],"ColsToProc=",11)==0) { 
      sscanf(argv[i],"ColsToProc=%ld",&ColsToProc);
    }    
    if (strncmp(argv[i],"PhaseErrorTgt=",14)==0) { 
      sscanf(argv[i],"PhaseErrorTgt=%f",&PhaseErrorTgt);
    }
    if (strncmp(argv[i],"MaxIterations=",11)==0) { 
      sscanf(argv[i],"MaxIterations=%ld",&MaxIterations);
    }  
    if (strncmp(argv[i],"RowSumDBFile=%s",13)==0) {
      sscanf(argv[i],"RowSumDBFile=%s",RowSumDBFileName);
    }
    if (strncmp(argv[i],"CycShiftFile=%s",13)==0) {
      sscanf(argv[i],"CycShiftFile=%s",CycShiftFileName);
    }
    if (strncmp(argv[i],"PhErrorFile=%s",12)==0) {
      sscanf(argv[i],"PhErrorFile=%s",PhErrorFileName);
    }    
  }  // end for i

  /* Some checks and assignments */
  if (RowsToProc==DEF_VAL_I) RowsToProc=Rows-StartRow;
  if (ColsToProc==DEF_VAL_I) ColsToProc=Cols-StartCol;
  
  /* Messages */
  fprintf(msg,"MESSAGES:\n");
  fprintf(msg,"Input file              : %s\n",argv[1]);
  fprintf(msg,"Output file             : %s\n",argv[2]);
  fprintf(msg,"Rows                    : %ld\n", Rows);
  fprintf(msg,"Columns                 : %ld\n",Cols);
  fprintf(msg,"Start row to process    : %ld\n",StartRow);
  fprintf(msg,"Rows to process         : %ld\n",RowsToProc);
  fprintf(msg,"Start column to process : %ld\n",StartCol);
  fprintf(msg,"Columns to process      : %ld\n",ColsToProc);  
  if (ProgressiveWin) {
    fprintf(msg,"Window threshold        : progressive\n");
  }
  else {
    fprintf(msg,"Window threshold        : %f dB\n",WinThresh);
  }
  fprintf(msg,"RMS phase error target  : %f deg\n",PhaseErrorTgt);
  if (MaxIterations != DEF_VAL_I) {
    fprintf(msg,"Max iterations          : %ld\n",MaxIterations);
  }
  else {
    fprintf(msg,"Max iterations          : unspecified\n");
  }
  if (strcmp(RowSumDBFileName,"") != 0)
    fprintf(msg,"Row sum dB file         : %s\n",RowSumDBFileName);
  if (strcmp(CycShiftFileName,"") != 0)
    fprintf(msg,"Cyclic shift file       : %s\n",CycShiftFileName);
  if (strcmp(PhErrorFileName,"") != 0)
    fprintf(msg,"Phase error file        : %s\n",PhErrorFileName);        
  fprintf(msg,"\n");

  /* Misc */
  TimeStart = time(NULL);
  BytesPerVal = 2*sizeof(float);  // size of complex value
  ColsToProcX2 = 2*ColsToProc;
  ColsToProcD2 = (Int4B)((float)ColsToProc/2.0);
  FFTSize = nextPowerOfTwo(ColsToProc); // use one size for all (simpler)
  FFTSizeX2 = 2*FFTSize;
  PhaseErrorTgt *= PI/180.0; // convert to radians
 
  /* Allocate array space */
  InData = (float **)malloc(sizeof(float *)*RowsToProc);
  InDataBlock = (float *)malloc(BytesPerVal*RowsToProc*ColsToProc);
  ShiftData = (float **)malloc(sizeof(float *)*RowsToProc);
  ShiftDataBlock = (float *)malloc(BytesPerVal*RowsToProc*ColsToProc);
  RowData = (float *)malloc(sizeof(float)*ColsToProcX2);
  RowPow = (float *)malloc(sizeof(float)*ColsToProc);
  FFTData = (float *)malloc(sizeof(float)*FFTSizeX2);
  ConjMultSum = (float *)malloc(sizeof(float)*FFTSizeX2);
  PhError = (float *)malloc(sizeof(float)*FFTSize);
 
  if (InDataBlock==NULL || InData==NULL || ShiftDataBlock==NULL ||
      ShiftData==NULL || RowData==NULL || RowPow==NULL ||
      FFTData==NULL || ConjMultSum==NULL || PhError==NULL)
    { fprintf(msg,"ERROR - in array mem allocation\n"); exit(1); }

  for (row=0; row<RowsToProc; row++) {    // Arrange 2-D array indexing
    InData[row] = &InDataBlock[row*ColsToProcX2];
    ShiftData[row] = &ShiftDataBlock[row*ColsToProcX2];
  }
 
  /* Open input file, move to start, read data, and close */
  InFile  = fopen(argv[1],"rb");
  if (InFile==NULL) { 
     printf("ERROR - Input file %s not opened/found\n",argv[1]);
     exit(1);
  }
  fseek(InFile,StartRow*BytesPerVal*Cols+StartCol*BytesPerVal,SEEK_SET); 

  BytesToSkip = (Cols-ColsToProc)*BytesPerVal;
  for (row=0;row<RowsToProc;row++) {
    fread(RowData,sizeof(float),ColsToProcX2,InFile);
    for (j=0;j<ColsToProcX2;j++) {
      InData[row][j] = RowData[j];  // assign data to 2-D array
    }  // end for indx
    if (row != RowsToProc-1) fseek(InFile,BytesToSkip,SEEK_CUR);
  }  // end for row

  fclose(InFile);


  /* ITERATE UNTIL PHASE ERROR TARGET REACHED (or MaxInterations reached) *
   * -------------------------------------------------------------------- */

  if (ProgressiveWin) WinColSize = ColsToProc;
  PhaseErrorRMS = 999999.9;  // set arbitrary large
  Iteration = 0;
  while (PhaseErrorRMS > PhaseErrorTgt) {

    if (MaxIterations != DEF_VAL_I && Iteration >= MaxIterations) {
      fprintf(msg,"MaxIterations exceeded - breaking!\n");
      break;
    }

    /* Centre Shift Data and store in new array. The max is shifted to 
     * column ColsToProcD2 */
    for (row=0;row<RowsToProc;row++) {  // repeat for all rows of image

      /* Find max power value for row */
      MaxVal = -999.9;
      indxi = 0; indxq = 1;
      for (col=0;col<ColsToProc;col++) {
        tmpi = InData[row][indxi++]; indxi++;
        tmpq = InData[row][indxq++]; indxq++;
        if ( (tmpi*tmpi+tmpq*tmpq) > MaxVal) {
          MaxVal = tmpi*tmpi+tmpq*tmpq;
          MaxCol = col;
        }
      }

      /* Store in cyclically shifted array */
      ShiftX2 = 2*(ColsToProcD2-MaxCol); // max moved to column ColsToProcD2
      for (i=0;i<ColsToProcX2;i++) {
        if (ShiftX2 < 0) {
          ShiftData[row][i] = InData[row][(i-ShiftX2)%ColsToProcX2];
        }
        else { // for +ve shift
          ShiftData[row][i] = InData[row][(ColsToProcX2+i-ShiftX2)%ColsToProcX2];
        } 
      }  // end for i

    } // end for row
     
    /* Calc window size, if not progressive window */
    if (ProgressiveWin) {
      WinLeftCol = ColsToProcD2 - (int)((float)WinColSize/2.0);
      WinRightCol = ColsToProcD2 + (int)((float)WinColSize/2.0);
      if (WinColSize%2==0) WinRightCol = WinRightCol-1; // to get correct win size
    }
    else {
      /* sum power values of shifted data over range */
      for (col=0;col<ColsToProc;col++) RowPow[col] = 0.0;  // init RowPow array
      for (row=0;row<RowsToProc;row++) {
         indxi=0; indxq=1;
         for (col=0;col<ColsToProc;col++) {
           tmpi = ShiftData[row][indxi++]; indxi++;
           tmpq = ShiftData[row][indxq++]; indxq++;
           RowPow[col] += tmpi*tmpi + tmpq*tmpq;       
         }  // end for col          
      } // end for row

      /* find window threshold points on left and right */
      WinLeftCol = ColsToProcD2;
      WinRightCol = ColsToProcD2;
      WinThreshPow = RowPow[ColsToProcD2]*(float)pow(10.0,WinThresh/10.0);
      for (col=ColsToProcD2;col>=0;col--) {
        if (RowPow[col] < WinThreshPow) {
          WinLeftCol = col;
          break;
        }
      }
      for (col=ColsToProcD2;col<ColsToProc;col++) {
        if (RowPow[col] < WinThreshPow) {
          WinRightCol = col;
          break;
        }
      }
      WinColSize = WinRightCol - WinLeftCol + 1;
    }  // end else not using a progressive window

    if (WinColSize < 3) {
      fprintf(msg,"Window size %ld less than 3 - breaking!\n",WinColSize);
      break;
    }

    /* Repeat for each row of data */
    WinColSizeX2 = 2*WinColSize;
    for (j=0;j<FFTSizeX2;j++) ConjMultSum[j] = 0.0; // init ConjMultSum array
    for (row=0;row<RowsToProc;row++) {

      /* FFT az line of windowed data (put max at FFT index zero) */
      indx = 0;
      for (col=ColsToProcD2;col<WinRightCol+1;col++) { // put max at zero for FFT
        FFTData[indx++] = ShiftData[row][2*col];
        FFTData[indx++] = ShiftData[row][2*col+1];
      }
      for (col=0;col<(FFTSize-WinColSize);col++) { // zero pad FFT
        FFTData[indx++] = 0.0;
        FFTData[indx++] = 0.0;
      }
      for (col=WinLeftCol;col<ColsToProcD2;col++) { // put max at zero for FFT
        FFTData[indx++] = ShiftData[row][2*col];
        FFTData[indx++] = ShiftData[row][2*col+1];
      }      
      CFFT_float(FFTData-1,FFTSize,-1); // in-place FFT

      /* Add to the sum over range of az conj multiplies */
      indxi = 2; indxq = 3;
      for (col=1;col<FFTSize;col++) {
        tmpi = FFTData[indxi]; tmpq = FFTData[indxq];
        tmpi2 = FFTData[indxi-2]; tmpq2 = -FFTData[indxq-2]; // take conj
        ConjMultSum[indxi++] += tmpi*tmpi2 - tmpq*tmpq2; // add to sum
        ConjMultSum[indxq++] += tmpi*tmpq2 + tmpq*tmpi2; // add to sum
        indxi++; indxq++;            
      }
      
    }  // end for row

    /* Calc the phase error by integrating the phase error difference for
     * each az freq. Also calc the RMS phase error. */
    PhError[0] = 0.0;  // assign zero integrated phase error for first sample
    PhaseErrorRMS = 0;
    indxi = 0; indxq =1;    
    for (col=1;col<FFTSize;col++) {
      PhError[col] = PhError[col-1]+
                     (float)atan2(ConjMultSum[indxq++],ConjMultSum[indxi++]);
      indxi++; indxq++;
      PhaseErrorRMS += PhError[col]*PhError[col]; 
    }
    PhaseErrorRMS = (float)sqrt(PhaseErrorRMS/(float)FFTSize);
    fprintf(msg,"Iteration = %ld, Window Size = %ld, RMS Phase Error = %f deg\n",Iteration,
            WinColSize,PhaseErrorRMS*180.0/PI);

    /* Break if phase error target reached, else apply phase correction */
    if (PhaseErrorRMS < PhaseErrorTgt) {
      fprintf(msg,"RMS phase error target reached - breaking!\n");
      break;
    }
    else { 

      /* Apply phase correction to image by multiply in freq domain */
      for (row=0;row<RowsToProc;row++) {
        for (j=0;j<ColsToProcX2;j++) FFTData[j] = InData[row][j];
        for (j=ColsToProcX2;j<FFTSizeX2;j++) FFTData[j] = 0.0;
        CFFT_float(FFTData-1,FFTSize,-1); // in-place FFT of az line of image      

        indxi=0; indxq=1;
        for (j=0;j<FFTSize;j++) {  // multiply with conj of phase error
          tmpi = FFTData[indxi]; tmpq = FFTData[indxq]; 
          tmpi2 = cos(PhError[j]); tmpq2 = -sin(PhError[j]);
          FFTData[indxi++] = tmpi*tmpi2 - tmpq*tmpq2;
          FFTData[indxq++] = tmpi*tmpq2 + tmpq*tmpi2;
          indxi++; indxq++;
        } // end for j
        CFFT_float(FFTData-1,FFTSize,1); // in-place IFFT of corrected az line
        for (j=0;j<ColsToProcX2;j++)
          InData[row][j] = FFTData[j]/(float)FFTSize; // write back to 2-D array
      } // end for row

    } // end else not within phase error target
    
    /* Update the window size, if progressive window */
    if (ProgressiveWin) WinColSize = (Int4B)((float)WinColSize/2.0);

    Iteration++; // increment the iteration counter
    
  }  /* ------ END ITERATE WHILE PhaseErrorRMS < PhaseErrorTgt ----- */


  /* Write focused image to output file *
   * ---------------------------------- */
    
  OutFile = fopen(argv[2],"wb");
  if (InFile==NULL) {
     printf("ERROR - Output file %s not opened/found\n",argv[2]);
     exit(1);
  }

  for (row=0;row<RowsToProc;row++) {
    for (j=0;j<ColsToProcX2;j++) {
      RowData[j] = InData[row][j];
    }
    fwrite(RowData,sizeof(float),ColsToProcX2,OutFile);
  }

  fclose(OutFile);


  /* Write out optional output files *
   * ------------------------------- */

  /* Write out cyclic shifted image to file, if requested */
  if (strcmp(CycShiftFileName,"")!=0) {
    TmpFile  = fopen(CycShiftFileName,"wb");
    if (TmpFile==NULL) { 
      printf("ERROR - cyclic shift file %s not opened for write\n",argv[1]);
      exit(1);
    }
    for (row=0;row<RowsToProc;row++) {
      for (j=0;j<ColsToProcX2;j++) {
        RowData[j] = ShiftData[row][j];
      }
      fwrite(RowData,sizeof(float),ColsToProcX2,TmpFile);
    }
    printf("Final cyclic shifted data written to file %s!\n",CycShiftFileName);
    fclose(TmpFile);
  }  // end if strcmp 

  /* Write out RowPow in dB to file, if requested */
  if (strcmp(RowSumDBFileName,"")!=0) {
    TmpFile  = fopen(RowSumDBFileName,"wt");
    if (TmpFile==NULL) { 
      printf("ERROR - row sum dB file %s not opened for write\n",argv[1]);
      exit(1);
    }
    for (col=0;col<ColsToProc;col++)
       fprintf(TmpFile,"%f\n",(float)(10.0*log10(RowPow[col]))); // dB values 
    printf("Final row sum dB data written to file %s!\n",RowSumDBFileName);
    fclose(TmpFile);
  }  // end if strcmp

  /* Write out phase error to file, if requested */
  if (strcmp(PhErrorFileName,"")!=0) {
    TmpFile  = fopen(PhErrorFileName,"wt");
    if (TmpFile==NULL) { 
      printf("ERROR - phase error file %s not opened for write\n",argv[1]);
      exit(1);
    }
    for (col=0;col<FFTSize;col++)
       fprintf(TmpFile,"%f\n",PhError[col]*180.0/PI); //  
    fclose(TmpFile);
    printf("Final phase written to file %s!\n",PhErrorFileName);
  }  // end if strcmp

  /*--------end optional output files-----------*/

  /* Tidy up */
  free(InData[0]);    // frees InDataBlock
  free(InData);
  free(ShiftData[0]); // frees ShiftDataBlock
  free(ShiftData);
  free(RowData);
  free(RowPow);
  free(FFTData);
  free(ConjMultSum);
  free(PhError);

  TimeEnd = time(NULL); 
  fprintf(msg,"PGA autofocus done in %ld secs (%.3f min)!\n",TimeEnd-TimeStart,
          (float)(TimeEnd-TimeStart)/60.0);

  return(0);

} /* END MAIN PROG */
