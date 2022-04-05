/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2parse.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Code to parse SASAR LBR file.

VERSION/AUTHOR/DATE : 1999-01-16 / Jasper Horrell / 1999-01-16
COMMENTS: 
First version adapted from G2CEOS program (CEOS prog ver. 1998-11-24).

VERSION/AUTHOR/DATE : 1999-01-20 / Jasper Horrell / 1999-01-20
COMMENTS:
Minor clean up. Added PRIIDIncr field to DataChannelStruct. Changed 
DataChannelStruct so that all units MKS (Hz, sec, etc) 

VERSION/AUTHOR/DATE : 0.2 /jmh / 1999-07-22
COMMENTS:
(IMUStruct) GPS_UTC_sec changed to double from Int4B.

VERSION/AUTHOR/DATE : 
COMMENTS:

=========================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#include"g2func.h"
#include"g2parse.h"

/* Misc defns limited to this file  */
#define PROG_VERSION "1999-01-20"
#define FEET2METRES       0.3048

/* Fixed fields, for now */
#define DEBUG             0


/*****************************************/
/* Converts from raw Msg 5 format to understandable */
double ConvMsg5Val(Int4B RawVal, Int4B MSB)
{
  return ( ((double)ConvertLong(RawVal)) / (pow(2,31)/(double)MSB) );
}  

/***************************************************************/
/* Function to read in binary data message from IMU/GPS
   FIN 3110 system of Marconi and convert to readable format */
void ReadIMUMsg5(FILE *LBRFile,struct IMUStruct *IMU)
{
  struct Msg5Struct Msg5;
  double tmp;
  
  /* read in raw binary IMU values */
  fread( &Msg5.GPS_UTC, sizeof(Msg5), 1, LBRFile);

  /* Convert to readable format and in units of hrs,min,sec,deg,metres,m/s */
  tmp = ConvMsg5Val(Msg5.GPS_UTC,131072 );
  IMU->GPS_UTC_hrs  = (Int4B) (tmp / 3600.0);
  tmp -=  IMU->GPS_UTC_hrs * 3600.0;
  IMU->GPS_UTC_min = (Int4B) (tmp / 60.0);
  tmp -= IMU->GPS_UTC_min * 60.0;
  IMU->GPS_UTC_sec = tmp;

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

/***************************************************************
  Function to write header line in ASCII IMU motion file
  (see WriteMotion function below) */
void WriteMotionHeader(FILE *OutF)
{
  fprintf(OutF,"LBR_PRI GPS_UTC(secs)");
  fprintf(OutF," InertLat(deg) InertLong(deg) BaroInertHgt(m)");
  fprintf(OutF," INTLat(deg) INTLong(deg) INTHgt(deg)");
  fprintf(OutF," GPSLat(deg) GPSLong(deg) GPSHgt(m)");
  fprintf(OutF,"\n");
}

/***************************************************************/
/* Function to write out IMU/GPS data to ASCII file (appends)*/
void WriteMotion(FILE *OutF, struct IMUStruct IMU,Int4B PRI)
{
  fprintf(OutF,"%ld\t%.6f",
           PRI,IMU.GPS_UTC_hrs*3600+IMU.GPS_UTC_min*60+IMU.GPS_UTC_sec);
  fprintf(OutF,"\t%.8f\t%.8f\t%.4f",
           IMU.InertLat,IMU.InertLong,IMU.BaroInertHgt);
//  fprintf(OutF,"\t%.4f\t%.4f\t%.4f",         
//           IMU.InertVelNorth,IMU.InertVelEast,IMU.InertVelVert);
  fprintf(OutF,"\t%.8f\t%.8f\t%.4f",         
           IMU.INTLat,IMU.INTLong,IMU.INTHgt);
//  fprintf(OutF,"\t%.4f\t%.4f\t%.4f",         
//           IMU.INTVelNorth,IMU.INTVelEast,IMU.INTVelVert);
  fprintf(OutF,"\t%.8f\t%.8f\t%.4f",
           IMU.GPSLat,IMU.GPSLong,IMU.GPSHgt);
//  fprintf(OutF,"\t%hd\t%.4f\t%.4f",         
//           IMU.GPSRxStatus,IMU.GPSEstHorizErr,IMU.GPSEstVertErr);
//  fprintf(OutF,"\t%.4f\t%.4f\t%.4f",         
//           IMU.GPSVelNorth,IMU.GPSVelEast,IMU.GPSVelVert);
//  fprintf(OutF,"\t%.4f\t%.4f\t%.4f",         
//           IMU.InertRollAngle,IMU.InertPitchAngle,IMU.InertHeadingTrue);
//  fprintf(OutF,"\t%.4f\t%.4f\t%.4f",         
//           IMU.INTRollAngle,IMU.INTPitchAngle,IMU.INTHeadingTrue);

  fprintf(OutF,"\n");

} /* end function WriteMotion() */


/*********************************************/
void ParseDataChannels(FILE *msg,
                       FILE *LBRFile,
                       struct CntrlStruct *Gen,
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
  for (ii=0;ii<Gen->NumDataChannels;ii++) {
    Warns = 0; 

    /* 
    fprintf(msg,"Parsing data channel %d config info...\n",ii);
     */
    
    /* Parse general configuration info */
    PrevWarns = Warns;
    MaxO = 1000; 
    if (SeekP(LBRFile,"// General Configuration","",-1,0))
      { fprintf(msg,
	  "ERROR - general config header not found!\n"); exit(1);}
    else
      { St = ftell(LBRFile); }    
    if (SeekP(LBRFile,"Sensor",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Sensor,STRING_SPACE)) Warns++; 
    if (SeekP(LBRFile,"Location",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Location,STRING_SPACE)) Warns++;  
    if (SeekP(LBRFile,"Mission Number",":",MaxO,St)) Warns++; 
    else if (ReadStr(LBRFile,Chan[ii].Mission,STRING_SPACE)) Warns++; 
    if (SeekP(LBRFile,"Track Number",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Track,STRING_SPACE)) Warns++;  
    if (SeekP(LBRFile,"Date",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Date,STRING_SPACE)) Warns++;
    if (SeekP(LBRFile,"Nominal Ground Speed",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].NomGroundSpeed)) Warns++;      
    if (SeekP(LBRFile,"Average Terrain Height",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].AveTerrainHeight))Warns++;

    Chan[ii].GenCommentLen = 0;
    sprintf(Chan[ii].GenComment,"%s","");
    if (SeekP(LBRFile,"Length of comments",":",MaxO,St))
      { Warns++; }
    else if (!fscanf(LBRFile,"%ld",&Chan[ii].GenCommentLen))
      { Warns++; }  
    else if (SeekP(LBRFile,"Comments",":",MaxO,St))
      { Warns++; }
    else
      {
      if (Chan[ii].GenCommentLen > COMMENT_SPACE)
	Chan[ii].GenCommentLen = COMMENT_SPACE;
      for (i=0;i<Chan[ii].GenCommentLen;i++)
	Chan[ii].GenComment[i]=fgetc(LBRFile);  /* read gen comment */
      Chan[ii].GenComment[Chan[ii].GenCommentLen] = '\0';
      }
    
    
    if (Warns > PrevWarns)
      fprintf(msg,
	 "WARNINGS - during gen config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Parse LBR configuration info */
    PrevWarns = Warns;
    MaxO = 500;
    if (SeekP(LBRFile,"// LBR Configuration","",-1,0))
      { fprintf(msg,
	  "ERROR - LBR configuration header not found!\n"); exit(1);}
    else
      { St = ftell(LBRFile); }
    if (SeekP(LBRFile,"Enable PRI",":",MaxO,St)) Warns++;
    else if (ReadChar(LBRFile, &Chan[ii].PRICorrectionFlg)) Warns++;

    if (Warns > PrevWarns)
      fprintf(msg,
	 "WARNINGS - during LBR config parse : %ld!!\n",Warns-PrevWarns);
    
    /* Parse TMC configuration info */
    PrevWarns = Warns;
    MaxO = 1000;
    if (SeekP(LBRFile,"// TMC Configuration","",-1,0)) { 
      fprintf(msg,"ERROR - TMC configuration header not found!\n"); 
      exit(1);
    }
    else { 
      St = ftell(LBRFile);
    }

    if (Chan[ii].Radar == 'A') {
      if (SeekP(LBRFile,"RADAR A Sampling Rate",":",MaxO,St))
	    Warns++;
      else if (!fscanf(LBRFile,"%lf",&Chan[ii].ADFreq))
	    Warns++;
	  Chan[ii].ADFreq *= 1.0e+06;  /* convert to Hz */  
      if (SeekP(LBRFile,"RADAR A Sampling Bins",":",MaxO,St))
	    Warns++;
      else if (!fscanf(LBRFile,"%ld",&Chan[ii].RngBins)) { 
        fprintf(msg,"ERROR - Num rng bins not read!\n"); 
        exit(1); 
      } 
    } /* end if Chan[ii].Radar == 'A' */
    else if (Chan[ii].Radar == 'B') {
      if (SeekP(LBRFile,"RADAR B Sampling Rate",":",MaxO,St))
	    Warns++;
      else if (!fscanf(LBRFile,"%lf",&Chan[ii].ADFreq))
	    Warns++; 
	  Chan[ii].ADFreq *= 1.0e+06;  /* convert to Hz */  	    
      if (SeekP(LBRFile,"RADAR B Sampling Bins",":",MaxO,St))
	    Warns++;
      else if (!fscanf(LBRFile,"%ld",&Chan[ii].RngBins)) { 
        fprintf(msg,"ERROR - Num rng bins not read!\n"); 
        exit(1);
      } 
    } /* end else if Chan[ii].Radar == 'B' */
    else { 
      fprintf(msg,"ERROR - Radar %c unknown!\n",Chan[ii].Radar);
      exit(1);
    } /* end else */	

    if (Chan[ii].RngBins < 115) {
      fprintf(msg,
	 "ERROR - current data format requires at least 115 range bins\n");
      fprintf(msg,
	 "Could be fixed in later software versions with zero padding\n");
      exit(1);
    }
    
    if (Warns > PrevWarns)
      fprintf(msg,"WARNINGS - during TMC config parse : %ld!!\n",
              Warns-PrevWarns);

    /* Parse HBR configuration info */    

    PrevWarns = Warns;
    Fatals = 0;
    MaxO = 1100;
    if (Chan[ii].HBR == '1') { 
      if (SeekP(LBRFile,"// HBR 1 Configuration","",-1,0))
	    Fatals++;
	}
    else if (Chan[ii].HBR == '2') { 
      if (SeekP(LBRFile,"// HBR 2 Configuration","",-1,0))
	    Fatals++; 
	}
    else if (Chan[ii].HBR == '3') { 
      if (SeekP(LBRFile,"// HBR 3 Configuration","",-1,0))
	    Fatals++; 
	}
    else if (Chan[ii].HBR == '4') { 
      if (SeekP(LBRFile,"// HBR 4 Configuration","",-1,0))
	    Fatals++; 
	}
    else if (Chan[ii].HBR == '5') { 
      if (SeekP(LBRFile,"// HBR 5 Configuration","",-1,0))
	    Fatals++; 
	}
    else if (Chan[ii].HBR == '6') { 
      if (SeekP(LBRFile,"// HBR 6 Configuration","",-1,0))
	    Fatals++; 
	}
    else { 
      Fatals++; 
    }

    if (Fatals != 0) { 
      fprintf(msg,"ERROR - HBR config header not found!\n");
      exit(1); 
    }
    else { 
      St = ftell(LBRFile); 
    }
   
    if (SeekP(LBRFile,"Bits per I and Q",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%hd",&Chan[ii].Bits)) Warns++; 
    if (SeekP(LBRFile,"Master PRF",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].MasterPRF)) Warns++;
    if (SeekP(LBRFile,"Carrier Frequency",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].CarrierFreq)) Warns++;
    Chan[ii].CarrierFreq *= 1.0e+06; /* convert to Hz */
    if (SeekP(LBRFile,"Polarization",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Polarization,3)) Warns++;
    if (SeekP(LBRFile,"Range Compression",":",MaxO,St)) Warns++;
    else if (ReadChar(LBRFile, &Chan[ii].RngCompressedFlg)) Warns++;
    if (SeekP(LBRFile,"Pulse Modulation",":",MaxO,St)) Warns++;
    else if (ReadStr(LBRFile,Chan[ii].Modulation,STRING_SPACE)) Warns++;
    if (SeekP(LBRFile,"Pulse Length",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].PulseLen)) Warns++;
    Chan[ii].PulseLen *= 1.0e-06; /* convert to secs */
    if (SeekP(LBRFile,"Pulse Bandwidth",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].PulseBandwidth)) Warns++;
    Chan[ii].PulseBandwidth *= 1.0e+06; /* convert to Hz */
    if (SeekP(LBRFile,"Pulse Power",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].PulsePower)) Warns++;
    Chan[ii].PulsePower *= 1.0e+03; /* convert to Watts */
    if (SeekP(LBRFile,"Antenna azimuth b",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].AzBeamwidth)) Warns++;
    if (SeekP(LBRFile,"Antenna elevation b",":",MaxO,St)) Warns++ ;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].ElevBeamwidth)) Warns++;
    if (SeekP(LBRFile,"Antenna depression a",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].DepressAngle)) Warns++;
    if (SeekP(LBRFile,"Antenna squint",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].SquintAngle)) Warns++;
    if (SeekP(LBRFile,"Calibration tone",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].CalToneLevel)) Warns++;
    if (SeekP(LBRFile,"Calibration frequency",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].CalFreq)) Warns++;
    if (SeekP(LBRFile,"System Losses",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].SysLosses)) Warns++;
    if (SeekP(LBRFile,"Measured Delay0",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].Delay0Meas)) Warns++;
    Chan[ii].Delay0Meas *= 1.0e-06; /* convert to secs */
    if (SeekP(LBRFile,"Calculated delay0",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].Delay0Calc)) Warns++;
    Chan[ii].Delay0Calc *= 1.0e-06; /* convert to secs */
    if (SeekP(LBRFile,"Calculated STC en",":",MaxO,St)) Warns++;
    else if (ReadChar(LBRFile,&Chan[ii].STCFlg)) Warns++;
    if (SeekP(LBRFile,"STC/DGC Assigned ch",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%hd",&Chan[ii].DGCChannel)) Warns++;
    if (SeekP(LBRFile,"Receiver noise figure",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%lf",&Chan[ii].RxNoiseFig)) Warns++;
    if (SeekP(LBRFile,"DC offset",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%hd",&Chan[ii].DCOffset)) Warns++;
    if (SeekP(LBRFile,"Presum Ratio",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%hd",&Chan[ii].PresumRatio)) Warns++;
    if (SeekP(LBRFile,"Digital Gain",":",MaxO,St)) Warns++;
    else if (!fscanf(LBRFile,"%ld",&Chan[ii].DigGain)) Warns++;
    
    if (Warns > PrevWarns)
      fprintf(msg,
	 "WARNINGS - during HBR config parse : %ld!!\n",Warns-PrevWarns);

    /* Find start PRIID (don't use value under TMC config as unreliable) */
    if (SeekP(LBRFile,"RADAR B Status bit definition","",
        SEEKP_NOLIM,SEEKP_SET) ){
      fprintf(msg,"ERROR - in search for 'RADAR B Status bit'..!\n");
      exit(0); 
    }
    St = ftell(LBRFile);

    if (Chan[ii].Radar == 'A') {
      if (SeekP(LBRFile,"RADAR A Start PRIID Number",":",SEEKP_NOLIM,St))
	    Warns++;
      else if (!fscanf(LBRFile,"%ld",&Chan[ii].StartPRI)) { 
        fprintf(msg,"ERROR - Start PRI ID not read!\n"); 
        exit(1); 
      } 
    }
    else if (Chan[ii].Radar == 'B') {
      if (SeekP(LBRFile,"RADAR B Start PRIID Number",":",SEEKP_NOLIM,St)) {
	    Warns++;
	  }  
      else if (!fscanf(LBRFile,"%ld",&Chan[ii].StartPRI)) { 

        fprintf(msg,"ERROR - Start PRI ID not read!\n");
        exit(1);
      } 
    }
  
    /* Find end PRIID */
    fseek(LBRFile,-500,SEEK_END); /* Move close to end of file */
    St = ftell(LBRFile);
    MaxO = 510;
    Fatals = 0;
    if (Chan[ii].Radar == 'A')
      { if (SeekP(LBRFile,"RADAR A Stop PRIID Number",":",MaxO,St))
	Fatals++; }
    else if (Chan[ii].Radar == 'B')
      { if (SeekP(LBRFile,"RADAR B Stop PRIID Number",":",MaxO,St))
	Fatals++; }   
    if (!fscanf( LBRFile,"%ld",&Chan[ii].EndPRI)) Fatals++;  

    if (Fatals > 0)
      { fprintf(msg,
	   "ERROR - Radar stop PRI ID not read!\n"); exit(1); }

    Chan[ii].PRIIDIncr = (Int4B)( (Chan[ii].RngBins*2) / BYTES_PER_SECTOR);
    Chan[ii].EndPRI = Chan[ii].EndPRI - Chan[ii].PRIIDIncr; /* note end PRI in LBR
                                                      file not inclusive */
                                                      
    Chan[ii].PRIs = (Int4B)( (Chan[ii].EndPRI-Chan[ii].StartPRI)/
                    Chan[ii].PRIIDIncr + 1);
    
    /* Find IMU data for PRIID close to scene centre, if exists */
    if (SeekP(LBRFile,"Sync IDA","",-1,0)) { /* no IMU data */
      fprintf(msg,"WARNING - no IMU/GPS data found!\n");
      Chan[ii].SceneCenUTC_hrs = DEFAULT_VAL_I;
      Chan[ii].SceneCenUTC_min = DEFAULT_VAL_I;
      Chan[ii].SceneCenUTC_sec = DEFAULT_VAL_I;;
      Chan[ii].SceneCenLat = DEFAULT_VAL_F;
      Chan[ii].SceneCenLong = DEFAULT_VAL_F;
      Chan[ii].SceneCenTrueHead = DEFAULT_VAL_F;
    }
    else { /* IMU data exists */
      PRICen = Chan[ii].StartPRI +
	           (Int4B)(( (double)Chan[ii].EndPRI-
		       (double)Chan[ii].StartPRI )/2.0);

      if (Chan[ii].Radar == 'B')
	    SeekP(LBRFile,"IDB","",50,-1); /* find 'IDB' for Radar B */
      fscanf(LBRFile,"%ld",&PRINext);  /* Read in PRI */          

      while ( PRINext < PRICen ) {
	    if (SeekP(LBRFile,"\nSync IDA","",-1,-1)) { 
	      break; 
	    } 
	    else {
          if (Chan[ii].Radar == 'B')
            SeekP(LBRFile,"IDB","",50,-1); /* find 'IDB' for Radar B */
          fscanf(LBRFile,"%ld",&PRINext);  /* Read in PRI */     
	    } /* end else */
	  } /* end while */

      /* Read in IMU data line (note fp not updated in while loop
	  if Sync not found). If Sync found, the line read here is the
	  one equal to or immediately after the centre PRI  */	
      SeekP(LBRFile,"\n","",-1,-1);
      ReadIMUMsg5(LBRFile,&IMUDat); /* Read and convert IMU record */

      /* Assign values */
      Chan[ii].SceneCenUTC_hrs  = IMUDat.GPS_UTC_hrs;
      Chan[ii].SceneCenUTC_min  = IMUDat.GPS_UTC_min;
      Chan[ii].SceneCenUTC_sec  = IMUDat.GPS_UTC_sec;
      Chan[ii].SceneCenLat      = IMUDat.INTLat;
      Chan[ii].SceneCenLong     = IMUDat.INTLong; 
      Chan[ii].SceneCenHgt      = IMUDat.INTHgt;
      Chan[ii].SceneCenTrueHead = IMUDat.INTHeadingTrue;

      /* 
      fprintf(msg,
	  "LBR scene centre GPS UTC (hh:min:sec): %02ld:%02ld:%02ld\n",
	  Chan[ii].SceneCenUTC_hrs,Chan[ii].SceneCenUTC_min,
	  Chan[ii].SceneCenUTC_sec);
	  fprintf(msg,
	  "LBR scene centre Lat / Long / Height : %f deg/ %f deg/ %f m\n",
	  Chan[ii].SceneCenLat,Chan[ii].SceneCenLong,
	  Chan[ii].SceneCenHgt); 
	  */
	      
    }  /* end else IMU data exists */
 
       
    /* Display total warnings */
    /* 
    fprintf(msg,
       "Warnings in parsing config info (channel %d) : %ld\n",ii,Warns);
    */
    
    } /* End for ii (DataChannel) loop */

} /* end ParseDataChannels function */

/**************************************/

void ParseLBRFile(FILE *msg,
                  FILE *LBRFile,
                  struct CntrlStruct *Gen,
	              struct DataChannelStruct Chan[])
{
  Int2B DataChannel=0;

  /* Read Radars and HBRs used */
  if (SeekP(LBRFile,"// General Configuration","",-1,0)) {
    fprintf(msg,"ERROR at start - General config header not found!\n");
    exit(1);
  }
 
  if (SeekP(LBRFile,"Systems enabled",":",200,-1)) {
    fprintf(msg,"ERROR - 'Systems enabled' parameter not found!\n");
    exit(1);
  }
    
  while (DataChannel < MAX_DATA_CHANNELS) {
    if (SeekP(LBRFile,"RADAR_","",100,-1)) {
      Chan[DataChannel].Radar = getc(LBRFile);      
      if (SeekP(LBRFile,"_HBR","",100,-1)) {
        Chan[DataChannel].HBR = getc(LBRFile);
        Gen->NumDataChannels++;
      }  /* end if SeekP _HBR */
      else {
        DataChannel = MAX_DATA_CHANNELS;
      }  
    }  /* end if SeekP RADAR_ */
    else {
      DataChannel = MAX_DATA_CHANNELS;
    }  
    DataChannel++;
  }  /* end while DataChannel */
    
  if (Gen->NumDataChannels==0) { 
    fprintf(msg,"ERROR - No valid data channels!\n"); 
    exit(1);
  }
  else {
    fprintf(msg,"Number of data channels = %d\n",Gen->NumDataChannels);
  }  

  /* Parse data channels */
  ParseDataChannels(msg,LBRFile,Gen,Chan); /* note Gen is pointer here hence no amplisand */

} /* end ParseLBRFile function */
