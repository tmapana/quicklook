/*==========================================
COPYRIGHT: (C) UCT Radar Remote Sensing Group 1999
FILE NAME: g2func.c
CODE CONTROLLER: Jasper Horrell
DESCRIPTION: 
Part of G2 SAR processor. Contains misc functions for parsing etc.

VERSION/AUTHOR/DATE : 1b / Jasper Horrell (mostly) / 1997-11-28
COMMENTS: 
Rename to as part of the new g2 code (post ATP) - no changes yet apart 
from different included header file name

VERSION/AUTHOR/DATE : G2Func / Jasper Horrell / 1997-12-05
COMMENTS:
Rename. Also rename include file.

VERSION/AUTHOR/DATE : 1999-01-27 / Richard Lord, Mark Gebhardt,
                      Jasper Horrell  / 1999-01-27
COMMENTS:
Add Richard's median function. add Mark's changes to GetString() 
(should be backward compatible).

VERSION/AUTHOR/DATE : 1.0 / Jasper Horrell / 1999-08-04
COMMENTS:
Add nextPowerOfTwo() function.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 1999-09-29
COMMENTS:
Changed functions so that sucessful completion returns zero, and
error returns 1.

VERSION/AUTHOR/DATE : 1.1 / Jasper Horrell / 2000-04-04
COMMENTS:
Added CFFT_float function.

=========================================*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#include"g2func.h"

/**********************************************************************/
/*  FUNCTION :  HostBigEndian() */
/* Returns 0 if big-endian host machine, else returns 1 */

Int2B HostBigEndian()
{
  Int4B tmpI=1;
  unsigned char *ptrI;
  
  ptrI = (unsigned char *)&tmpI;
  if (  (*ptrI == 0) && (*(ptrI+3)== 1) )
    return(1);

  return(0);
}  
    


/**********************************************************************/
/* FUNCTION : ReadStr */
/* Read a string which may include white space up to an end of 
   line character or specified length (whichever comes first). 
   Remove the end of line char from end of string. Ignore first
   leading white space, if exists. Returns 1 for error, 0 for success. */

Int2B ReadStr(FILE *fp, char Str[], Int2B MaxLen)
{
  if (fgetc(fp)!=' ') fseek(fp,-1,1);
  fgets(Str,MaxLen,fp);
 
  if (strlen(Str)==0) 
    { return(1); }
  else if (strlen(Str)>1)
    {
    if (Str[strlen(Str)-2]=='\r' && Str[strlen(Str)-1]=='\n') /*DOS*/
      { Str[strlen(Str)-2]='\0'; }
    else if (Str[strlen(Str)-1]=='\n') 
      { Str[strlen(Str)-1]='\0'; }
    }
  else /* if strlen 1 */
    { if (Str[0]=='\n') Str[0]='\0'; }
 
  if (strlen(Str)==0) return(1);  /* for zero string (recheck) */

  return(0); /* for non zero string successfully read */
}

/**********************************************************************/
/* FUNCTION : ReadStr2 */
/* Read a string which may include white space up to an end of 
   line character or specified length (whichever comes first). 
   Remove the end of line char from end of string. Ignore first
   leading white space, if exists. Remove any trailing white space. */

Int2B ReadStr2(FILE *fp, char Str[], Int2B MaxLen)
{
  Int2B i;  

  if (fgetc(fp)!=' ') fseek(fp,-1,1);
  fgets(Str,MaxLen,fp);

  if (strlen(Str)==0) 
    { return(1); }
  else if (strlen(Str)>1)
    {
    if (Str[strlen(Str)-2]=='\r' && Str[strlen(Str)-1]=='\n') /*DOS*/
      { Str[strlen(Str)-2]='\0'; }
    else if (Str[strlen(Str)-1]=='\n') 
      { Str[strlen(Str)-1]='\0'; }
    }
  else /* if strlen 1 */
    { if (Str[0]=='\n') Str[0]='\0'; }

  if (strlen(Str)==0) return(1); 

  i = strlen(Str)-1;
  while (Str[i]==' ' && i>=0)
    { Str[i--] = '\0'; }  /* remove trailing white space */
   
  if (strlen(Str)==0) return(1);  /* for zero string (recheck) */

  return(0); /* for non zero string successfully read */
}

/*********************************************************************/
/* FUNCTION : ReadChar */
/* Read a char from a stream, but ignore the first leading
   white space, if exists */

Int2B ReadChar(FILE *fp, char *Ch)
{
  *Ch=fgetc(fp);
  if (*Ch==' ') *Ch = fgetc(fp);
  if (*Ch==EOF) { *Ch=' '; return(1); }

  return(0); /* if char successfully read */ 
}


/***********************************************************************/
/* FUNCTION : FindStr */
/* Read from current file pointer posn until a specified
   string is found (assumes file opened as binary) */

Int2B FindStr(FILE *fp, char SearchStr[])
{
  char ReadStr[80]="";
  Int2B i=0,StrSize,StrFound=0;
  Int4B FileSize, PtrPosn;

  /* Find size of file */
  PtrPosn = ftell(fp); /* note current pointer posn */
  fseek(fp,0,2);       /* move to end of file */
  FileSize = ftell(fp);
  fseek(fp,PtrPosn,0); /* move pointer back to orig posn */

  /* find size of string */
  StrSize = strlen(SearchStr);

  while (i<StrSize) {
    if (ftell(fp) >= FileSize) break;
    ReadStr[i] = fgetc(fp);
    if (ReadStr[i] == SearchStr[i])
      i++;
    else
      i=0;
    if (i==StrSize) StrFound=1;
  }

  return(StrFound);

} /* end function FindStr */


/***********************************************************************/
/* FUNCTION : SeekP (generic seek for parameters in file) */
/* Read from current or specified  posn in file until a specified string is 
   found. Then search for a second string. If either string 1 or string
   2 is not found, move file pointer back to posn at start of
   search. String 2 may be null, if required.
   Search only up to MaxOffset bytes from start posn. 
   (NB: Written to work with binary files, but also seems to work with
   text files.). */

Int2B SeekP(FILE *fp,     /* Pointer of file to search */
               char Str1[],  /* First string to find */
               char Str2[],  /* Second string to find ("" is OK) */
               Int4B MaxOffset, /* Extent of search in bytes from
				   start posn (<0 for no limit). */
               Int4B Whence) /* Start posn of search in bytes from
				start of file (-1 for search from
				current posn, 0 for search from
				file start). */
{

  char ReadStr[80]="",SearchStr[80];
  Int2B i=0,StrSize,Str1Found=0,Str2Found=0;
  Int4B FileSize,PtrPosn,Offset=0;

  /* Find size of file and perform check */
  PtrPosn = ftell(fp);
  fseek(fp,0,2);       /* move to end of file */
  FileSize = ftell(fp);
  if (Whence > FileSize-1) /* exit function if, non-valid call */
    return (1);        

  if (Whence != -1 )   /* if not searching from current posn */
    PtrPosn = Whence; 
  fseek(fp,PtrPosn,0); /* move pointer to start search posn */

  
  /* find size of string */
  strcpy(SearchStr,Str1); /* needs new string to work! */
  StrSize = strlen(SearchStr);
  
  /* Find string1, then string2 */
  while (i<StrSize) {
    if (ftell(fp) >= FileSize) 
      break;
    if ((MaxOffset >= 0) && (Offset++>MaxOffset)) 
      break;
    ReadStr[i] = fgetc(fp);
    if (ReadStr[i] == SearchStr[i]) { 
      i++; 
    }
    else { 
      i=0; 
    }
    if (i==StrSize && !Str1Found) {
      Str1Found=1;
      strcpy(SearchStr,Str2);
      StrSize=strlen(SearchStr);
      if (StrSize==0) {  /* if Str2 is a null string */
        Str2Found=1;
        break;    
	  }
      else {
        i=0; 
      }
    } /* end if i */
    else if (i==StrSize && Str1Found) { 
      Str2Found = 1; 
    }
  }  /* end while i < StrSize */

  /* if either string not found, move fp to posn at search start */
  if (!Str2Found) fseek(fp,PtrPosn,0);

  if (Str2Found) 
    return(0);  /* on success */
  else 
    return(1); 

} /* end function SeekP */


/***********************************************************************/
/* FUNCTION : FindParam */
/* Read from begin of file until a specified param is found. Then 
   search for short arrow (=>) (assumes file opened as binary).
   This function now largely superceded by the SeekP function.  */

Int2B FindParam(FILE *fp, char ParamStr[])
{
  char ReadStr[80]="",SearchStr[80];
  Int2B i=0,StrSize,StrFound=0,ArrowFound=0;
  Int4B FileSize;

  /* Find size of file */
  fseek(fp,0,2);       /* move to end of file */
  FileSize = ftell(fp);
  fseek(fp,0,0); /* move pointer to start of file */

  /* find size of string */
  strcpy(SearchStr,ParamStr); /* needs new string to work! */
  StrSize = strlen(SearchStr);
  
  /* Find param, then arrow */
  while (i<StrSize) {
    if (ftell(fp) >= FileSize) break;
    ReadStr[i] = fgetc(fp);
    if (ReadStr[i] == SearchStr[i])
      i++;
    else
      i=0;
    if (i==StrSize && !StrFound) {
      strcpy(SearchStr,"=>");
      StrFound=1;
      i=0;
      StrSize = strlen("=>");
    }
    else if (i==StrSize && StrFound)
      ArrowFound = 1;
  }  /* end while */

  if (ArrowFound)
    return(0); /* on success */
  else
    return(1);

} /* end function FindParam */


/*************************************************************************/
/* NAME: seekarrow()							
   DESCRIPTION: Reads characters from a file until an
                         arrow (==>) found	
  PARAMETERS:							
        fp           local     R    Pointer to file to read */		 


Int2B seekarrow(FILE *fp)

{
 char arr1,arr2,arr3;

 arr1 = getc(fp); arr2 = getc(fp); arr3 = getc(fp);
 while ( (arr1!='=')||(arr2!='=')||(arr3!='>') )
      { arr1 = arr2; arr2 = arr3; arr3 = getc(fp); }
 return(0);
}

/*************************************************************************/
/* FUNCTION : skipchar()
    DESCRIPTION: Reads characters from a file until a certain chosen 
              character ch is reached					 
    PARAMETERS:								 
        fp           local     R    Pointer to file to read		 
        ch           local     R    Character to read until   */         

Int2B skipchar(FILE *fp,char ch)
{
 char tch;       /* Store of character read in from the file */
 do
  {
   tch = getc(fp);     /* Reads character from file fp */
  }
 while ( (tch!= ch) && (tch != EOF) );  /* Reads while character is not */
					/* the end character or the EOF */
 return(0);
}  /* end function skipchar */


/*************************************************************************/
/* FUNCTION: GetString()					
   DESCRIPTION: Reads string from a file    				 
   PARAMETERS:		  
          fp           local     R    Pointer to file to read	       
        string       local     R    Char pointer to start of string    */

Int2B GetString(FILE *fp,char *string)
{
 Int2B n=-1;		        /* Counter of string postion */

 /* Get initial spaces */
 do
 	{
  	string[0] = getc(fp);      /* Reads character from file fp */
  }
 /* Reads while character is not a space*/
 while (string[n] == ' ');
 fseek(fp,-1,SEEK_CUR);

 do
  {
   string[++n] = getc(fp);      /* Reads character from file fp */
  }
 /* Reads while character is not a space, end of line, or end-of-file*/
 while ( (string[n] != ' ') && (string[n] != '\n') && (string[n] != EOF));
 string[n] = '\0';		/* Indicates end of string */
				/* the end character or the EOF */
 return(0);
}  /* end function GetString */




/*************************************************************/
/* FUNCTION: CFFT.C
    DISCRIPTION: Performs an complex FFT on a row of data
    INPUT: data[]pointer to the address of start of the array - 1
        nn -> number of columns in a row (must be a power of 2) 
        isign -> 1 forward fft 
                -1 inverse fft  */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void CFFT(double data[],Int4B nn,Int4B isign)

{
  Int4B n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2)
    {
    if (j > i)
      {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
      }  /* end if j */
    m=n >> 1;
    while (m >= 2 && j > m)
      {
      j -= m;
      m >>= 1;
      }  /* end while m loop */
    j += m;
    }  /* end for i loop */
  mmax=2;
  while (n > mmax)
    {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2)
      {
      for (i=m;i<=n;i+=istep)
	{
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
	} /* end for i loop */
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      }  /* end for m loop */
    mmax=istep;
    }  /* end while n>mmax */
}  /* end function CFFT */

#undef SWAP


/*************************************************************/
/* FUNCTION: CFFT_float.C (floating point version of CFFT)
    DISCRIPTION: Performs an complex FFT on a row of data
    INPUT: data[]pointer to the address of start of the array - 1
        nn -> number of columns in a row (must be a power of 2) 
        isign -> 1 forward fft 
                -1 inverse fft  */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void CFFT_float(float data[],Int4B nn,Int4B isign)

{
  Int4B n,mmax,m,j,istep,i;
  float wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2)
    {
    if (j > i)
      {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
      }  /* end if j */
    m=n >> 1;
    while (m >= 2 && j > m)
      {
      j -= m;
      m >>= 1;
      }  /* end while m loop */
    j += m;
    }  /* end for i loop */
  mmax=2;
  while (n > mmax)
    {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=(float)sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=(float)sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2)
      {
      for (i=m;i<=n;i+=istep)
	{
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
	} /* end for i loop */
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      }  /* end for m loop */
    mmax=istep;
    }  /* end while n>mmax */
}  /* end function CFFT */

#undef SWAP



/*****************************************************************************/
/* NAME: Cmult()      
   DESCRIPTION: Performs complex multiplication          		
   PARAMETERS:								
        z1           local           1st complex number 
        z2           local           2nd complex number 
        z                            Pointer to answer z   */          

Int2B Cmult(double z1[2],double z2[2], double *z)

{
 z[0] = z1[0]*z2[0] - z1[1]*z2[1];
 z[1] = z1[0]*z2[1] + z1[1]*z2[0];

 return(0);
}



/**************************************************************************/
/* NAME: IQ2power()                     
    DESCRIPTION: Calculates power from complex array 
    PARAMETERS: Note 'size' is the no. of complex pts. */

Int2B IQ2power(double pow_array[],double iq_array[], Int4B size)

{
  Int4B i,indxi=0,indxq=1;

  for (i=0;i<size;i++)
     {
     pow_array[i]=iq_array[indxi]*iq_array[indxi]+
              iq_array[indxq]*iq_array[indxq];
     indxi++;indxi++;indxq++;indxq++; 
     }   
   return(0);
 
 }  /* end function IQ2Power */


/**********************************************************************/
/* NAME: CorrelMaxIndex()
   DESCRIPTION: Finds the index shift corresponding to the maximum 
   value of the cyclic correlation between two +ve, real, and 
   equi-sized arrays. The returned index shift at the max is +ve 
   where +1 indicates a shift by one index of array 2 (i.e.
   array2[0] becomes  array2[1], ..., array2[size-1] --> array2[0]) */

Int4B CorrelMaxIndex(double array1[],double array2[],Int4B size)

{
  Int4B max_index=-999,indx1,indx2,shift=0;
  double sum,max_sum=0;
 
  while (shift<size)  /* Repeat for each possible shift */
    {
    sum = 0;
    indx1 = 0;
    indx2 = size-shift;
    while (indx1<size) /* Repeat for all points in array 1 */
      {
      if (indx2 == size) indx2=0; 
      sum += array1[indx1++]*array2[indx2++];
      }
    if (sum>max_sum) { max_sum=sum; max_index=shift; }
    shift++;
    }   
  
  return (max_index);
}  /* end function CorrelMaxIndex */


/***********************************************************************/
/* FUNCTION ConvertLong */
/* function to reverse byte order of 4 byte int */

Int4B ConvertLong(Int4B data)
{
  Int4B temp;
  temp = ( ((data>>24)& 0x000000ff) | ((data>>8) & 0x0000ff00) |
	   ((data<<8) & 0x00ff0000) | ((data<<24) & 0xff000000) );
  return (temp);
} /* end function ConvertLong */



/***********************************************************************/
/* FUNCTION FindDateTime */
/* Finds current date and time from system clock in format
    YYYYMMDD and hhmmssdd. Need to include <time.h> somewhere  */

void FindDateTime(char DateNow[], char TimeNow[])
{
   Int4B SecCount=0,SecNow,incr;
   Int2B year=1970, month=1, day=1,
	 hour=0, min=0, sec=0, dd=0;

   SecNow = (Int4B)time(NULL);

   /* find year */
   while (SecCount < SecNow)
     {
     if (year%4 == 0 ) incr = 366*24*3600;
     else incr = 365*24*3600;

     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; year++; }
     }

   /* find month */
   while (SecCount < SecNow)
     {
     if (month==1 || month==3 || month==5 || month==7
	 || month==8 || month==10 || month==12)
       incr = 31*24*3600;
     else if (month == 2)
       {
       if (year%4 == 0) incr = 29*24*3600;
       else incr = 28*24*3600;
       }
     else
       incr = 30*24*3600;

     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; month++; }
     }  /* end while loop for month */

   /* find day */
   while (SecCount < SecNow)
     {
     incr = 24*3600;
     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; day++; }
     }

   /* find hour */
    while (SecCount < SecNow)
     {
     incr = 3600;
     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; hour++; }
     }

   /* find min */
   while (SecCount < SecNow)
     {
     incr = 60;
     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; min++; }
     }

   /* find sec */
   while (SecCount < SecNow)
     {
     incr = 1;
     if (SecCount+incr >= SecNow) break;
     else { SecCount += incr; sec++; }
     }

   sprintf(DateNow,"%04d%02d%02d",year,month,day);
   sprintf(TimeNow,"%02d%02d%02d%02d",hour,min,sec,dd);

   return;

}  /* end FindDateTime function */


/**********************************************************************/
Int4B BigEndWriteInt4B(FILE *fp,Int4B value)
{
  int i;
  unsigned char unchar;

  for (i=3;i>=0;i--) {
    unchar = (unsigned char)(( value>>(8*i)) & 0x000000FF );
    fprintf(fp,"%c",unchar); }

  return 4;  /* bytes written */
}


/**********************************************************************/
Int2B BigEndWriteInt2B(FILE *fp,Int2B value)
{
  int i;
  unsigned char unchar;

  for (i=1;i>=0;i--) {
    unchar = (unsigned char)(( value>>(8*i)) & 0x00FF );
    fprintf(fp,"%c",unchar); }

  return 2;  /* bytes written */
}

/**********************************************************************/
Int4B BigEndReadInt4B(FILE *fp)
{
  unsigned char InByte[4];
  Int4B Value;

  fread(InByte,sizeof(unsigned char),4,fp);
  Value = (( (Int4B)InByte[0]<<24 ) & 0xff000000) |
          (( (Int4B)InByte[1]<<16 ) & 0x00ff0000) |
          (( (Int4B)InByte[2]<<8 ) & 0x0000ff00) |
          ( (Int4B)InByte[3] & 0x000000ff );

  return Value;
}

/**********************************************************************/
Int2B BigEndReadInt2B(FILE *fp)
{
  unsigned char InByte[2];
  Int2B Value;

  fread(InByte,sizeof(unsigned char),2,fp);
  Value = (( (Int2B)InByte[0]<<8 ) & 0xff00) |
          ( (Int2B)InByte[1] & 0x00ff );

  return Value;
}

/***********************************************************************/
/* FUNCTION median */
/* function to sort the input array and to return the median value */
/* n is the number of elements of array */

double median(double array[], Int4B n)
{
  Int4B i, j;
  double a;
  for (j=1; j<n; j++)  {
     a = array[j];
     i = j - 1;
     while ((i >= 0) && (array[i] > a)) {
        array[i + 1] = array[i];
        i -= 1;
     }
     array[i + 1] = a;
   }    
  return (array[(n - 1) / 2]);
}

/***********************************************************************/
/* FUNCTION nextPowerOfTwo */
/* Returns next power of two >= min value supplied */

Int4B nextPowerOfTwo(Int4B MinSize)
{
  Int4B nextPower = 1;
  while (nextPower < MinSize) {
    nextPower *= 2;
  }
  return (nextPower);
}

