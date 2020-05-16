/*
 * This routine now supports GRIB definition 1. According to an informal talk
 * with ECMWF this output is generated as subcenter 255 from ECMWF. This is
 * different to the former version. Additional to prior versions Local tables
 * are now supported. It is assumed that the standard table has the number 128.
 * Additional tables can be build in. As well the parameter century is supported.
 * There is no necessity to start runs by year 1. True dates can be used as long as
 * the years are > 0.
 * All subroutines connected to the GRIB output can be found in MODULE mo_grib. 
 *
 * Authors:
 *
 * Arno Hellbach, DKRZ, April 1996, original version (last update)
 * Luis Kornblueh, MPI, November 1998, extension for GRIB edition 1
 * Uwe Schulzweida, MPI, March 2000, EMOS compatible results
 *
 */

#include "util_fortran.h" 

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <float.h>

#if defined(FORTRANCAPS)
#define codegb5_ CODEGB5
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define codegb5_ codegb5__
#elif !defined(FORTRANUNDERSCORE)
#define codegb5_ codegb5
#endif

int GRB_Debug = 0;

typedef unsigned char CHAR;

#define NINT(x) ((x) < 0 ? (int)((x)-.5) : (int)((x)+.5))

#define TRUE  1
#define FALSE 0

#define SIGN 0x800000

#define Put1Byte(Value)  (IGrib[z++] = (Value))
#define Put2Byte(Value) ((IGrib[z++] = (Value) >>  8),(IGrib[z++] = (Value)))
#define Put3Byte(Value) ((IGrib[z++] = (Value) >> 16),(IGrib[z++] = (Value) >> 8),(IGrib[z++] = (Value)))

long util_pack(LONG *up, char *cp, long bc, long tc);

static INT    LocalTable     ;
static INT    CenterID       ;
static INT    ModelID        ;
static INT    GridDefinition ;
static INT    Block2Included ;
static INT    Code           ;
static INT    LevelType      ;
static INT    Level1         ;
static INT    Level2         ;
static INT    Year           ;
static INT    Month          ;
static INT    Day            ;
static INT    Hour           ;
static INT    Minute         ;
static INT    TimeUnit       ;
static INT    Time1          ;
static INT    Time2          ;
static INT    TimeRangeFlag  ;
static INT    AverageIncluded;
static INT    AverageMissing ;
static INT    Century        ;
static INT    SubCenter      ;
static INT    DecimalScaling ;

static INT    Representation ;
static INT    Truncation     ;

static INT    Longitudes     ;
static INT    Latitudes      ;
static INT    lon0, lon1     ;
static INT    lat0, lat1     ;
static INT    Increment      ;

static INT    EditionNumber  ;

static int z;
static int grib_length;

int    One = 1;
int    BitsPerInt = sizeof(int) * 8;
int    Debug = 0;

static LONG *IGrib;

static INT *ivct = NULL;



#include <stdio.h>
#include <stdarg.h>

int _ExitOnError   = 1;	/* If set to 1, exit on error       */
int _Verbose = 1;	/* If set to 1, errors are reported */
int _Debug   = 0;       /* If set to 1, debugging           */


void
SysError(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  perror("system error message ");
	
  exit(1);
}

void
Error(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "Error (%s) : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);

  if ( _ExitOnError ) exit(1);
}

void
Warning(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

  if ( _Verbose )
    {
       fprintf(stderr, "Warning (%s) : ", caller);
      vfprintf(stderr, fmt, args);
       fprintf(stderr, "\n");
    }

  va_end(args);
}

void
Message(const char *caller, const char *fmt, ...)
{
  va_list args;
	
  va_start(args, fmt);

   fprintf(stderr, "%-18s : ", caller);
  vfprintf(stderr, fmt, args);
   fprintf(stderr, "\n");

  va_end(args);
}

double
decfp2(int kexp, int kmant)
{
  /*

    Purpose:
    --------

    Convert GRIB representation of a floating point
    number to machine representation.

    Input Parameters:
    -----------------

    kexp    - 8 Bit signed exponent.
    kmant   - 24 Bit mantissa.

    Output Parameters:
    ------------------

    Return value   - Floating point number represented
                     by kexp and kmant.

    Method:
    -------

    Floating point number represented as 8 bit exponent
    and 24 bit mantissa in integer values converted to
    machine floating point format.

    Externals:
    ----------

    None.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Rewritten from DECFP, to conform to programming standards.
    Sign bit on 0 value now ignored, if present.
    If using 32 bit reals, check power of 16 is not so small as to
    cause overflows (underflows!); this causes warning to be given
    on Fujitsus.

    Author:
    -------

    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130
     
  */

  static char func[] = "decfp2";
  double pval;
  static int iexp, isign;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Convert value of 0.0. Ignore sign bit.              */
  /* ----------------------------------------------------------------- */

  if ( GRB_Debug ) Message(func, "KEXP = %d  KMANT = %d", kexp, kmant);
  /*
  if ( (kexp == 128 || kexp == 0) && kmant == 0 )
  */
  if ( (kexp == 128) || (kexp == 0) || (kexp == 255) )
    {
      pval = 0.0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  /*  Sign of value. */

  iexp = kexp;
  isign = 1;

  if ( iexp >= 128 )
    {
      iexp += -128;
      isign = -1;
    }

  /*  Decode value. */

  pval = isign * pow(2.0, -24.0) * kmant * pow(16.0, (double)(iexp - 64));

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return to calling routine.                           */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug ) Message(func, "Returned value = %d", pval);

  return (pval);
} /* decfp2 */

void
confp3(double pval, int *kexp, int *kmant, int kbits, int kround)
{
  /*

    Purpose:
    --------

    Convert floating point number from machine
    representation to GRIB representation.

    Input Parameters:
    -----------------

       pval    - Floating point number to be converted.
       kbits   - Number of bits in computer word.
       kround  - Conversion type.
                 0 , Closest number in GRIB format less than
                     original number.
                 1 , Closest number in GRIB format to the
                     original number (equal to, greater than or
                     less than original number).

    Output Parameters:
    ------------------

       kexp    - 8 Bit signed exponent.
       kmant   - 24 Bit mantissa.

    Method:
    -------

    Floating point number represented as 8 bit signed
    exponent and 24 bit mantissa in integer values.

    Externals.
    ----------

    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    WMO Manual on Codes re GRIB representation.

    Comments:
    ---------

    Routine aborts if an invalid conversion type parameter
    is used or if a 24 bit mantissa is not produced.

    Author:
    -------
     
    John Hennessy   ECMWF   18.06.91

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

  */

  static char func[] = "confp3";
  double zval, rpowref;
  static double zref, zeps;
  static int iexp, isign;
  static int iround;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1 . Initialise                                          */
  /* ----------------------------------------------------------------- */

  /*  Check conversion type parameter. */

  iround = kround;
  if ( iround != 0 && iround != 1 )
    {
      Error(func,"Invalid conversion type = %d", iround);

      /*  If not aborting, arbitrarily set rounding to 'up'. */
     iround = 1;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 2 . Convert value of zero.                              */
  /* ----------------------------------------------------------------- */

  if ( ! fabs(pval) > 0.0 )
    {
      *kexp  = 0;
      *kmant = 0;
      iexp   = 0;
      isign  = 0;
      goto LABEL900;
    }

  /* ----------------------------------------------------------------- */
  /*   Section 3 . Convert other values.                               */
  /* ----------------------------------------------------------------- */

  zeps = 1.0e-12;
  if ( kbits == 32 ) zeps = 1.0e-8;
  zref = pval;

  /*  Sign of value. */

  isign = 0;
  if ( zref < 0.0 )
    {
      isign = 128;
      zref  = - zref;
    }

  /*  Exponent. */

  iexp = (int) (log(zref) * (1.0 / log(16.0)) + 64.0 + 1.0 + zeps);

  if ( iexp < 0   ) iexp = 0;
  if ( iexp > 127 ) iexp = 127;

  rpowref  = zref / pow(16.0, (double)(iexp - 70));

  /*  Mantissa. */

  if ( iround == 0 )
    {

      /*  Closest number in GRIB format less than original number. */
      /*    Truncate for positive numbers. */
      /*    Round up for negative numbers. */

      if ( isign == 0 )
	*kmant = (int) rpowref;
      else
	*kmant = NINT(rpowref + 0.5);
    }
  else
    {

      /*  Closest number in GRIB format to the original number   */
      /*  (equal to, greater than or less than original number). */

      *kmant = NINT(rpowref);
    }

  /*  Check that mantissa value does not exceed 24 bits. */
  /*  If it does, adjust the exponent upwards and recalculate */
  /*  the mantissa. */
  /*  16777215 = 2**24 - 1 */

  if ( *kmant > 16777215 )
    {

    LABEL350:

      ++iexp;

      /*  Check for exponent overflow during adjustment  */

      if ( iexp > 127 )
	{
          Message(func, "Exponent overflow");
          Message(func, "Original number = %30.20f", pval);
          Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
		  isign, iexp, kmant);

	  Error(func, "Exponent overflow");

	  /*  If not aborting, arbitrarily set value to zero  */

          Message(func, "Value arbitrarily set to zero.");
          kexp  = 0;
          kmant = 0;
          iexp  = 0;
          isign = 0;
          goto LABEL900;
	}

      rpowref  = zref / pow(16.0, (double)(iexp - 70));

      if ( iround == 0 )
	{

	  /*  Closest number in GRIB format less than original number. */
	  /*  Truncate for positive numbers. */
	  /*  Round up for negative numbers. */

	  if ( isign == 0 )
	    *kmant = (int) rpowref;
	  else
	    *kmant = NINT(rpowref + 0.5);
	}
      else
	{

	  /*  Closest number in GRIB format to the original number */
	  /*  (equal to, greater or less than original number). */

	  *kmant = NINT(rpowref);
	}

      /*  Repeat calculation (with modified exponent) if still have */
      /*  mantissa overflow. */

      if ( *kmant > 16777215 ) goto LABEL350;
    }

  /*  Add sign bit to exponent. */

  *kexp = iexp + isign;

  /* ----------------------------------------------------------------- */
  /*   Section 9. Return                                               */
  /* ----------------------------------------------------------------- */

LABEL900:

  if ( GRB_Debug )
    {
      Message(func, "Conversion type parameter = %4d", kround);
      Message(func, "Original number = %30.20f", pval);

      zval = decfp2(*kexp, *kmant);

      Message(func, "Converted to      %30.20f", zval);
      Message(func, "Sign = %3d, Exponent = %3d, Mantissa = %12d",
	      isign, iexp, *kmant);
    }

  return;
} /* confp3 */

void
ref2ibm(double *pref, int kbits)
{
  /*

    Purpose:
    --------

    Code and check reference value in IBM format

    Input Parameters:
    -----------------

    pref       - Reference value
    kbits      - Number of bits per computer word.

    Output Parameters:
    ------------------

    pref       - Reference value

    Method:
    -------

    Codes in IBM format, then decides to ensure that reference 
    value used for packing is not different from that stored
    because of packing differences.

    Externals.
    ----------

    confp3    - Encode into IBM floating point format.
    decfp2    - Decode from IBM floating point format.

    Reference:
    ----------

    None.

    Comments:
    --------

    None.

    Author:
    -------

    J.D.Chambers     ECMWF      17:05:94

    Modifications:
    --------------

    Uwe Schulzweida   MPIfM   01/04/2001

    Convert to C from EMOS library version 130

  */

  static char func[] = "ref2ibm";
  static int itrnd;
  static int kexp, kmant;
  static double ztemp, zdumm;
  extern int GRB_Debug;

  /* ----------------------------------------------------------------- */
  /*   Section 1. Convert to and from IBM format.                      */
  /* ----------------------------------------------------------------- */

  /*  Convert floating point reference value to IBM representation. */

  itrnd = 1;
  zdumm = ztemp = *pref;
  confp3(zdumm, &kexp, &kmant, kbits, itrnd);

  if ( kexp == 0 && kmant == 0 ) return;

  /*  Set reference value to that actually stored in the GRIB code. */

  *pref = decfp2(kexp, kmant);

  /*  If the nearest number which can be represented in */
  /*  GRIB format is greater than the reference value,  */
  /*  find the nearest number in GRIB format lower      */
  /*  than the reference value.                         */

  if ( ztemp < *pref )
    {
      /*  Convert floating point to GRIB representation */
      /*  using truncation to ensure that the converted */
      /*  number is smaller than the original one.      */

      itrnd = 0;
      zdumm = *pref = ztemp;
      confp3(zdumm, &kexp, &kmant, kbits, itrnd);

      /*  Set reference value to that stored in the GRIB code. */

      *pref = decfp2(kexp, kmant);

      if ( ztemp < *pref )
	{
	  if ( GRB_Debug )
	    {
	      Message(func, "Reference value error.");
	      Message(func, "Notify Met.Applications Section.");
	      Message(func, "ZTEMP = ", ztemp);
	      Message(func, "PREF = ", pref);
	    }
	  *pref = ztemp;
	}
    }

  return;
} /* ref2ibm */


void split(int Exponent,int Mantissa,INT *i1,INT *i2,INT *i3,INT *i4)
{
   *i1 =  Exponent;
   *i2 =  Mantissa >> 16;
   *i3 = (Mantissa >>  8) & 255;
   *i4 =  Mantissa        & 255;
}

void PutnZero(int n)
{
   int i;
  
   for (i = z; i < z+n; i++) IGrib[i] = 0;
   z += n;
}

void Put1Real(REAL Value)
{
   int Exponent, Mantissa;
  
   confp3(Value, &Exponent, &Mantissa, BitsPerInt, One);
   Put1Byte (Exponent);
   Put3Byte (Mantissa);
}

/* GRIB block 0 - indicator block */

void PutGribBlock0(void)
{
  IGrib[0] = 'G';
  IGrib[1] = 'R';
  IGrib[2] = 'I';
  IGrib[3] = 'B';

  /* 
   * IGrib[4]-IGrib[6] contains full length of grib record. 
   * included before finished CODEGB
   */

  z = 7;   
  Put1Byte(EditionNumber); 
  z = 8;
}

/* GRIB block 5 - end block */

void PutGribBlock5(void)
{
  IGrib[z++] = '7';
  IGrib[z++] = '7';
  IGrib[z++] = '7';
  IGrib[z++] = '7';
  grib_length = z;
  while (z & 7) IGrib[z++] = 0;
}

/* GRIB block 1 - product definition block. */

void PutGribBlock1()
{
  Put3Byte(            28); /*  0 Length of Block 1        */
  Put1Byte(LocalTable    ); /*  3 Local table number       */
  Put1Byte(CenterID      ); /*  4 Identification of Centre */
  Put1Byte(ModelID       ); /*  5 Identification of model  */
  Put1Byte(GridDefinition); /*  6 Grid definition          */
  Put1Byte(Block2Included); /*  7 Block 2 included         */
  Put1Byte(Code          ); /*  8 Parameter Code           */
  Put1Byte(LevelType     ); /*  9 Type of level            */
  Put1Byte(Level1        ); /* 10 Level1                   */
  Put1Byte(Level2        ); /* 11 Level2                   */
  Put1Byte(Year          ); /* 12 Year of Century          */
  Put1Byte(Month         ); /* 13 Month                    */
  Put1Byte(Day           ); /* 14 Day                      */
  Put1Byte(Hour          ); /* 15 Hour                     */
  Put1Byte(Minute        ); /* 16 Minute                   */
  Put1Byte(TimeUnit      ); /* 17 Time unit = day          */
  Put1Byte(Time1         ); /* 18 Time 1                   */
  Put1Byte(Time2         ); /* 19 Time 2                   */
  Put1Byte(TimeRangeFlag ); /* 20 Timerange flag           */
  Put2Byte(             0); /* 21 Average 1/2              */
  Put1Byte(             0); /* 23 Missing in average       */
  Put1Byte(Century       ); /* 24 Century                  */
  Put1Byte(SubCenter     ); /* 25 Subcenter                */
  Put2Byte(DecimalScaling); /* 26 Units decimal scale      */
}

/* GRIB block 2 - grid description block */

void PutGribBlock2(REAL vcts[], INT vctlen)
{
  int Exponent, Mantissa;
  INT i, nv, pv;
  INT BlockLength = 32;
  
  if (LevelType == 109) {
    nv = vctlen;
    pv = BlockLength+1;
    BlockLength += vctlen << 2;
  } else {
    nv = 0;
    pv = 255;
  }

  Put3Byte(BlockLength    ); /*  0- 2 Length of Block 2 Byte 0                */
  Put1Byte(nv             ); /*  3    Number of vertical coordinate parameter */
  Put1Byte(pv             ); /*  4    octet number of location of VCT 
                                   or location of the list of numbers per row 
			           or 255 if neither is present */
  Put1Byte(Representation ); /*  5    Gauss=4 ,Spectral = 50   */
  
  if (Representation == 50)
    {
      Put2Byte(Truncation    ); /*  6- 7 Pentagonal resolution J  */
      Put2Byte(Truncation    ); /*  8- 9 Pentagonal resolution K  */
      Put2Byte(Truncation    ); /* 10-11 Pentagonal resolution M  */
      Put1Byte(             1); /* 12    Representation type      */
      Put1Byte(             1); /* 13    Representation mode      */
      PutnZero(            18); /* 14-31 reserved                 */
    }
  else
    {
      if (lat0 < 0) {lat0 = -lat0; lat0 = (lat0 ^ SIGN);}    
      if (lon0 < 0) {lon0 = -lon0; lon0 = (lon0 ^ SIGN);}     
      if (lat1 < 0) {lat1 = -lat1; lat1 = (lat1 ^ SIGN);}  
      if (lon1 < 0) {lon1 = -lon1; lon1 = (lon1 ^ SIGN);}     
      
      Put2Byte(Longitudes    ); /*  6- 7 Longitudes               */
      Put2Byte(Latitudes     ); /*  8- 9 Latitudes                */
      Put3Byte(lat0          ); /* 10-12 Latitude  of Origin      */
      Put3Byte(lon0          ); /* 13-15 Longitude of Origin      */
      Put1Byte(           128); /* 16    Resolution flag          */
      Put3Byte(lat1          ); /* 17-19 Latitude  of Extreme     */
      Put3Byte(lon1          ); /* 20-22 Longitude of Extreme     */
      Put2Byte(Increment     ); /* 23-24 i - direction increment  */
      Put2Byte((Latitudes/2) ); /* 25-26 Latitudes Pole->Equator  */
      Put1Byte(             0); /* 27    Scanning mode            */
      PutnZero(             4); /* 28-31 reserved                 */
    }
  
  if (LevelType == 109)
    {
      if (ivct == NULL)
	{
	  ivct = (INT *) malloc((size_t) vctlen * 4 * sizeof(INT));
	  for (i = 0; i < vctlen; ++i)
	    {
              confp3(vcts[i], &Exponent, &Mantissa, BitsPerInt, One);
	      split(Exponent,Mantissa,&ivct[4*i],&ivct[4*i+1],&ivct[4*i+2],&ivct[4*i+3]);
	    }
	}
#pragma _CRI ivdep
#pragma vdir nodep 
      for (i = 0 ; i < vctlen * 4 ; ++i) IGrib[z+i] = ivct[i];
      z += vctlen * 4;
    }
}

/************************************/
/* GRIB BLOCK 4 - BINARY DATA BLOCK */
/************************************/

void PutGribBlock4(REAL data[], INT Dim, INT IBits)
{
  INT i,BlockLength,PackStart,Flag;
  INT Scale;
  LONG pval;
  INT byte_per_value;
  REAL Factor,Fmin,Fmax,ZScale;
  
  byte_per_value = IBits >> 3;

  if (Representation == 50)
    {
      PackStart   = 1;
      Flag        = 128 + 8;
      BlockLength = 16 + byte_per_value * (Dim - 1);
    }
  else
    {
      PackStart   = 0;
      Flag        = 8;
      BlockLength = 12 + byte_per_value * Dim;
    }
  
  Fmin = Fmax = data[PackStart];
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
  for (i = PackStart+1 ; i < Dim ; ++i)
    {
      if (Fmin > data[i]) Fmin = data[i];
      if (Fmax < data[i]) Fmax = data[i];
    }

  ZScale = (Fmax - Fmin) / (pow(2.0,(double)(IBits + 1)) - 1);
  if (ZScale >= 10e-100) ZScale = log(ZScale) / log(2.0) + 2.0;
  else                   ZScale = 0.0;
  
  Scale = floor(ZScale);
  
  Factor = 1;

  if (Fmax > Fmin)
    {
      if (Scale < 0) Factor = pow(2.0,(double)(-Scale));
      else           Factor = 1.0 / pow(2.0,(double)(Scale));
      if (Scale < 0) Scale = 32768 - Scale;
    }
  
  ref2ibm(&Fmin, BitsPerInt);

  Put3Byte(BlockLength);    /*  0-2 Length of Block 4        */
  Put1Byte(Flag       );    /*  3   Flag & Unused bits       */
  Put2Byte(Scale      );    /*  4-5 Scale factor             */
  Put1Real(Fmin       );    /*  6-9 Reference value          */
  Put1Byte(IBits      );    /* 10   Packing size             */
  
  if (PackStart) Put1Real(data[0]);
  
  if      (IBits ==  8)
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for (i = PackStart ; i < Dim; ++i)
	{
          IGrib[z] = (data[i] - Fmin) * Factor + 0.5;
          ++z;
	}
    }
  else if (IBits == 16)
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for (i = PackStart ; i < Dim; ++i)
	{
          pval = (data[i] - Fmin) * Factor + 0.5;
          IGrib[z  ] = pval >>  8;
          IGrib[z+1] = pval;
          z += 2;
	}
    }
  else if (IBits == 24)
    {
#if   defined (CRAY)
#pragma _CRI ivdep
#elif defined (SX)
#pragma vdir nodep
#elif defined (__uxp__)
#pragma loop novrec
#endif
      for (i = PackStart ; i < Dim ; ++i)
	{
          pval = (data[i] - Fmin) * Factor + 0.5;
          IGrib[z  ] =  pval >> 16;
          IGrib[z+1] =  pval >>  8;
          IGrib[z+2] =  pval;
          z += 3;
	}
    }
  else
    {
      printf(" Unimplemented packing factor %d\n",(int)IBits);
      exit(1);
    }

  Put1Byte(0);              /*  Fillbyte                     */
}


FORTRAN_CALL
void codegb5_(REAL *Pdata, INT *datalen,
	      INT *gbits,  INT *wbits,
	      INT *Pisec1, INT *Pisec2,
	      REAL *Pvcts, INT *vctslen,
	      INT *ogrib,  INT *griblen, INT *gribwords,
	      INT *round,  INT *ierr)
{
  INT zall;
  INT *isec1=Pisec1, *isec2=Pisec2;
  REAL *vcts=Pvcts, *data=Pdata;
  
  IGrib = (LONG *) malloc((size_t) (*datalen * 4 + 500) * sizeof(LONG));

  EditionNumber  = 1; 
  
  LocalTable     = *isec1++; 
  CenterID       = *isec1++;
  ModelID        = *isec1++;
  GridDefinition = *isec1++;
  Block2Included = *isec1++;
  Code           = *isec1++;
  LevelType      = *isec1++;
  /*  Level2         = *isec1++;       
      Level1         = *isec1++; 
      Level1         = Level1 << 8;
      Level2         = Level2 & 255; */
  /*
  Level1         = *isec1++;
  Level2         = *isec1++;
  */
  Level2         = *isec1++;
  Level1         = *isec1++;
  if (Level2 > 255) {
      Level1         = Level2 >> 8;
      Level2         = Level2 & 255;
  }
  else {
     Level1 = 0;
  }
  Year           = *isec1++;
  Month          = *isec1++;
  Day            = *isec1++;
  Hour           = *isec1++;
  Minute         = *isec1++;
  TimeUnit       = *isec1++;
  Time1          = *isec1++;
  Time2          = *isec1++;
  TimeRangeFlag  = *isec1++;
  isec1 += 2;
  Century        = *isec1++;
  SubCenter      = *isec1++;
  DecimalScaling = *isec1++;
  
  Representation = *isec2++;
  if (Representation == 50)
    {
      Truncation = *isec2++;
    }
  else
    {
      Longitudes = *isec2++;
      Latitudes  = *isec2++;
      lat0       = *isec2++;
      lon0       = *isec2++;
      *isec2++;
      lat1       = *isec2++;
      lon1       = *isec2++;
      Increment  = *isec2++;
    }
  
  PutGribBlock0();
  PutGribBlock1();
  PutGribBlock2(vcts,*vctslen);
  PutGribBlock4(data,*datalen,*gbits);
  PutGribBlock5();
  
  *griblen   = z;
  *gribwords = (z << 3) / *wbits;

  if (*gribwords > *datalen+1000) {
     fprintf(stderr,"codegb5 error: griblen = %d\n", *gribwords);
     fprintf(stderr,"codegb5 error: datalen = %d\n", *datalen+1000);
  }

  zall = z;

  z = 4;

  Put3Byte(grib_length);

  z = zall;

#if defined (CRAY)
      _pack(IGrib, (char *)ogrib, z, -1L);
#else
  util_pack(IGrib, (char *)ogrib, z, -1L);
#endif
  
  free(IGrib);
}

/* pack 8-bit bytes from 64-bit words to a packed buffer */

long util_pack(LONG *up, char *cp, long bc, long tc)
{
  char *cp0;
  LONG upi, *up0, *ip0, *ip1, *ip2, *ip3, *ip4, *ip5, *ip6, *ip7;
  long head, trail, inner, i, j;
  
  /* Bytes until first word boundary in destination buffer */

  head = ( (long ) cp ) & 7;
  if ( head != 0 ) head = 8 - head;

  inner = bc - head;

  /* Trailing bytes which do not make a full word */

  trail = inner & 7;

  /* Number of bytes/words to be processed in fast loop */

  inner = inner - trail ; inner = inner/8;

  ip0 = up + head;
  ip1 = ip0 + 1;
  ip2 = ip0 + 2;
  ip3 = ip0 + 3;
  ip4 = ip0 + 4;
  ip5 = ip0 + 5;
  ip6 = ip0 + 6;
  ip7 = ip0 + 7;
  up0 = (LONG *) (cp + head);

  /* Here we should process any bytes until the first word boundary 
   * of our destination buffer 
   * That code is missing so far  because our output buffer is 
   * word aligned by FORTRAN 
   */

  j = 0;

#ifdef WORDS_BIGENDIAN
#pragma _CRI ivdep
#pragma vdir nodep
#ifdef __uxpch__
#pragma loop novrec
#endif
  for ( i = 0 ; i < inner ; i++ ) {
    upi =             (   ip0[j]         << 56 ) 
                   |  ( ( ip1[j] & 255 ) << 48 )
                   |  ( ( ip2[j] & 255 ) << 40 )
                   |  ( ( ip3[j] & 255 ) << 32 )
                   |  ( ( ip4[j] & 255 ) << 24 ) ;
    up0[i] = upi   |  ( ( ip5[j] & 255 ) << 16 )
                   |  ( ( ip6[j] & 255 ) <<  8 )
                   |    ( ip7[j] & 255 ) ;
    j = j + 8;
  }
#else
  for ( i = 0 ; i < inner ; i++ ) {
    upi =             (   ip7[j]         << 56 ) 
                   |  ( ( ip6[j] & 255 ) << 48 )
                   |  ( ( ip5[j] & 255 ) << 40 )
                   |  ( ( ip4[j] & 255 ) << 32 )
                   |  ( ( ip3[j] & 255 ) << 24 ) ;
    up0[i] = upi   |  ( ( ip2[j] & 255 ) << 16 )
                   |  ( ( ip1[j] & 255 ) <<  8 )
                   |    ( ip0[j] & 255 ) ;
    j = j + 8;
  }
#endif

  cp0 = (char *) ( up0 + inner );
  if ( trail > 0 ) {
    up0[inner] = 0;
    for ( i = 0 ; i < trail ; i ++ ) {
      *cp0 = (char ) ip0[8*inner+i];
      cp0++;
    }
  }

  if ( tc != -1 ) {
    bc++;
    *cp0 = (char) tc;
  }

  return(bc);
}
