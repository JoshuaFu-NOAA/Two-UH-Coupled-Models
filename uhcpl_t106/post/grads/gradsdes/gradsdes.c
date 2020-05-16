#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gradsdes.h"
/* #include <console.h> */
/* #include <sioux.h> */

/*  Routine to map gridded GRIB files.
    The size and scaling of the grids is essentially hardcoded */

/*  hipnt info:  1 - number of times in file
                 2 - number of records per time
                 0 - version number (1)
                 3 - number of records per grid
                 (used for verification on open)

    hfpnt info:  None

    intpnt info (for each record) :
                 0 - position of start of data in file
                 1 - position of start of bit map in file
                 2 - number of bits per data element

    fltpnt info :
                 0 - decimal scale factor for this record
                 1 - binary scale factor
                 2 - reference value

*/

#define gfile intape

/* Date/time structure.         */

struct dt {
  int yr;
  int mo;
  int dy;
  int hr;
  int mn;
};

struct gaindx {
  int type;                 /* Indexing file type             */
  int hinum;                /* Number of ints in header       */
  int hfnum;                /* Number of floats in header     */
  int intnum;               /* Number of index ints (long)    */
  int fltnum;               /* Number of index floats         */
  int   *hipnt;             /* Pointer to header int values   */
  float *hfpnt;             /* Pointer to header float values */
  int   *intpnt;            /* Pointer to int index values    */
  float *fltpnt;            /* Pointer to float index values  */
} ;
struct gaindx *pindx;       /* Index Strucure if indexed file */

#define REP_REGULAR      0
#define REP_GAUSS        4
#define REP_ROTREG      10
#define REP_SPECTRAL    50

int RepGrib, ReducedGrid;

struct grhdr {
  int vers;
  int len;
  int pdslen,gdslen,bmslen,bdslen;
  int id;
  int gdsflg,bmsflg;
  int parm;
  int ltyp;
  int level;
  int l1,l2;
  struct dt dtim;
  int ftu,p1,p2,tri;
  int fcstu,fcstt;
  int cent;
  float dsf;
  int gtyp,gicnt,gjcnt,gsf1,gsf2,gsf3;
  int bnumr;
  int bpos;
  int iflg;
  float bsf;
  float ref;
  int bnum;
  int dpos;
};

struct grhdr ghdr;

int gribhdr(struct grhdr *);

int lat1,lon1,lat2,lon2,latdi,londi;

int fpos;  /* File pointer into GRIB file */
int flen;
int scanflg;
int verboseflg=0;

int IGribRec;
int parm[MAX_CODES];
int pcode[MAX_GREC];
int ltyp[MAX_CODES];
int hinum[4];
int intnum[MAX_GREC][3];
float fltnum[MAX_GREC][3];
int gtype;  /* GRIB type */
int ghinum;
int ghfnum ;
int gintnum;
int gfltnum;
int *ghipnt;
float *ghfpnt;
int *gintpnt;
float *gfltpnt;

int GribRec = 0;
INT LYREV=0;
INT GAUSS=0;
INT QUIET=0;
int Debug=0;
INT LTYP=0;
char *pTime;
char *pIncr;
char Time[30], Incr[10], *IncrKey[] = {"mn","hr","dy","mo","yr"};

#define VERSION 3.40

void
Version(void)
{
   printf("gradsdes version %3.2f\n", VERSION);
/*
   v1.0  17 May 1996: initial version
   v1.1  11 Jun 1996: use lat1 lon1 lat2 lon2
                      add lev0 to use level 0
   v1.2  24 Jul 1996: new calculation of time increment
   v1.3  22 Aug 1996: use bit map data
                      set start record from 100 to 500 bytes
   v1.4   6 Nov 1996: add option -Codesfile
                      add option -mapfile
   v1.5  18 Dec 1996: reduce ifilename only when extention .grb
   v1.6  15 Jan 1997: used ADD from Codesfile
   v1.7  20 Jan 1997: bug in gribhdr: set imax variable!!!
   v1.8  20 Feb 1997: use gaindx (bug on DEC!)
   v1.9  27 Mar 1997: levels in Pascal
   v2.0  11 Apr 1997: add ^ for filenames
   v2.1  26 Aug 1997: set MAX_GREC 99999
                      if NLAT or NLON eq 1 set DX or DY to 1
   v2.2   1 Dec 1997: xydef file                    
   v2.3   2 Dec 1997: turn on/off yrev
   v2.4   3 Dec 1997: change format -> options
   v2.5   8 Jul 1998: GDS rec from 256 to 512-3
   v2.6  17 May 1999: zdef file
   v2.7  28 Sep 1999: T319 support
   v3.0  18 Sep 2000: change output information
                      add info in descriptor file
                      use minutes (mn) for time calculation
		      skip spectral codes
   v3.1  26 Sep 2000: optimize GRIB scan
                      skip data on reduce grids
   v3.2  28 Sep 2000: Rotated regular grid support
   v3.3  30 Nov 2000: use gauss grid also for zonal means
   v3.31 21 Feb 2001: Error if GRIB not found only for len > 200
   v3.32 23 Mar 2001: Calculate delta x and y if it is 0.
   v3.33 17 Apr 2001: gaussgrid for gaussgrids!
   v3.34 14 May 2001: write DATA and INDEX filename without path
   v3.35 12 Jul 2001: correct xdef and ydef for 1x1 grid
   v3.36 02 Jun 2001: correct ydef for non gaussian grids
   v3.37 19 Jun 2002: remove default yrev on
   v3.40 ?? Sep 2002: Calculate delta x and y from last and first x, y if 0
                      skip codes with different leveltype
		      bug fix layer level -> level1
   v3.x               set default Codesfile to "codes"
*/
}

void
Help(void)
{
   Version();

   printf("usage: gradsdes [-options] [ifile]\n");
   printf("with:\n");
   printf("   -c[ontrolfile] [cfile]  descriptor filename\n");
   printf("                           default: -c ifile.ctl\n");
   printf("   -C[odesfile] [Cfile]    variables abbreviation and description\n");
   printf("                           see file: /pf/m/m214003/bin/codes\n");
/*   printf("                           default: -C codes\n");*/
   printf("   -g[auss] on/off         turn on/off gaussian grid\n");
   printf("                           default: -g on  for T21,T30,T42,T63,T106,T213 and T319 grids\n");
   printf("                                    -g off for all other grids\n");
   printf("   -h[elp]                 give this help\n");
   printf("   -i[nputfile] ifile      gribfile\n");
   printf("   -I[ncrement] vvkk       time increment\n");
   printf("   -m[apfile] mfile        grib map filename\n");
   printf("   -T[ime] start-time      the starting date/time value\n");
   printf("   -q[uiet]                quiet mode\n");
   printf("   -V[ersion]              display version number\n");
   printf("   -xydef file             xydef file\n");
   printf("   -zdef file              zdef file\n");
   printf("   -y[rev] on/off          turn on/off yrev\n");
}

void Abort(char *errtext)
{
   fprintf (stderr, errtext);
    printf (        errtext);
    exit (1);
}

INT get_codes(INT code, char C_code[], char C_CODE[])
{
    INT level;
    char line[256],ch;
    char *endp[256];
    
    ADD=0.0;
    MULT=1.0;
    NTIME=0;
    rewind(CodesTape);
    while (fgets(line,256-1,CodesTape))
    {
       if (code == strtol(line,endp,10))
       {
          level=strtol(*endp, endp,10);
          while (**endp == ' ') (*endp)++;
          while ((ch=**endp) != ' ') 
          {
             *C_code=ch;
             C_code++;
             (*endp)++;
          }
          *C_code=0;
          ADD=strtod(*endp, endp);
          MULT=strtod(*endp, endp);
 /*         NTIME=strtol(*endp, endp,10);*/
          while (**endp == ' ') (*endp)++;
 /*         printf("%s",*endp);*/
          while (*C_CODE=**endp)
          {
             if (**endp == '\n')
             {
                *C_CODE=0;
                break;
             }
             C_CODE++;
             (*endp)++;
          }
          return(0);
       }
    }
}

INT AddLevel(INT ilev, INT nlev)
{
    if (Debug) fprintf(stderr,"ilev=%d nlev%d\n",ilev,nlev);
    LEVELS[nlev-1] = ilev;
    return (1);
}


INT SetLevel(INT icode, INT nc, INT ilev, INT nlev)
{
    INT i=0;

    while (i < MAX_CODES || S_Codes[i].code == 0)
    {
       if (S_Codes[i].code == icode)
       {
          if (Debug) fprintf(stderr,"%d %d %d\n", icode, i+1, nc);
          if (i+1 != nc) 
          {
            fprintf(stderr,"code = %d leveltype = %d  %d %d\n",
		    icode, S_Codes[i].leveltype, i+1,nc);
	    if ( icode == 0 )
	      Abort("SetLevel: code is 0\n");
	    else
	      Abort("SetLevel: code allready defined on other levels\n");
          }
          if (nlev > S_Codes[nc-1].level) AddLevel(ilev, nlev);
          S_Codes[nc-1].level = nlev; 
          if (S_Codes[nc-1].level == 1 ) S_Codes[nc-1].lev0 = ilev; 
       }
       i++;
    }

    return(0);
}

INT AddCode(INT icode, INT nc, int leveltype)
{
    static INT NC;
 
    NC++;
    if (NC != nc) printf("AddCode: NC=%d ne nc=%d\n", NC, nc);
    if (nc > MAX_CODES) Abort("AddCode: Code > MAX_CODES\n");

    S_Codes[nc-1].code = icode;
    S_Codes[nc-1].level = 0;
    S_Codes[nc-1].leveltype = leveltype;
    parm[nc-1]=ghdr.parm;
    ltyp[nc-1]=ghdr.ltyp;

    return(0);
}


void gribctl()
{
   INT I_Code, I_Date, I_Time, I_Lat, I_Lon, I_Res;
   INT I_Level;
   INT I_Year, I_Month, I_Day, I_Hour, I_Minute, I_Dt;
   INT I_Head[8], I_Dim, nwrite, i;
   INT GTIME;
   REAL GDIV;
   REAL alat[MAX_Latitudes], w[MAX_Latitudes], degpi, coa[MAX_Latitudes];
   INT ii,jj;
   INT lyrev,lxdim=0,lydim=0;
   char line[256];
   int leveltype = 0;
   
   I_Head[0] = ghdr.parm;
   leveltype = ghdr.ltyp;
   if (ghdr.ltyp == 100)
     {
       I_Head[1] = 100*ghdr.level;
       LTYP=1;
     }
   else if (ghdr.ltyp == 99)
     {
       I_Head[1] = ghdr.level;
       LTYP=1;
       leveltype = 100;
     }
   else if ( ghdr.ltyp == 110 )
     {
       I_Head[1] = ghdr.l1;
     }
   else if ( ghdr.ltyp == 112 )
     {
       I_Head[1] = ghdr.l1;
     }
   else
     {
       I_Head[1] = ghdr.level;
     }
   I_Head[2] = 10000*ghdr.dtim.yr + 100*ghdr.dtim.mo + ghdr.dtim.dy;
   I_Head[3] = 100*ghdr.dtim.hr + ghdr.dtim.mn;
   I_Head[4] = ghdr.gicnt;
   I_Head[5] = ghdr.gjcnt;
   I_Head[6] = 0;
   I_Head[7] = 0;
   
/*   GTIME=nblck1[15]; */

   I_Code  = I_Head[0]; 
   I_Level = I_Head[1];
   I_Date  = I_Head[2];
   I_Time  = I_Head[3];
   I_Lon   = I_Head[4];
   I_Lat   = I_Head[5];
   /*
   if (!QUIET)
      fprintf(stderr,"%d %d %d %d %d %d\n",I_Head[0],I_Head[1],I_Head[2],I_Head[3],I_Head[4],I_Head[5]);
      */
   I_Dim = I_Lon*I_Lat;

   I_Res   = 0;

   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 32  ) I_Res = 21;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 48  ) I_Res = 31;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 64  ) I_Res = 42;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 96  ) I_Res = 63;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 128 ) I_Res = 85;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 160 ) I_Res = 106;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 320 ) I_Res = 213;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 480 ) I_Res = 319;
   if ((I_Lon == 1 || I_Lon == I_Lat*2) && I_Lat == 512 ) I_Res = 341;

   if ( I_Res && GAUSS == 0 )      GAUSS = 1;
   if ( (RepGrib == REP_GAUSS) && GAUSS == 0 )  GAUSS = 1;
   if ( I_Lat == 1 ) GAUSS = 0;

   if (I_ires++)
   {
      if (O_Lat != I_Lat || O_Lon != I_Lon)
      {
         fprintf (stderr, "Resolution must not change!\n");
         fprintf (stderr, "LAT: first=%d  now=%d\n", O_Lat, I_Lat);
         fprintf (stderr, "LON: first=%d  now=%d\n", O_Lon, I_Lon);
         exit(1);
      }
/*printf(" %d %d %d %d\n" ,S_Codes[0].code,I_Code,S_Codes[0].lev0,I_Level);*/

      if (S_Codes[0].code == I_Code && S_Codes[0].lev0 == I_Level ) 
      {
	if (Debug) {
	  fprintf(stderr,">>> %d %d %d %d %d %d\n",
		  I_Head[0],I_Head[1],I_Head[2],I_Head[3],I_Head[4],I_Head[5]),
	    fprintf(stderr,">>> %d %d %d %d\n",I_Level,O_Level,I_Code,S_Codes[0].code);
	}
/*         if (O_CCode != O_ccode)
         {        
            fprintf(stderr,"Nummer of Codes not the same as at first!\n");
            fprintf(stderr,"NCODE: first=%d  now=%d\n", O_CCode, O_ccode);
            exit(1);
         }
         else
         {*/
            O_ccode = 1;
            DATE[O_ctime] = I_Date;
            TIME[O_ctime] = I_Time;

            if (O_ctime == 1 && !QUIET)
	      {
		fprintf (stdout, "\n  Date and time of other timesteps :\n");
		fprintf (stdout, "  TIMESTEP  YYYYMMDD  HHMM\n");
	      }

            O_ctime++;

	    if (!QUIET)
	      fprintf (stdout, "  %8d  %8d  %04d\n", O_ctime, I_Date, I_Time);

            O_Time  = I_Time;
            O_Date  = I_Date;
            O_Level = I_Level;

   /*      } */
      }
      else
      {
         O_Time = I_Time;
         O_Date = I_Date;

         if (O_ctime == 1)
         {
	   if (!QUIET)
	     fprintf (stderr, "  %4d %9d\n", I_Code, I_Level);

           if (O_Code != I_Code)
           {
              O_ccode++;
              O_Code = I_Code;
              O_CCode = O_ccode;
              AddCode(O_Code, O_ccode, leveltype);
              if (O_Level && !I_Level)
              {
                 O_clevel = 0;
                 O_Level = I_Level;
              }
              else if (I_Level)
              {
                 O_clevel = 1;
                 O_Level = I_Level;
                 SetLevel(O_Code, O_ccode, O_Level, O_clevel);
              }               
           }
           else
           {
              O_clevel++;
              O_Level = I_Level;
              if (!MaxLevel) 
              {
                 O_CLevel = O_clevel;
                 SetLevel(O_Code, O_ccode, O_Level, O_clevel);
              }        
           }
         }
      }
   }
   else
   {
      O_Lat = I_Lat;
      O_Lon = I_Lon;

      if (I_Res)
	{
	  I_dx = 360./O_Lon;
	  I_dy = 180./O_Lat;
	}
      else
	{
          if (londi == 0)
	    {
	      if ( (lon2 - lon1) != 0 )
		{
		  I_dx = (lon2 - lon1)/(1000.*(O_Lon-1));
		}
	      else
		I_dx = 360./O_Lon;
	    }
	  else
	    {
	      I_dx = londi / 1000.;
	      if (I_dx > 360./O_Lon) I_dx = 360./O_Lon;
	    }

          if (latdi == 0)
	    {
	      if ( (lat2 - lat1) != 0 )
		{
		  I_dy = (lat2 - lat1)/(1000.*(O_Lat-1));
		}
	      else
		I_dy = 180./O_Lat;
	    }
	  else
	    {
	      I_dy = fabs(latdi / 1000.);
	      if (I_dy > 180./(O_Lat-1)) I_dy = 180./(O_Lat-1);
	    }
	}

      I_north = 90.-I_dy/2.;
        
      O_Level = I_Level;
      O_Code = I_Code;
      O_Time = I_Time;
      O_Date = I_Date;

      O_ccode = 1;
      O_CCode = 1;
      O_clevel = 0;
/*      if (O_Level)*/ O_clevel = 1;
      O_CLevel = O_clevel;
      O_ctime = 1;

      if (!QUIET)
	{
	  fprintf (stdout, " \n");
	  fprintf (stdout, "  Startdate    : %d\n",   I_Date);
	  fprintf (stdout, "  Starttime    : %04d\n", I_Time);
	  fprintf (stdout, "  Nlon         : %4d\n",  I_Lon);
	  fprintf (stdout, "  Nlat         : %4d\n",  I_Lat);
          if (RepGrib == REP_REGULAR)
	    fprintf (stdout, "  Grid         : Regular\n");
          else if (RepGrib == REP_GAUSS)
	    fprintf (stdout, "  Grid         : Gaussian\n");
          else if (RepGrib == REP_ROTREG)
	    fprintf (stdout, "  Grid         : Rotated regular\n");
          else
	    fprintf (stdout, "  Grid         : Unsupported\n");

	  fprintf (stdout, " \n");

	  fprintf (stdout, "  Code     Level\n");
	  fprintf (stdout, "  %4d %9d\n", I_Code, I_Level);
	}

      AddCode(O_Code, O_ccode, leveltype);
      SetLevel(O_Code, O_ccode, O_Level, O_clevel);
      DATE[0] = I_Date;
      TIME[0] = I_Time;

      fprintf(cntape,"* Generated automatically by gradsdes version %3.2f\n",VERSION);
      fprintf(cntape,"*\n");

      if (ifile[0] == '/' )
	fprintf(cntape,"DSET  %s\n",ifile);
      else
	{
	  char *filename;
	  filename = strrchr(ifile, '/');
	  if ( filename == 0 ) filename = ifile;
	  else                 filename++;	  
	  fprintf(cntape, "DSET  ^%s\n", filename);
	}

      lyrev = 0;
      
      if ((lat1 != 0 || lat2 != 0) && (lon1 != 0 || lon2 != 0)) {
         if (lat2 < lat1) lyrev = 1;
         lxdim = 1;
         lydim = 1;
      }

      if (LYREV) lyrev = 1;
    
      if (lyrev) fprintf(cntape,"OPTIONS YREV\n");

      fprintf(cntape,"UNDEF  -9.9E33\n");
      if (mfile) {
        fprintf(cntape,"DTYPE GRIB\n");
        if (ifile[0] == '/' )
          fprintf(cntape,"INDEX  %s\n",mfile);
        else
	  {
	    char *filename;
	    filename = strrchr(mfile, '/');
	    if ( filename == 0 ) filename = mfile;
	    else                 filename++;	  
	    fprintf(cntape, "INDEX  ^%s\n", filename);
	  }
      }
      if (I_Res)
         fprintf(cntape,"TITLE T%d DATA\n",I_Res);
      else
         fprintf(cntape,"TITLE %dx%d GRID DATA\n",O_Lon,O_Lat);

     if (xyfile) {
        while (fgets(line,256-1,XYTape))  {
           fputs(line,cntape);
        }
     }
     else {
      if (lxdim)
        fprintf(cntape,"XDEF %d LINEAR %f %f\n",O_Lon,(double)lon1/1000.,I_dx);
      else {
        if (O_Lon == 1) I_dx=1.;
        fprintf(cntape,"XDEF %d LINEAR 0.0 %f\n",O_Lon,I_dx);
      }

      if (GAUSS == 1) {
        gaussg_(O_Lat/2, coa, w);
        degpi = 180. / acos(-1.);
        for (ii = 1; ii <= O_Lat/2; ++ii) {
          alat[ii - 1] = -asin(coa[ii - 1]) * degpi;
          alat[O_Lat - ii] = -alat[ii - 1];
        }
        fprintf(cntape,"YDEF %d LEVELS ",O_Lat);
        jj=0;
        for (ii=1; ii <= O_Lat; ii++)  {
           fprintf(cntape,"%7.3f ",alat[ii-1]); 
           jj++;
           if (jj == 6)  {
             fprintf(cntape,"\n");
             jj=0;
             if (ii != O_Lat) fprintf(cntape,"               ");
           }
        }
        if (jj) fprintf(cntape,"\n");
      }
      else  {
        if (lydim) {
	  if (lat1 < lat2)
            fprintf (cntape, "YDEF %d LINEAR %f %f\n", O_Lat, (double)lat1/1000., I_dy);
	  else
            fprintf (cntape, "YDEF %d LINEAR %f %f\n", O_Lat, (double)lat2/1000., I_dy);
	}
        else {
          if (O_Lat == 1) I_dy = 1.;
          fprintf (cntape, "YDEF %d LINEAR %f %f\n", O_Lat, -I_north, I_dy);
        }
      }
     }
   }
/*     
   if (ontape) {
     ADD=0.0;
     MULT=1.0;
     NTIME=0;
     GDIV=1.0;
     if (CodesTape)
         get_codes(I_Code, C_code, C_CODE);
        
     if (NTIME && GTIME) GDIV = (double) GTIME / (double) NTIME ;
      for (i=0; i< I_Dim; i++)
         Field[i] = field[i] * MULT / GDIV + ADD;
      nwrite = fwrite(Field, sizeof(float), I_Dim, ontape);
      if (nwrite != I_Dim)  {
         fprintf(stderr, "ERROR:\n");
         fprintf(stderr, "can't write array with %ld elements\n", I_Dim);
         exit(1);
      }
   }
*/ 
}

void endctl()
{
     INT i,jj, maxlev=0, level=0;
     INT dlev=0,dlevold;
     INT dt=1, dt0=1, iik=0, iik0=0;
     INT imn, imns, ihh, ihhs, iyy, iyys, imm, imms, idd, idds;
     INT idmn, idhh, idmm, idyy, iddd;
     char line[256];
     char cextra[99];
     char *cmons[]={"jan","feb","mar","apr","may","jun",
                    "jul","aug","sep","oct","nov","dec"};

     for (i=0; i<O_CCode; i++)
        if (S_Codes[i].level > maxlev) maxlev = S_Codes[i].level ;

     for (i=0; i<maxlev-1; i++)
     {
        dlev = LEVELS[i+1] - LEVELS[i];
        if (i >= 1) 
        {
           if (dlev != dlevold)
           {
              level=1;
              break;
           }
        }
        dlevold=dlev;
     }

/*     if (maxlev == 1) level = 1; */

     if (zfile) {
        while (fgets(line,256-1,ZTape))  {
           fputs(line,cntape);
        }
     }
     else {
       if (level) 
       {
          fprintf(cntape,"ZDEF %d LEVELS ",maxlev);
          jj=0;
          for (i=0; i<maxlev ; i++)  {
             if (LTYP) fprintf(cntape,"%g ", (double) LEVELS[i]/100.);
             else fprintf(cntape,"%d ", LEVELS[i]);
             jj++;
             if (jj == 10)  {
               fprintf(cntape,"\n");
               jj=0;
               if (i != (maxlev-1)) fprintf(cntape,"               ");
             }
          }
          if (jj) fprintf(cntape,"\n");
       }
       else
       {  
          if (maxlev==1) { LEVELS[0]=1; dlev=1; }
          fprintf(cntape,"ZDEF %d LINEAR %d %d\n",maxlev,LEVELS[0],dlev);
       }
     }
/* 
     if (!QUIET)
        for (i=0; i<O_ctime; i++)
           fprintf(stderr,"%d %d\n", DATE[i], TIME[i]); 
*/
     ihhs =  TIME[0] / 100;
     imns =  TIME[0] - ihhs*100;
     iyys =  DATE[0] / 10000;
     imms = (DATE[0] - iyys*10000) / 100;
     idds =  DATE[0] - iyys*10000 - imms*100;
     
     if (pTime == NULL)
     {
        if (imms < 1 || imms > 12)  imms=1;
        pTime = Time;
        sprintf (Time, "%02d:%02dZ%02d%s%04d", ihhs, imns, idds, cmons[imms-1], iyys);
        
     }

     for (i=1; i<O_ctime; i++)
     {
        ihhs =  TIME[i-1] / 100;
        imns =  TIME[i-1] - ihhs*100;
        iyys =  DATE[i-1] / 10000;
        imms = (DATE[i-1] - iyys*10000) / 100;
        idds =  DATE[i-1] - iyys*10000 - imms*100;

        ihh =  TIME[i] / 100;
        imn =  TIME[i] - ihh*100;
        iyy =  DATE[i] / 10000;
        imm = (DATE[i] - iyy*10000) / 100;
        idd =  DATE[i] - iyy*10000 - imm*100;

        idmn = imn - imns;
        idhh = ihh - ihhs;
        iddd = idd - idds;
        idmm = imm - imms;
        idyy = iyy - iyys;

        if (idmn != 0)  {
           dt = idmn + (idhh + (iddd + (idmm*30 + idyy*12)*30)*24)*60;
        }
        else if (idhh != 0)  {
           dt = idhh + (iddd + (idmm + idyy*12)*30)*24;
           iik = 1;
        }
        else if (iddd != 0)  {
           dt = iddd + (idmm + idyy*12)*30;
           iik = 2;
        }
        else if (idmm != 0)  {
           dt = idmm + idyy*12;
           iik = 3;
        }
        else if (idyy != 0)  {
           dt = idyy;
           iik = 4;
        }

	if (Debug)
	  fprintf (stderr, "%d %d %d %d %d %d %d\n", ihh, imn, iyy, imm, idd, dt, iik);
   
        if (i == 1)  {
          dt0  = dt;
          iik0 = iik;
        }
        else  {
          if (((dt0 - dt) != 0) || ((iik0 - iik) != 0)) 
            if (!QUIET) fprintf (stdout, "Warning: inconsistent time increment at timestep %d\n", i+1);
        }
     }
     
     if (dt0 <= 0) dt0 = 1;
     if (pIncr == NULL)
     {
        pIncr = Incr;
        sprintf (Incr, "%d%s", dt0, IncrKey[iik0]);
     }

     if (Debug)
       fprintf (stderr, "TDEF %d LINEAR %s %s\n", O_ctime, pTime, pIncr);

     fprintf (cntape, "TDEF %d LINEAR %s %s\n", O_ctime, pTime, pIncr);

     fprintf (cntape, "VARS %d\n", O_CCode);

     for (i=0; i<O_CCode; i++)
     {
        strcpy(C_code,"");

        if (CodesTape)
           get_codes(S_Codes[i].code, C_code, C_CODE);
        
        if (!strlen(C_code))
        {
           strcpy(C_code, "code");
           sprintf(&C_code[4],"%d",S_Codes[i].code);
           strcpy(C_CODE, "CODE");
           sprintf(&C_CODE[4],"%d",S_Codes[i].code);
        }
        sprintf(cextra,"99");
        if (mfile) 
           sprintf(cextra,"%d,%d,0",parm[i],ltyp[i]);
        if (S_Codes[i].level == 1 && S_Codes[i].lev0 == 0) S_Codes[i].level = 0;
/*        printf("lev0=%d\n",S_Codes[i].lev0);*/
        fprintf(cntape,"%-8s %2d %s %s\n",  C_code, S_Codes[i].level,cextra, C_CODE);
     }
     fprintf(cntape,"ENDVARS\n");


  gtype = 1;  /* GRIB type */
  ghinum = 4;
  ghfnum = 0;
  gintnum = 3 * IGribRec;
  gfltnum = 3 * IGribRec;
  ghipnt = NULL;
  ghfpnt = NULL;
  gintpnt = NULL;
  gfltpnt = NULL;

  pindx = (struct gaindx *)malloc(sizeof(struct gaindx));
  /* Set index stuff */

  pindx->type = 1;  /* GRIB type */
  pindx->hinum = 4;
  pindx->hfnum = 0;
  pindx->intnum = gintnum;
  pindx->fltnum = gfltnum;
  pindx->hipnt = NULL;
  pindx->hfpnt = NULL;
  pindx->intpnt = NULL;
  pindx->fltpnt = NULL;

  fwrite (pindx,sizeof(struct gaindx),1,mntape);

     hinum[0]=1;
     hinum[1]=1;
     hinum[2]=IGribRec;
     hinum[3]=255;
     fwrite(&hinum[0],sizeof(int),1,mntape);
     fwrite(&hinum[1],sizeof(int),1,mntape);
     fwrite(&hinum[2],sizeof(int),1,mntape);
     fwrite(&hinum[3],sizeof(int),1,mntape);

     for (i=0; i<IGribRec; i++)
     {
        fwrite(&intnum[i][0],sizeof(int),1,mntape);
        fwrite(&intnum[i][1],sizeof(int),1,mntape);
        fwrite(&intnum[i][2],sizeof(int),1,mntape);
     }

     for (i=0; i<IGribRec; i++)
     {
        ADD=0.0;
        MULT=1.0;
        if (CodesTape)
           get_codes(pcode[i], C_code, C_CODE); 
 
        fltnum[i][2] = fltnum[i][2] + ADD;
        fltnum[i][0] = fltnum[i][0] / MULT;
        
        fwrite(&fltnum[i][0],sizeof(float),1,mntape);
        fwrite(&fltnum[i][1],sizeof(float),1,mntape);
        fwrite(&fltnum[i][2],sizeof(float),1,mntape);
     }
}

int grbcntl(void)
{
  int cnt,rc,i,flg,iarg;
  char rec[512], *ch;
  int len, skip;
  int inum;
  int skipcode0 = 0, skiplevel0 = 0;
  int stopprint = 0;
  int skiprecord;
  int leveltype;
  int ic;

  skip    = 0;
  scanflg = 1;
  /* Get file size */

  fseek (gfile,0L,2);
  flen = ftell(gfile);

  /* Set up to skip appropriate amount.  */

  if (skip > -1)
    fpos = skip;
  else {
    fseek (gfile,0,0);
    rc = fread (rec,1,100,gfile);
    if (rc < 100) {
      printf ("I/O Error reading header\n");
      return (1);
    }
    len  = gagby(rec,88,4);
    fpos = len*2 + 100;
  }

  /* We are positioned.  Go read the first GRIB header */

  IGribRec = 0;
  while (1) {
    rc = gribhdr(&ghdr);
    if (rc) break;
    GribRec++;
    if (Debug) {
      printf("%d %d %d %d %d %d %d %d %d\n",  ghdr.parm, ghdr.level, 
                       ghdr.dtim.yr, ghdr.dtim.mo,
                       ghdr.dtim.dy,ghdr.dtim.hr,ghdr.dtim.mn,
                       ghdr.gicnt,ghdr.gjcnt);
      printf("code:%d  ltyp:%d \n",ghdr.parm,ghdr.ltyp);
      printf("intnum: %d %d %d\n",ghdr.dpos,ghdr.bpos,ghdr.bnum);
      printf("fltnum: %f %f %f\n",ghdr.dsf,ghdr.bsf,ghdr.ref);
    }

    skiprecord = 0;
    leveltype = ghdr.ltyp;
    if ( ghdr.ltyp == 99 ) leveltype = 100;
    /*
    ic = 0;
    for (ic = 0; ic < MAX_CODES; ic++)
      {
	if ( S_Codes[ic].code == 0 ) break;
	printf(">> %d ", S_Codes[ic].code);
      }
    printf("\n");
    */
    ic = 0;
    for (ic = 0; ic < MAX_CODES; ic++)
      {
	if ( S_Codes[ic].code == 0 ) break;
	if ( S_Codes[ic].code == ghdr.parm )
	  {
	    if ( S_Codes[ic].leveltype != leveltype )
	      {
		skiprecord = 1;
		/*
		printf("%d %d %d\n", S_Codes[ic].code,  S_Codes[ic].leveltype, ghdr.parm);
		*/
		break;
	      }
	  }
      }


    if (ghdr.iflg || ReducedGrid || skiprecord)
      {
	if (O_ctime <= 1)
	  {
	    int code, level;

	    code = ghdr.parm;
	    if (ghdr.ltyp == 100)
	      level = 100*ghdr.level;
	    else
	      level = ghdr.level;

	    if (skipcode0 == code && skiplevel0 == level) stopprint = 1;

            if (!stopprint)
	      {
		fprintf (stderr, "  %4d %9d", code, level);
		if (ReducedGrid)
		  fprintf (stderr, "    reduced grid skiped !!!\n");
		else if (ghdr.iflg)
		  fprintf (stderr, "    spectral data skiped !!!\n");
		else
		  fprintf (stderr, "    already defined skiped !!!\n");
	      }

	    if (!skipcode0)
	      {
		skipcode0 = code;
		skiplevel0 = level;
	      }
	  }
      }
    else
      {
	intnum[IGribRec][0] = ghdr.dpos;
	intnum[IGribRec][1] = -999;
	if (ghdr.bmsflg) intnum[IGribRec][1] = ghdr.bpos;
	intnum[IGribRec][2] = ghdr.bnum;
	fltnum[IGribRec][0] = ghdr.dsf;
	fltnum[IGribRec][1] = ghdr.bsf;
	fltnum[IGribRec][2] = ghdr.ref;
	pcode[IGribRec] = ghdr.parm;
	IGribRec++;

	gribctl();
      }

  }
  if (rc==99) {
     exit(1);
  }
  if (rc==50) {
    printf ("I/O Error reading GRIB file\n");
    printf ("Possible cause: premature EOF\n");
    return(1);
  }
  if (rc>1) {
    printf ("GRIB file format error \n");
    return (1);
  }

/*  printf ("Reached EOF\n"); */
/*  fclose (gfile); */
}
/* Read a GRIB header, and get needed info out of it.  */

int gribhdr (struct grhdr *ghdr) {
#define MAXREC 512*100
  char rec[MAXREC];
  int i,isc,len,rc,sign,mant;
  int cpos;
  int la1,la2,lo1,lo2,ladi,lodi;
  int imax;
  int gribrecfound;

  if (fpos+10 >= flen) return(1);

  /* Position at start of next record */

  rc = fseek(gfile,fpos,0);
  if (rc) return(50);

  /* Read start of record -- length, version number */

  imax = MAXREC - 12;
  if (fpos+imax>flen) imax = flen - fpos;

  if (verboseflg)
    printf ("recmax:%d\n", imax);

  rc = fread (rec, 1, imax, gfile);
  if (rc < imax) return(50);

  gribrecfound = 0;

  if (*rec!='G' || *(rec+1)!='R' || *(rec+2)!='I' || *(rec+3)!='B') {
    i = 0;
    if (scanflg) {
      i = 1;
      while ((i < MAXREC-17) && !gribrecfound) {
        isc = i;

	if (verboseflg)
	  printf ("%5d ->%c%c%c%c%<-\n", i, *(rec+i),*(rec+i+1),*(rec+i+2),*(rec+i+3));

        if (*(rec+i)=='G' && *(rec+i+1)=='R' && *(rec+i+2)=='I' && *(rec+i+3)=='B')
	{
          fpos = fpos + i;
          rc = fseek(gfile,fpos,0);
          if (rc) return (50);
          rc = fread(rec,1,8,gfile);
	  /*          if (rc < 8) return (50);*/
          if (rc < 8) return (-1);
          gribrecfound = 1;
        } else i++;
      }
    }
    if (verboseflg)
      printf ("scanflg:%d %d %d\n", isc, fpos, flen);

    if ((flen-fpos) < 200) return (-1);

    if (!gribrecfound) {
      printf ("GRIB header not found.  File apparantly not GRIB data\n");
      printf ("->%c%c%c%c%c%c<-\n",*rec,*(rec+1),*(rec+2),*(rec+3),*(rec+4),*(rec+5));
      return (99);
    }
  }

  cpos = fpos;
  ghdr->vers = gagby(rec,7,1);
  if (ghdr->vers > 1) {
    printf ("File is not GRIB version 0 or 1, 0 or 1 is required. \n");
    printf (" Version number is %i\n",ghdr->vers);
    return (99);
  }
  if (ghdr->vers == 0) {
    cpos += 4;
    rc = fseek(gfile,cpos,0);
    if (rc) return(50);
  } else {
    ghdr->len = gagby(rec,4,3);
    cpos = cpos + 8;
    rc = fseek(gfile,cpos,0);
    if (rc) return(50);
  }

  /* Get PDS length, read rest of PDS */

  rc = fread(rec,1,3,gfile);
  if (rc<3) return(50);
  len = gagby(rec,0,3);
  ghdr->pdslen = len;
  cpos = cpos + len;
  rc = fread(rec+3,1,len-3,gfile);
  if (rc<len-3) return(50);

  /* Get info from PDS */

  ghdr->id      = gagby(rec,6,1);
  ghdr->gdsflg  = gagbb(rec+7,0,1);
  ghdr->bmsflg  = gagbb(rec+7,1,1);
  ghdr->parm    = gagby(rec,8,1);
  ghdr->ltyp    = gagby(rec,9,1);
  ghdr->level   = gagby(rec,10,2);
  ghdr->l1      = gagby(rec,10,1);
  ghdr->l2      = gagby(rec,11,1);
  ghdr->dtim.yr = gagby(rec,12,1);
  ghdr->dtim.mo = gagby(rec,13,1);
  ghdr->dtim.dy = gagby(rec,14,1);
  ghdr->dtim.hr = gagby(rec,15,1);
  ghdr->dtim.mn = gagby(rec,16,1);
  ghdr->ftu     = gagby(rec,17,1);
  ghdr->tri     = gagby(rec,20,1);

  if (ghdr->tri == 10) {
    ghdr->p1    = gagby(rec,18,2);
    ghdr->p2    = 0;
  } else {
    ghdr->p1    = gagby(rec,18,1);
    ghdr->p2    = gagby(rec,19,1);
  }

  if (len > 24) {
    ghdr->cent = gagby(rec,24,1);
    ghdr->dtim.yr = ghdr->dtim.yr + (ghdr->cent-1)*100;
  } else {
    ghdr->cent = -999;
    if (ghdr->dtim.yr > 49) ghdr->dtim.yr += 1900;
    else ghdr->dtim.yr += 2000;
  }
  if (len > 25) {
    ghdr->dsf = (float)gagbb(rec+26,1,15);
    i = gagbb(rec+26,0,1);
    if (i) ghdr->dsf = -1.0*ghdr->dsf;
    ghdr->dsf = pow(10.0,ghdr->dsf);
  } else ghdr->dsf = 1.0;

  if (verboseflg) {
    printf ("GRIB %i record at %i length %i\n",ghdr->vers,fpos,ghdr->len);
    printf (" id,gdsflg,bmsflg,parm,ltyp,level = %i %i %i %i %i %i %i %i\n",
	    ghdr->id,ghdr->gdsflg,ghdr->bmsflg,ghdr->parm,ghdr->ltyp,
	    ghdr->level,ghdr->l1,ghdr->l2);
    printf (" Cent = %i\n",ghdr->cent);
    printf (" date = %i %i %i %i %i %i\n",ghdr->cent,ghdr->dtim.yr,
	    ghdr->dtim.mo,ghdr->dtim.dy,ghdr->dtim.hr,ghdr->dtim.mn);
    printf (" ftunit,trind,p1,p2 = %i %i %i %i\n",ghdr->ftu,
	    ghdr->tri,ghdr->p1,ghdr->p2);
    printf (" dsf = %g\n",ghdr->dsf);
  }

  /* If it is there, get info from GDS */

  if (ghdr->gdsflg) {
    rc = fread(rec,1,3,gfile);
    if (rc<3) return(50);
    len = gagby(rec,0,3);
    cpos = cpos + len;
    /*    if (len>256) {*/
    if (len>MAXREC-3) {
      printf ("GDS error: bad length %i\n",len);
      return (50);
    }
    ghdr->gdslen = len;
    rc = fread(rec+3,1,len-3,gfile);
    if (rc<len-3) return(50);
    ghdr->gtyp = gagby(rec,4,1);
    ReducedGrid = 0;
    if (ghdr->gtyp != 255 && ghdr->gtyp != 0)
      {
        ReducedGrid = (ghdr->gdslen - 32 - 4 * gagby(rec,3,1));
        if (gagby(rec,6,2) != 65535) ReducedGrid = 0;
      }
    RepGrib = gagby(rec,5,1);
    if (RepGrib == REP_REGULAR || RepGrib == REP_GAUSS || RepGrib == REP_ROTREG) {
      if (ReducedGrid) {
	ghdr->gjcnt = gagby(rec,8,2);
	ghdr->gicnt = 2 * ghdr->gjcnt;
      } else {
	ghdr->gicnt = gagby(rec,6,2);
	ghdr->gjcnt = gagby(rec,8,2);
      }
    }
    else {
      if (RepGrib != REP_SPECTRAL) {
	fprintf (stderr, "warning: unsupported GRIB representation %d\n", RepGrib);
      }
    }

    la1 = gagbb(rec+10,1,23);
    i = gagbb(rec+10,0,1);
    if (i) la1 = -1*la1;
    lo1 = gagbb(rec+13,1,23);
    i = gagbb(rec+13,0,1);
    if (i) lo1 = -1*lo1;
    la2 = gagbb(rec+17,1,23);
    i = gagbb(rec+17,0,1);
    if (i) la2 = -1*la2;
    lo2 = gagbb(rec+20,1,23);
    i = gagbb(rec+20,0,1);
    if (i) lo2 = -1*lo2;
    /*
    lodi = gagby(rec,23,2);
    ladi = gagby(rec,25,2);
    */
    lodi = gagbb(rec+23,1,15);
    i = gagbb(rec+23,0,1);
    if (i) lodi = -1*lodi;
    ladi = gagbb(rec+25,1,15);
    i = gagbb(rec+25,0,1);
    if (i) ladi = -1*ladi;

    ghdr->gsf1 = gagbb(rec+27,0,1);
    ghdr->gsf2 = gagbb(rec+27,1,1);
    ghdr->gsf3 = gagbb(rec+27,2,1);
    lat1=la1;lat2=la2;lon1=lo1;lon2=lo2;latdi=ladi;londi=lodi;
    if (verboseflg) {
      printf (" gds len,typ,i,j,flags = %i %i %i %i %i %i %i\n",
	      len,ghdr->gtyp,ghdr->gicnt,ghdr->gjcnt,ghdr->gsf1,ghdr->gsf2,ghdr->gsf3);
      printf ("  gds la1,lo1,la2,lo2,ladi,lodi = %i %i %i %i %i %i\n",
	      la1,lo1,la2,lo2,ladi,lodi);
    }
  } else ghdr->gdslen = 0;

  /* Get necessary info about BMS if it is there */

  if (ghdr->bmsflg) {
    rc = fread(rec,1,6,gfile);
    if (rc<6) return(50);
    len = gagby(rec,0,3);
    ghdr->bmslen = len;
    ghdr->bnumr = gagby(rec,4,2);
    ghdr->bpos = cpos+6;
    cpos = cpos + len;
    rc = fseek(gfile,cpos,0);
if (verboseflg) 
    printf (" Bit Map Section: pos = %i\n",ghdr->bpos);
  } else ghdr->bmslen = 0;

  /* Get necessary info from data header */

  rc = fread(rec,1,11,gfile);
  if (rc<11) return(50);
  len = gagby(rec,0,3);
  ghdr->bdslen = len;
  ghdr->iflg = gagbb(rec+3,0,2);
  /*
  if (ghdr->iflg) {
    printf ("GRIB data bds flag is set!\n");
  }
  */
  i = gagby(rec,4,2);
  if (i>32767) i = 32768-i;
  ghdr->bsf = pow(2.0,(float)i);

  i = gagby(rec,6,1);
  sign = 0;
  if (i>127) {
    sign = 1;
    i = i - 128;
  }
  mant = gagby(rec,7,3);
  if (sign) mant = -mant;
  ghdr->ref = (float)mant * pow(16.0,(float)(i-70));

  ghdr->bnum = gagby(rec,10,1);
  ghdr->dpos = cpos+11;
if (verboseflg) {
  printf (" data len,flg,bsf,ref,bnum,dpos %i %i %g %g %i %i\n",len,ghdr->iflg,
    ghdr->bsf,ghdr->ref,ghdr->bnum,ghdr->dpos);
}

  if (ghdr->vers==0) {
    fpos = fpos + 8 + ghdr->pdslen + ghdr->gdslen +
                      ghdr->bmslen + ghdr->bdslen;
if (verboseflg) {
    printf ("Lengths: pds,gds,bms,bds = %i %i %i %i \n",
      ghdr->pdslen,ghdr->gdslen,ghdr->bmslen,ghdr->bdslen);
}
  } else fpos = fpos + ghdr->len;

/*  free (ghdr); */
  return(0);

}

int main(argc,argv)
int       argc;
char          *argv[];
{
   INT code,i, ic;
   unsigned char *cp;
   INT c, len,ifilelen;
   char  *Progname,*cstring;
   char *gauss_arg,*yrev_arg;

/*     SIOUXSettings.asktosaveonclose = false; */
/*     argc = ccommand(&argv); */
    Progname = argv[0];

    while (--argc && (*++argv)[0] == '-')  {
        c = *++argv[0];
        cstring = argv[0];
        len = strlen(cstring);
        switch (c)  {
        case 'c':
            if (!strncmp(cstring, "controlfile", len))  {
                cfile=(++argv)[0];
                --argc;
                break;
            }
            if (!strncmp(cstring, "codesfile", 5))  {
                codesfile=(++argv)[0];
                --argc;
                break;
            }
        case 'C':
            if (!strncmp(cstring, "Codesfile", len))  {
                codesfile=(++argv)[0];
                --argc;
                break;
            }
        case 'g':
            if (!strncmp(cstring, "gauss", len))  { 
                gauss_arg=(++argv)[0];
                if (!strncmp(gauss_arg,"on", 2)) GAUSS=1;
                if (!strncmp(gauss_arg,"off", 3)) GAUSS=0;
                --argc;
                break;
            }
        case 'h':
            if (!strncmp(cstring, "help", len))  { 
                Help();
                return(0);
                break;
            }
        case 'i':
            if (!strncmp(cstring, "inputfile", len))  { 
                ifile=(++argv)[0];
                --argc;
                break;
            }
        case 'I':
            if (!strncmp(cstring, "Increment", len))  { 
                pIncr=(++argv)[0];
                --argc;
                break;
            }
        case 'm':
            if (!strncmp(cstring, "mapfile", len))  {
                mfile=(++argv)[0];
                --argc;
                break;
            }
        case 'o':
            if (!strncmp(cstring, "outputfile", len))  {
                ofile=(++argv)[0];
                --argc;
                break;
            }
        case 'q':
            if (!strncmp(cstring, "quiet", len))  {
                QUIET=1;
                break;
            }
        case 'd':
            if (!strncmp(cstring, "debug", len))  {
                Debug=1;
                break;
            }
        case 'v':
            if (!strncmp(cstring, "verbose", len))  {
                verboseflg=1;
                break;
            }
        case 'V':
            if (!strncmp(cstring, "Version", len))  {
                Version();
                return(0);
                break;
            }
        case 'T':
            if (!strncmp(cstring, "Time", len))  { 
                pTime=(++argv)[0];
                --argc;
                break;
            }
        case 'x':
            if (!strncmp(cstring, "xydef", len))  {
                xyfile=(++argv)[0];
                --argc;
                break;
            }
        case 'y':
            if (!strncmp(cstring, "yrev", len))  { 
                yrev_arg=(++argv)[0];
                if (!strncmp(yrev_arg,"on", 2)) LYREV=1;
                if (!strncmp(yrev_arg,"off", 3)) LYREV=0;
                --argc;
                break;
            }
        case 'z':
            if (!strncmp(cstring, "zdef", len))  {
                zfile=(++argv)[0];
                --argc;
                break;
            }
        default:
/*            usage(Progname); */
/*            error("illegal option %s\n", cstring);*/
            Help();
            printf("illegal option %s\n", cstring);
            exit(1);
            break;
        }
    }
    if (argc && !ifile)  {
        ifile = argv[0];
        argc--;
    }
    if (argc && !ofile) {
        ofile = (++argv)[0];
        argc--;
    }
    if (argc)  {
/*        usage(Progname);*/
/*        error("too many arguments");*/
    }

    /*******************/
    /* open input file */
    /*******************/

    if (!ifile) {
      Help();
      return(0);
    }
    else
      {
	if (!QUIET)
	  fprintf (stdout, "  Input File   : %s\n", ifile);
      }

    intape = fopen(ifile,"rb");
    if (intape == 0)
      {
	fprintf (stderr, "could not open input file %s\n", ifile);
	exit(1);
      }

    filename = strrchr(ifile,'/');
    if (filename == NULL) filename = ifile;
    else               filename++ ;

    fext = strrchr(ifile,'.');
    ifilelen=strlen(ifile);
    if (ifilelen >=4)  {
      if (strcmp(&ifile[ifilelen-3],"grb"))
	fext=0;
    }
    else {
      fext=0;
    }

    if (fext == NULL) strcpy (iname,ifile);
    else strncpy(iname,ifile,fext-ifile);

    /*********************/
    /* open control file */
    /*********************/

    if (!cfile) {
      cfile=cname;
      strcpy(cname,iname);
      strcat(cname,".ctl");
    }
    if (!QUIET)
      fprintf (stdout, "  Control File : %s\n", cfile);

    cntape = fopen(cfile,"w");
    if (cntape == NULL)
      {
	fprintf (stderr, "could not open control file %s\n", cfile);
	exit(1);
      }
    /* printf("base : %d\n",&cntape->_base);*/

    /*********************/
    /* open map file     */
    /*********************/

    if (!mfile) {
      mfile=mname;
      strcpy(mname,iname);
      strcat(mname,".gmp");
    }
    if (!QUIET)
      fprintf (stdout, "  Map File     : %s\n", mfile);

    mntape = fopen(mfile,"wb");
    if (mntape == NULL)
      {
	fprintf (stderr, "could not open map file %s\n", mfile);
	exit(1);
      }

    /*********************/
    /* open codes file   */
    /*********************/

    if (codesfile)  {
      if (!QUIET)
	fprintf (stdout, "  Code File : %s\n", codesfile);
      CodesTape = fopen(codesfile,"r");
      if (CodesTape == NULL)
	{
	  fprintf (stderr, "could not open codes file %s\n", codesfile);
	  exit(1);
	}
    }

    /*********************/
    /* open xydef file   */
    /*********************/

    if (xyfile)  {
      if (!QUIET)
	fprintf (stdout, "  xydef File : %s\n", xyfile);
      XYTape = fopen(xyfile,"r");
      if (XYTape == NULL)
	{
	  fprintf (stderr, "could not open xydef file %s\n", xyfile);
	  exit(1);
	}
    }

    /*********************/
    /* open zdef file   */
    /*********************/

    if (zfile)  {
      if (!QUIET)
	fprintf (stdout, "  zdef File : %s\n", zfile);
      ZTape = fopen(zfile,"r");
      if (ZTape == NULL)
	{
	  fprintf (stderr, "could not open zdef file %s\n", zfile);
	  exit(1);
	}
    }

    for ( ic = 0; ic < MAX_CODES; ic++ )
      S_Codes[ic].code = 0;

    grbcntl();
    if (intape) fclose(intape);
    if (O_CCode) endctl();
    if (cntape) fclose(cntape);

    if (!QUIET)
      fprintf (stdout, "\n  Processed %d of %d records\n", IGribRec, GribRec);

    return(0);
}
