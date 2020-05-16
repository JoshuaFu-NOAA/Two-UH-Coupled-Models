#define REAL double
#define INT  long
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

#define SEEK_CUR 1

#define MAX_TIME     19999
#define MAX_LEVELS      99
#define MAX_CODES      199
/*#define MAX_GREC     MAX_CODES*MAX_LEVELS*/
/*#define MAX_GREC     299999*/
#define MAX_GREC     499999
#define MAX_Latitudes  2048


char    *filename,*fext;
char    *ifile;
char    *ofile;
char    *cfile;
char    *mfile;
char    *codesfile;
char    *xyfile;
char    *zfile;
char    iname[256];
char    oname[256];
char    cname[256];
char    mname[256];


REAL ADD,MULT;
INT NTIME;

FILE  *intape;
FILE  *cntape;
FILE  *ontape;
FILE  *mntape;
FILE  *CodesTape;
FILE  *XYTape;
FILE  *ZTape;

INT I_code, I_level, I_clev, I_time, I_ires, I_ilev;
INT O_Date, O_Time, O_Lat, O_Lon;
INT O_CCode, O_ccode, O_ctime, O_CLevel, O_clevel;
INT O_Year, O_Month, O_Day, O_Hour, O_Minute, O_Dt; 
INT O_Code, O_Level, MaxLevel;
REAL I_dx, I_dy, I_north;

struct CODES {
  INT code;
  REAL add;
  REAL mult;
  INT level;
  INT lev0;
  int leveltype;
};

char C_code[20], C_CODE[60];
INT LEVELS[MAX_LEVELS];
INT DATE[MAX_TIME];
INT TIME[MAX_TIME];

struct CODES S_Codes[MAX_CODES];

void Abort(char *errtext);

INT get_codes(INT code, char C_code[], char C_CODE[]);

INT AddLevel(INT ilev, INT nlev);

INT SetLevel(INT icode, INT nc, INT ilev, INT nlev);

INT AddCode(INT icode, INT nc, int leveltype);

void gribctl(void);
int grbcntl(void);

void endctl(void);
int gagby (char *, int, int);
int gagbb (char *, int, int);
extern int gaussg_(int, double *, double *);
