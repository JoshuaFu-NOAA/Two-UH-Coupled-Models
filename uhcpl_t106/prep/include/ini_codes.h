#ifndef INI_CODES_H
#define INI_CODES_H

#define MAXCODE_PTsurfini 28

typedef struct
{
  int code;
  char *name, *longname, *unit;
}
Table;

typedef struct
{
  int ncodes;
  int nlon;
  int nlat;
  int varid[MAXCODE_PTsurfini];
  int code[MAXCODE_PTsurfini];
  double *var[MAXCODE_PTsurfini];
}
SINITIAL;
/*
int PT_getidx_code(Table *ptable, int code);
*/
#endif
