#ifndef PP_H
#define PP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "config.h"
#include "netcdf.h"

typedef struct
{
  int nc_file_id;		/* netCDF file id */
  int nc_access_mode;		/* access mode for that file */
  int nc_dims[NC_MAX_VAR_DIMS];	/* dimension vector */
  int climate;			/* climate data =1 or (monthly means =0) */
  char nc_creation_program[NC_MAX_NAME];	/* name of this program */
  char nc_creation_user[NC_MAX_NAME];	/* who has run this program */
  char nc_creation_date[NC_MAX_NAME];	/* date of creation of netCDF file */
  char nc_binary_source[NC_MAX_NAME];	/* original binary data type (CRAY/IEEE) */
  char nc_file_type[NC_MAX_NAME];	/* initital or restart file ... */
  char nc_conventions[NC_MAX_NAME];	/* convention */
  char nc_model[NC_MAX_NAME];	/* model */
  char nc_institution[NC_MAX_NAME];	/* institution */
  char nc_data_source[NC_MAX_NAME];	/* source of data  set by -a  */
}
netCDF_file;

typedef struct
{
  int nm, nn, nk;		/* spherical truncation */
  int nmp1, nnp1, nkp1;		/* spherical truncation + 1 */
  int nsp;			/* number of wavenumbers per latitude and level */
  int ngl, nlon, nlev;		/* dimensions in gridspace */
  int nhgl;			/* ngl/2 */
  int nlevp1;			/* levels + 1 */
  int nlp2, n2;			/* nlon + 2 */
  int nhtrac;			/* number of tracer in history file - gl buffer */
  int ntrac;			/* number of tracer */
  int nstep;			/* Saved timestep */
  int nrd;			/* number of switches ... */
  double dt;			/* timestep */
  double *switches;		/* pointer to switch vector */
}
DIMENSIONS;

typedef struct
{
  char label[8][81];
}
LABELS;

typedef struct
{
  int fdate;			/* time of the initial data */
  int ftime;
  int vdate;			/* verifying time of the initial data */
  int vtime;
}
TIMES;

typedef struct
{
  char def[5];
  UINT8 ddr1_init;
  UINT8 *ddr1, *ddr2, *ddr3, *ddr4, *ddr5;
  UINT8 nrecl4;
  size_t bddr1, bddr2, bddr3, bddr4, bddr5;
}
HEADER;

typedef struct
{
  size_t nstep;			/* number of time step */
  size_t nlat, nlon;		/* number of Gaussian latitudes and longitudes */
  size_t nlev;
  double *step;
  double *lev;
  double *lat, *lon;		/* Gaussian latitudes and longitudes */
  double *vct;			/* vertical transformation coefficients (all in one vector) */
  double *vct_a;		/* vertical transformation coefficients set A */
  double *vct_b;		/* vertical transformation coefficients set B */
  int nvclev;			/* number of elements per vector (for vct_a and vct_b) */

}
GRID;

typedef struct
{
  INT8 icode, ilevel, idate, itime, nlon, nlat, idisp1, idisp2;
}
SERVICE;

typedef struct
{
  double *geosp, *ts, *ws, *wl, *sn, *slm, *az0, *alb, *varp, *varor, *forest,
    *vgrat, *vlt, *wsmx, *fao;
}
INITIAL;

typedef struct
{
  int icode;
  char *name;
  char *longname;
  char *unit;
  double *wert;
  int da;
  int varid;
}
VAR;


/* end of header.h */
void set_Progname (char *progarg);
void set_Commandline (int argc, char *argv[]);

void gaussian_latitudes (double *lat, int nlat);
void tgrid_lat (double *lat, int nlat);
void tgrid_lon (double *lon, int nlon);
int get_pplat (int ilat, int nlat);
int get_offset (int it, int iz, int iy, int ix, int nt, int nz, int ny,
		int nx);

void cray2ieee (UINT8 * crayf, UINT8 * ieeef, int *nf);
void util_gwunpk (double *a, double *b, double *c, double *d, double *p,
		  int *klen);

int read_fblk (FILE * fp, char *fortran_block);
int read_i8array_fblk (FILE * fp, INT8 * i8array, size_t dim);
int read_r8array_fblk (FILE * fp, REAL8 * r8array, size_t dim);
int read_i8array_bin (FILE * fp, INT8 * i8array, size_t dim);
int read_r8array_bin (FILE * fp, REAL8 * r8array, size_t dim);

int write_fblk (FILE * fp, char *fortran_block);
int write_i8array_fblk (FILE * fp, INT8 * i8array, size_t dim);
int write_r8array_fblk (FILE * fp, REAL8 * r8array, size_t dim);
int write_i8array_bin (FILE * fp, INT8 * i8array, size_t dim);
int write_r8array_bin (FILE * fp, REAL8 * r8array, size_t dim);

void read_srv8_fblk (FILE * ostream, double *data, GRID * grid,
		     SERVICE * serv);
void write_srv8_fblk (FILE * ostream, double *data, GRID * grid,
		      SERVICE * serv);

double min_darray (double *parray, size_t dim);
double max_darray (double *parray, size_t dim);
double mean_darray (double *parray, size_t dim);

#endif
