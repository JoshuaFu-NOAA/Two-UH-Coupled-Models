#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "pp.h"

int Debug = 0;
int Verbose = 0;
char *Progname = NULL;
char *Commandline = NULL;


int
get_pplat (int ilat, int nlat)
{
  return (((ilat % 2) * (nlat - (ilat - 1) / 2 - 1) +
	   ((ilat + 1) % 2) * ilat / 2));
}

int
get_offset (int it, int iz, int iy, int ix, int nt, int nz, int ny, int nx)
{
  return (((it * nz + iz) * ny + iy) * nx + ix);
}



/*
 *  Convert Floating point, Cray to IEEE 64-bit
 *
 *
 *      input : cf      Cray Floating point numbers
 *              nf      Number of elements in cf
 *      output: ieeef   IEEE Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */

void
cray2ieee (UINT8 * crayf, UINT8 * ieeef, int *nf)
{

  const UINT8 cray_expbias = 0x4000;	/* Cray exponent bias */
  const UINT8 ieee_expbias = 0x3ff;	/* IEEE exponent bias */

  const UINT8 mask1 = 0x8000000000000000;
  const UINT8 mask2 = 0x7fff;
  const UINT8 mask3 = 0xfffffffffffff;
  const UINT8 mask4 = 0x7ff;

  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < *nf; i++)
    {
      *(ieeef + i) = (*(crayf + i) & mask1)
	|
	((((*(crayf + i) >> 48) & mask2) - cray_expbias + ieee_expbias - 1) <<
	 52) | ((*(crayf + i) << 5) & mask3);
    }
  /*
   * handle 0.0, underflow and overflow :
   *
   * if the exponent field goes negative then the Cray number was
   * either 0.0 or too small to represent in IEEE, in either case
   * set the IEEE result to all 0's which is the IEEE representation for 0.0.
   *
   * if the exponent field is too large for the IEEE field, set
   * the result to the indefinite value.
   */

  /* Cray's internal conversion seems to skip this underflow part and instead
     of indef sets 0 .... */

  for (i = 0; i < *nf; i++)
    {
      if ((((*(crayf + i) >> 48) & mask2) - cray_expbias + ieee_expbias - 1) <
	  0)
	{
	  *(ieeef + i) = 0;
	}
      else
	if ((((*(crayf + i) >> 48) & mask2) - cray_expbias + ieee_expbias - 1)
	    > mask4)
	{
	  *(ieeef + i) = 0;
	}
    }

  return;
}

/*
 *  Convert Floating point, IEEE 64-bit, to Cray Floating point
 *
 *	input : ieeef	IEEE Floating point numbers (double precision)
 *		nf	Number of elements in ieeef
 *	output: crayf	Cray Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa	unused
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */
void
ieee2cray (UINT8 *ieeef, UINT8 *crayf, int *nf)
{
  const UINT8 cray_expbias = 0x4000;      /* Cray exponent bias */
  const UINT8 ieee_expbias = 0x3ff;       /* IEEE exponent bias */
  const UINT8 implied = 0x10000000000000; /* implied bit in IEEE mantissa */

  const UINT8 mask1 = 0x8000000000000000;
  const UINT8 mask3 = 0xfffffffffffff;
  const UINT8 mask4 = 0x7ff;

  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < *nf; i++) {
    if (*(ieeef+i) == 0) {
      *(crayf+i) = 0;
    } else {
      *(crayf+i) = (*(ieeef+i) & mask1)
	 |
 	  ((((*(ieeef+i) >> 52) & mask4)-ieee_expbias+cray_expbias+1) << 48)
 	 |
 	  (((*(ieeef+i) & mask3)+implied) >> 5);
    }
  }

  return;
}

/*  this is the new corrected gaussian_latitudes function */

void
gaussian_latitudes (double *pa, int nlat)
{
  /*
   * Compute Gaussian latitudes.  On return pa contains the
   * sine of the latitudes starting closest to the north pole and going
   * toward the south
   *
   */

  double pmod;			/* last modification */
  const int itemax = 20;

  int isym, iter, ins2, jn, j;
  double za, zw, zan;
  double z, zk, zkm1, zkm2, zx, zxn, zldn, zmod;

  /*
   * Perform the Newton loop
   * Find 0 of Legendre polynomial with Newton loop
   */

  ins2 = nlat / 2 + nlat % 2;

  for (j = 0; j < ins2; j++)
    {
      z = (double) (4 * (j + 1) - 1) * M_PI / (double) (4 * nlat + 2);
      pa[j] = cos (z + 1.0 / (tan (z) * (double) (8 * nlat * nlat)));
    }

  for (j = 0; j < ins2; j++)
    {

      za = pa[j];

      iter = 0;
      do
	{
	  iter++;
	  zk = 0.0;

	  /* Newton iteration step */

	  zkm2 = 1.0;
	  zkm1 = za;
	  zx = za;
	  for (jn = 2; jn <= nlat; jn++)
	    {
	      zk =
		((double) (2 * jn - 1) * zx * zkm1 -
		 (double) (jn - 1) * zkm2) / (double) (jn);
	      zkm2 = zkm1;
	      zkm1 = zk;
	    }
	  zkm1 = zkm2;
	  zldn = ((double) (nlat) * (zkm1 - zx * zk)) / (1. - zx * zx);
	  zmod = -zk / zldn;
	  zxn = zx + zmod;
	  zan = zxn;

	  /* computes weight */

	  zkm2 = 1.0;
	  zkm1 = zxn;
	  zx = zxn;
	  for (jn = 2; jn <= nlat; jn++)
	    {
	      zk =
		((double) (2 * jn - 1) * zx * zkm1 -
		 (double) (jn - 1) * zkm2) / (double) (jn);
	      zkm2 = zkm1;
	      zkm1 = zk;
	    }
	  zkm1 = zkm2;
	  zw = (1.0 - zx * zx) / ((double) (nlat * nlat) * zkm1 * zkm1);
	  za = zan;
	}
      while (iter <= itemax && fabs (zmod) >= DBL_EPSILON);

      pmod = zmod;
      pa[j] = zan;
      /*
         pw[j] = 2.0*zw;
       */
    }

  for (j = 0; j < nlat / 2; j++)
    {
      isym = nlat - (j + 1);
      pa[isym] = -pa[j];
      /*
         pw[isym] = pw[j];
       */
    }

  return;
}

void
tgrid_lat (double *lat, int nlat)
{
  int i;

  gaussian_latitudes (lat, nlat);

  for (i = 0; i < nlat; i++)
    {
      lat[i] = asin (lat[i]) / M_PI * 180.0;
    }
}

void
tgrid_lon (double *lon, int nlon)
{
  int i;
  double lon0, dlon;

  lon0 = 0.0;
  dlon = 360.0 / nlon;
  for (i = 0; i < nlon; i++)
    {
      lon[i] = lon0 + i * dlon;
    }
}


/*
 * util_gwunpk - unpacks one input array into four output arrays
 * each packed field is held in 16 bits
 * it is shifted by 7 bits and floated
 *
 * p is the packed input array of length klen words
 * a, b, c, d are gwunpked output arrays of length klen words
 */

void
util_gwunpk (double *a, double *b, double *c, double *d, double *p, int *klen)
{
  UINT8 ia, ib, ic, id, ip;
  UINT8 imask;

  const UINT8 mask16 = 0xffff000000000000;

  const int itrunc = 7;

  int jl;

  imask = mask16 >> (48 - itrunc);

  for (jl = 0; jl < *klen; jl++)
    {

      memcpy (&ip, p + jl, sizeof (UINT8));

      ia = ip >> (48 - itrunc) & imask;
      ib = ip >> (32 - itrunc) & imask;
      ic = ip >> (16 - itrunc) & imask;
      id = ip << itrunc & imask;

      *(a + jl) = (double) ia;
      *(b + jl) = (double) ib;
      *(c + jl) = (double) ic;
      *(d + jl) = (double) id;
    }

  return;
}

void
check_blocklen (int len1, int len2)
{
  if (len1 != len2)
    {
      fprintf (stderr, "error: unexpected block length\n");
      fprintf (stderr, "       expected block length: %d\n", len1);
      fprintf (stderr, "           read block length: %d\n", len2);
    }
}

int
read_i8array_fblk (FILE * fp, INT8 * i8array, size_t dim)
{
  char block[4];
  int blocklen;

  blocklen = read_fblk (fp, block);
  if (!feof (fp))
    check_blocklen (dim * sizeof (INT8), blocklen);
  else
    return (-1);

  read_i8array_bin (fp, i8array, dim);

  blocklen = read_fblk (fp, block);
  if (!feof (fp))
    check_blocklen (dim * sizeof (INT8), blocklen);
  else
    return (-1);

  return (0);
}

int
read_i8array_bin (FILE * fp, INT8 * i8array, size_t dim)
{
  fread (i8array, sizeof (INT8), dim, fp);

  return (0);
}

int
read_r8array_fblk (FILE * fp, REAL8 * r8array, size_t dim)
{
  char block[4];
  int blocklen;

  blocklen = read_fblk (fp, block);
  if (!feof (fp))
    check_blocklen (dim * sizeof (REAL8), blocklen);
  else
    return (-1);

  read_r8array_bin (fp, r8array, dim);

  blocklen = read_fblk (fp, block);
  if (!feof (fp))
    check_blocklen (dim * sizeof (REAL8), blocklen);
  else
    return (-1);

  return (0);
}

int
read_r8array_bin (FILE * fp, REAL8 * r8array, size_t dim)
{
  fread (r8array, sizeof (REAL8), dim, fp);

  return (0);
}

int
read_fblk (FILE * fp, char *fortran_block)
{
  int blocklen;

  fread (fortran_block, sizeof (char), 4, fp);
  memcpy (&blocklen, fortran_block, 4);

  return (blocklen);
}

int
write_i8array_fblk (FILE * fp, INT8 * i8array, size_t dim)
{
  char block[4];
  int blocklen;

  blocklen = dim * sizeof (INT8);
  memcpy (&block, &blocklen, 4);

  write_fblk (fp, block);

  write_i8array_bin (fp, i8array, dim);

  write_fblk (fp, block);

  return (0);
}

int
write_i8array_bin (FILE * fp, INT8 * i8array, size_t dim)
{
  fwrite (i8array, sizeof (INT8), dim, fp);

  return (0);
}

int
write_r8array_fblk (FILE * fp, REAL8 * r8array, size_t dim)
{
  char block[4];
  int blocklen;

  blocklen = dim * sizeof (REAL8);
  memcpy (&block, &blocklen, 4);

  write_fblk (fp, block);

  write_r8array_bin (fp, r8array, dim);

  write_fblk (fp, block);

  return (0);
}

int
write_r8array_bin (FILE * fp, REAL8 * r8array, size_t dim)
{
  fwrite (r8array, sizeof (REAL8), dim, fp);

  return (0);
}

int
write_fblk (FILE * fp, char *fortran_block)
{

  fwrite (fortran_block, sizeof (char), 4, fp);

  return (0);
}

double
min_darray (double *parray, size_t dim)
{
  size_t irun = dim;
  double xmin;

  xmin = DBL_MAX * .9;
  while (irun--)
    {
      if (*parray < xmin)
	xmin = *parray;
      parray++;
    }

  return (xmin);
}

double
max_darray (double *parray, size_t dim)
{
  size_t irun = dim;
  double xmax;

  xmax = - DBL_MAX * .9;
  while (irun--)
    {
      if (*parray > xmax)
	xmax = *parray;
      parray++;
    }

  return (xmax);
}

double
mean_darray (double *parray, size_t dim)
{
  size_t irun = dim;
  double xmean;

  xmean = 0.0;
  while (irun--)
    {
      xmean += *parray++;
    }
  xmean /= dim;

  return (xmean);
}

void
set_Progname (char *progarg)
{
  extern char *Progname;

  Progname = strrchr (progarg, '/');
  if (Progname)
    Progname++;
  else
    Progname = progarg;
}

void
set_Commandline (int argc, char *argv[])
{
  extern char *Progname;
  extern char *Commandline;
  static char line[1024];
  int i;

  line[0] = '\0';
  if (Progname)
    strcat (line, Progname);
  else
    strcat (line, argv[0]);

  for (i = 1; i < argc; i++)
    {
      strcat (line, " ");
      strcat (line, argv[i]);
    }
  Commandline = line;
}

void
write_srv8_fblk (FILE * ostream, double *data, GRID * grid, SERVICE * serv)
{
  int dlen, doffset;
  int istep, ilev;
  int nstep, nlev, nlat, nlon;
  INT8 header[8];
  static int irec = 0;
  static int printheader = 0;
  extern int Verbose;

  nlon = grid->nlon;
  nlat = grid->nlat;
  nlev = grid->nlev;
  nstep = grid->nstep;


  dlen = nlon * nlat;

  header[0] = serv->icode;
  header[3] = serv->itime;
  header[4] = serv->nlon;
  header[5] = serv->nlat;
  header[6] = serv->idisp1;
  header[7] = serv->idisp2;

  if (Verbose && !printheader)
    {
      printheader = 1;
      fprintf (stderr, " Rec :     Date   Time Code  Level  Lon  Lat :"
	       "     Minimum        Mean     Maximum\n");
    }

  irec++;
  for (istep = 0; istep < nstep; istep++)
    {
      if (grid->step == NULL)
	header[2] = istep;
      else
	header[2] = grid->step[istep];

      for (ilev = 0; ilev < nlev; ilev++)
	{
	  if (grid->lev == NULL)
	    {
	      if (nlev <= 1)
		header[1] = 0;
	      else
		header[1] = ilev + 1;
	    }
	  else
	    header[1] = (int) grid->lev[ilev];

	  if (Verbose)
	    fprintf (stderr, "%4d : %8d %6d %4d%7d %4d %4d :", irec,
		     (int) header[2], (int) header[3], (int) header[0],
		     (int) header[1], (int) header[4], (int) header[5]);
	  /* write header */
	  write_i8array_fblk (ostream, header, 8);

	  doffset = get_offset (istep, ilev, 0, 0, nstep, nlev, nlat, nlon);

	  /* write data */
	  write_r8array_fblk (ostream, &data[doffset], dlen);

	  if (Verbose)
	    fprintf (stderr, "%12.6g%12.6g%12.6g\n",
		     min_darray (&data[doffset], dlen),
		     mean_darray (&data[doffset], dlen),
		     max_darray (&data[doffset], dlen));
	  irec++;
	}
    }
  irec--;
}

void
read_srv8_fblk (FILE * istream, double *data, GRID * grid, SERVICE * serv)
{
  int irec;
  int dlen, doffset;
  int istep, ilev;
  int nstep, nlev, nlat, nlon;
  INT8 header[8];
  static int printheader = 0;
  extern int Verbose;

  nlon = grid->nlon;
  nlat = grid->nlat;
  nlev = grid->nlev;
  nstep = grid->nstep;

  dlen = nlon * nlat;

  if (Verbose && !printheader)
    {
      printheader = 1;
      fprintf (stderr, " Rec :     Date   Time Code  Level  Lon  Lat :"
	       "     Minimum        Mean     Maximum\n");
    }

  irec = 1;
  for (istep = 0; istep < nstep; istep++)
    {
      for (ilev = 0; ilev < nlev; ilev++)
	{
	  /* read header */
	  read_i8array_fblk (istream, header, 8);
	  if (feof (istream))
	    return;
	  if (Verbose)
	    fprintf (stderr, "%4d : %8d %6d %4d%7d %4d %4d :", irec,
		     (int) header[2], (int) header[3], (int) header[0],
		     (int) header[1], (int) header[4], (int) header[5]);

	  if (irec == 1)
	    serv->icode = header[0];
	  if (istep == 0)
	    grid->lev[ilev] = header[1];

	  doffset = get_offset (istep, ilev, 0, 0, nstep, nlev, nlat, nlon);

	  /* read data */
	  read_r8array_fblk (istream, &data[doffset], dlen);

	  if (Verbose)
	    fprintf (stderr, "%12.6g%12.6g%12.6g\n",
		     min_darray (&data[doffset], dlen),
		     mean_darray (&data[doffset], dlen),
		     max_darray (&data[doffset], dlen));
	  irec++;
	}
    }
}
