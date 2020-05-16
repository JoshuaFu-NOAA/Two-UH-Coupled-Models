#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prep.h"


double *
read_year3_cray (FILE * istream, GRID * grid, size_t * rb)
{
  static double *f;
  int n, npos;
  long noffset;
  size_t lrb;
  UINT8 *in, *out;
  int no, nlat, nlon, ilen;

  *rb = 0;
  nlat = grid->nlat;
  nlon = grid->nlon;
  ilen = nlat * nlon * 12;

  /* now  go to begin of data section */

  rewind (istream);

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));

  for (no = 0; no < 12; no++)
    {
      noffset = no * nlon * nlat;
      for (n = 0; n < nlat; n++)
	{
	  npos = ((n % 2) * (nlat - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
	  
	  fseek(istream, nlon*sizeof(UINT8), SEEK_CUR);
	  lrb = fread (in + npos * nlon + noffset, sizeof (UINT8), nlon, istream);
	  fseek(istream, nlon*sizeof(UINT8), SEEK_CUR);
	  *rb += 3*lrb;
	}
    }
  out = (UINT8 *) f;

#ifndef CRAY
  cray2ieee (in, out, &ilen);
#else
  memcpy (f, in, ilen * sizeof (INT8));
#endif

  free (in);

  return (f);
}

double *
read_year_cray (FILE * istream, GRID * grid, size_t * rb)
{
  static double *f;
  int n, npos;
  long noffset;
  size_t lrb;
  UINT8 *in, *out;
  int no, nlat, nlon, ilen;

  *rb = 0;
  nlat = grid->nlat;
  nlon = grid->nlon;
  ilen = nlat * nlon * 12;

  /* now  go to begin of data section */

  rewind (istream);

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));

  for (no = 0; no < 12; no++)
    {
      noffset = no * nlon * nlat;
      for (n = 0; n < nlat; n++)
	{
	  npos = ((n % 2) * (nlat - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
	  lrb =
	    fread (in + npos * nlon + noffset, sizeof (UINT8), nlon, istream);
	  *rb += lrb;
	}
    }
  out = (UINT8 *) f;

#ifndef CRAY
  cray2ieee (in, out, &ilen);
#else
  memcpy (f, in, ilen * sizeof (INT8));
#endif

  free (in);

  return (f);
}

void
write_year_cray (FILE * stream, GRID * grid, double *data)
{
  static double *f;
  int n, npos;
  long noffset;
  size_t lrb;
  UINT8 *in, *out;
  int no, nlat, nlon, ilen;

  nlat = grid->nlat;
  nlon = grid->nlon;
  ilen = nlat * nlon * 12;

  in = (UINT8 *) data;

  out = (UINT8 *) malloc (ilen * sizeof (UINT8));


#ifndef CRAY
  ieee2cray (in, out, &ilen);
#else
  memcpy (out, data, ilen * sizeof (INT8));
#endif

  for (no = 0; no < 12; no++)
    {
      noffset = no * nlon * nlat;
      for (n = 0; n < nlat; n++)
	{
	  npos = ((n % 2) * (nlat - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
	  lrb =
	    fwrite (out + npos * nlon + noffset, sizeof (UINT8), nlon, stream);
	}
    }

  free (in);
}
