#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "prep.h"

double *
read_sp (FILE * file, int ilen, HEADER * info, size_t * rb)
{
  static double *f;
  unsigned long long *in, *out;

  in = (unsigned long long *) malloc (ilen * sizeof (unsigned long long));
  f = (double *) malloc (ilen * sizeof (double));
  *rb = fread (in, sizeof (unsigned long long), ilen, file);
  out = (unsigned long long *) f;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_gl (FILE * file, int ilen, HEADER * info, DIMENSIONS * model,
	 size_t * rb)
{
  static double *f;
  int n, npos, ilenpp;
  size_t lrb;
  unsigned long long *in, *out;

  *rb = 0;
  ilenpp = ilen / model->ngl;

  in = (unsigned long long *) malloc (ilen * sizeof (unsigned long long));
  f = (double *) malloc (ilen * sizeof (double));
  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb =
	fread (in + npos * ilenpp, sizeof (unsigned long long), ilenpp, file);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
#endif

  free (in);

  return (f);
}

double *
read_g3 (int no, long data_offset, FILE * file, int ilen, HEADER * info,
	 DIMENSIONS * model, size_t * rb)
{
  static double *f;
  int n, npos, ilenpp;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  *rb = 0;
  ilenpp = ilen / model->ngl;

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) ((no - 1) * ilenpp * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  noffset = (long) (14 * ilenpp * sizeof (unsigned long long));

  in = (unsigned long long *) malloc (ilen * sizeof (unsigned long long));
  f = (double *) malloc (ilen * sizeof (double));
  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb =
	fread (in + npos * ilenpp, sizeof (unsigned long long), ilenpp, file);
      if (n < model->ngl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

  if (no == 9)
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
  else
    {
#ifndef CRAY
      if (!strncmp (info->def, "CRAY", 4))
	{
	  cray2ieee (in, out, &ilen);
	}
      else
	{
	  memcpy (f, in, ilen * sizeof (unsigned long long));
	}
#endif
    }

  free (in);

  return (f);
}


double *
read_sgl (int no, long data_offset, FILE * file, HEADER * info,
	  DIMENSIONS * model, size_t * rb)
{
  static double *f;
  int n, npos, ilengl, ilen, ilenpp, ilenra, ilenrb, ilengl2;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  *rb = 0;

  /* Has to be defined from model dimensions in previous header ... */

  ilengl2 = model->nlon * model->nlev * model->ngl;
  ilengl =
    (2 * model->nlp2 * model->nlev +
     model->nlon * model->nlev * model->nhtrac + 1) * model->ngl;
  ilenpp = ilengl / model->ngl;
  ilenra = model->nlp2 * model->nlev;
  ilenrb = model->nlon * model->nlev * model->nhtrac;

  if (no == 3 && model->nhtrac == 0)
    {
      fprintf (stderr, "Request for tracer useless because no saved.\n");
      return (NULL);
    }

  printf ("\ndata sets    : %d\n", 2 + model->nhtrac);
  printf ("total length : %d\n", ilengl);
  printf ("record length : %d\n", ilenpp);

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  if (no == 1)
    {
      noffset = 0L;
      ilen = ilenra;
    }
  else if (no == 2)
    {
      noffset = (long) ilenra *sizeof (unsigned long long);
      ilen = ilenra;
    }
  else if (no == 3)
    {
      noffset = (long) (2 * ilenra) * sizeof (unsigned long long);
      ilen = ilenrb;
    }
  else
    {
      ilen = ilengl2;
      noffset =
	(long) (ilengl + (no - 4) * ilen) * sizeof (unsigned long long);
    }

  fseek (file, noffset, SEEK_CUR);	/* begin .... */

  if (no > 3)
    {
      in = (unsigned long long *) malloc (ilen * sizeof (unsigned long long));
      f = (double *) malloc (ilen * sizeof (double));
      lrb = fread (in, sizeof (unsigned long long), ilen, file);
      *rb += lrb;
    }
  else
    {
      noffset = (long) ((ilenpp - ilen) * sizeof (unsigned long long));	/* offset between latitudes */

      in =
	(unsigned long long *) malloc (ilen * model->ngl *
				       sizeof (unsigned long long));
      f = (double *) malloc (ilen * model->ngl * sizeof (double));

      for (n = 0; n < model->ngl; n++)
	{
	  npos =
	    ((n % 2) * (model->ngl - (n - 1) / 2 - 1) +
	     ((n + 1) % 2) * n / 2);
	  lrb =
	    fread (in + npos * ilen, sizeof (unsigned long long), ilen, file);
	  if (n < model->ngl - 1)
	    fseek (file, noffset, SEEK_CUR);
	  *rb += lrb;
	}
      ilen = ilen * model->ngl;
    }
  out = (unsigned long long *) f;


#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
#endif

  free (in);

  return (f);
}

double *
read_f1 (int no, long data_offset, FILE * file, DIMENSIONS * model,
	 HEADER * info, size_t * rb)
{
  static double *f;
  int n, ilenf1, ilenr, ilenpp;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  /* Has to be defined from model dimensions in header ... */

  ilenf1 = 6 * (model->nmp1 * 2 * model->nlev) * model->nhgl;
  ilenpp = 6 * (model->nmp1 * 2 * model->nlev);
  ilenr = model->nmp1 * 2 * model->nlev;

  if (no > 6)
    {
      fprintf (stderr, "requested field not available.\n");
      return (NULL);
    }
  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) ((no - 1) * ilenr * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  /* set skip space between single records */

  noffset = (long) ((ilenpp - ilenr) * sizeof (unsigned long long));

  in =
    (unsigned long long *) malloc (ilenr * model->nhgl *
				   sizeof (unsigned long long));
  f = (double *) malloc (ilenr * model->nhgl * sizeof (double));

  for (n = 0; n < model->nhgl; n++)
    {
      lrb = fread (in + n * ilenr, sizeof (unsigned long long), ilenr, file);
      if (n < model->nhgl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

  ilenr = ilenr * model->nhgl;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilenr);
    }
  else
    {
      memcpy (f, in, ilenr * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_f3 (int no, long data_offset, FILE * file, DIMENSIONS * model,
	 HEADER * info, size_t * rb)
{
  static double *f;
  int n, ilenf3, ilenr, ilenpp;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  /* Has to be defined from model dimensions in previous header ... */

  ilenf3 = 2 * (model->nmp1 * 2 * model->nlev) * model->nhgl;
  ilenpp = 2 * (model->nmp1 * 2 * model->nlev);
  ilenr = model->nmp1 * 2 * model->nlev;

  if (no > 2)
    {
      fprintf (stderr, "requested field not available.\n");
      return (NULL);
    }
  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) ((no - 1) * ilenr * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  /* set skip space between single records */

  noffset = (long) ((ilenpp - ilenr) * sizeof (unsigned long long));

  in =
    (unsigned long long *) malloc (ilenr * model->nhgl *
				   sizeof (unsigned long long));
  f = (double *) malloc (ilenr * model->nhgl * sizeof (double));

  for (n = 0; n < model->nhgl; n++)
    {
      lrb = fread (in + n * ilenr, sizeof (unsigned long long), ilenr, file);
      if (n < model->nhgl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

  ilenr = ilenr * model->nhgl;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilenr);
    }
  else
    {
      memcpy (f, in, ilenr * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_f4 (int no, long data_offset, FILE * file, DIMENSIONS * model,
	 HEADER * info, size_t * rb)
{
  static double *f;
  int n, ilenf4, ilen, ilenra, ilenrb, ilenpp;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  /* Has to be defined from model dimensions in previous header ... */

  ilenf4 =
    (4 * (model->nmp1 * 2 * model->nlevp1) + 4 * model->nlev) * model->nhgl;
  ilenpp = (4 * (model->nmp1 * 2 * model->nlevp1) + 4 * model->nlev);
  ilenra = model->nmp1 * 2 * model->nlevp1;
  ilenrb = model->nlev;

  if (no > 8)
    {
      fprintf (stderr, "requested field not available.\n");
      return (NULL);
    }

  if (no <= 4)
    {
      noffset = (long) ((no - 1) * ilenra * sizeof (unsigned long long));
      ilen = ilenra;
    }
  else
    {
      noffset =
	(long) ((4 * ilenra + (no - 5) * ilenrb) *
		sizeof (unsigned long long));
      ilen = ilenrb;
    }
  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  fseek (file, noffset, SEEK_CUR);

  /* set skip space between single records */

  noffset = (long) ((ilenpp - ilen) * sizeof (unsigned long long));

  in =
    (unsigned long long *) malloc (ilen * model->nhgl *
				   sizeof (unsigned long long));
  f = (double *) malloc (ilen * model->nhgl * sizeof (double));

  for (n = 0; n < model->nhgl; n++)
    {
      lrb = fread (in + n * ilen, sizeof (unsigned long long), ilen, file);
      if (n < model->nhgl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

  ilen = ilen * model->nhgl;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_g1r (long data_offset, FILE * file, int ilenpp, int ilenr, int start,
	  HEADER * info, DIMENSIONS * model, size_t * rb)
{
  static double *f;
  static int count = 0;
  int n, npos;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  count++;
  *rb = 0;

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) (start * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  noffset = (long) ((ilenpp - ilenr) * sizeof (unsigned long long));

  in =
    (unsigned long long *) malloc (ilenr * model->ngl *
				   sizeof (unsigned long long));
  f = (double *) malloc (ilenr * model->ngl * sizeof (double));

  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb =
	fread (in + npos * ilenr, sizeof (unsigned long long), ilenr, file);
      if (n < model->ngl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }

  printf ("fort.35: field %3d - %10d bytes\n", count, *rb);

  out = (unsigned long long *) f;

  ilenr = ilenr * model->ngl;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilenr);
    }
  else
    {
      memcpy (f, in, ilenr * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_g2r (int no, long data_offset, FILE * file, int ilen, HEADER * info,
	  DIMENSIONS * model, size_t * rb)
{
  static double *f;
  int n, npos, ilenpp;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  *rb = 0;
  ilenpp = ilen / model->ngl;

  if (no > 2)
    {
      fprintf (stderr, "Requested field does not exist.\n");
      return (NULL);
    }

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) ((no - 1) * ilenpp * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  noffset = (long) (ilenpp * sizeof (unsigned long long));

  in = (unsigned long long *) malloc (ilen * sizeof (unsigned long long));
  f = (double *) malloc (ilen * sizeof (double));
  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb =
	fread (in + npos * ilenpp, sizeof (unsigned long long), ilenpp, file);
      if (n < model->ngl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (unsigned long long));
    }
#endif
  free (in);

  return (f);
}

double *
read_g3r (long data_offset, FILE * file, int ilenpp, int ilenr, int start,
	  HEADER * info, DIMENSIONS * model, size_t * rb, int convtoieee)
{
  static double *f;
  static int count = 0;
  int n, npos;
  long noffset;
  size_t lrb;
  unsigned long long *in, *out;

  *rb = 0;
  count++;

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) (start * sizeof (unsigned long long));
  fseek (file, noffset, SEEK_CUR);

  noffset = (long) ((ilenpp - ilenr) * sizeof (unsigned long long));

  in =
    (unsigned long long *) malloc (ilenr * model->ngl *
				   sizeof (unsigned long long));
  f = (double *) malloc (ilenr * model->ngl * sizeof (double));

  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb =
	fread (in + npos * ilenr, sizeof (unsigned long long), ilenr, file);
      if (n < model->ngl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (unsigned long long *) f;

  ilenr = ilenr * model->ngl;

if (convtoieee == 0)
  {
    memcpy (f, in, ilenr * sizeof (UINT8));
  }
else
  {
#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilenr);
    }
  else
    {
      memcpy (f, in, ilenr * sizeof (UINT8));
    }
#else
      memcpy (f, in, ilenr * sizeof (UINT8));
#endif
  }
  free (in);

  return (f);
}
