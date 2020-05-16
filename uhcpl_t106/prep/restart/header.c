#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "prep.h"

int
header_info (FILE * file, HEADER * info)
{
  char *type;
  int nrecl1;
  unsigned long long ddr1_init;

  fread (&ddr1_init, sizeof (unsigned long long), 1, file);
  type = (char *) &ddr1_init;
  strncpy (info->def, type, 4);
  if (strncmp (info->def, "CRAY", 4) && strncmp (info->def, "IEEE", 4))
    {
      strcpy (info->def, "CRAY");
    }
  info->def[4] = '\0';
  printf ("Binary type           : %s\n\n", info->def);

  nrecl1 = ddr1_init & 0xffffffff;

  rewind (file);

  return (nrecl1);
}

int
read_header_records (FILE * file, HEADER * info, int nrecl1)
{
  int nrecl, nrecs;

  info->bddr1 = 0;
  info->bddr2 = 0;
  info->bddr3 = 0;
  info->bddr4 = 0;
  info->bddr5 = 0;

  /* extract first header record */

  info->ddr1 =
    (unsigned long long *) malloc (nrecl1 * sizeof (unsigned long long));
  info->bddr1 = fread (info->ddr1, sizeof (unsigned long long), nrecl1, file);
  nrecs = (int) info->ddr1[14];
  printf ("Header record 1 length: %d\n", nrecl1);

  /* extract second header record */

  if (nrecs > 1)
    {
      nrecl = (int) info->ddr1[2];
      printf ("Header record 2 length: %d\n", nrecl);
      info->ddr2 =
	(unsigned long long *) malloc (nrecl * sizeof (unsigned long long));
      info->bddr2 =
	fread (info->ddr2, sizeof (unsigned long long), nrecl, file);
    }

  /* extract third header record */

  if (nrecs > 2)
    {
      nrecl = (int) info->ddr2[2];
      printf ("Header record 3 length: %d\n", nrecl);
      info->ddr3 =
	(unsigned long long *) malloc (nrecl * sizeof (unsigned long long));
      info->bddr3 =
	fread (info->ddr3, sizeof (unsigned long long), nrecl, file);
    }

  /* extract fourth header record */

  if (nrecs > 3)
    {
      nrecl = (int) info->ddr3[2];
      printf ("Header record 4 length: %d\n", nrecl);
      info->ddr4 =
	(unsigned long long *) malloc (nrecl * sizeof (unsigned long long));
      info->bddr4 =
	fread (info->ddr4, sizeof (unsigned long long), nrecl, file);
      info->nrecl4 = nrecl;
    }

  /* extract fifth header record */

  if (nrecs > 4)
    {
      nrecl = (int) info->ddr4[2];
      printf ("Header record 5 length: %d\n", nrecl);
      info->ddr5 =
	(unsigned long long *) malloc (nrecl * sizeof (unsigned long long));
      info->bddr5 =
	fread (info->ddr5, sizeof (unsigned long long), nrecl, file);
    }

  return (nrecs);
}

int
read_special_header (FILE * file, HEADER * info, DIMENSIONS * model,
		     int nrecl1)
{
  int nrecs;

  info->bddr1 = 0;
  info->bddr2 = 0;
  info->bddr3 = 0;
  info->bddr4 = 0;
  info->bddr5 = 0;

  /* extract first header record */

  info->ddr1 =
    (unsigned long long *) malloc (nrecl1 * sizeof (unsigned long long));
  info->bddr1 = fread (info->ddr1, sizeof (unsigned long long), nrecl1, file);
  nrecs = 1;
  printf ("Header record 1 length: %d\n", nrecl1);

  /* extract information */

  model->nhtrac = (int) info->ddr1[1];
  model->nlon = (int) info->ddr1[2];
  model->nlev = (int) info->ddr1[3];

  printf ("Number of tracer          : %3d\n", model->nhtrac);
  printf ("Number of longitudes      : %3d\n", model->nlon);
  printf ("Number of vertical levels : %3d\n", model->nlev);

  return (nrecs);
}

void
extract_block_1 (unsigned long long *ddr1, TIMES * all, LABELS * run)
{

  /* extract information from record 1 */

  int i, start;
  size_t end;
  char *current_label, *tmp;

  all->fdate = ddr1[9];
  all->ftime = ddr1[8];
  all->vdate = ddr1[11];
  all->vtime = ddr1[10];

  current_label = (char *) (ddr1 + ddr1[1]);
  for (i = 0; i < 8; i++)
    {
      *(current_label + 79) = '\0';
      start = strspn (current_label, " ");
      if (start > 78)
	{
	  start = 0;
	  *current_label = '\0';
	}

      end = strlen (current_label);
      tmp = current_label + end - 1;
      while (*tmp == ' ')
	tmp--;
      *(++tmp) = '\0';

      if (strlen (current_label + start))
	{
	  printf ("%s\n", current_label + start);
	}
      strcpy (run->label[i], current_label + start);
      current_label += 80;
    }
  printf ("\n");

  return;
}

void
extract_block_3 (unsigned long long *ddr3, DIMENSIONS * model, GRID * grid,
		 char *def)
{
  int isp0, isp1, isp2;
  int i, nv, ivct, nc;
  unsigned long long *in, *out;

  model->nm = (int) ddr3[30];
  model->nn = (int) ddr3[31];
  model->nk = (int) ddr3[32];
  model->nkp1 = model->nk + 1;
  model->nmp1 = model->nm + 1;
  model->nnp1 = model->nn + 1;
  isp0 = model->nm + model->nn - model->nk;
  isp1 = model->nnp1 * model->nmp1;
  isp2 = isp0 * (isp0 + 1) / 2;
  model->nsp = isp1 - isp2;
  model->ngl = (int) ddr3[11];
  model->nhgl = model->ngl / 2;
  model->nlon = (int) ddr3[ddr3[14]];
  model->nlev = (int) ddr3[36];
  model->nlevp1 = model->nlev + 1;
  model->nlp2 = model->nlon + 2;
  grid->nvclev = (int) ddr3[38];
  nv = (int) ddr3[39];
  ivct = (int) ddr3[41];
  nc = nv * grid->nvclev;

  grid->lat = (double *) malloc (model->ngl * sizeof (double));
  grid->lon = (double *) malloc (model->nlon * sizeof (double));

  tgrid_lat (grid->lat, model->ngl);
  tgrid_lon (grid->lon, model->nlon);

  grid->vct = (double *) malloc (nc * sizeof (double));
  grid->vct_a = (double *) malloc (grid->nvclev * sizeof (double));
  grid->vct_b = (double *) malloc (grid->nvclev * sizeof (double));

  in = (unsigned long long *) ddr3 + ivct - 1;
  out = (unsigned long long *) grid->vct;

#ifndef CRAY
  if (!strncmp (def, "CRAY", 4))
    {
      cray2ieee (in, out, &nc);
    }
  else
    {
      memcpy (grid->vct, in, nc * sizeof (unsigned long long));
    }
#endif
  for (i = 0; i < grid->nvclev; i++)
    {
      grid->vct_a[i] = grid->vct[i];
      grid->vct_b[i] = grid->vct[grid->nvclev + i];
    }

  return;
}

void
extract_block_4 (unsigned long long *ddr4, DIMENSIONS * model, char *def)
{
  int istart;
  unsigned long long *in, *out;

  istart = (int) ddr4[30];

  in = (unsigned long long *) ddr4 + ddr4[7] - 1;
  out = (unsigned long long *) &model->dt;
#ifndef CRAY
  if (!strncmp (def, "CRAY", 4))
    {
      int one = 1;
      cray2ieee (in, out, &one);
    }
  else
    {
      memcpy (&model->dt, in, sizeof (unsigned long long));
    }
#endif

  model->nrd = (int) ddr4[31];
  model->nstep = (int) ddr4[ddr4[10] + 6];
  model->switches = (double *) malloc (model->nrd * sizeof (double));

  in = (unsigned long long *) ddr4 + istart - 1;
  out = (unsigned long long *) model->switches;

#ifndef CRAY
  if (!strncmp (def, "CRAY", 4))
    {
      cray2ieee (in, out, &(model->nrd));
    }
  else
    {
      memcpy (model->switches, in, model->nrd * sizeof (unsigned long long));
    }
#endif

  return;
}

void
compare (DIMENSIONS * d1, DIMENSIONS * d2, TIMES * t1, TIMES * t2, GRID * g1,
	 GRID * g2)
{
  printf ("\nComparison of base times and model dimensions:\n");
  printf
    ("------------------------------------------------------------------\n");
  printf
    ("                                            (fort.23) (fort.24)\n");
  printf
    ("------------------------------------------------------------------\n");
  printf ("Date (f)                                 :   %7d   %7d\n",
	  t1->fdate, t2->fdate);
  printf ("Time (f)                                 :   %7d   %7d\n",
	  t1->ftime, t2->ftime);
  printf ("Date (v)                                 :   %7d   %7d\n",
	  t1->vdate, t2->vdate);
  printf ("Time (v)                                 :   %7d   %7d\n",
	  t1->vtime, t2->vtime);
  printf
    ("------------------------------------------------------------------\n");
  printf ("Number of Gaussian latitudes             :   %5d   %7d\n", d1->ngl,
	  d2->ngl);
  printf ("Number of longitudes                     :   %5d   %7d\n",
	  d1->nlon, d2->nlon);
  printf ("Number of vertical levels                :   %5d   %7d\n",
	  d1->nlev, d2->nlev);
  printf ("Number of vertical transformation levels :   %5d   %7d\n",
	  g1->nvclev, g2->nvclev);
  printf
    ("------------------------------------------------------------------\n");
  printf
    ("Spectral truncation (fort.23)            : m = %d, n = %d, k = %d\n",
     d1->nm, d1->nn, d1->nk);
  printf
    ("Spectral truncation (fort.24)            : m = %d, n = %d, k = %d\n",
     d2->nm, d2->nn, d2->nk);
  printf
    ("------------------------------------------------------------------\n\n");

  return;
}
