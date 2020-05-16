#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "prep.h"

void
swap_u_long_long_8 (unsigned long long *tni8)
/**************************************************************************
swap_u_long_8		swap an unsigned long long integer
***************************************************************************/
{
  *tni8 = (((*tni8 >> 56) & 0xff) | ((*tni8 & 0xff) << 56) |
	   ((*tni8 >> 40) & 0xff00) | ((*tni8 & 0xff00) << 40) |
	   ((*tni8 >> 24) & 0xff0000) | ((*tni8 & 0xff0000) << 24) |
	   ((*tni8 >> 8) & 0xff000000) | ((*tni8 & 0xff000000) << 8));
}

int
read_header_records (FILE * file, HEADER * info, int nrecl1)
{
  int nrecl, nrecs, i;
  extern int Debug;

  info->bddr1 = 0;
  info->bddr2 = 0;
  info->bddr3 = 0;
  info->bddr4 = 0;
  info->bddr5 = 0;

  /* extract first header record */

  info->ddr1 = (UINT8 *) malloc (nrecl1 * sizeof (UINT8));
  info->bddr1 = fread (info->ddr1, sizeof (UINT8), nrecl1, file);

  if (Debug)
    fprintf (stderr, "Header record 1 length: %d\n", nrecl1);

#ifdef i386
  if (Debug)
    fprintf (stderr, "Byte-swaping..\n");

  for (i = 0; i < nrecl1; i++)
    swap_u_long_long_8 ((unsigned long long *) (info->ddr1 + i));
#endif

  nrecs = (int) info->ddr1[14];


  /* extract second header record */

  if (nrecs > 1)
    {
      nrecl = (int) info->ddr1[2];
      if (Debug)
	fprintf (stderr, "Header record 2 length: %d\n", nrecl);
      info->ddr2 = (UINT8 *) malloc (nrecl * sizeof (UINT8));
      info->bddr2 = fread (info->ddr2, sizeof (UINT8), nrecl, file);

#ifdef i386
      if (Debug)
	fprintf (stderr, "Byte-swaping..\n");

      for (i = 0; i < nrecl; i++)
	swap_u_long_long_8 ((unsigned long long *) (info->ddr2 + i));
#endif
    }
  /* extract third header record */

  if (nrecs > 2)
    {
      nrecl = (int) info->ddr2[2];
      if (Debug)
	fprintf (stderr, "Header record 3 length: %d\n", nrecl);
      info->ddr3 = (UINT8 *) malloc (nrecl * sizeof (UINT8));
      info->bddr3 = fread (info->ddr3, sizeof (UINT8), nrecl, file);

#ifdef i386
      if (Debug)
	fprintf (stderr, "Byte-swaping..\n");

      for (i = 0; i < nrecl; i++)
	swap_u_long_long_8 ((unsigned long long *) (info->ddr3 + i));
#endif
    }

  /* extract fourth header record */

  if (nrecs > 3)
    {
      nrecl = (int) info->ddr3[2];
      if (Debug)
	fprintf (stderr, "Header record 4 length: %d\n", nrecl);
      info->ddr4 = (UINT8 *) malloc (nrecl * sizeof (UINT8));
      info->bddr4 = fread (info->ddr4, sizeof (UINT8), nrecl, file);
      info->nrecl4 = nrecl;
#ifdef i386
      if (Debug)
	fprintf (stderr, "Byte-swaping..\n");
      for (i = 0; i < nrecl; i++)
	swap_u_long_long_8 ((unsigned long long *) (info->ddr4 + i));
#endif
    }

  /* extract fifth header record */

  if (nrecs > 4)
    {
      nrecl = (int) info->ddr4[2];
      if (Debug)
	fprintf (stderr, "Header record 5 length: %d\n", nrecl);
      info->ddr5 = (UINT8 *) malloc (nrecl * sizeof (UINT8));
      info->bddr5 = fread (info->ddr5, sizeof (UINT8), nrecl, file);

#ifdef i386
      if (Debug)
	fprintf (stderr, "Byte-swaping..\n");
      for (i = 0; i < nrecl; i++)
	swap_u_long_long_8 ((unsigned long long *) (info->ddr5 + i));
#endif
    }

  return (nrecs);
}

void
extract_block_1 (UINT8 * ddr1, TIMES * all, LABELS * run)
{

  /* extract information from record 1 */

  int i, j, start;
  size_t end;
  char *current_label, *tmp;
  extern int Debug;

  all->fdate = ddr1[9];
  all->ftime = ddr1[8];
  all->vdate = ddr1[11];
  all->vtime = ddr1[10];

#ifdef i386
  /* swap_u_long_long_8(ddr1 + ddr1[1]); */
#endif

  current_label = (char *) (ddr1 + ddr1[1]);
  for (i = 0; i < 8; i++)
    {
#ifdef i386
      for (j = 0; j < 10; j++)
	{
	  swap_u_long_long_8 ((unsigned long long *) current_label + j);
	}
#endif
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
	  if (Debug)
	    fprintf (stderr, "%s\n", current_label + start);
	}
      strcpy (run->label[i], current_label + start);
      current_label += 80;
    }
  if (Debug)
    fprintf (stderr, "\n");

  return;
}

void
extract_block_3 (UINT8 * ddr3, DIMENSIONS * model, GRID * grid, char *def)
{
  int isp0, isp1, isp2;
  int i, nv, ivct, nc;
  UINT8 *in, *out;

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

  in = (UINT8 *) ddr3 + ivct - 1;
  out = (UINT8 *) grid->vct;

#ifndef CRAY
  if (!strncmp (def, "CRAY", 4))
    {
      cray2ieee (in, out, &nc);
    }
  else
    {
      memcpy (grid->vct, in, nc * sizeof (UINT8 *));
    }
#endif
  for (i = 0; i < grid->nvclev; i++)
    {
      grid->vct_a[i] = grid->vct[i];
      grid->vct_b[i] = grid->vct[grid->nvclev + i];
    }

  return;
}


int
header_info (FILE * file, HEADER * info)
{
  char *type;
  int nrecl1;
  UINT8 ddr1_init, *s;
  extern int Debug;


  fread (&ddr1_init, sizeof (UINT8), 1, file);
  type = (char *) &ddr1_init;

  strncpy (info->def, type, 4);
  if (strncmp (info->def, "CRAY", 4) && strncmp (info->def, "IEEE", 4))
    {
      strcpy (info->def, "CRAY");
    }
  if (Debug)
    fprintf (stderr, "Binary type           : %s\n\n", info->def);

  /*=================== new test =======================*/
#ifdef i386
  s = malloc (sizeof (ddr1_init));
  memcpy (s, &ddr1_init, sizeof (ddr1_init));

  swap_u_long_long_8 (s);

  ddr1_init = *s;

  if (Debug)
    fprintf (stderr,
	     "Host is little endian type, ==> Byte swaping is needed  \n");

#endif
  nrecl1 = ddr1_init & 0xffffffff;

  rewind (file);
  return (nrecl1);
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
