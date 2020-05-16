#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"

void
usage (void)
{
  extern char *Progname;

  fprintf (stderr,
	   "usage: %s  [-H  cray_unit23-file ] [-f file-format] cray_unit24-file  snc-file \n",
	   Progname);
  fprintf (stderr,
	   "       convert ECHAM4 surface initial data from CRAY to netCDF \n");
  fprintf (stderr, "               -H : read header from cray_unit23-file\n");

  return;
}

void
correct_snow (double *sn, double *slm, DIMENSIONS * d)
{
  extern int Debug;
  int ngl, nlp2;
  int jl, jrow, index, ipo;

  ipo = 0;
  ngl = d->ngl;
  nlp2 = d->nlp2;

  if (ngl == 32)
    {
      if (Debug)
	fprintf (stderr, "\nMinor corrections to initial snow field !!! ");
      for (jl = 51; jl <= 61; jl++)
	{
	  jrow = 2;
	  index = (jrow - 1) * nlp2 + (jl - 1);
	  if (*(sn + index) != 10.)
	    {
	      *(sn + index) = 10.;
	      ipo++;
	    }
	}
      if (Debug)
	{
	  if (ipo)
	    fprintf (stderr, "\n  set %d snow points to 10.\n\n", ipo);
	  else
	    fprintf (stderr, " ==>  no changes\n\n");
	}
    }
  else if (ngl == 64)
    {
      if (Debug)
	fprintf (stderr, "\nMinor corrections to initial snow field !!! ");
      for (jl = 101; jl <= 120; jl++)
	{
	  jrow = 3;
	  index = (jrow - 1) * nlp2 + (jl - 1);
	  if (*(slm + index) >= 0.5)
	    {
	      if (*(sn + index) != 10.)
		{
		  *(sn + index) = 10.;
		  ipo++;
		}
	    }
	}
      for (jl = 101; jl <= 102; jl++)
	{
	  jrow = 4;
	  index = (jrow - 1) * nlp2 + (jl - 1);
	  if (*(sn + index) != 10.)
	    {
	      *(sn + index) = 10.;
	      ipo++;
	    }
	}
      jl = 25;
      jrow = 57;
      index = (jrow - 1) * nlp2 + (jl - 1);
      if (*(sn + index) != 10.)
	{
	  *(sn + index) = 10.;
	  ipo++;
	}
      if (Debug)
	{
	  if (ipo)
	    fprintf (stderr, "\n  set %d snow points to 10.\n\n", ipo);
	  else
	    fprintf (stderr, "==>   no changes\n\n");
	}
    }
}
int
main (int argc, char *argv[])
{
  HEADER header_23, header_24;
  DIMENSIONS model_23, model_24;
  TIMES all_23, all_24;
  LABELS run_23, run_24;
  GRID grid_23, grid_24;
  SINITIAL vars;
  int c, ivar;
  double *tmp;
  int dlen, ilat, ilev, ioff, ioff2;
  int nlev, nlat, nlon, nlp2;

  size_t bvar, tmpbvar;
  int fileformat = 1;
  int ndim, skip;
  int lrecl, nrecl1, nrecs;
  extern Table Table_surfini[MAXCODE_PTsurfini];
  extern int Debug;
  extern int Verbose;

/* global file pointer for read data */
  FILE *file_23=0, *file_24;
  char *ifile23 = NULL, *ifile24 = NULL, *ofile = NULL;


  set_Progname (argv[0]);
  set_Commandline (argc, argv);

  while ((c = getopt (argc, argv, "i:o:H:f:dhv")) != EOF)
    switch (c)
      {
      case 'd':
	Debug = 1;
	break;
      case 'v':
	Verbose = 1;
	break;
      case 'i':
	ifile24 = optarg;
	break;
      case 'o':
	ofile = optarg;
	break;
      case 'H':
	ifile23 = optarg;
	break;
      case 'f':
	fileformat = atoi (optarg);
	break;
      case 'h':
	usage ();
	exit (0);
      }

  if (optind < argc)
    while (optind < argc)
      {
	ifile24 = argv[optind++];
	ofile = argv[optind++];
      }

  if (ifile24 == NULL || ofile == NULL)
    {
      usage ();
      exit (1);
    }


  if (Verbose)
    {
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile24);
      fprintf (stderr, "  > %s\n\n", ofile);
      if (ifile23 != NULL)
	fprintf (stderr, "  header %s\n\n", ifile23);
    }


  if (ifile23)
    if ((file_23 = fopen (ifile23, "r")) == NULL)
      {
	fprintf (stderr, "Can't open file:  %s ! \n ", ifile23);
	exit (1);
      }

  if ((file_24 = fopen (ifile24, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file:  %s ! \n ", ifile24);
      exit (1);
    }

  /* check output fileformat for old or new ECHAM5 format */
  if (fileformat < 1 || fileformat > 2)
    {
      fprintf (stderr, "Unsupported output fileformat: %d\n", fileformat);
      exit (1);
    }

  /* file fort.24 should be the right one use with fort.23 !!!!!!! */
  if (Debug)
    fprintf (stderr, "\nExtract information from %s (fort.24):\n\n", ifile24);

  /* check type of file and total header length */

  nrecl1 = header_info (file_24, &header_24);

  /* extract header records and print information */

  nrecs = read_header_records (file_24, &header_24, nrecl1);

  nrecs = (int) header_24.ddr1[14];
  if (Debug)
    fprintf (stderr, "\nNumber of records     : %d\n", nrecs);
  lrecl = (int) header_24.ddr1[17];
  if (Debug)
    fprintf (stderr, "Total header length   : %ld\n\n",
	     (long) lrecl * sizeof (UINT8));

  /* extract information from record 1 */

  extract_block_1 (header_24.ddr1, &all_24, &run_24);

  /* extract information from record 3 */

  extract_block_3 (header_24.ddr3, &model_24, &grid_24, header_24.def);

  if (ifile23)
    {
      if (Debug)
	fprintf (stderr,
		 "\nExtract header information from %s (fort.23):\n\n",
		 ifile23);

      /* check type of file and total header length */

      nrecl1 = header_info (file_23, &header_23);

      /* extract header records and print information */

      nrecs = read_header_records (file_23, &header_23, nrecl1);

      nrecs = (int) header_23.ddr1[14];
      if (Debug)
	fprintf (stderr, "\nNumber of records     : %d\n", nrecs);
      lrecl = (int) header_23.ddr1[17];
      if (Debug)
	fprintf (stderr, "Total header length   : %ld\n\n",
		 (long) lrecl * sizeof (UINT8));

      /* extract information from record 1 */

      extract_block_1 (header_23.ddr1, &all_23, &run_23);

      /* extract information from record 3 */

      extract_block_3 (header_23.ddr3, &model_23, &grid_23, header_23.def);

      /* compare 23 & 24 ... */

      compare (&model_23, &model_24, &all_23, &all_24, &grid_23, &grid_24);
      memcpy (&model_24, &model_23, sizeof (model_23));
      memcpy (&all_24, &all_23, sizeof (all_23));
      memcpy (&grid_24, &grid_23, sizeof (grid_23));
      compare (&model_23, &model_24, &all_23, &all_24, &grid_23, &grid_24);
    }


  /* read and convert fort.24 fields */

  {
    int ilen1;
    long data_offset;

    ilen1 = model_24.nlp2 * model_24.ngl;

    data_offset = ftell (file_24);

    /* ============================  new  ========================= */
    tmpbvar = 0;
    bvar = 0;
    skip = 0;
    for (ivar = 0; ivar < 19; ivar++)
      {
	vars.code[ivar] = Table_surfini[ivar].code;

	tmp =
	  read_g3 (ivar + 1 - skip, data_offset, file_24, ilen1, &header_24,
		   &model_24, &tmpbvar);

	/* removing zeroes from nlp2  */
	nlev = model_24.nlev;
	nlat = model_24.ngl;
	nlon = model_24.nlon;
	nlp2 = model_24.nlp2;

	dlen = nlev * nlat * nlon;
	vars.var[ivar] = (double *) malloc (dlen * sizeof (double));
	for (ilev = 0; ilev < nlev; ilev++)
	  {
	    for (ilat = 0; ilat < nlat; ilat++)
	      {
		ioff = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlon);
		ioff2 = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlp2);
		memcpy (&vars.var[ivar][ioff], &tmp[ioff2],
			nlon * sizeof (double));
	      }
	  }

	bvar += tmpbvar;
	if (Debug)
	  printf (" reading var: %d  with %d bytes from file %s at pos:%d \n",
		  ivar + 1, tmpbvar, ifile24, ivar + 1 - skip);

	if (ivar == 8)
	  {
	    /* unpack directional orographic variances */

	    ndim = model_24.ngl * model_24.nlp2;

	    vars.var[ivar + 1] =
	      (double *) malloc ((size_t) ndim * sizeof (double));
	    vars.var[ivar + 2] =
	      (double *) malloc ((size_t) ndim * sizeof (double));
	    vars.var[ivar + 3] =
	      (double *) malloc ((size_t) ndim * sizeof (double));
	    vars.var[ivar + 4] =
	      (double *) malloc ((size_t) ndim * sizeof (double));

	    util_gwunpk (vars.var[ivar + 1], vars.var[ivar + 2],
			 vars.var[ivar + 3], vars.var[ivar + 4],
			 vars.var[ivar], &ndim);

	    vars.code[ivar + 1] = Table_surfini[ivar + 1].code;
	    vars.code[ivar + 2] = Table_surfini[ivar + 2].code;
	    vars.code[ivar + 3] = Table_surfini[ivar + 3].code;
	    vars.code[ivar + 4] = Table_surfini[ivar + 4].code;
	    if (Debug)
	      fprintf (stderr, "unpack directional orographic variances\n");
	    ivar += 4;
	    skip = 4;
	  }

      }

    vars.ncodes = ivar;

    /* test, if whole file is processed ... */
    {
      size_t read_bytes;
      off_t file_bytes;

      struct stat ifile_info;
      stat (ifile24, &ifile_info);
      file_bytes = ifile_info.st_size;

      read_bytes = (header_24.bddr1 + header_24.bddr2 + header_24.bddr3 +
		    header_24.bddr4 + header_24.bddr5 +
		    bvar) * sizeof (UINT8);
      if (Debug)
	fprintf (stderr, "Bytes read : %ld of %ld in file %s\n",
		 (long) read_bytes, (long) file_bytes, ifile24);
      if (read_bytes != file_bytes)
	{
	  fprintf (stderr, "Did not read all bytes from file %s ==> Exit \n",
		   ifile24);
	  exit (1);
	}
    }

  }

  /* Minor corrections to initial snow field */
  /* correct_snow (sn, slm, &model_24);  */
  correct_snow (vars.var[4], vars.var[5], &model_24);


  /* correct the 2-digit year to four digits */
  if (all_24.fdate < 1000000)
    {
      all_24.fdate = 19000000 + all_24.fdate;
      if (Debug)
	fprintf (stderr, "set year to four digits: %i ! \n", all_24.fdate);
    }
  if (all_24.vdate < 1000000)
    {
      all_24.vdate = 19000000 + all_24.vdate;
    }
  /* create netcdf file and setup global attributes */


  write_snc (&grid_24, NULL, ofile, &vars, fileformat, &model_24, &run_24,
	     &all_24);

  return (0);
}
