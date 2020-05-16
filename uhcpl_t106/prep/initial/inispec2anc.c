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
	   "usage: %s  [-f file-format] cray_unit23-file  anc-file \n",
	   Progname);
  fprintf (stderr,
	   "       convert ECHAM4 atmosphere initial data from CRAY to netCDF \n");

  return;
}

int
main (int argc, char *argv[])
{
  HEADER header_23;
  DIMENSIONS model_23;
  TIMES all_23;
  LABELS run_23;
  GRID grid_23;
  int fileformat = 1;
  int c;
  extern int Debug;
  extern int Verbose;

/* global file pointer for read data */
  FILE *file_23;
/* global accessible netCDF file information */
  netCDF_file initial;
/* pointer to variables read */
  double *svo, *sd, *stp;	/* global spectral fields   */
  double *q;			/* gl1 buffer fields */
  double *dtmp;
/* check out number of bytes read during processing, compare later with file size */

  size_t bsvo, bsd, bstp;
  size_t bq;

  int nrecl1, nrecs;
  char *ifile23 = NULL, *ofile = NULL;


  set_Progname (argv[0]);
  set_Commandline (argc, argv);

  while ((c = getopt (argc, argv, "i:o:f:dhv")) != EOF)
    switch (c)
      {
      case 'd':
	Debug = 1;
	break;
      case 'v':
	Verbose = 1;
	break;
      case 'i':
	ifile23 = optarg;
	break;
      case 'o':
	ofile = optarg;
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
	ifile23 = argv[optind++];
	ofile = argv[optind++];
      }

  if (ifile23 == NULL || ofile == NULL)
    {
      usage ();
      exit (1);
    }

  if (Verbose)
    {
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile23);
      fprintf (stderr, "  > %s\n\n", ofile);
    }

  if ((file_23 = fopen (ifile23, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file:  %s ! \n ", ifile23);
      exit (1);
    }

  /* check fileformat for old or new ECHAM5 format */
  if (fileformat < 1 || fileformat > 2)
    {
      fprintf (stderr, "Unsupported output fileformat: %d\n", fileformat);
      exit (1);
    }

  /* file fort.24 should be the right one use with fort.23 !!!!!!! */
  if (Debug)
    fprintf (stderr, "\nExtract information from %s (fort.23):\n\n", ifile23);

  /* check type of file and total header length */

  nrecl1 = header_info (file_23, &header_23);

  /* extract header records and print information */

  nrecs = read_header_records (file_23, &header_23, nrecl1);

  nrecs = (int) header_23.ddr1[14];
  if (Debug)
    fprintf (stderr, "\nNumber of records     : %d\n", nrecs);

  /* extract information from record 1 */

  extract_block_1 (header_23.ddr1, &all_23, &run_23);

  /* extract information from record 3 */

  extract_block_3 (header_23.ddr3, &model_23, &grid_23, header_23.def);

  /* read and convert fort.23 fields */

  {
    int ilen1, ilen2, ilen3;
    int dlen, ilat, ilev, ioff, ioff2;
    int nlev, nlat, nlon, nlp2;

    ilen1 = model_23.nlev * 2 * model_23.nsp;
    ilen2 = model_23.nlevp1 * 2 * model_23.nsp;
    ilen3 = model_23.nlp2 * model_23.nlev * model_23.ngl;

    svo = read_sp (file_23, ilen1, &header_23, &bsvo);
    sd = read_sp (file_23, ilen1, &header_23, &bsd);
    stp = read_sp (file_23, ilen2, &header_23, &bstp);

    /* skip spectral humidity field */

    fseek (file_23, (long) ilen1 * sizeof (UINT8), SEEK_CUR);

    /* read latitude by latitude humidity (could be slow ?) or all at once ? */

    dtmp = read_gl (file_23, ilen3, &header_23, &model_23, &bq);


    nlev = model_23.nlev;
    nlat = model_23.ngl;
    nlon = model_23.nlon;
    nlp2 = model_23.nlp2;
    /*
    if (nlat == 48 && nlon == 96 && model_23.nn == 30)
      {
	model_23.nn = 31;
	model_23.nm = 31;
	model_23.nk = 31;
      }
      */
    if (Debug)
      fprintf (stderr, "removing zeroes at nlp2 for variable Q \n\n");

    /*  ============   deleting zeroes from Q at nlp2  ======================  */

    dlen = nlev * nlat * nlon;
    q = (double *) malloc (dlen * sizeof (double));
    for (ilev = 0; ilev < nlev; ilev++)
      {
	for (ilat = 0; ilat < nlat; ilat++)
	  {
	    ioff = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlon);
	    ioff2 = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlp2);
	    memcpy (&q[ioff], &dtmp[ioff2], nlon * sizeof (double));
	  }
      }

    /* test, if whole file is processed ... */
    {
      size_t read_bytes;
      off_t file_bytes;

      struct stat ifile_info;
      stat (ifile23, &ifile_info);
      file_bytes = ifile_info.st_size;

      read_bytes = (header_23.bddr1 + header_23.bddr2 + header_23.bddr3 +
		    header_23.bddr4 + header_23.bddr5 + bsvo + bsd + bstp +
		    ilen1 + bq) * sizeof (UINT8);
      if (Debug)
	fprintf (stderr, "Bytes read : %ld of %ld in file %s\n",
		 (long) read_bytes, (long) file_bytes, ifile23);
      if (read_bytes != file_bytes)
	{
	  fprintf (stderr,
		   "Bytes read not equal to data size !: ==> EXIT \n");
	  exit (1);
	}
    }
  }

  /* correct the 2-digit year to four digits */
  if (all_23.fdate < 1000000)
    {
      all_23.fdate = 19000000 + all_23.fdate;
      if (Debug)
	fprintf (stderr, "set year to four digits: %i ! \n", all_23.fdate);
    }
  if (all_23.vdate < 1000000)
    {
      all_23.vdate = 19000000 + all_23.vdate;
    }

  /* create netcdf file and setup global attributes */

  write_anc (&initial, &header_23, &model_23, &grid_23, &run_23,
	     &all_23, ofile, svo, sd, stp, q, fileformat);
  return (0);
}
