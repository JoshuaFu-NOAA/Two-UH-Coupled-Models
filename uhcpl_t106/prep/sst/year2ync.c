#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"

#define  T21YEAR     196608
#define  T30YEAR     442368
#define  T42YEAR     786432
#define  T63YEAR    1769472
#define  T106YEAR   4915200

#define  T21LAT          32
#define  T30LAT          48
#define  T42LAT          64
#define  T63LAT          96
#define  T106LAT        160


void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s  [-c code] [-y year]  year-file  ync-file \n",
	   Progname);
  fprintf (stderr,
	   "convert annual cycle monthly mean data from CRAY to netCDF \n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  struct stat ifile_info;
  off_t file_bytes;
  int nlon, ngl;
  int ocode, oyear;
  FILE *istream;		/* echam4 unit 20 */
  double *year;			/* pointer to variables read, global field   */
  size_t byear;
  char *ifile = NULL, *ofile = NULL;
  int c;
  extern int Debug;
  extern int Verbose;
  int echam3format = 0;


  set_Progname (argv[0]);
  set_Commandline (argc, argv);

  ocode = 0;
  oyear = 0;

  while ((c = getopt (argc, argv, "c:y:i:o:dhv")) != EOF)
    switch (c)
      {
      case 'c':
	ocode = atoi (optarg);
	break;
      case 'y':
	oyear = atoi (optarg);
	break;
      case 'd':
	Debug = 1;
	break;
      case 'v':
	Verbose = 1;
	break;
      case 'i':
	ifile = optarg;
	break;
      case 'o':
	ofile = optarg;
	break;
      case 'h':
	usage ();
	exit (0);
      }

  if (optind < argc)
    while (optind < argc)
      {
	ifile = argv[optind++];
	ofile = argv[optind++];
      }

  if (ifile == NULL || ofile == NULL)
    {
      usage ();
      exit (1);
    }

  if (Verbose)
    {
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile);
      fprintf (stderr, "  > %s\n\n", ofile);
    }

  if (ocode <= 0)
    ocode = 77;

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s ! \n ", ifile);
      exit (1);
    }

  stat (ifile, &ifile_info);
  file_bytes = ifile_info.st_size;

  if (file_bytes == T21YEAR)
    ngl = 32;
  else if (file_bytes == T30YEAR)
    ngl = 48;
  else if (file_bytes == T42YEAR)
    ngl = 64;
  else if (file_bytes == T63YEAR)
    ngl = 96;
  else if (file_bytes == T106YEAR)
    ngl = 160;
  else
    {
      if (file_bytes == 3*T21YEAR)
	ngl = 32;
      else if (file_bytes == 3*T30YEAR)
	ngl = 48;
      else if (file_bytes == 3*T42YEAR)
	ngl = 64;
      else if (file_bytes == 3*T63YEAR)
	ngl = 96;
      else if (file_bytes == 3*T106YEAR)
	ngl = 160;
      else
	{
	  fprintf (stderr, "size of inputfile is unexpected\n");
	  exit (1);
	}
      echam3format = 1;
    }

  nlon = 2 * ngl;

  grid.nlat = ngl;
  grid.nlon = nlon;

  grid.lat = (double *) malloc (ngl * sizeof (double));
  grid.lon = (double *) malloc (nlon * sizeof (double));

  tgrid_lat (grid.lat, ngl);
  tgrid_lon (grid.lon, nlon);

  /* read and convert unit.20 fields */

  if ( echam3format )
    year = read_year3_cray (istream, &grid, &byear);
  else
    year = read_year_cray (istream, &grid, &byear);


  {
    size_t read_bytes;

    read_bytes = byear * sizeof (UINT8);
    if (read_bytes != file_bytes)
      {
	fprintf (stderr, "Bytes read : %ld of %ld\n", (long) read_bytes,
		 (long) file_bytes);
	exit (1);
      }
  }

  /* create year netcdf file (ync) and setup global attributes */

  write_ync (&grid, ofile, year, oyear, ocode, 1);

  return (0);
}
