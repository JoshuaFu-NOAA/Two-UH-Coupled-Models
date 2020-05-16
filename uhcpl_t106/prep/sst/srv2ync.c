#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"


void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s  [-c code] [-y year]  srv-file  ync-file \n",
	   Progname);
  fprintf (stderr,
	   "convert annual cycle monthly mean data from SERVICE to netCDF\n");

  return;
}

int
main (int argc, char *argv[])
{
  SERVICE serv;
  GRID grid;
  FILE *istream;		/* echam4 unit 20 */
  double *data;			/* pointer to variables read */
  char *ifile = NULL, *ofile = NULL;
  int date;
  int ocode, oyear;
  int nlon, nlat, nlev, nstep;
  INT8 header[8];
  int c;
  extern int Debug;
  extern int Verbose;


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

  /* ---------  read year data from srv8 file  ------ */

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for reading ! \n", ifile);
      exit (1);
    }

  rewind (istream);
  read_i8array_fblk (istream, header, 8);

  if (ocode <= 0)
    ocode = header[0];

  date = header[2];

  nlon = header[4];
  nlat = header[5];
  nlev = 1;
  nstep = 12;

  grid.nlon = nlon;
  grid.nlat = nlat;
  grid.nlev = nlev;
  grid.nstep = nstep;

  grid.lon = (double *) malloc (nlon * sizeof (double));
  grid.lat = (double *) malloc (nlat * sizeof (double));
  grid.lev = (double *) malloc (nlev * sizeof (double));
  grid.step = (double *) malloc (nstep * sizeof (double));

  tgrid_lon (grid.lon, nlon);
  tgrid_lat (grid.lat, nlat);

  data = (double *) malloc (nstep * nlev * nlat * nlon * sizeof (double));

  rewind (istream);
  read_srv8_fblk (istream, data, &grid, &serv);

  if (oyear <= 0)
    oyear = date / 10000;

  /* create year netcdf file (ync) and setup global attributes */

  write_ync (&grid, ofile, data, oyear, ocode, 0);

  return (0);
}
