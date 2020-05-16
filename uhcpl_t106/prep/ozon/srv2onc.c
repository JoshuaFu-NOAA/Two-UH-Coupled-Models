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

  fprintf (stderr,
	   "usage: %s  [-c code] [-a \"data source\"] srv-file  nc-file \n",
	   Progname);
  fprintf (stderr,
	   "       convert ECHAM4 Ozone data from SERVICE to netCDF\n");

  return;
}

int
main (int argc, char *argv[])
{
  SERVICE serv;
  GRID grid;
  FILE *istream;
  char *ifile = NULL, *ofile = NULL;
  char dsource[NC_MAX_NAME];
  int nstep, nlev, nlat, nlon;
  int irec, ocode=0;
  int dlen, doffset;
  int tmp;
  double *ozon;
  INT8 header[8];
  int c;
  extern int Debug;
  extern int Verbose;

  set_Progname (argv[0]);
  set_Commandline (argc, argv);
  strcpy (dsource, "unknown");

  while ((c = getopt (argc, argv, "c:a:i:o:dhv")) != EOF)
    switch (c)
      {
      case 'c':
	ocode = atoi (optarg);
	break;
      case 'a':
	strcpy (dsource, optarg);
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
      fprintf (stdout, "\n");
      fprintf (stdout, "  < %s  \n", ifile);
      fprintf (stdout, "  > %s\n\n", ofile);
    }


  dlen = 0;
  doffset = 0;
  tmp = 0;

  /* ---------  read ozon data from srv8 file  ------ */

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for reading ! \n", ifile);
      exit (1);
    }

  /*------------------ get nlev --------------------------------*/

  rewind (istream);
  irec = 0;
  for (irec = 0; irec < 99; irec++)
    {
      /* read header */
      read_i8array_fblk (istream, header, 8);

      if (irec == 0)
	{
	  serv.icode = header[0];
	  serv.ilevel = header[1];
	  serv.idate = header[2];
	  serv.itime = header[3];
	  serv.nlon = header[4];
	  serv.nlat = header[5];
	  serv.idisp1 = header[6];
	  serv.idisp2 = header[7];
	  dlen = serv.nlon * serv.nlat;
	  doffset = sizeof (double) * dlen + 4 * 2;
	}
      else
	{
	  if (header[2] != serv.idate)
	    break;
	}

      fseek (istream, doffset, SEEK_CUR);
    }

  if (Debug)
    fprintf (stderr, " Found %d levels\n", irec);

  if (ocode <= 0)
    ocode = serv.icode;

  nlon = serv.nlon;
  nlat = serv.nlat;
  nlev = irec;
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

  ozon = (double *) malloc (nstep * nlev * nlat * nlon * sizeof (double));

  rewind (istream);
  read_srv8_fblk (istream, ozon, &grid, &serv);

  /* -----  create netcdf file and setup global attributes  ----- */

  write_onc (argv, &grid, ofile, ozon, ocode, dsource);

  return (0);
}
