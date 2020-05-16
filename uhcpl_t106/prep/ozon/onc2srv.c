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

  fprintf (stderr, "usage: %s  netcdf-file  srv-file\n", Progname);
  fprintf (stderr,
	   "       convert ECHAM4 ozone data from netCDF to SERVICE \n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  SERVICE serv;
  FILE *ostream;		/* global file pointer for writing data */
  double *ozon;			/* pointer to variables read */
  int istep;
  int nstep, nlev, nlat, nlon;
  int iyear;
  int icode, ierr;
  char *ifile = NULL, *ofile = NULL, date_unit[20];
  int c;
  extern int Debug;
  extern int Verbose;


  set_Progname (argv[0]);
  set_Commandline (argc, argv);

  while ((c = getopt (argc, argv, "i:o:dhv")) != EOF)
    switch (c)
      {
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

  /* read netcdf file fields */

  ozon = read_onc (&grid, ifile, date_unit, &icode);

  nlon = grid.nlon;
  nlat = grid.nlat;
  nlev = grid.nlev;
  nstep = grid.nstep;

  if (Debug)
    {
      fprintf (stderr, "nlon:  %d\n", nlon);
      fprintf (stderr, "nlat:  %d\n", nlat);
      fprintf (stderr, "nlev:  %d\n", nlev);
      fprintf (stderr, "nstep: %d\n", nstep);
    }

  iyear = 0;

  ierr = sscanf (date_unit, "days since %4d-1-1", &iyear);

  if (ierr == EOF || iyear <= 0 || iyear > 9999)
    {
      fprintf (stderr, "unsupported date unit: %s\n", date_unit);
    }

  grid.step = (double *) malloc (nstep * sizeof (double));
  for (istep = 0; istep < nstep; istep++)
    grid.step[istep] = iyear * 10000 + (istep + 1) * 100 + 1;

  serv.icode = icode;
  serv.idate = grid.step[0];

  if (Debug)
    {
      if (iyear == 1)
	fprintf (stderr,
		 "  file contains climate data, date is set to:  %d \n",
		 (int) serv.idate);
      else
	fprintf (stderr,
		 "  file contains monthly means, date is set to:  %d \n",
		 (int) serv.idate);
    }

  if ((ostream = fopen (ofile, "w")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for writing ! \n", ofile);
      exit (1);
    }

  serv.itime = 0;
  serv.nlon = grid.nlon;
  serv.nlat = grid.nlat;
  serv.idisp1 = 0;
  serv.idisp2 = 0;

  write_srv8_fblk (ostream, ozon, &grid, &serv);

  return (0);
}
