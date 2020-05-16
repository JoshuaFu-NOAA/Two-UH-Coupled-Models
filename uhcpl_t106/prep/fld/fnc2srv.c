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

  fprintf (stderr, "usage: %s  netcdf-file  srvfile \n", Progname);
  fprintf (stderr, "convert optional global fields from netCDF to SERVICE\n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  SERVICE serv;
  FILE *ostream;		/* global file pointer for writing data */
  double *fld;			/* pointer to variables read */
  int iyear, month, day;
  int icode, nfield;
  char *ifile = NULL, *ofile = NULL;
  int c;
  extern int Debug;
  extern int Verbose;

  set_Progname (argv[0]);

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

  fld = read_fnc (&grid, ifile, &icode);

  nfield = grid.nlev;
  grid.lev = NULL;
  grid.nstep = 1;
  if (Debug)
    {
      if (nfield > 1)
	fprintf (stderr, " Found %d fields\n", nfield);
      fprintf (stderr, "nlon:  %d\n", grid.nlon);
      fprintf (stderr, "nlat:  %d\n", grid.nlat);
      fprintf (stderr, "nlev:  %d\n", grid.nlev);
      fprintf (stderr, "nstep: %d\n", grid.nstep);
    }


  iyear = 0;
  month = 1;
  day = 1;

  serv.icode = icode;
  serv.idate = iyear * 10000 + month * 100 + day;

  if (Debug)
    {
      if (iyear == 1)
	fprintf (stderr,
		 "  file contains climate data, date is set to:  %d \n",
		 (int) serv.idate);
      else
	fprintf (stderr,
		 "  file contains global field, date is set to:  %d \n",
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
  serv.ilevel = 0;
  serv.idisp1 = 0;
  serv.idisp2 = 0;

  write_srv8_fblk (ostream, fld, &grid, &serv);


  return (0);
}
