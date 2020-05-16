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

  fprintf (stderr, "usage: %s  snc-file  srv-file\n", Progname);
  fprintf (stderr,
	   "       convert ECHAM surface initial data from netCDF to SERVICE \n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  TIMES t;
  SERVICE serv;
  SINITIAL vars;
  FILE *ostream;		/* global file pointer for writing data */
  char *ifile = NULL, *ofile = NULL;
  int ivar;
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
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile);
      fprintf (stderr, "  > %s\n\n", ofile);
    }

  /* read netcdf file fields */

  read_snc (&grid, &t, ifile, &vars);

  /* ++++++    now check the date and write srv-file +++++++   */
  if (t.fdate < 1000000)
    {
      serv.idate = 19000000 + t.fdate;
      if (Verbose)
	fprintf (stderr, " Set year to four digits: %d\n\n",
		 (int) serv.idate);
    }
  else
    serv.idate = t.fdate;

  serv.itime = t.ftime / 100;
  serv.nlon = grid.nlon;
  serv.nlat = grid.nlat;
  serv.ilevel = 0;
  serv.idisp1 = 0;
  serv.idisp2 = 0;

  if ((ostream = fopen (ofile, "w")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for writing !\n", ofile);
      exit (1);
    }

  grid.nlev = 1;
  grid.nstep = 1;
  grid.lev = (double *) malloc (grid.nlev * sizeof (double));
  grid.lev[0] = 0;
  grid.step = (double *) malloc (grid.nstep * sizeof (double));
  grid.step[0] = serv.idate;

  for (ivar = 0; ivar < vars.ncodes; ivar++)
    {
      serv.icode = vars.code[ivar];
      write_srv8_fblk (ostream, vars.var[ivar], &grid, &serv);
    }

  return (0);
}
