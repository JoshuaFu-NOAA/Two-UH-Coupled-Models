#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"

#define MAXFLD 10000
void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s   srv-file  ync-file \n", Progname);
  fprintf (stderr, "convert optional global fields from SERVICE to netCDF\n");

  return;
}

int
main (int argc, char *argv[])
{
  SERVICE serv;
  GRID grid;
  FILE *istream;
  double *fld;
  char *ifile = NULL, *ofile = NULL;
  int nstep, nlev, nlat, nlon;
  int icode;
  INT8 header[8];
  int c;
  int nfields;
  int dlen, doffset=0;
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
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile);
      fprintf (stderr, "  > %s\n\n", ofile);
    }

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for reading ! \n", ifile);
      exit (1);
    }
  /*--------------------------------------------------------*/

  rewind (istream);


  for (nfields = 0; nfields < MAXFLD; nfields++)
    {
      /* read header : set the date to 99999 to break also if date =0 */
      header[2] = 99999;

      read_i8array_fblk (istream, header, 8);

      if (nfields == 0)
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
	  doffset = sizeof (double) * dlen + 8;
	}
      else
	{
	  if (header[2] != serv.idate || header[2] == 99999)
	    break;
	}

      fseek (istream, doffset, SEEK_CUR);
    }

  if (Debug)
    fprintf (stderr, " Found %d fields\n", nfields);

  icode = serv.icode;
  nlon = serv.nlon;
  nlat = serv.nlat;
  /* we use the levels for the different fields */
  nlev = nfields;
  nstep = 1;
  grid.nlon = nlon;
  grid.nlat = nlat;
  grid.nlev = nlev;
  grid.nstep = nstep;


  grid.lon = (double *) malloc (nlon * sizeof (double));
  grid.lat = (double *) malloc (nlat * sizeof (double));
  grid.lev = (double *) malloc (nlev * sizeof (double));

  tgrid_lon (grid.lon, nlon);
  tgrid_lat (grid.lat, nlat);

  fld = (double *) malloc (nstep * nfields * nlat * nlon * sizeof (double));

  rewind (istream);
  read_srv8_fblk (istream, fld, &grid, &serv);

  /* we have to set the # of the field */
  for (c = 0; c < nfields; c++)
    grid.lev[c] = (double) (c + 1);

  /* -----  create netcdf file and setup global attributes  ----- */

  write_fnc (&grid, ofile, fld, icode);


  return (0);
}
