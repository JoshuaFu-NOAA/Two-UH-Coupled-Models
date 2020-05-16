#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"

#define  T21FLD      16384
#define  T30FLD      36864
#define  T42FLD      65536
#define  T63FLD     147456
#define  T106FLD    409600

#define  T21LAT         32
#define  T30LAT         48
#define  T42LAT         64
#define  T63LAT         96
#define  T106LAT       160

int nfield = 1;

void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s  fldfile  srvfile \n", Progname);
  fprintf (stderr, "convert optional global fields from CRAY to SERVICE\n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  SERVICE serv;
  struct stat ifile_info;
  off_t file_bytes;
  int nlon, ngl;
  FILE *istream;		/* echam4 unit */
  FILE *ostream;
  double *fld;			/* pointer to variables read, global field   */
  size_t bfld;
  char *ifile = NULL, *ofile = NULL;
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

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s ! \n ", ifile);
      exit (1);
    }

  stat (ifile, &ifile_info);
  file_bytes = ifile_info.st_size;

  if (file_bytes == T21FLD)
    ngl = 32;
  else if (file_bytes == T30FLD)
    ngl = 48;
  else if (file_bytes == T42FLD)
    ngl = 64;
  else if (file_bytes == T63FLD)
    ngl = 96;
  else if (file_bytes == T106FLD)
    ngl = 160;
  else
    {
      fprintf (stderr, "size of inputfile is unexpected");
      exit (1);
    }

  nlon = 2 * ngl;

  grid.nlat = ngl;
  grid.nlon = nlon;
  grid.nlev = 1;
  grid.nstep = 1;


  grid.lat = (double *) malloc (ngl * sizeof (double));
  grid.lon = (double *) malloc (nlon * sizeof (double));

  tgrid_lat (grid.lat, ngl);
  tgrid_lon (grid.lon, nlon);

  if (Debug)
    {
      fprintf (stderr, "nlon:  %d\n", grid.nlon);
      fprintf (stderr, "nlat:  %d\n", grid.nlat);
      fprintf (stderr, "nlev:  %d\n", grid.nlev);
      fprintf (stderr, "nstep: %d\n", grid.nstep);
    }

  /* read and convert unit.21 fields */

  fld = read_fld_cray (istream, &grid, &bfld);
  {
    size_t read_bytes;

    read_bytes = bfld * sizeof (UINT8);
    if (read_bytes != file_bytes)
      {
	fprintf (stderr, "Bytes read : %d of %d\n", read_bytes,
		 (int) file_bytes);
	exit (1);
      }
  }

  /* create srv file */

  if ((ostream = fopen (ofile, "w")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for writing ! \n", ofile);
      exit (1);
    }

  serv.icode = 99;
  serv.ilevel = 0;
  serv.nlon = grid.nlon;
  serv.nlat = grid.nlat;
  serv.idate = 0;
  serv.itime = 0;
  serv.idisp1 = 0;
  serv.idisp2 = 0;

  write_srv8_fblk (ostream, fld, &grid, &serv);


  return (0);
}
