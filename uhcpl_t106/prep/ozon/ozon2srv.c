#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"


#define  T21OZON     107520
#define  T30OZON     161280
#define  T42OZON     215040
#define  T63OZON     322560
#define  T106OZON    537600

#define  T21LAT          32
#define  T30LAT          48
#define  T42LAT          64
#define  T63LAT          96
#define  T106LAT        160

#define MAXLEV 35

int level[MAXLEV] = {
  18, 23, 31, 41, 55, 74, 98,
  131, 175, 233, 311, 414, 552, 737,
  982, 1310, 1750, 2330, 3100, 4140, 5520,
  7360, 9810, 13100, 17400, 23300, 31000, 41300,
  55100, 73500, 79000, 84900, 91200, 98000, 101100
};

void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s ozonfile  srvfile \n", Progname);
  fprintf (stderr,
	   "       convert ECHAM4 Ozone data from CRAY to SERVICE \n");

  return;
}

int
main (int argc, char *argv[])
{
  GRID grid;
  SERVICE serv;
  struct stat ifile_info;
  off_t file_bytes;
  int istep, ilev;
  int nstep, nlev, nlat, nlon;
  double *ozon;			/* pointer to variables read, global field   */
  FILE *istream;		/* echam4 unit 21 */
  FILE *ostream;
  size_t bozon;
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
      fprintf (stderr, "Can't open file %s for reading ! \n", ifile);
      exit (1);
    }

  stat (ifile, &ifile_info);
  file_bytes = ifile_info.st_size;

  if (file_bytes == T21OZON)
    nlat = 32;
  else if (file_bytes == T30OZON)
    nlat = 48;
  else if (file_bytes == T42OZON)
    nlat = 64;
  else if (file_bytes == T63OZON)
    nlat = 96;
  else if (file_bytes == T106OZON)
    nlat = 160;
  else
    {
      fprintf (stderr, "size of inputfile is unexpected");
      exit (1);
    }

  nlon = 1;
  nlev = MAXLEV;
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
  for (ilev = 0; ilev < nlev; ilev++)
    grid.lev[ilev] = level[ilev];
  for (istep = 0; istep < nstep; istep++)
    grid.step[istep] = 10000 + (istep + 1) * 100 + 1;

  /* read and convert unit.21 fields */

  ozon = read_ozon_cray (istream, &grid, &bozon);
  {
    size_t read_bytes;

    read_bytes = bozon * sizeof (UINT8);
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
  serv.itime = 0;
  serv.nlon = nlon;
  serv.nlat = nlat;
  serv.idisp1 = 0;
  serv.idisp2 = 0;

  write_srv8_fblk (ostream, ozon, &grid, &serv);

  return (0);
}
