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
	   "usage: %s  -t snc-template-file [-f file-format]  srv-file  snc-file \n",
	   Progname);
  fprintf (stderr,
	   "       convert ECHAM surface initial data from SERVICE to netCDF \n");

  return;
}


int calres(int nlat)
{
  int res = 0;

  res = (nlat*2 - 1) / 3;

  return (res);
}


int
main (int argc, char *argv[])
{
  DIMENSIONS d;
  LABELS run;
  TIMES t;
  SERVICE serv;
  GRID grid;
  SINITIAL vars;		/* field with structs of vars { icode,name,wert,da } */
  FILE *istream, *tstream;
  char *ifile = NULL, *tfile = NULL, *ofile = NULL;
  int nlat, nlon, nlev, nstep, irec, ivar;
  int dlen, tmp;
  double *dtmp;
  INT8 header[8];
  int fileformat = 1;
  int c;
  extern int Debug;
  extern int Verbose;


  set_Progname (argv[0]);
  set_Commandline (argc, argv);

  while ((c = getopt (argc, argv, "i:o:t:f:dhv")) != EOF)
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
      case 't':
	tfile = optarg;
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
	ifile = argv[optind++];
	ofile = argv[optind++];
      }

  if (ifile == NULL || ofile == NULL)
    {
      usage ();
      exit (1);
    }

  if (tfile == NULL && fileformat == 1)
    {
      usage ();
      fprintf (stderr,
	       "          you have to define a template-file to be able to write file - format:1 \n");
      exit (1);
    }

  if (Verbose)
    {
      fprintf (stderr, "\n");
      fprintf (stderr, "  < %s  \n", ifile);
      if (tfile)
	fprintf (stderr, "      <===  using %s as template \n", tfile);
      fprintf (stderr, "  > %s\n\n", ofile);
    }

  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s for reading ! \n", ifile);
      exit (1);
    }
  if (tfile)
    {
      if ((tstream = fopen (tfile, "r")) == NULL)
	{
	  fprintf (stderr, "Can't open file %s for reading ! \n", tfile);
	  exit (1);
	}
      fclose (tstream);
    }
  /* check fileformat */
  if (fileformat < 1 || fileformat > 2)
    {
      fprintf (stderr, "Unsupported file-format: %d\n", fileformat);
      exit (1);
    }

  tmp = 0;
  dlen = 0;
  irec = 0;
  rewind (istream);
  read_i8array_fblk (istream, header, 8);
  serv.icode = header[0];
  serv.ilevel = header[1];
  serv.idate = header[2];
  serv.itime = header[3];
  serv.nlon = header[4];
  serv.nlat = header[5];
  serv.idisp1 = header[6];
  serv.idisp2 = header[7];
  dlen = serv.nlon * serv.nlat;
  nlon = serv.nlon;
  nlat = serv.nlat;
  nlev = 1;
  nstep = 1;
  /* we need this for write_snc  */
  grid.nlon = nlon;
  grid.nlat = nlat;
  grid.nlev = nlev;
  grid.nstep = nstep;
  d.nlon = nlon;
  d.ngl = nlat;
  t.fdate = serv.idate;
  t.ftime = serv.itime;
  t.vdate = serv.idate;
  t.vtime = serv.itime;
  /* label[0] =9999 indikates not to write labels in netCDF file */
  /*  strcpy (run.label[0], "9999"); */
  if (!tfile)
    {
      strcpy (run.label[0], "Orographic data derived from GTOPO30 (USGS EROS data center)");
      strcpy (run.label[1], "Land surface data derived from AVHRR (USGS EROS data center)");
      strcpy (run.label[2], "FAO data from the soil map of the world");
      strcpy (run.label[3], "SNOW data from ECMWF Reanalysis");
      strcpy (run.label[4], " ");
      strcpy (run.label[5], " ");
      strcpy (run.label[6], " ");
      strcpy (run.label[7], " ");
    }

  d.nn = calres(nlat);
  d.nm = d.nn;
  d.nk = d.nn;

  grid.lon = (double *) malloc (nlon * sizeof (double));
  grid.lat = (double *) malloc (nlat * sizeof (double));
  grid.lev = (double *) malloc (nlev * sizeof (double));
  grid.step = (double *) malloc (nstep * sizeof (double));

  if (!tfile)
    {
      tgrid_lon (grid.lon, nlon);
      tgrid_lat (grid.lat, nlat);
    }

  dtmp = (double *) malloc (dlen * sizeof (double));
  rewind (istream);
  if (Verbose)
    {
      fprintf (stderr, " Used data from SERVICE file\n");
      fprintf (stderr, " ===========================\n");
    }
  for (ivar = 0; ivar < MAXCODE_PTsurfini; ivar++)
    {
      read_srv8_fblk (istream, dtmp, &grid, &serv);
      if (feof (istream))
	break;
      else
	{
	  vars.var[ivar] = (double *) malloc (dlen * sizeof (double));
	  vars.code[ivar] = serv.icode;
	  memcpy (vars.var[ivar], dtmp, dlen * sizeof (double));
	}

      if (ivar == 0)
	{
	  tmp = serv.idate;
	}
      else
	{
	  if (tmp != serv.idate)
	    {
	      fprintf (stderr, " dates of vars in file %s differ \n", ifile);
	      exit (1);
	    }
	}
    }
  vars.ncodes = ivar;
  free (dtmp);

  /* -----  create netcdf file and setup global attributes  ----- */

  write_snc (&grid, tfile, ofile, &vars, fileformat, &d, &run, &t);

  return (0);
}
