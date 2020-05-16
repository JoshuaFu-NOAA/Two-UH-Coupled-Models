#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "pp.h"


void
usage (void)
{
  extern char *Progname;

  fprintf (stderr, "usage: %s  srvfile \n", Progname);

  return;
}

/*
typedef struct
{
  INT8 header[8];
  INT8 zcode;
  INT8 *code = &zcode;
  INT8 icode, ilevel, idate, itime, nlon, nlat, idisp1, idisp2;
}
ServiceHeader;

ServiceHeader srvhead;
*/            
int
main (int argc, char *argv[])
{
  FILE *istream;
  char *ifile = NULL;
  int irec, dlen;
  double *fld;
  INT8 header[8];
  int c;
  extern int Debug;


  set_Progname (argv[0]);

  while ((c = getopt(argc, argv, "i:o:dh")) != EOF)
    switch (c) {
    case 'd':
      Debug = 1;
      break;
    case 'i':
      ifile = optarg;
      break;
    case 'h':
      usage ();
      exit (0);
    }

  if (optind < argc)
      while (optind < argc)
	{
	  ifile = argv[optind++];
	}

  if (ifile == NULL)
    {
      usage ();
      exit (1);
    }

  if (Debug)
    {
      fprintf (stdout, "\n");
      fprintf (stdout, "  < %s \n\n", ifile);
    }


  if ((istream = fopen (ifile, "r")) == NULL)
    {
      fprintf (stderr, "Can't open file %s ! \n ", ifile);
      exit (1);
    }


  if (Debug)
    fprintf(stderr, " Rec :     Date   Time Code  Level Nlon Nlat :"
                    "     Minimum        Mean     Maximum\n");

  dlen = 0;
  irec = 1;
  while (1)
    {
      /* read header */
      read_i8array_fblk (istream, header, 8);

      if (feof(istream)) break;

      if (Debug)
	fprintf(stderr, "%4d : %8d %6d %4d%7d %4d %4d :",
		irec, (int)header[2], (int)header[3], (int)header[0],
                      (int)header[1], (int)header[4], (int)header[5]);

     if (irec == 1)
	{
	  dlen = header[4] * header[5];
	  fld  = (double *) malloc (dlen * sizeof (double));
	}
      irec++;

      /* read data */
      read_r8array_fblk (istream, fld, dlen);

      if (Debug)
	fprintf(stderr, "%12.6g%12.6g%12.6g\n",
	 	 min_darray (fld, dlen),
		mean_darray (fld, dlen),
		 max_darray (fld, dlen));
    }

  return (0);
}
