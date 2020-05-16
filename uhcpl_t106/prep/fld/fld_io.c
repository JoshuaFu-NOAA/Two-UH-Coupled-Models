#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <unistd.h>
#include <pwd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "prep.h"
#include "io_nc.h"

void
write_fnc (GRID * grid, char *ofile, double *fld, int icode)
{
  int nlatid, nlonid, nfldid;
  int nc_glid, nc_lonid;
  int io_var_id;
  int ngl, nlon, nfld;
  int IO_file_id;
  int ivar;
  time_t current_time;
  struct passwd *user;
  char ttype[] = "Mean";
  netCDF_file ncheader;		/* global accessible netCDF file information */
  char *varname = NULL, *varlname = NULL, *varunit = NULL;
  extern char *Commandline;
  extern int Debug;
  extern int MaxCode_fld;
  extern Table Table_fld[];
  extern Table Table_dim[];
  extern char institute[];

  ngl = grid->nlat;
  nlon = grid->nlon;
  nfld = grid->nlev;
  ncheader.climate = 1;

  strcpy (ncheader.nc_conventions, "COARDS");
  strcpy (ncheader.nc_model, "ECHAM");
  strcpy (ncheader.nc_institution, institute);
  strcpy (ncheader.nc_creation_program, Commandline);
  current_time = time (NULL);
  strcpy (ncheader.nc_creation_date, ctime (&current_time));
  ncheader.nc_creation_date[strlen (ncheader.nc_creation_date) - 1] = '\0';
  user = getpwuid (getuid ());
  strcpy (ncheader.nc_creation_user, user->pw_gecos);

#ifndef CRAY
  strcpy (ncheader.nc_binary_source, "IEEE");
#else
  strcpy (ncheader.nc_binary_source, "CRAY");
#endif

  if (ncheader.climate == 1)
    sprintf (ttype, "Climate data");
  else
    sprintf (ttype, "Monthly means");

  if (Debug)
    fprintf (stderr, " Type of data is :   %s \n", ttype);


  ncheader.nc_access_mode = NC_WRITE;
  IO_create (ofile, ncheader.nc_access_mode, &(ncheader.nc_file_id));


  IO_file_id = ncheader.nc_file_id;


  for (ivar = 0; ivar < MaxCode_fld; ivar++)
    {
      if (Debug)
	fprintf (stderr, "check code %d (%d) from codetable\n",
		 Table_fld[ivar].code, ivar);
      if (icode == Table_fld[ivar].code)
	{
	  varname = Table_fld[ivar].name;
	  varlname = Table_fld[ivar].longname;
	  varunit = Table_fld[ivar].unit;
	  strcpy (ncheader.nc_file_type, "FIELD file");
	  break;
	}
    }

  if (varname == NULL)
    {
      if (Debug)
	fprintf (stderr, "no code found in codetable ==> using code %d   \n",
		 Table_fld[0].code);
      icode = Table_fld[ivar].code;
      varname = Table_fld[0].name;
      varlname = Table_fld[0].longname;
      varunit = Table_fld[0].unit;
      strcpy (ncheader.nc_file_type, "FIELD file");
    }

  IO_def_dim (IO_file_id, Table_dim[0].name, nlon, &nlonid);
  IO_def_dim (IO_file_id, Table_dim[1].name, ngl, &nlatid);
  IO_def_dim (IO_file_id, "nfield", nfld, &nfldid);

  IO_put_att_text (IO_file_id, NC_GLOBAL, "Conventions",
		   strlen (ncheader.nc_conventions), ncheader.nc_conventions);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "model",
		   strlen (ncheader.nc_model), ncheader.nc_model);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "file_type",
		   strlen (ncheader.nc_file_type), ncheader.nc_file_type);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "source_type",
		   strlen (ncheader.nc_binary_source),
		   ncheader.nc_binary_source);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "title",
		   strlen (ncheader.nc_file_type), ncheader.nc_file_type);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "institution",
		   strlen (ncheader.nc_institution), ncheader.nc_institution);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "user",
		   strlen (ncheader.nc_creation_user),
		   ncheader.nc_creation_user);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "history",
		   strlen (ncheader.nc_creation_program),
		   ncheader.nc_creation_program);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "created",
		   strlen (ncheader.nc_creation_date),
		   ncheader.nc_creation_date);


  /* Coordinates, but not following the netCDF convention, */
  /* do it an other way ? */


  ncheader.nc_dims[0] = nlonid;
  IO_def_var (IO_file_id, Table_dim[0].name, NC_DOUBLE, 1, ncheader.nc_dims,
	      &nc_lonid);
  IO_put_att_text (IO_file_id, nc_lonid, "units", strlen (Table_dim[0].unit),
		   Table_dim[0].unit);
  IO_put_att_text (IO_file_id, nc_lonid, "long_name",
		   strlen (Table_dim[0].longname), Table_dim[0].longname);

  ncheader.nc_dims[0] = nlatid;
  IO_def_var (IO_file_id, Table_dim[1].name, NC_DOUBLE, 1, ncheader.nc_dims,
	      &nc_glid);
  IO_put_att_text (IO_file_id, nc_glid, "units", strlen (Table_dim[1].unit),
		   Table_dim[1].unit);
  IO_put_att_text (IO_file_id, nc_glid, "long_name",
		   strlen (Table_dim[1].longname), Table_dim[1].longname);


  /*  test the nfields , do we need this , COARDS ?  */
  ncheader.nc_dims[0] = nfldid;
  IO_def_var (IO_file_id, "nfield", NC_INT, 1, ncheader.nc_dims, &nfldid);
  IO_put_att_text (IO_file_id, nfldid, "long_name", 16, "Number of Fields");
  /* So long no more real data fields are defined ... */


  ncheader.nc_dims[0] = nfldid;
  ncheader.nc_dims[1] = nlatid;
  ncheader.nc_dims[2] = nlonid;


  IO_def_var (IO_file_id, varname, NC_DOUBLE, 3, ncheader.nc_dims,
	      &io_var_id);

  IO_put_att_text (IO_file_id, io_var_id, "long_name", strlen (varlname),
		   varlname);
  IO_put_att_text (IO_file_id, io_var_id, "units", strlen (varunit), varunit);

  IO_enddef (IO_file_id);

  IO_put_var_double (IO_file_id, nc_glid, grid->lat);
  IO_put_var_double (IO_file_id, nc_lonid, grid->lon);
  IO_put_var_double (IO_file_id, nfldid, grid->lev);

  IO_put_var_double (IO_file_id, io_var_id, fld);

  IO_close (IO_file_id);

  return;
}


/* #################################################################### */

double *
read_fnc (GRID * grid, char *ifilename, int *icode)
{
  int nlatid, nlonid, nfldid;
  int nc_glid, nc_lonid, nc_levid, nc_varid;
  int dlen;
  int IO_file_id, status;
  double *fld=0;
  int varcode, ivar;
  char *varname=0;
  netCDF_file ncheader;		/* global accessible netCDF file information */
  extern int Debug;
  extern int MaxCode_fld;
  extern Table Table_fld[];
  extern Table Table_dim[];

  ncheader.nc_access_mode = NC_NOWRITE;

  IO_open (ifilename, ncheader.nc_access_mode, &(ncheader.nc_file_id));

  IO_file_id = ncheader.nc_file_id;


  IO_inq_dimid (IO_file_id, Table_dim[0].name, &nlonid);
  IO_inq_dimid (IO_file_id, Table_dim[1].name, &nlatid);
  IO_inq_dimid (IO_file_id, "nfield", &nfldid);

  IO_inq_dimlen (IO_file_id, nlatid, &grid->nlat);
  IO_inq_dimlen (IO_file_id, nlonid, &grid->nlon);
  IO_inq_dimlen (IO_file_id, nfldid, &grid->nlev);

  IO_inq_varid (IO_file_id, Table_dim[0].name, &nc_lonid);
  IO_inq_varid (IO_file_id, Table_dim[1].name, &nc_glid);
  IO_inq_varid (IO_file_id, "nfield", &nc_levid);

  *icode = 0;
  dlen = (grid->nlat) * (grid->nlon) * (grid->nlev);
  /* no code defined => varname is FIELD */
  for (ivar = 0; ivar < MaxCode_fld; ivar++)
    {
      varcode = Table_fld[ivar].code;
      varname = Table_fld[ivar].name;

      if (Debug)
	fprintf (stderr, "  Try to read code:   %d (%s) \n\n", varcode,
		 varname);
      status = nc_inq_varid (IO_file_id, varname, &nc_varid);
      if (status == NC_NOERR)
	{
	  *icode = varcode;
	  fld = (double *) malloc (dlen * sizeof (double));
	  IO_get_var_double (IO_file_id, nc_varid, fld);
	  break;
	}
    }

  IO_close (IO_file_id);


  if (*icode == 0)
    {
      fprintf (stderr, "no code found\n");
      exit (-1);
    }

  if (Debug)
    fprintf (stderr, "  Code:   %d (%s) \n\n", *icode, varname);

  return (fld);
}

double *
read_fld_cray (FILE * istream, GRID * grid, size_t * rb)
{
  static double *f;
  size_t lrb;
  UINT8 *in, *out;
  int nlat, nlon, ilen;

  *rb = 0;
  nlat = grid->nlat;
  nlon = grid->nlon;
  ilen = nlat * nlon;

  /* now  go to begin of data section */

  rewind (istream);

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));

  lrb = fread (in, sizeof (UINT8), ilen, istream);
  *rb += lrb;

  out = (UINT8 *) f;

#ifndef CRAY
  cray2ieee (in, out, &ilen);
#else
  memcpy (f, in, ilen * sizeof (INT8));
#endif

  free (in);

  return (f);
}
