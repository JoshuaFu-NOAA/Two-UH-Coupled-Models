#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prep.h"

#include <time.h>
#include <unistd.h>
#include <pwd.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "io_nc.h"

void
write_ync (GRID * grid, char *ofile, double *year, int oyear, int ocode,
	   int fromcray)
{
  int nlatid, nlonid;
  int nc_glid, nc_lonid, ntimeid, nc_timeid;
  int nc_yearid;
  int ivar;
  int ngl, nlon;
  int IO_file_id;
  time_t current_time;
  struct passwd *user;
  char ttype[99];
  char date_unit[20];
  int itime[12] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
  netCDF_file ncheader;		/* global accessible netCDF file information */
  char *varname = NULL, *varlname = NULL, *varunit = NULL;
  extern char *Commandline;
  extern int Debug;
  extern int MaxCode_year;
  extern Table Table_year[];
  /* extern int MaxCode_dim;  not used yet */
  extern Table Table_dim[];
  extern char institute[];


  ngl = grid->nlat;
  nlon = grid->nlon;

  if (oyear <= 0)
    {
      ncheader.climate = 1;
      oyear = 1;
    }
  else
    {
      ncheader.climate = 0;
    }
  sprintf (date_unit, "days since %4d-1-1", oyear);
  if (Debug)
    fprintf (stderr, "\n  Set start date as : %s \n ", date_unit);

  /* when reading from SRV check if this is climatedata */
  if (oyear == 1)
    ncheader.climate = 1;
  else
    ncheader.climate = 0;
  /* setup of global attributes for netCDF file (creation information) */

  strcpy (ncheader.nc_conventions, "COARDS");
  strcpy (ncheader.nc_model, "ECHAM");
  strcpy (ncheader.nc_institution, institute);
  strcpy (ncheader.nc_creation_program, Commandline);
  current_time = time (NULL);
  strcpy (ncheader.nc_creation_date, ctime (&current_time));
  ncheader.nc_creation_date[strlen (ncheader.nc_creation_date) - 1] = '\0';
  user = getpwuid (getuid ());
  strcpy (ncheader.nc_creation_user, user->pw_gecos);
  /* year2ync sets 1 for fromcray and srv2ync 0  */
  if (fromcray)
    strcpy (ncheader.nc_binary_source, "CRAY");
  else
    strcpy (ncheader.nc_binary_source, "IEEE");

  if (ncheader.climate == 1)
    sprintf (ttype, "Climate data");
  else
    sprintf (ttype, "Monthly means");

  if (Debug)
    fprintf (stderr, " Type of data is :   %s \n", ttype);

  ncheader.nc_access_mode = NC_WRITE;
  IO_create (ofile, ncheader.nc_access_mode, &(ncheader.nc_file_id));


  IO_file_id = ncheader.nc_file_id;

  for (ivar = 0; ivar < MaxCode_year; ivar++)
    {
      if (Debug)
	fprintf (stderr, "check code %d (%d) from codetable\n",
		 Table_year[ivar].code, ivar);
      if (ocode == Table_year[ivar].code)
	{
	  varname = Table_year[ivar].name;
	  varlname = Table_year[ivar].longname;
	  varunit = Table_year[ivar].unit;
	  strcpy (ncheader.nc_file_type, "Annual cycle file");
	  break;
	}
    }
  if (ivar == MaxCode_year)
    {
      fprintf (stderr, "can't find code %d in codetable\n", ocode);
      return;
    }

  if (Debug)
    fprintf (stderr, "  Code:   %d (%s) \n\n", ocode, varname);

  IO_def_dim (IO_file_id, Table_dim[0].name, nlon, &nlonid);
  IO_def_dim (IO_file_id, Table_dim[1].name, ngl, &nlatid);
  IO_def_dim (IO_file_id, "time", 12, &ntimeid);

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

  ncheader.nc_dims[0] = ntimeid;
  IO_def_var (IO_file_id, "time", NC_INT, 1, ncheader.nc_dims, &nc_timeid);
  IO_put_att_text (IO_file_id, nc_timeid, "units", strlen (date_unit),
		   date_unit);

  IO_put_att_text (IO_file_id, nc_timeid, "long_name", strlen (ttype), ttype);


  /* So long no more real data fields are defined ... */

  ncheader.nc_dims[2] = nlonid;
  ncheader.nc_dims[1] = nlatid;
  ncheader.nc_dims[0] = ntimeid;

  IO_def_var (IO_file_id, varname, NC_DOUBLE, 3, ncheader.nc_dims,
	      &nc_yearid);
  IO_put_att_text (IO_file_id, nc_yearid, "long_name", strlen (varlname),
		   varlname);
  IO_put_att_text (IO_file_id, nc_yearid, "units", strlen (varunit), varunit);

  IO_enddef (IO_file_id);

  IO_put_var_int (IO_file_id, nc_timeid, itime);

  IO_put_var_double (IO_file_id, nc_glid, grid->lat);
  IO_put_var_double (IO_file_id, nc_lonid, grid->lon);

  IO_put_var_double (IO_file_id, nc_yearid, year);

  IO_close (IO_file_id);

  return;
}

/* #################################################################### */

double *
read_ync (GRID * grid, char *ifilename, char *date_unit, int *icode)
{
  int nlatid, nlonid, ntimeid;
  int nc_glid, nc_lonid, nc_timeid;
  int nc_varid, dlen, ivar;
  int IO_file_id, status;
  double *year = NULL;
  int varcode;
  char *varname = NULL;
  netCDF_file ncheader;		/* global accessible netCDF file information */
  extern int Debug;
  extern int MaxCode_year;
  extern Table Table_year[];


  ncheader.nc_access_mode = NC_NOWRITE;

  IO_open (ifilename, ncheader.nc_access_mode, &(ncheader.nc_file_id));

  IO_file_id = ncheader.nc_file_id;

  IO_inq_dimid (IO_file_id, "lat", &nlatid);	/* (status != NC_NO... */
  IO_inq_dimid (IO_file_id, "lon", &nlonid);
  IO_inq_dimid (IO_file_id, "time", &ntimeid);

  IO_inq_dimlen (IO_file_id, nlatid, &grid->nlat);
  IO_inq_dimlen (IO_file_id, nlonid, &grid->nlon);
  IO_inq_dimlen (IO_file_id, ntimeid, &grid->nstep);

  IO_inq_varid (IO_file_id, "lat", &nc_glid);
  IO_inq_varid (IO_file_id, "lon", &nc_lonid);
  IO_inq_varid (IO_file_id, "time", &nc_timeid);

  IO_get_att_text (IO_file_id, nc_timeid, "units", date_unit);

  *icode = 0;
  dlen = (grid->nlat) * (grid->nlon) * (grid->nstep);

  for (ivar = 0; ivar < MaxCode_year; ivar++)
    {
      varcode = Table_year[ivar].code;
      varname = Table_year[ivar].name;
      if (Debug)
	fprintf (stderr, "  Try to read code:   %d (%s) \n\n", varcode,
		 varname);
      status = nc_inq_varid (IO_file_id, varname, &nc_varid);
      if (status == NC_NOERR)
	{
	  *icode = varcode;
	  year = (double *) malloc (dlen * sizeof (double));
	  IO_get_var_double (IO_file_id, nc_varid, year);
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

  return (year);
}
