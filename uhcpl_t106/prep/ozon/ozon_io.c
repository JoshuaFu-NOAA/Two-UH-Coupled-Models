#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/stat.h>

#include "prep.h"
#include "io_nc.h"

void
write_onc (char *argv[], GRID * grid, char *ofile, double *ozon, int ocode,
	   char *dsource)
{
  int nlatid, nlonid, nlevid;
  int nc_glid, nc_lonid, nc_levid, ntimeid, nc_timeid;
  int io_var_id;
  int ngl, nlon, nlev;
  int IO_file_id;
  int ivar;
  time_t current_time;
  struct passwd *user;
  char ttype[] = "Mean";
  char date_unit[20];
  int itime[12] = { 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330 };
  netCDF_file ncheader;		/* global accessible netCDF file information */
  char *varname=0, *varlname=0, *varunit=0;
  extern char *Commandline;
  extern int Debug;
  extern int MaxCode_ozon;
  extern Table Table_ozon[];
  extern Table Table_dim[];
  extern char institute[];

  ngl = grid->nlat;
  nlon = grid->nlon;
  nlev = grid->nlev;

  ncheader.climate = 1;

  sprintf (date_unit, "days since 1-1-1");
  if (Debug)
    fprintf (stderr, "\n  Set start date as : %s \n ", date_unit);

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

  for (ivar = 0; ivar < MaxCode_ozon; ivar++)
    {
      if (Debug)
	fprintf (stderr, "check code %d (%d) from codetable\n",
		 Table_ozon[ivar].code, ivar);
      if (ocode == Table_ozon[ivar].code)
	{
	  varname = Table_ozon[ivar].name;
	  varlname = Table_ozon[ivar].longname;
	  varunit = Table_ozon[ivar].unit;
	  strcpy (ncheader.nc_file_type, "OZON file");
	  break;
	}
    }
  if (ivar == MaxCode_ozon)
    {
      fprintf (stderr, "can't find code %d in codetable\n", ocode);
      return;
    }

  if (Debug)
    fprintf (stderr, "  Code:   %d (%s) \n\n", ocode, varname);

  IO_def_dim (IO_file_id, Table_dim[0].name, nlon, &nlonid);
  IO_def_dim (IO_file_id, Table_dim[1].name, ngl, &nlatid);
  IO_def_dim (IO_file_id, "level", nlev, &nlevid);
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
  IO_put_att_text (IO_file_id, NC_GLOBAL, "data_source", strlen (dsource),
		   dsource);
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



  /*
   * Coordinates, but not following the netCDF convention,
   * do it an other way ?
   *
   */

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

  ncheader.nc_dims[0] = nlevid;
  IO_def_var (IO_file_id, "level", NC_DOUBLE, 1, ncheader.nc_dims, &nc_levid);
  IO_put_att_text (IO_file_id, nc_levid, "units", 2, "Pa");
  IO_put_att_text (IO_file_id, nc_levid, "long_name", 10, "Ozon Level");

  ncheader.nc_dims[0] = ntimeid;
  IO_def_var (IO_file_id, "time", NC_INT, 1, ncheader.nc_dims, &nc_timeid);
  IO_put_att_text (IO_file_id, nc_timeid, "units", strlen (date_unit),
		   date_unit);
  IO_put_att_text (IO_file_id, nc_timeid, "long_name", strlen (ttype), ttype);


  /* So long no more real data fields are defined ... */

  ncheader.nc_dims[0] = ntimeid;
  ncheader.nc_dims[1] = nlevid;
  ncheader.nc_dims[2] = nlatid;
  ncheader.nc_dims[3] = nlonid;

  IO_def_var (IO_file_id, varname, NC_DOUBLE, 4, ncheader.nc_dims,
	      &io_var_id);
  IO_put_att_text (IO_file_id, io_var_id, "long_name", strlen (varlname),
		   varlname);
  IO_put_att_text (IO_file_id, io_var_id, "units", strlen (varunit), varunit);

  IO_enddef (IO_file_id);

  IO_put_var_int (IO_file_id, nc_timeid, itime);

  IO_put_var_double (IO_file_id, nc_glid, grid->lat);
  IO_put_var_double (IO_file_id, nc_lonid, grid->lon);
  IO_put_var_double (IO_file_id, nc_levid, grid->lev);

  IO_put_var_double (IO_file_id, io_var_id, ozon);

  IO_close (IO_file_id);

  return;
}

/* #################################################################### */

double *
read_onc (GRID * grid, char *ifilename, char *date_unit, int *icode)
{
  int nlatid, nlonid, nlevid, ntimeid;
  int nc_glid, nc_lonid, nc_levid, nc_timeid;
  int nc_ozonid, ilen, ivar;
  int IO_file_id, status;
  double *ozon=0;
  int varcode;
  char *varname=0;
  netCDF_file ncheader;		/* global accessible netCDF file information */
  extern int Debug;
  extern int MaxCode_ozon;
  extern Table Table_ozon[];



  ncheader.nc_access_mode = NC_NOWRITE;

  IO_open (ifilename, ncheader.nc_access_mode, &(ncheader.nc_file_id));

  IO_file_id = ncheader.nc_file_id;

  IO_inq_dimid (IO_file_id, "lat", &nlatid);	/* (status != NC_NO... */
  IO_inq_dimid (IO_file_id, "lon", &nlonid);
  IO_inq_dimid (IO_file_id, "level", &nlevid);
  IO_inq_dimid (IO_file_id, "time", &ntimeid);

  IO_inq_dimlen (IO_file_id, nlatid, &grid->nlat);
  IO_inq_dimlen (IO_file_id, nlonid, &grid->nlon);
  IO_inq_dimlen (IO_file_id, nlevid, &grid->nlev);
  IO_inq_dimlen (IO_file_id, ntimeid, &grid->nstep);

  IO_inq_varid (IO_file_id, "lat", &nc_glid);
  IO_inq_varid (IO_file_id, "lon", &nc_lonid);
  IO_inq_varid (IO_file_id, "level", &nc_levid);
  IO_inq_varid (IO_file_id, "time", &nc_timeid);

  grid->lev = (double *) malloc (grid->nlev * sizeof (double));
  IO_get_var_double (IO_file_id, nc_levid, grid->lev);

  IO_get_att_text (IO_file_id, nc_timeid, "units", date_unit);


  for (ivar = 0; ivar < MaxCode_ozon; ivar++)
    {
      varcode = Table_ozon[ivar].code;
      varname = Table_ozon[ivar].name;
      if (Debug)
	fprintf (stderr, "  Try to read code:   %d (%s) \n\n", varcode,
		 varname);
      status = nc_inq_varid (IO_file_id, varname, &nc_ozonid);
      if (status == NC_NOERR)
	{
	  *icode = varcode;
	  ilen = (grid->nlev) * (grid->nlat) * (grid->nlon) * (grid->nstep);
	  ozon = (double *) malloc (ilen * sizeof (double));
	  IO_get_var_double (IO_file_id, nc_ozonid, ozon);
	  break;
	}
    }

  if (Debug)
    fprintf (stderr, "  Code:   %d (%s) \n\n", *icode, varname);

  IO_close (IO_file_id);

  return (ozon);
}

double *
read_ozon_cray (FILE * istream, GRID * grid, size_t * rb)
{
  static double *f;
  long doffset;
  size_t lrb;
  int npos;
  int istep, ilev, ilat, nlat, nlev, ilen, nstep, nlon;
  UINT8 *in, *out;

  *rb = 0;

  nlon = grid->nlon;
  nlat = grid->nlat;
  nlev = grid->nlev;
  nstep = grid->nstep;

  ilen = nstep * nlev * nlat * nlon;

  /* now  go to begin of data section */

  rewind (istream);

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));

  for (ilat = 0; ilat < nlat; ilat++)
    {
      npos = get_pplat (ilat, nlat);

      for (istep = 0; istep < nstep; istep++)
	{
	  for (ilev = 0; ilev < nlev; ilev++)
	    {
	      doffset =
		get_offset (istep, ilev, npos, 0, nstep, nlev, nlat, nlon);
	      lrb = fread (in + doffset, sizeof (UINT8), 1, istream);
	      *rb += lrb;
	    }
	}
    }
  out = (UINT8 *) f;

#ifndef CRAY
  cray2ieee (in, out, &ilen);
#else
  memcpy (f, in, ilen * sizeof (INT8));
#endif

  free (in);

  return (f);
}
