#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <netcdf.h>


extern int Debug;


void
IO_create (const char *path, int cmode, int *ncidp)
{
  int status;

  status = nc_create (path, cmode, ncidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_create : %s %d %d\n", path, cmode, *ncidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_create : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_open (const char *path, int omode, int *ncidp)
{
  int status;

  status = nc_open (path, omode, ncidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_open : %s %d %d\n", path, omode, *ncidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_open : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_redef (int ncid)
{
  int status;

  status = nc_redef (ncid);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_redef : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_enddef (int ncid)
{
  int status;

  status = nc_enddef (ncid);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_enddef : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_close (int ncid)
{
  int status;

  status = nc_close (ncid);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_close : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_def_dim (int ncid, const char *name, size_t len, int *dimidp)
{
  int status;

  status = nc_def_dim (ncid, name, len, dimidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_def_dim : %d %s %d %d\n", ncid, name, len, *dimidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_def_dim : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_dimid (int ncid, const char *name, int *dimidp)
{
  int status;

  status = nc_inq_dimid (ncid, name, dimidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_dimid : %d %s %d\n", ncid, name, *dimidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_dimid : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_dim (int ncid, int dimid, char *name, size_t * lengthp)
{
  int status;

  status = nc_inq_dim (ncid, dimid, name, lengthp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_dim : %d %d %s %d\n", ncid, dimid, name,
	     *lengthp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_dim : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_dimname (int ncid, int dimid, char *name)
{
  int status;

  status = nc_inq_dimname (ncid, dimid, name);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_dimname : %d %d %s\n", ncid, dimid, name);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_dimname : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_dimlen (int ncid, int dimid, size_t * lengthp)
{
  int status;

  status = nc_inq_dimlen (ncid, dimid, lengthp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_dimlen : %d %d %d\n", ncid, dimid, *lengthp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_dimlen : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_rename_dim (int ncid, int dimid, const char *name)
{
  int status;

  status = nc_rename_dim (ncid, dimid, name);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_rename_dim : %d %d %s\n", ncid, dimid, name);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_rename_dim : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_def_var (int ncid, const char *name, nc_type xtype, int ndims,
	    const int dimids[], int *varidp)
{
  int status;

  status = nc_def_var (ncid, name, xtype, ndims, dimids, varidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_def_var : %d %s %d %d %d\n", ncid, name, xtype,
	     ndims, *varidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_def_var : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_varid (int ncid, const char *name, int *varidp)
{

  int status;

  status = nc_inq_varid (ncid, name, varidp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_varid : %d %s %d\n", ncid, name, *varidp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_varid: %s\n", nc_strerror (status));
      exit (7);
    }
}
void
IO_inq_nvars (int ncid, int *nvarsp)
{
  int status;

  status = nc_inq_nvars (ncid, nvarsp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_nvars : %d %d \n", ncid, *nvarsp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_nvars : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_var (int ncid, int varid, char *name, nc_type * xtypep, int *ndimsp,
	    int dimids[], int *nattsp)
{
  int status;

  status = nc_inq_var (ncid, varid, name, xtypep, ndimsp, dimids, nattsp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_var : %d %d %s  \n", ncid, varid, name);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_var : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_varname (int ncid, int varid, char *name)
{
  int status;

  status = nc_inq_varname (ncid, varid, name);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_varname : %d %d %s\n", ncid, varid, name);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_varname : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_vartype (int ncid, int varid, nc_type * xtypep)
{
  int status;

  status = nc_inq_vartype (ncid, varid, xtypep);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_vartype : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_vartype : %s\n", nc_strerror (status));
      exit (7);
    }
}
void
IO_inq_varndims (int ncid, int varid, int *ndimsp)
{
  int status;

  status = nc_inq_varndims (ncid, varid, ndimsp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_varndims : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_varndims : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_vardimid (int ncid, int varid, int dimids[])
{
  int status;

  status = nc_inq_vardimid (ncid, varid, dimids);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_vardimid : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_vardimid: %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_varnatts (int ncid, int varid, int *nattsp)
{
  int status;

  status = nc_inq_varnatts (ncid, varid, nattsp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_varnatts : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_varnatts : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_text (int ncid, int varid, const char *tp)
{
  int status;

  status = nc_put_var_text (ncid, varid, tp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_text : %d %d %s \n", ncid, varid, tp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_text : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_uchar (int ncid, int varid, const unsigned char *up)
{
  int status;

  status = nc_put_var_uchar (ncid, varid, up);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_uchar : %d %d %s \n", ncid, varid, up);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_uchar : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_schar (int ncid, int varid, const signed char *cp)
{
  int status;

  status = nc_put_var_schar (ncid, varid, cp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_schar : %d %d %s \n", ncid, varid, cp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_schar : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_short (int ncid, int varid, const short *sp)
{
  int status;

  status = nc_put_var_short (ncid, varid, sp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_short : %d %d %hd \n", ncid, varid, *sp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_short : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_int (int ncid, int varid, const int *ip)
{
  int status;

  status = nc_put_var_int (ncid, varid, ip);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_int : %d %d %d \n", ncid, varid, *ip);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_int : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_long (int ncid, int varid, const long *lp)
{
  int status;

  status = nc_put_var_long (ncid, varid, lp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_long : %d %d %ld \n", ncid, varid, *lp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_long : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_float (int ncid, int varid, const float *fp)
{
  int status;

  status = nc_put_var_float (ncid, varid, fp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_float : %d %d %f \n", ncid, varid, *fp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_float : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_var_double (int ncid, int varid, const double *dp)
{
  int status;

  status = nc_put_var_double (ncid, varid, dp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_double : %d %d %f \n", ncid, varid, *dp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_var_double : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_text (int ncid, int varid, char *tp)
{
  int status;

  status = nc_get_var_text (ncid, varid, tp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_text : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_text : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_uchar (int ncid, int varid, unsigned char *up)
{
  int status;

  status = nc_get_var_uchar (ncid, varid, up);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_uchar : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_uchar : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_schar (int ncid, int varid, signed char *cp)
{
  int status;

  status = nc_get_var_schar (ncid, varid, cp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_schar : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_schar : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_short (int ncid, int varid, short *sp)
{
  int status;

  status = nc_get_var_short (ncid, varid, sp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_short : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_short : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_int (int ncid, int varid, int *ip)
{
  int status;

  status = nc_get_var_int (ncid, varid, ip);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_int : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_int : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_long (int ncid, int varid, long *lp)
{
  int status;

  status = nc_get_var_long (ncid, varid, lp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_long : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_long : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_float (int ncid, int varid, float *fp)
{
  int status;

  status = nc_get_var_float (ncid, varid, fp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_var_float : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_float : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_var_double (int ncid, int varid, double *dp)
{
  int status;

  status = nc_get_var_double (ncid, varid, dp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_var_double : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_var_double : %s\n", nc_strerror (status));
      exit (7);
    }
}
void
IO_copy_att (int ncid_in, int varid_in, const char *name, int ncid_out,
	     int varid_out)
{
  int status;

  status = nc_copy_att (ncid_in, varid_in, name, ncid_out, varid_out);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_copy_att : %d %d %s %d %d\n", ncid_in, varid_out,
	     name, ncid_out, varid_out);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_copy_att : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_att_text (int ncid, int varid, const char *name, size_t len,
		 const char *tp)
{
  int status;

  status = nc_put_att_text (ncid, varid, name, len, tp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_att_text : %d %d %s \n", ncid, varid, tp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_att_text : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_att_int (int ncid, int varid, const char *name, nc_type xtype,
		size_t len, const int *ip)
{
  int status;

  status = nc_put_att_int (ncid, varid, name, xtype, len, ip);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_att_int : %d %d %d \n", ncid, varid, *ip);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_att_int : %s \n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_att_float (int ncid, int varid, const char *name, nc_type xtype,
		  size_t len, const float *fp)
{
  int status;

  status = nc_put_att_float (ncid, varid, name, xtype, len, fp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_att_float : %d %d %f \n", ncid, varid, *fp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_att_float : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_put_att_double (int ncid, int varid, const char *name, nc_type xtype,
		   size_t len, const double *dp)
{
  int status;

  status = nc_put_att_double (ncid, varid, name, xtype, len, dp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_put_att_double : %d %d %f \n", ncid, varid, *dp);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_put_att_double : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_att_text (int ncid, int varid, char *name, char *tp)
{
  int status;

  status = nc_get_att_text (ncid, varid, name, tp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_att_text : %d %d \n", ncid, varid);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_att_text : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_att_int (int ncid, int varid, char *name, int *ip)
{
  int status;

  status = nc_get_att_int (ncid, varid, name, ip);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_att_int : %d %d \n", ncid, varid);
  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_att_int : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_att_float (int ncid, int varid, char *name, float *fp)
{
  int status;

  status = nc_get_att_float (ncid, varid, name, fp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_att_float : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_att_float : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_get_att_double (int ncid, int varid, char *name, double *dp)
{
  int status;

  status = nc_get_att_double (ncid, varid, name, dp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_get_att_double : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_get_att_double : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_att (int ncid, int varid, const char *name, nc_type * xtypep,
	    size_t * lenp)
{
  int status;

  status = nc_inq_att (ncid, varid, name, xtypep, lenp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_att : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_att : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_atttype (int ncid, int varid, const char *name, nc_type * xtypep)
{
  int status;

  status = nc_inq_atttype (ncid, varid, name, xtypep);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_atttype : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_atttype : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_attlen (int ncid, int varid, const char *name, size_t * lenp)
{
  int status;

  status = nc_inq_attlen (ncid, varid, name, lenp);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_attlen : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_attlen : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_attname (int ncid, int varid, int attnum, char *name)
{
  int status;

  status = nc_inq_attname (ncid, varid, attnum, name);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_attname : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_attname : %s\n", nc_strerror (status));
      exit (7);
    }
}

void
IO_inq_attid (int ncid, int varid, const char *name, int *attnump)
{
  int status;

  status = nc_inq_attid (ncid, varid, name, attnump);

  if (Debug || status != NC_NOERR)
    fprintf (stderr, "IO_inq_attid : %d %d \n", ncid, varid);

  if (status != NC_NOERR)
    {
      fprintf (stderr, "IO_inq_attid : %s\n", nc_strerror (status));
      exit (7);
    }
}
