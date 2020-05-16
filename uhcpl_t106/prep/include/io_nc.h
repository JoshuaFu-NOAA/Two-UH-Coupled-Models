
#include <netcdf.h>

void IO_create (const char *path, int cmode, int *idp);
void IO_open (const char *path, int omode, int *idp);
void IO_redef (int ncid);
void IO_enddef (int ncid);
void IO_close (int ncid);

void IO_def_dim (int ncid, const char *name, size_t len, int *idp);
void IO_inq_dimid (int ncid, const char *name, int *dimidp);
void IO_inq_dim (int ncid, int dimid, char *name, size_t * lengthp);
void IO_inq_dimname (int ncid, int dimid, char *name);
void IO_inq_dimlen (int ncid, int dimid, size_t * lengthp);
void IO_rename_dim (int ncid, int dimid, const char *name);
void IO_def_var (int ncid, const char *name, nc_type xtype, int ndims,
		 const int dimids[], int *varidp);
void IO_inq_varid (int ncid, const char *name, int *varidp);
void IO_inq_nvars (int ncid, int *nvarsp);
void IO_inq_var (int ncid, int varid, char *name, nc_type * xtypep,
		 int *ndimsp, int dimids[], int *nattsp);
void IO_inq_varname (int ncid, int varid, char *name);
void IO_inq_vartype (int ncid, int varid, nc_type * xtypep);
void IO_inq_varndims (int ncid, int varid, int *ndimsp);
void IO_inq_vardimid (int ncid, int varid, int dimids[]);
void IO_inq_varnatts (int ncid, int varid, int *nattsp);

/* void IO_put_var1_type() not yet defined */
/* void IO_put_vara_type() not yet defined */
/* void IO_put_vars_type() not yet defined */
/* void IO_put_varm_type() not yet defined */

void IO_copy_att (int ncid_in, int varid_in, const char *name, int ncid_out,
		  int varid_out);
void IO_put_var_text (int ncid, int varid, const char *tp);
void IO_put_var_uchar (int ncid, int varid, const unsigned char *up);
void IO_put_var_schar (int ncid, int varid, const signed char *cp);
void IO_put_var_short (int ncid, int varid, const short *sp);
void IO_put_var_int (int ncid, int varid, const int *ip);
void IO_put_var_long (int ncid, int varid, const long *lp);
void IO_put_var_float (int ncid, int varid, const float *fp);
void IO_put_var_double (int ncid, int varid, const double *dp);

void IO_get_var_text (int ncid, int varid, char *tp);
void IO_get_var_uchar (int ncid, int varid, unsigned char *up);
void IO_get_var_schar (int ncid, int varid, signed char *cp);
void IO_get_var_short (int ncid, int varid, short *sp);
void IO_get_var_int (int ncid, int varid, int *ip);
void IO_get_var_long (int ncid, int varid, long *lp);
void IO_get_var_float (int ncid, int varid, float *fp);
void IO_get_var_double (int ncid, int varid, double *dp);

void IO_put_att_text (int ncid, int varid, const char *name, size_t len,
		      const char *tp);
void IO_put_att_int (int ncid, int varid, const char *name, nc_type xtype,
		     size_t len, const int *ip);
void IO_put_att_float (int ncid, int varid, const char *name, nc_type xtype,
		       size_t len, const float *fp);
void IO_put_att_double (int ncid, int varid, const char *name, nc_type xtype,
			size_t len, const double *dp);

void IO_get_att_text (int ncid, int varid, char *name, char *tp);
void IO_get_att_int (int ncid, int varid, char *name, int *ip);
void IO_get_att_float (int ncid, int varid, char *name, float *fp);
void IO_get_att_double (int ncid, int varid, char *name, double *dp);

void IO_inq_att (int ncid, int varid, const char *name, nc_type * xtypep,
		 size_t * lenp);
void IO_inq_atttype (int ncid, int varid, const char *name, nc_type * xtypep);
void IO_inq_attlen (int ncid, int varid, const char *name, size_t * lenp);
void IO_inq_attname (int ncid, int varid, int attnum, char *name);
void IO_inq_attid (int ncid, int varid, const char *name, int *attnump);
