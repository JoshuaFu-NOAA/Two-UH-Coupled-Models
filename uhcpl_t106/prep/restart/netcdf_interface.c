#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/stat.h>

#include "prep.h"
#include "io_nc.h"


void ncheader (char *argv[], netCDF_file * restart, HEADER * info,
	       DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t);

/*  make_netcdf_31 writes data from file 31, 33, 34 to file 31 */
void
make_netcdf_31 (char *argv[], netCDF_file * restart, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *fsvo, double *favo, double *fsu,
		double *fau, double *fsv, double *fav, double *fsd,
		double *fad, double *fstp, double *fatp, double *fstpm,
		double *fatpm, double *fsu0, double *fau0, double *fsdu0,
		double *fadu0)
{
  int nglid, nhglid, nlonid, nlevid, nlevp1id, nspid, nvclevid;
  int n2id, nlp2id, nmp1id, nhtracid, nswitchesid;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid, nc_switchesid;
  int nc_fsvoid, nc_favoid, nc_fsuid, nc_fauid, nc_fsvid, nc_favid;
  int nc_fsdid, nc_fadid;
  int nc_fstpid, nc_fatpid, nc_fstpmid, nc_fatpmid, nc_fsu0id, nc_fau0id,
    nc_fsdu0id, nc_fadu0id;

  /* setup of global attributes for netCDF file (creation information) */

  restart->nc_access_mode = NC_WRITE;
  IO_create (ofile, restart->nc_access_mode, &(restart->nc_file_id));

  IO_def_dim (restart->nc_file_id, "ngl", d->ngl, &nglid);
  IO_def_dim (restart->nc_file_id, "nhgl", d->nhgl, &nhglid);
  IO_def_dim (restart->nc_file_id, "nlon", d->nlon, &nlonid);
  IO_def_dim (restart->nc_file_id, "nlp2", d->nlp2, &nlp2id);
  IO_def_dim (restart->nc_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (restart->nc_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (restart->nc_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (restart->nc_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (restart->nc_file_id, "n2", 2, &n2id);
  IO_def_dim (restart->nc_file_id, "nmp1", d->nmp1, &nmp1id);
  IO_def_dim (restart->nc_file_id, "nhtrac", d->nhtrac, &nhtracid);
  IO_def_dim (restart->nc_file_id, "nswitches", d->nrd, &nswitchesid);

  ncheader (argv, restart, info, d, g, run, t);

  /*
   * Coordinates, but not following the netCDF convention,
   * do it another way ?
   *
   */

  restart->nc_dims[0] = nglid;
  IO_def_var (restart->nc_file_id, "lat", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_glid);
  IO_put_att_text (restart->nc_file_id, nc_glid, "units", 9, "degrees_N");
  IO_put_att_text (restart->nc_file_id, nc_glid, "long_name", 17,
		   "Gaussian latitude");

  restart->nc_dims[0] = nlonid;
  IO_def_var (restart->nc_file_id, "lon", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_lonid);
  IO_put_att_text (restart->nc_file_id, nc_lonid, "units", 9, "degrees_E");
  IO_put_att_text (restart->nc_file_id, nc_lonid, "long_name", 9,
		   "longitude");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_a", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_b", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  restart->nc_dims[0] = nswitchesid;
  IO_def_var (restart->nc_file_id, "switches", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_switchesid);
  IO_put_att_text (restart->nc_file_id, nc_switchesid, "long_name", 19,
		   "diagnostic switches");


  /* gl1-buffer q, x, and xt fields */

  restart->nc_dims[0] = nhglid;
  restart->nc_dims[1] = nmp1id;
  restart->nc_dims[2] = n2id;
  restart->nc_dims[3] = nlevid;


  IO_def_var (restart->nc_file_id, "FSVO", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fsvoid);
  IO_put_att_text (restart->nc_file_id, nc_fsvoid, "long_name", 52,
		   "vorticity (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsvoid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FAVO", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_favoid);
  IO_put_att_text (restart->nc_file_id, nc_favoid, "long_name", 56,
		   "vorticity (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_favoid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FSU", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fsuid);
  IO_put_att_text (restart->nc_file_id, nc_fsuid, "long_name", 44,
		   "u (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsuid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FAU", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fauid);
  IO_put_att_text (restart->nc_file_id, nc_fauid, "long_name", 48,
		   "u (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fauid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FSV", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fsvid);
  IO_put_att_text (restart->nc_file_id, nc_fsvid, "long_name", 44,
		   "v (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsvid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FAV", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_favid);
  IO_put_att_text (restart->nc_file_id, nc_favid, "long_name", 48,
		   "v (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_favid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FSD", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fsdid);
  IO_put_att_text (restart->nc_file_id, nc_fsdid, "long_name", 52,
		   "divergence (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsdid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FAD", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fadid);
  IO_put_att_text (restart->nc_file_id, nc_fadid, "long_name", 56,
		   "divergence (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fadid, "units", 1, "");

  restart->nc_dims[0] = nhglid;
  restart->nc_dims[1] = nmp1id;
  restart->nc_dims[2] = n2id;
  restart->nc_dims[3] = nlevp1id;

  IO_def_var (restart->nc_file_id, "FSTP", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fstpid);
  IO_put_att_text (restart->nc_file_id, nc_fstpid, "long_name", 53,
		   "temperature (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fstpid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FATP", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fatpid);
  IO_put_att_text (restart->nc_file_id, nc_fatpid, "long_name", 57,
		   "temperature (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fatpid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FSTPM", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fstpmid);
  IO_put_att_text (restart->nc_file_id, nc_fstpmid, "long_name", 53,
		   "temperature (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fstpmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FATPM", NC_DOUBLE, 4, restart->nc_dims,
	      &nc_fatpmid);
  IO_put_att_text (restart->nc_file_id, nc_fatpmid, "long_name", 57,
		   "temperature (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fatpmid, "units", 1, "");

  restart->nc_dims[0] = nhglid;
  restart->nc_dims[1] = nlevid;

  IO_def_var (restart->nc_file_id, "FSU0", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_fsu0id);
  IO_put_att_text (restart->nc_file_id, nc_fsu0id, "long_name", 44,
		   "u0 (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsu0id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FAU0", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_fau0id);
  IO_put_att_text (restart->nc_file_id, nc_fau0id, "long_name", 48,
		   "u0 (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fau0id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FSDU0", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_fsdu0id);
  IO_put_att_text (restart->nc_file_id, nc_fsdu0id, "long_name", 44,
		   "? (symmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fsdu0id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "FADU0", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_fadu0id);
  IO_put_att_text (restart->nc_file_id, nc_fadu0id, "long_name", 48,
		   "? (antisymmetric part of Fourier coefficients)");
  IO_put_att_text (restart->nc_file_id, nc_fadu0id, "units", 1, "");


  IO_enddef (restart->nc_file_id);

  /* from file 31 */
  IO_put_var_double (restart->nc_file_id, nc_glid, g->lat);
  IO_put_var_double (restart->nc_file_id, nc_lonid, g->lon);
  IO_put_var_double (restart->nc_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (restart->nc_file_id, nc_vctbid, g->vct_b);
  IO_put_var_double (restart->nc_file_id, nc_switchesid, d->switches);

  IO_put_var_double (restart->nc_file_id, nc_fsvoid, fsvo);
  IO_put_var_double (restart->nc_file_id, nc_favoid, favo);

  IO_put_var_double (restart->nc_file_id, nc_fsuid, fsu);
  IO_put_var_double (restart->nc_file_id, nc_fauid, fau);

  IO_put_var_double (restart->nc_file_id, nc_fsvid, fsv);
  IO_put_var_double (restart->nc_file_id, nc_favid, fav);

  /* from file 33 */
  IO_put_var_double (restart->nc_file_id, nc_fsdid, fsd);
  IO_put_var_double (restart->nc_file_id, nc_fadid, fad);


  /* from file 34 */
  IO_put_var_double (restart->nc_file_id, nc_fstpid, fstp);
  IO_put_var_double (restart->nc_file_id, nc_fatpid, fatp);

  IO_put_var_double (restart->nc_file_id, nc_fstpmid, fstpm);
  IO_put_var_double (restart->nc_file_id, nc_fatpmid, fatpm);

  IO_put_var_double (restart->nc_file_id, nc_fsu0id, fsu0);
  IO_put_var_double (restart->nc_file_id, nc_fau0id, fau0);

  IO_put_var_double (restart->nc_file_id, nc_fsdu0id, fsdu0);
  IO_put_var_double (restart->nc_file_id, nc_fadu0id, fadu0);

  IO_close (restart->nc_file_id);

  return;
}

/*  ============================================================================== */

void
make_netcdf_32 (char *argv[], netCDF_file * restart, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *q, double *x, double *xt, double *lammp,
		double *phimp, double *sigmp)
{
  int nglid, nhglid, nlonid, nlevid, nlevp1id, nspid, nvclevid;
  int n2id, nlp2id, nmp1id, nhtracid, nswitchesid;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid, nc_switchesid;
  int nc_qid, nc_xid, nc_xtid;
  int nc_lammpid, nc_phimpid, nc_sigmpid;

  /* setup of global attributes for netCDF file (creation information) */

  restart->nc_access_mode = NC_WRITE;
  IO_create (ofile, restart->nc_access_mode, &(restart->nc_file_id));


  IO_def_dim (restart->nc_file_id, "ngl", d->ngl, &nglid);
  IO_def_dim (restart->nc_file_id, "nhgl", d->nhgl, &nhglid);
  IO_def_dim (restart->nc_file_id, "nlon", d->nlon, &nlonid);
  IO_def_dim (restart->nc_file_id, "nlp2", d->nlp2, &nlp2id);
  IO_def_dim (restart->nc_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (restart->nc_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (restart->nc_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (restart->nc_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (restart->nc_file_id, "n2", 2, &n2id);
  IO_def_dim (restart->nc_file_id, "nmp1", d->nmp1, &nmp1id);
  IO_def_dim (restart->nc_file_id, "nhtrac", d->nhtrac, &nhtracid);
  IO_def_dim (restart->nc_file_id, "nswitches", d->nrd, &nswitchesid);

  ncheader (argv, restart, info, d, g, run, t);

  /*
   * Coordinates, but not following the netCDF convention,
   * do it another way ?
   *
   */

  restart->nc_dims[0] = nglid;
  IO_def_var (restart->nc_file_id, "lat", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_glid);
  IO_put_att_text (restart->nc_file_id, nc_glid, "units", 9, "degrees_N");
  IO_put_att_text (restart->nc_file_id, nc_glid, "long_name", 17,
		   "Gaussian latitude");

  restart->nc_dims[0] = nlonid;
  IO_def_var (restart->nc_file_id, "lon", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_lonid);
  IO_put_att_text (restart->nc_file_id, nc_lonid, "units", 9, "degrees_E");
  IO_put_att_text (restart->nc_file_id, nc_lonid, "long_name", 9,
		   "longitude");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_a", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_b", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  restart->nc_dims[0] = nswitchesid;
  IO_def_var (restart->nc_file_id, "switches", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_switchesid);
  IO_put_att_text (restart->nc_file_id, nc_switchesid, "long_name", 19,
		   "diagnostic switches");


  /* gl1-buffer q, x and xt fields */

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "Q", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_qid);
  IO_put_att_text (restart->nc_file_id, nc_qid, "long_name", 17,
		   "specific humidity");
  IO_put_att_text (restart->nc_file_id, nc_qid, "units", 5, "kg/kg");

  IO_def_var (restart->nc_file_id, "X", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_xid);
  IO_put_att_text (restart->nc_file_id, nc_xid, "long_name", 20,
		   "liquid water content");
  IO_put_att_text (restart->nc_file_id, nc_xid, "units", 5, "kg/kg");

  if (d->nhtrac > 0)
    {
      restart->nc_dims[0] = nglid;
      restart->nc_dims[1] = nhtracid;
      restart->nc_dims[2] = nlevid;
      restart->nc_dims[3] = nlonid;

      IO_def_var (restart->nc_file_id, "xt", NC_DOUBLE, 4, restart->nc_dims,
		  &nc_xtid);
      IO_put_att_text (restart->nc_file_id, nc_qid, "long_name", 6, "tracer");
      IO_put_att_text (restart->nc_file_id, nc_qid, "units", 1, "");
    }

  /* gl1-buffer lammp, phimp and sigmp fields */

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlonid;

  IO_def_var (restart->nc_file_id, "LAMMP", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_lammpid);
  IO_put_att_text (restart->nc_file_id, nc_lammpid, "long_name", 5, "lammp");
  IO_put_att_text (restart->nc_file_id, nc_lammpid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "PHIMP", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_phimpid);
  IO_put_att_text (restart->nc_file_id, nc_phimpid, "long_name", 5, "phimp");
  IO_put_att_text (restart->nc_file_id, nc_phimpid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "SIGMP", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_sigmpid);
  IO_put_att_text (restart->nc_file_id, nc_sigmpid, "long_name", 5, "sigmp");
  IO_put_att_text (restart->nc_file_id, nc_sigmpid, "units", 1, "");

  IO_enddef (restart->nc_file_id);

  IO_put_var_double (restart->nc_file_id, nc_glid, g->lat);
  IO_put_var_double (restart->nc_file_id, nc_lonid, g->lon);
  IO_put_var_double (restart->nc_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (restart->nc_file_id, nc_vctbid, g->vct_b);
  IO_put_var_double (restart->nc_file_id, nc_switchesid, d->switches);

  IO_put_var_double (restart->nc_file_id, nc_qid, q);
  IO_put_var_double (restart->nc_file_id, nc_xid, x);
  if (d->nhtrac > 0)
    {
      IO_put_var_double (restart->nc_file_id, nc_xtid, xt);
    }
  IO_put_var_double (restart->nc_file_id, nc_lammpid, lammp);
  IO_put_var_double (restart->nc_file_id, nc_phimpid, phimp);
  IO_put_var_double (restart->nc_file_id, nc_sigmpid, sigmp);

  IO_close (restart->nc_file_id);

  return;
}

/*  ============================================================================== */

void
make_netcdf_35 (char *argv[], netCDF_file * restart, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *vom1, double *dm1, double *tm1,
		double *alpsm1, double *qm1, double *xm1, double *xtm1)
{
  int nglid, nhglid, nlonid, nlevid, nlevp1id, nspid, nvclevid;
  int n2id, nlp2id, nmp1id, nhtracid, nswitchesid;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid, nc_switchesid;
  int nc_vom1id, nc_dm1id, nc_tm1id, nc_alpsm1id, nc_qm1id, nc_xm1id,
    nc_xtm1id;

  /* setup of global attributes for netCDF file (creation information) */

  restart->nc_access_mode = NC_WRITE;
  IO_create (ofile, restart->nc_access_mode, &(restart->nc_file_id));

  IO_def_dim (restart->nc_file_id, "ngl", d->ngl, &nglid);
  IO_def_dim (restart->nc_file_id, "nhgl", d->nhgl, &nhglid);
  IO_def_dim (restart->nc_file_id, "nlon", d->nlon, &nlonid);
  IO_def_dim (restart->nc_file_id, "nlp2", d->nlp2, &nlp2id);
  IO_def_dim (restart->nc_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (restart->nc_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (restart->nc_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (restart->nc_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (restart->nc_file_id, "n2", 2, &n2id);
  IO_def_dim (restart->nc_file_id, "nmp1", d->nmp1, &nmp1id);
  IO_def_dim (restart->nc_file_id, "nhtrac", d->nhtrac, &nhtracid);
  IO_def_dim (restart->nc_file_id, "nswitches", d->nrd, &nswitchesid);

  ncheader (argv, restart, info, d, g, run, t);

  /*
   * Coordinates, but not following the netCDF convention,
   * do it another way ?
   *
   */

  restart->nc_dims[0] = nglid;
  IO_def_var (restart->nc_file_id, "lat", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_glid);
  IO_put_att_text (restart->nc_file_id, nc_glid, "units", 9, "degrees_N");
  IO_put_att_text (restart->nc_file_id, nc_glid, "long_name", 17,
		   "Gaussian latitude");

  restart->nc_dims[0] = nlonid;
  IO_def_var (restart->nc_file_id, "lon", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_lonid);
  IO_put_att_text (restart->nc_file_id, nc_lonid, "units", 9, "degrees_E");
  IO_put_att_text (restart->nc_file_id, nc_lonid, "long_name", 9,
		   "longitude");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_a", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_b", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  restart->nc_dims[0] = nswitchesid;
  IO_def_var (restart->nc_file_id, "switches", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_switchesid);
  IO_put_att_text (restart->nc_file_id, nc_switchesid, "long_name", 19,
		   "diagnostic switches");


  /* g1-buffer  fields vom1, dm1, tm1, alpsm1, qm1, xm1, and xtm1 */

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "VOM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_vom1id);
  IO_put_att_text (restart->nc_file_id, nc_vom1id, "long_name", 4, "vom1");
  IO_put_att_text (restart->nc_file_id, nc_vom1id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "DM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_dm1id);
  IO_put_att_text (restart->nc_file_id, nc_dm1id, "long_name", 3, "dm1");
  IO_put_att_text (restart->nc_file_id, nc_dm1id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "TM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_tm1id);
  IO_put_att_text (restart->nc_file_id, nc_tm1id, "long_name", 3, "tm1");
  IO_put_att_text (restart->nc_file_id, nc_tm1id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "QM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_qm1id);
  IO_put_att_text (restart->nc_file_id, nc_qm1id, "long_name", 3, "qm1");
  IO_put_att_text (restart->nc_file_id, nc_qm1id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "XM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_xm1id);
  IO_put_att_text (restart->nc_file_id, nc_xm1id, "long_name", 3, "xm1");
  IO_put_att_text (restart->nc_file_id, nc_xm1id, "units", 1, "");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "ALPSM1", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_alpsm1id);
  IO_put_att_text (restart->nc_file_id, nc_alpsm1id, "long_name", 6,
		   "alpsm1");
  IO_put_att_text (restart->nc_file_id, nc_alpsm1id, "units", 1, "");

  if (d->ntrac > 0)
    {
      restart->nc_dims[0] = nglid;
      restart->nc_dims[1] = nhtracid;
      restart->nc_dims[2] = nlevid;
      restart->nc_dims[3] = nlonid;

      IO_def_var (restart->nc_file_id, "XTM1", NC_DOUBLE, 4, restart->nc_dims,
		  &nc_xtm1id);
      IO_put_att_text (restart->nc_file_id, nc_xtm1id, "long_name", 4,
		       "xtm1");
      IO_put_att_text (restart->nc_file_id, nc_xtm1id, "units", 1, "");
    }

  IO_enddef (restart->nc_file_id);

  IO_put_var_double (restart->nc_file_id, nc_glid, g->lat);
  IO_put_var_double (restart->nc_file_id, nc_lonid, g->lon);
  IO_put_var_double (restart->nc_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (restart->nc_file_id, nc_vctbid, g->vct_b);
  IO_put_var_double (restart->nc_file_id, nc_switchesid, d->switches);

  IO_put_var_double (restart->nc_file_id, nc_vom1id, vom1);
  IO_put_var_double (restart->nc_file_id, nc_dm1id, dm1);
  IO_put_var_double (restart->nc_file_id, nc_tm1id, tm1);
  IO_put_var_double (restart->nc_file_id, nc_alpsm1id, alpsm1);
  IO_put_var_double (restart->nc_file_id, nc_qm1id, qm1);
  IO_put_var_double (restart->nc_file_id, nc_xm1id, xm1);
  if (d->ntrac > 0)
    {
      IO_put_var_double (restart->nc_file_id, nc_xtm1id, xtm1);
    }

  IO_close (restart->nc_file_id);

  return;
}

/*  ============================================================================== */

void
make_netcdf_36 (char *argv[], netCDF_file * restart, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *um1, double *vm1)
{
  int nglid, nhglid, nlonid, nlevid, nlevp1id, nspid, nvclevid;
  int n2id, nlp2id, nmp1id, nhtracid, nswitchesid;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid, nc_switchesid;
  int nc_um1id, nc_vm1id;

  /* setup of global attributes for netCDF file (creation information) */

  restart->nc_access_mode = NC_WRITE;
  IO_create (ofile, restart->nc_access_mode, &(restart->nc_file_id));


  IO_def_dim (restart->nc_file_id, "ngl", d->ngl, &nglid);
  IO_def_dim (restart->nc_file_id, "nhgl", d->nhgl, &nhglid);
  IO_def_dim (restart->nc_file_id, "nlon", d->nlon, &nlonid);
  IO_def_dim (restart->nc_file_id, "nlp2", d->nlp2, &nlp2id);
  IO_def_dim (restart->nc_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (restart->nc_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (restart->nc_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (restart->nc_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (restart->nc_file_id, "n2", 2, &n2id);
  IO_def_dim (restart->nc_file_id, "nmp1", d->nmp1, &nmp1id);
  IO_def_dim (restart->nc_file_id, "nhtrac", d->nhtrac, &nhtracid);
  IO_def_dim (restart->nc_file_id, "nswitches", d->nrd, &nswitchesid);

  ncheader (argv, restart, info, d, g, run, t);

  /*
   * Coordinates, but not following the netCDF convention,
   * do it another way ?
   *
   */

  restart->nc_dims[0] = nglid;
  IO_def_var (restart->nc_file_id, "lat", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_glid);
  IO_put_att_text (restart->nc_file_id, nc_glid, "units", 9, "degrees_N");
  IO_put_att_text (restart->nc_file_id, nc_glid, "long_name", 17,
		   "Gaussian latitude");

  restart->nc_dims[0] = nlonid;
  IO_def_var (restart->nc_file_id, "lon", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_lonid);
  IO_put_att_text (restart->nc_file_id, nc_lonid, "units", 9, "degrees_E");
  IO_put_att_text (restart->nc_file_id, nc_lonid, "long_name", 9,
		   "longitude");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_a", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_b", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  restart->nc_dims[0] = nswitchesid;
  IO_def_var (restart->nc_file_id, "switches", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_switchesid);
  IO_put_att_text (restart->nc_file_id, nc_switchesid, "long_name", 19,
		   "diagnostic switches");


  /* gl1-buffer q, x, and xt fields */

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "UM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_um1id);
  IO_put_att_text (restart->nc_file_id, nc_um1id, "long_name", 1, "um1");
  IO_put_att_text (restart->nc_file_id, nc_um1id, "units", 1, "");

  IO_def_var (restart->nc_file_id, "VM1", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_vm1id);
  IO_put_att_text (restart->nc_file_id, nc_vm1id, "long_name", 1, "vm1");
  IO_put_att_text (restart->nc_file_id, nc_vm1id, "units", 1, "");

  IO_enddef (restart->nc_file_id);

  IO_put_var_double (restart->nc_file_id, nc_glid, g->lat);
  IO_put_var_double (restart->nc_file_id, nc_lonid, g->lon);
  IO_put_var_double (restart->nc_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (restart->nc_file_id, nc_vctbid, g->vct_b);
  IO_put_var_double (restart->nc_file_id, nc_switchesid, d->switches);

  IO_put_var_double (restart->nc_file_id, nc_um1id, um1);
  IO_put_var_double (restart->nc_file_id, nc_vm1id, vm1);

  IO_close (restart->nc_file_id);

  return;
}

/*  ============================================================================== */

void
make_netcdf_37 (char *argv[], netCDF_file * restart, HEADER * info,
		DIMENSIONS * d, GRID * g, LABELS * run, TIMES * t,
		char *ofile, double *geospm, double *tsm, double *wsm,
		double *wlm, double *snm, double *slmm, double *az0m,
		double *albm, double *varpm, double *varorm, double *forestm,
		double *vgratm, double *vltm, double *wsmxm, double *faom,
		double *tdm, double *tdclm, double *apsm, double *aprlm,
		double *aprcm, double *aprsm, double *ustrgwm,
		double *vstrgwm, double *vdisgwm, double *aclcovm,
		double *temp2m, double *dew2m, double *wind10m, double *u10m,
		double *v10m, double *sradsm, double *tradsm, double *srad0m,
		double *trad0m, double *vdism, double *ustrm, double *vstrm,
		double *ahfsm, double *evapm, double *ahflm, double *wlm1m,
		double *tsm1m, double *tdm1m, double *wsm1m, double *wdm1m,
		double *snm1m, double *emterm, double *trsolm,
		double *runoffm, double *srad0um, double *sradsum,
		double *tradsum, double *albedom, double *tsurfm,
		double *tsnm, double *td3m, double *td4m, double *td5m,
		double *tsnm1m, double *td3m1m, double *td4m1m,
		double *td5m1m, double *tdclm1m, double *seaicem,
		double *sicedm, double *ustar3m, double *teffm, double *glacm,
		double *aclcm, double *aclcacm, double *snmelm,
		double *runtocm, double *tslinm, double *dsnacm,
		double *t2maxm, double *t2minm, double *tsmaxm,
		double *tsminm, double *wimaxm, double *topmaxm,
		double *aclcvm, double *qvim, double *alwcvim,
		double *runlndm, double *rgcgnm, double *sodifm,
		double *srafsm, double *trafsm, double *sraf0m,
		double *traf0m, double *emtefm, double *trsofm,
		double *tclfsm, double *sclfsm, double *tclf0m,
		double *sclf0m, double *tkem, double *tkem1m, double *auxil1m,
		double *auxil2m)
{
  int nglid, nhglid, nlonid, nlevid, nlevp1id, nspid, nvclevid;
  int n2id, nlp2id, nmp1id, nhtracid, nswitchesid;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid, nc_switchesid;

  int nc_dummyid;

  int nc_geospmid, nc_tsmid, nc_wsmid, nc_wlmid, nc_snmid, nc_slmmid,
    nc_az0mid, nc_albmid, nc_varpmid;
  int nc_varormid, nc_forestmid, nc_vgratmid, nc_vltmid, nc_wsmxmid,
    nc_faomid, nc_tdmid, nc_tdclmid;
  int nc_apsmid, nc_aprlmid, nc_aprcmid, nc_aprsmid, nc_ustrgwmid,
    nc_vstrgwmid, nc_vdisgwmid;
  int nc_aclcovmid, nc_temp2mid, nc_dew2mid, nc_wind10mid, nc_u10mid,
    nc_v10mid, nc_sradsmid;
  int nc_tradsmid, nc_srad0mid, nc_trad0mid, nc_vdismid, nc_ustrmid,
    nc_vstrmid, nc_ahfsmid;
  int nc_evapmid, nc_ahflmid, nc_wlm1mid, nc_tsm1mid, nc_tdm1mid, nc_wsm1mid,
    nc_wdm1mid, nc_snm1mid;
  int nc_emtermid, nc_trsolmid, nc_runoffmid, nc_srad0umid, nc_sradsumid,
    nc_tradsumid;
  int nc_albedomid, nc_tsurfmid, nc_tsnmid, nc_td3mid, nc_td4mid, nc_td5mid,
    nc_tsnm1mid;
  int nc_td3m1mid, nc_td4m1mid, nc_td5m1mid, nc_tdclm1mid, nc_seaicemid,
    nc_sicedmid;
  int nc_ustar3mid, nc_teffmid, nc_glacmid, nc_aclcmid, nc_aclcacmid,
    nc_snmelmid, nc_runtocmid;
  int nc_tslinmid, nc_dsnacmid, nc_t2maxmid, nc_t2minmid, nc_tsmaxmid,
    nc_tsminmid, nc_wimaxmid;
  int nc_topmaxmid, nc_aclcvmid, nc_qvimid, nc_alwcvimid, nc_runlndmid,
    nc_rgcgnmid, nc_sodifmid;
  int nc_srafsmid, nc_trafsmid, nc_sraf0mid, nc_traf0mid, nc_emtefmid,
    nc_trsofmid, nc_tclfsmid;
  int nc_sclfsmid, nc_tclf0mid, nc_sclf0mid, nc_tkemid, nc_tkem1mid,
    nc_auxil1mid, nc_auxil2mid;


  /* setup of global attributes for netCDF file (creation information) */

  restart->nc_access_mode = NC_WRITE;
  IO_create (ofile, restart->nc_access_mode, &(restart->nc_file_id));

  IO_def_dim (restart->nc_file_id, "ngl", d->ngl, &nglid);
  IO_def_dim (restart->nc_file_id, "nhgl", d->nhgl, &nhglid);
  IO_def_dim (restart->nc_file_id, "nlon", d->nlon, &nlonid);
  IO_def_dim (restart->nc_file_id, "nlp2", d->nlp2, &nlp2id);
  IO_def_dim (restart->nc_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (restart->nc_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (restart->nc_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (restart->nc_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (restart->nc_file_id, "n2", 2, &n2id);
  IO_def_dim (restart->nc_file_id, "nmp1", d->nmp1, &nmp1id);
  IO_def_dim (restart->nc_file_id, "nhtrac", d->nhtrac, &nhtracid);
  IO_def_dim (restart->nc_file_id, "nswitches", d->nrd, &nswitchesid);

  ncheader (argv, restart, info, d, g, run, t);

  /*
   * Coordinates, but not following the netCDF convention,
   * do it another way ?
   *
   */

  restart->nc_dims[0] = nglid;
  IO_def_var (restart->nc_file_id, "lat", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_glid);
  IO_put_att_text (restart->nc_file_id, nc_glid, "units", 9, "degrees_N");
  IO_put_att_text (restart->nc_file_id, nc_glid, "long_name", 17,
		   "Gaussian latitude");

  restart->nc_dims[0] = nlonid;
  IO_def_var (restart->nc_file_id, "lon", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_lonid);
  IO_put_att_text (restart->nc_file_id, nc_lonid, "units", 9, "degrees_E");
  IO_put_att_text (restart->nc_file_id, nc_lonid, "long_name", 9,
		   "longitude");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_a", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  restart->nc_dims[0] = nvclevid;
  IO_def_var (restart->nc_file_id, "vct_b", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (restart->nc_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  restart->nc_dims[0] = nswitchesid;
  IO_def_var (restart->nc_file_id, "switches", NC_DOUBLE, 1, restart->nc_dims,
	      &nc_switchesid);
  IO_put_att_text (restart->nc_file_id, nc_switchesid, "long_name", 19,
		   "diagnostic switches");


  /* g3-buffer */

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "GEOSPM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_geospmid);
  IO_put_att_text (restart->nc_file_id, nc_geospmid, "long_name", 32,
		   "surface geopotential (orography)");
  IO_put_att_text (restart->nc_file_id, nc_geospmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "TSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsmid);
  IO_put_att_text (restart->nc_file_id, nc_tsmid, "long_name", 19,
		   "surface temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsmid, "units", 1, "K");

  IO_def_var (restart->nc_file_id, "WSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wsmid);
  IO_put_att_text (restart->nc_file_id, nc_wsmid, "long_name", 12,
		   "soil wetness");
  IO_put_att_text (restart->nc_file_id, nc_wsmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "WLM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wlmid);
  IO_put_att_text (restart->nc_file_id, nc_wlmid, "long_name", 22,
		   "skin reservoir content");
  IO_put_att_text (restart->nc_file_id, nc_wlmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "SNM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_snmid);
  IO_put_att_text (restart->nc_file_id, nc_snmid, "long_name", 10,
		   "snow depth");
  IO_put_att_text (restart->nc_file_id, nc_snmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "SLMM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_slmmid);
  IO_put_att_text (restart->nc_file_id, nc_slmmid, "long_name", 13,
		   "land sea mask");
  IO_put_att_text (restart->nc_file_id, nc_slmmid, "units", 15,
		   "0: sea, 1: land");

  IO_def_var (restart->nc_file_id, "AZ0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_az0mid);
  IO_put_att_text (restart->nc_file_id, nc_az0mid, "long_name", 24,
		   "surface roughness length");
  IO_put_att_text (restart->nc_file_id, nc_az0mid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "ALBM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_albmid);
  IO_put_att_text (restart->nc_file_id, nc_albmid, "long_name", 25,
		   "surface background albedo");
  IO_put_att_text (restart->nc_file_id, nc_albmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "VARPM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_varpmid);
  IO_put_att_text (restart->nc_file_id, nc_varpmid, "long_name", 41,
		   "directional orographic variance (packed)");
  IO_put_att_text (restart->nc_file_id, nc_varpmid, "units", 3, "m^2");

  IO_def_var (restart->nc_file_id, "VARORM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_varormid);
  IO_put_att_text (restart->nc_file_id, nc_varormid, "long_name", 19,
		   "orographic variance");
  IO_put_att_text (restart->nc_file_id, nc_varormid, "units", 3, "m^2");

  IO_def_var (restart->nc_file_id, "FORESTM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_forestmid);
  IO_put_att_text (restart->nc_file_id, nc_forestmid, "long_name", 15,
		   "vegetation type");
  IO_put_att_text (restart->nc_file_id, nc_forestmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "VGRATM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vgratmid);
  IO_put_att_text (restart->nc_file_id, nc_vgratmid, "long_name", 16,
		   "vegetation ratio");
  IO_put_att_text (restart->nc_file_id, nc_vgratmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "VLTM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vltmid);
  IO_put_att_text (restart->nc_file_id, nc_vltmid, "long_name", 15,
		   "leaf area index");
  IO_put_att_text (restart->nc_file_id, nc_vltmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "WSMXM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wsmxmid);
  IO_put_att_text (restart->nc_file_id, nc_wsmxmid, "long_name", 22,
		   "field capacity of soil");
  IO_put_att_text (restart->nc_file_id, nc_wsmxmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "FAOM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_faomid);
  IO_put_att_text (restart->nc_file_id, nc_faomid, "long_name", 30,
		   "FAO data set (soil data flags)");
  IO_put_att_text (restart->nc_file_id, nc_faomid, "units", 5, "0...5");

  IO_def_var (restart->nc_file_id, "TDM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tdmid);
  IO_put_att_text (restart->nc_file_id, nc_tdmid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_tdmid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TDCLM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tdclmid);
  IO_put_att_text (restart->nc_file_id, nc_tdclmid, "long_name", 16,
		   "soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_tdclmid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "APSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_apsmid);
  IO_put_att_text (restart->nc_file_id, nc_apsmid, "long_name", 16,
		   "Surface pressure ");
  IO_put_att_text (restart->nc_file_id, nc_apsmid, "units", 3, "hPa");

  IO_def_var (restart->nc_file_id, "APRLM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_aprlmid);
  IO_put_att_text (restart->nc_file_id, nc_aprlmid, "long_name", 25,
		   "large scale precipitation");
  IO_put_att_text (restart->nc_file_id, nc_aprlmid, "units", 4, "mm/d");

  IO_def_var (restart->nc_file_id, "APRCM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_aprcmid);
  IO_put_att_text (restart->nc_file_id, nc_aprcmid, "long_name", 25,
		   "convective  precipitation");
  IO_put_att_text (restart->nc_file_id, nc_aprcmid, "units", 4, "mm/d");

  IO_def_var (restart->nc_file_id, "APRSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_aprsmid);
  IO_put_att_text (restart->nc_file_id, nc_aprsmid, "long_name", 9,
		   "snow fall");
  IO_put_att_text (restart->nc_file_id, nc_aprsmid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "USTRGWM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_ustrgwmid);
  IO_put_att_text (restart->nc_file_id, nc_ustrgwmid, "long_name", 22,
		   "u-gravity wave stress");
  IO_put_att_text (restart->nc_file_id, nc_ustrgwmid, "units", 9,
		   "mN/m`a2`n");

  IO_def_var (restart->nc_file_id, "VSTRGWM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vstrgwmid);
  IO_put_att_text (restart->nc_file_id, nc_vstrgwmid, "long_name", 22,
		   "v-gravity wave stress");
  IO_put_att_text (restart->nc_file_id, nc_vstrgwmid, "units", 9,
		   "mN/m`a2`n");

  IO_def_var (restart->nc_file_id, "VDISGWM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vdisgwmid);
  IO_put_att_text (restart->nc_file_id, nc_vdisgwmid, "long_name", 25,
		   "gravity wave dissipation");
  IO_put_att_text (restart->nc_file_id, nc_vdisgwmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "ACLCOVM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_aclcovmid);
  IO_put_att_text (restart->nc_file_id, nc_aclcovmid, "long_name", 17,
		   "total cloud cover");
  IO_put_att_text (restart->nc_file_id, nc_aclcovmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "TEMP2M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_temp2mid);
  IO_put_att_text (restart->nc_file_id, nc_temp2mid, "long_name", 14,
		   "2m temperature");
  IO_put_att_text (restart->nc_file_id, nc_temp2mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "DEW2M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_dew2mid);
  IO_put_att_text (restart->nc_file_id, nc_dew2mid, "long_name", 24,
		   "2m dew point temperature");
  IO_put_att_text (restart->nc_file_id, nc_dew2mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "WIND10M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wind10mid);
  IO_put_att_text (restart->nc_file_id, nc_wind10mid, "long_name", 13,
		   "10m windspeed");
  IO_put_att_text (restart->nc_file_id, nc_wind10mid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "U10M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_u10mid);
  IO_put_att_text (restart->nc_file_id, nc_u10mid, "long_name", 16,
		   "10m u-velocity");
  IO_put_att_text (restart->nc_file_id, nc_u10mid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "V10M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_v10mid);
  IO_put_att_text (restart->nc_file_id, nc_v10mid, "long_name", 16,
		   "10m v-velocity");
  IO_put_att_text (restart->nc_file_id, nc_v10mid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "SRADSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sradsmid);
  IO_put_att_text (restart->nc_file_id, nc_sradsmid, "long_name", 27,
		   "net surface solar radiation");
  IO_put_att_text (restart->nc_file_id, nc_sradsmid, "units", 8, "W/m`a2`n");
  IO_def_var (restart->nc_file_id, "TRADSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tradsmid);
  IO_put_att_text (restart->nc_file_id, nc_tradsmid, "long_name", 29,
		   "net surface thermal radiation");
  IO_put_att_text (restart->nc_file_id, nc_tradsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "SRAD0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_srad0mid);
  IO_put_att_text (restart->nc_file_id, nc_srad0mid, "long_name", 23,
		   "net top solar radiation");
  IO_put_att_text (restart->nc_file_id, nc_srad0mid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "TRAD0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_trad0mid);
  IO_put_att_text (restart->nc_file_id, nc_trad0mid, "long_name", 27,
		   "top thermal radiation (OLR)");
  IO_put_att_text (restart->nc_file_id, nc_trad0mid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "VDISM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vdismid);
  IO_put_att_text (restart->nc_file_id, nc_vdismid, "long_name", 26,
		   "boundary layer dissipation");
  IO_put_att_text (restart->nc_file_id, nc_vdismid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "USTRM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_ustrmid);
  IO_put_att_text (restart->nc_file_id, nc_ustrmid, "long_name", 16,
		   "surface u-stress");
  IO_put_att_text (restart->nc_file_id, nc_ustrmid, "units", 9, "mN/m`a2`n");

  IO_def_var (restart->nc_file_id, "VSTRM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_vstrmid);
  IO_put_att_text (restart->nc_file_id, nc_vstrmid, "long_name", 16,
		   "surface v-stress");
  IO_put_att_text (restart->nc_file_id, nc_vstrmid, "units", 9, "mN/m`a2`n");

  IO_def_var (restart->nc_file_id, "AHFSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_ahfsmid);
  IO_put_att_text (restart->nc_file_id, nc_ahfsmid, "long_name", 26,
		   "surface sensible heat flux");
  IO_put_att_text (restart->nc_file_id, nc_ahfsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "EVAPM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_evapmid);
  IO_put_att_text (restart->nc_file_id, nc_evapmid, "long_name", 19,
		   "surface evaporation");
  IO_put_att_text (restart->nc_file_id, nc_evapmid, "units", 4, "mm/d");

  IO_def_var (restart->nc_file_id, "AHFLM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_ahflmid);
  IO_put_att_text (restart->nc_file_id, nc_ahflmid, "long_name", 24,
		   "surface latent heat flux");
  IO_put_att_text (restart->nc_file_id, nc_ahflmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "WLM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wlm1mid);
  IO_put_att_text (restart->nc_file_id, nc_wlm1mid, "long_name", 32,
		   "skin reservoir content of plants");
  IO_put_att_text (restart->nc_file_id, nc_wlm1mid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "TSM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsm1mid);
  IO_put_att_text (restart->nc_file_id, nc_tsm1mid, "long_name", 19,
		   "surface temperature ");
  IO_put_att_text (restart->nc_file_id, nc_tsm1mid, "units", 1, "K");

  IO_def_var (restart->nc_file_id, "TDM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tdm1mid);
  IO_put_att_text (restart->nc_file_id, nc_tdm1mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_tdm1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "WSM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wsm1mid);
  IO_put_att_text (restart->nc_file_id, nc_wsm1mid, "long_name", 12,
		   "soil wetness");
  IO_put_att_text (restart->nc_file_id, nc_wsm1mid, "units", 2, "cm");

  IO_def_var (restart->nc_file_id, "WDM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wdm1mid);
  IO_put_att_text (restart->nc_file_id, nc_wdm1mid, "long_name", 1, "wdm1m");
  IO_put_att_text (restart->nc_file_id, nc_wdm1mid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "SNM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_snm1mid);
  IO_put_att_text (restart->nc_file_id, nc_snm1mid, "long_name", 10,
		   "snow depth");
  IO_put_att_text (restart->nc_file_id, nc_snm1mid, "units", 1, "m");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevp1id;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "EMTERM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_emtermid);
  IO_put_att_text (restart->nc_file_id, nc_emtermid, "long_name", 1,
		   "emterm");
  IO_put_att_text (restart->nc_file_id, nc_emtermid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "TRSOLM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_trsolmid);
  IO_put_att_text (restart->nc_file_id, nc_trsolmid, "long_name", 1,
		   "trsolm");
  IO_put_att_text (restart->nc_file_id, nc_trsolmid, "units", 1, "");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "RUNOFFM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_runoffmid);
  IO_put_att_text (restart->nc_file_id, nc_runoffmid, "long_name", 14,
		   "surface runoff");
  IO_put_att_text (restart->nc_file_id, nc_runoffmid, "units", 4, "mm/d");

  IO_def_var (restart->nc_file_id, "SRAD0UM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_srad0umid);
  IO_put_att_text (restart->nc_file_id, nc_srad0umid, "long_name", 26,
		   "top solar radiation upward");
  IO_put_att_text (restart->nc_file_id, nc_srad0umid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "SRADSUM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sradsumid);
  IO_put_att_text (restart->nc_file_id, nc_sradsumid, "long_name", 30,
		   "surface solar radiation upward");
  IO_put_att_text (restart->nc_file_id, nc_sradsumid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "TRADSUM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tradsumid);
  IO_put_att_text (restart->nc_file_id, nc_tradsumid, "long_name", 32,
		   "surface thermal radiation upward");
  IO_put_att_text (restart->nc_file_id, nc_tradsumid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "ALBEDOM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_albedomid);
  IO_put_att_text (restart->nc_file_id, nc_albedomid, "long_name", 14,
		   "surface albedo");
  IO_put_att_text (restart->nc_file_id, nc_albedomid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "TSURFM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsurfmid);
  IO_put_att_text (restart->nc_file_id, nc_tsurfmid, "long_name", 20,
		   "surface temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsurfmid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TSNM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsnmid);
  IO_put_att_text (restart->nc_file_id, nc_tsnmid, "long_name", 17,
		   "snow temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsnmid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD3M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td3mid);
  IO_put_att_text (restart->nc_file_id, nc_td3mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td3mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD4M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td4mid);
  IO_put_att_text (restart->nc_file_id, nc_td4mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td4mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD5M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td5mid);
  IO_put_att_text (restart->nc_file_id, nc_td5mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td5mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TSNM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsnm1mid);
  IO_put_att_text (restart->nc_file_id, nc_tsnm1mid, "long_name", 16,
		   "snow temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsnm1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD3M1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td3m1mid);
  IO_put_att_text (restart->nc_file_id, nc_td3m1mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td3m1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD4M1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td4m1mid);
  IO_put_att_text (restart->nc_file_id, nc_td4m1mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td4m1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TD5M1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_td5m1mid);
  IO_put_att_text (restart->nc_file_id, nc_td5m1mid, "long_name", 22,
		   "deep soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_td5m1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "TDCLM1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tdclm1mid);
  IO_put_att_text (restart->nc_file_id, nc_tdclm1mid, "long_name", 17,
		   "soil temperature");
  IO_put_att_text (restart->nc_file_id, nc_tdclm1mid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "SEAICEM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_seaicemid);
  IO_put_att_text (restart->nc_file_id, nc_seaicemid, "long_name", 13,
		   "sea ice cover");
  IO_put_att_text (restart->nc_file_id, nc_seaicemid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "SICEDM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sicedmid);
  IO_put_att_text (restart->nc_file_id, nc_sicedmid, "long_name", 13,
		   "sea ice dept");
  IO_put_att_text (restart->nc_file_id, nc_sicedmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "USTAR3M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_ustar3mid);
  IO_put_att_text (restart->nc_file_id, nc_ustar3mid, "long_name", 10,
		   "ustar`a3`n");
  IO_put_att_text (restart->nc_file_id, nc_ustar3mid, "units", 13,
		   "m`a3`n/s`a3`n");

  IO_def_var (restart->nc_file_id, "TEFFM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_teffmid);
  IO_put_att_text (restart->nc_file_id, nc_teffmid, "long_name", 36,
		   "(effective) sea-ice skin temperature");
  IO_put_att_text (restart->nc_file_id, nc_teffmid, "units", 6, "`ao`nC");

  IO_def_var (restart->nc_file_id, "GLACM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_glacmid);
  IO_put_att_text (restart->nc_file_id, nc_glacmid, "long_name", 12,
		   "glacier mask");
  IO_put_att_text (restart->nc_file_id, nc_glacmid, "units", 13,
		   "0: no, 1: yes");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "ACLCM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_aclcmid);
  IO_put_att_text (restart->nc_file_id, nc_aclcmid, "long_name", 11,
		   "cloud cover");
  IO_put_att_text (restart->nc_file_id, nc_aclcmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "ACLCACM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_aclcacmid);
  IO_put_att_text (restart->nc_file_id, nc_aclcacmid, "long_name", 20,
		   "liquid water content");
  IO_put_att_text (restart->nc_file_id, nc_aclcacmid, "units", 5, "kg/kg");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "SNMELM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_snmelmid);
  IO_put_att_text (restart->nc_file_id, nc_snmelmid, "long_name", 9,
		   "snow melt");
  IO_put_att_text (restart->nc_file_id, nc_snmelmid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "RUNTOCM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_runtocmid);
  IO_put_att_text (restart->nc_file_id, nc_runtocmid, "long_name", 25,
		   "surface runoff into ocean");
  IO_put_att_text (restart->nc_file_id, nc_runtocmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "TSLINM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tslinmid);
  IO_put_att_text (restart->nc_file_id, nc_tslinmid, "long_name", 34,
		   "land: residual surface heat budget");
  IO_put_att_text (restart->nc_file_id, nc_tslinmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "DSNACM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_dsnacmid);
  IO_put_att_text (restart->nc_file_id, nc_dsnacmid, "long_name", 17,
		   "snow depth change");
  IO_put_att_text (restart->nc_file_id, nc_dsnacmid, "units", 4, "mm/d");

  IO_def_var (restart->nc_file_id, "T2MAXM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_t2maxmid);
  IO_put_att_text (restart->nc_file_id, nc_t2maxmid, "long_name", 22,
		   "maximum 2m-temperature");
  IO_put_att_text (restart->nc_file_id, nc_t2maxmid, "units", 6, "ao`nC");

  IO_def_var (restart->nc_file_id, "T2MINM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_t2minmid);
  IO_put_att_text (restart->nc_file_id, nc_t2minmid, "long_name", 22,
		   "minimum 2m-temperature");
  IO_put_att_text (restart->nc_file_id, nc_t2minmid, "units", 6, "ao`nC");

  IO_def_var (restart->nc_file_id, "TSMAXM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsmaxmid);
  IO_put_att_text (restart->nc_file_id, nc_tsmaxmid, "long_name", 27,
		   "maximum surface temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsmaxmid, "units", 6, "ao`nC");

  IO_def_var (restart->nc_file_id, "TSMINM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tsminmid);
  IO_put_att_text (restart->nc_file_id, nc_tsminmid, "long_name", 27,
		   "minimum surface temperature");
  IO_put_att_text (restart->nc_file_id, nc_tsminmid, "units", 6, "ao`nC");

  IO_def_var (restart->nc_file_id, "WIMAXM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_wimaxmid);
  IO_put_att_text (restart->nc_file_id, nc_wimaxmid, "long_name", 22,
		   "maximum 10m-wind speed");
  IO_put_att_text (restart->nc_file_id, nc_wimaxmid, "units", 3, "m/s");

  IO_def_var (restart->nc_file_id, "TOPMAXM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_topmaxmid);
  IO_put_att_text (restart->nc_file_id, nc_topmaxmid, "long_name", 29,
		   "maximum convective cloud tops");
  IO_put_att_text (restart->nc_file_id, nc_topmaxmid, "units", 1, "m");

  IO_def_var (restart->nc_file_id, "ACLCVM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_aclcvmid);
  IO_put_att_text (restart->nc_file_id, nc_aclcvmid, "long_name", 17,
		   "total cloud cover");
  IO_put_att_text (restart->nc_file_id, nc_aclcvmid, "units", 1, "%");

  IO_def_var (restart->nc_file_id, "QVIM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_qvimid);
  IO_put_att_text (restart->nc_file_id, nc_qvimid, "long_name", 40,
		   "vertically integrated specific humidity");
  IO_put_att_text (restart->nc_file_id, nc_qvimid, "units", 9, "kg/m`a2`n");

  IO_def_var (restart->nc_file_id, "ALWCVIM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_alwcvimid);
  IO_put_att_text (restart->nc_file_id, nc_alwcvimid, "long_name", 41,
		   "vertically integrated liquid water cont.");
  IO_put_att_text (restart->nc_file_id, nc_alwcvimid, "units", 10,
		   "kg/m`a2`n");

  IO_def_var (restart->nc_file_id, "RUNLNDM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_runlndmid);
  IO_put_att_text (restart->nc_file_id, nc_runlndmid, "long_name", 37,
		   "surface runoff not running into ocean");
  IO_put_att_text (restart->nc_file_id, nc_runlndmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "RGCGNM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_rgcgnmid);
  IO_put_att_text (restart->nc_file_id, nc_rgcgnmid, "long_name", 21,
		   "heat capacity of soil");
  IO_put_att_text (restart->nc_file_id, nc_rgcgnmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "SODIFM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sodifmid);
  IO_put_att_text (restart->nc_file_id, nc_sodifmid, "long_name", 16,
		   "soil diffusivity");
  IO_put_att_text (restart->nc_file_id, nc_sodifmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "SRAFSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_srafsmid);
  IO_put_att_text (restart->nc_file_id, nc_srafsmid, "long_name", 40,
		   "net surf. solar radiation  (clear sky)");
  IO_put_att_text (restart->nc_file_id, nc_srafsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "TRAFSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_trafsmid);
  IO_put_att_text (restart->nc_file_id, nc_trafsmid, "long_name", 40,
		   "net surf. thermal radiation (clear sky)");
  IO_put_att_text (restart->nc_file_id, nc_trafsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "SRAF0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sraf0mid);
  IO_put_att_text (restart->nc_file_id, nc_sraf0mid, "long_name", 40,
		   "net top solar radiation    (clear sky)");
  IO_put_att_text (restart->nc_file_id, nc_sraf0mid, "units", 9, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "TRAF0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_traf0mid);
  IO_put_att_text (restart->nc_file_id, nc_traf0mid, "long_name", 40,
		   "net top thermal radiation  (clear sky)");
  IO_put_att_text (restart->nc_file_id, nc_traf0mid, "units", 8, "W/m`a2`n");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = n2id;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "EMTEFM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_emtefmid);
  IO_put_att_text (restart->nc_file_id, nc_emtefmid, "long_name", 6,
		   "emtefm");
  IO_put_att_text (restart->nc_file_id, nc_emtefmid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "TRSOFM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_trsofmid);
  IO_put_att_text (restart->nc_file_id, nc_trsofmid, "long_name", 6,
		   "trsofm");
  IO_put_att_text (restart->nc_file_id, nc_trsofmid, "units", 1, "");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "TCLFSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tclfsmid);
  IO_put_att_text (restart->nc_file_id, nc_tclfsmid, "long_name", 32,
		   "surface thermal cloud forcing");
  IO_put_att_text (restart->nc_file_id, nc_tclfsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "SCLFSM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sclfsmid);
  IO_put_att_text (restart->nc_file_id, nc_sclfsmid, "long_name", 28,
		   "surface solar cloud forcing");
  IO_put_att_text (restart->nc_file_id, nc_sclfsmid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "TCLF0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_tclf0mid);
  IO_put_att_text (restart->nc_file_id, nc_tclf0mid, "long_name", 26,
		   "top thermal cloud forcing");
  IO_put_att_text (restart->nc_file_id, nc_tclf0mid, "units", 8, "W/m`a2`n");

  IO_def_var (restart->nc_file_id, "SCLF0M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_sclf0mid);
  IO_put_att_text (restart->nc_file_id, nc_sclf0mid, "long_name", 24,
		   "top solar cloud forcing");
  IO_put_att_text (restart->nc_file_id, nc_sclf0mid, "units", 8, "W/m`a2`n");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlevid;
  restart->nc_dims[2] = nlp2id;

  IO_def_var (restart->nc_file_id, "TKEM", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_tkemid);
  IO_put_att_text (restart->nc_file_id, nc_tkemid, "long_name", 24,
		   "turbulent kinetic energy");
  IO_put_att_text (restart->nc_file_id, nc_tkemid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "TKEM1M", NC_DOUBLE, 3, restart->nc_dims,
	      &nc_tkem1mid);
  IO_put_att_text (restart->nc_file_id, nc_tkem1mid, "long_name", 30,
		   "turbulent kinetic energy (t-1)");
  IO_put_att_text (restart->nc_file_id, nc_tkem1mid, "units", 1, "");

  restart->nc_dims[0] = nglid;
  restart->nc_dims[1] = nlp2id;

  IO_def_var (restart->nc_file_id, "AUXIL1M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_auxil1mid);
  IO_put_att_text (restart->nc_file_id, nc_auxil1mid, "long_name", 1,
		   "auxil1m");
  IO_put_att_text (restart->nc_file_id, nc_auxil1mid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "AUXIL2M", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_auxil2mid);
  IO_put_att_text (restart->nc_file_id, nc_auxil2mid, "long_name", 1,
		   "auxil2m");
  IO_put_att_text (restart->nc_file_id, nc_auxil2mid, "units", 1, "");

  IO_def_var (restart->nc_file_id, "DRAINM", NC_DOUBLE, 2, restart->nc_dims,
	      &nc_dummyid);

  IO_enddef (restart->nc_file_id);

  IO_put_var_double (restart->nc_file_id, nc_glid, g->lat);
  IO_put_var_double (restart->nc_file_id, nc_lonid, g->lon);
  IO_put_var_double (restart->nc_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (restart->nc_file_id, nc_vctbid, g->vct_b);
  IO_put_var_double (restart->nc_file_id, nc_switchesid, d->switches);

  IO_put_var_double (restart->nc_file_id, nc_geospmid, geospm);
  IO_put_var_double (restart->nc_file_id, nc_tsmid, tsm);
  IO_put_var_double (restart->nc_file_id, nc_wsmid, wsm);
  IO_put_var_double (restart->nc_file_id, nc_wlmid, wlm);
  IO_put_var_double (restart->nc_file_id, nc_snmid, snm);
  IO_put_var_double (restart->nc_file_id, nc_slmmid, slmm);
  IO_put_var_double (restart->nc_file_id, nc_az0mid, az0m);
  IO_put_var_double (restart->nc_file_id, nc_albmid, albm);
  IO_put_var_double (restart->nc_file_id, nc_varpmid, varpm);
  IO_put_var_double (restart->nc_file_id, nc_varormid, varorm);
  IO_put_var_double (restart->nc_file_id, nc_forestmid, forestm);
  IO_put_var_double (restart->nc_file_id, nc_vgratmid, vgratm);
  IO_put_var_double (restart->nc_file_id, nc_vltmid, vltm);
  IO_put_var_double (restart->nc_file_id, nc_wsmxmid, wsmxm);
  IO_put_var_double (restart->nc_file_id, nc_faomid, faom);
  IO_put_var_double (restart->nc_file_id, nc_tdmid, tdm);
  IO_put_var_double (restart->nc_file_id, nc_tdclmid, tdclm);
  IO_put_var_double (restart->nc_file_id, nc_apsmid, apsm);
  IO_put_var_double (restart->nc_file_id, nc_aprlmid, aprlm);
  IO_put_var_double (restart->nc_file_id, nc_aprcmid, aprcm);
  IO_put_var_double (restart->nc_file_id, nc_aprsmid, aprsm);
  IO_put_var_double (restart->nc_file_id, nc_ustrgwmid, ustrgwm);
  IO_put_var_double (restart->nc_file_id, nc_vstrgwmid, vstrgwm);
  IO_put_var_double (restart->nc_file_id, nc_vdisgwmid, vdisgwm);
  IO_put_var_double (restart->nc_file_id, nc_aclcovmid, aclcovm);
  IO_put_var_double (restart->nc_file_id, nc_temp2mid, temp2m);
  IO_put_var_double (restart->nc_file_id, nc_dew2mid, dew2m);
  IO_put_var_double (restart->nc_file_id, nc_wind10mid, wind10m);
  IO_put_var_double (restart->nc_file_id, nc_u10mid, u10m);
  IO_put_var_double (restart->nc_file_id, nc_v10mid, v10m);
  IO_put_var_double (restart->nc_file_id, nc_sradsmid, sradsm);
  IO_put_var_double (restart->nc_file_id, nc_tradsmid, tradsm);
  IO_put_var_double (restart->nc_file_id, nc_srad0mid, srad0m);
  IO_put_var_double (restart->nc_file_id, nc_trad0mid, trad0m);
  IO_put_var_double (restart->nc_file_id, nc_vdismid, vdism);
  IO_put_var_double (restart->nc_file_id, nc_ustrmid, ustrm);
  IO_put_var_double (restart->nc_file_id, nc_vstrmid, vstrm);
  IO_put_var_double (restart->nc_file_id, nc_ahfsmid, ahfsm);
  IO_put_var_double (restart->nc_file_id, nc_evapmid, evapm);
  IO_put_var_double (restart->nc_file_id, nc_ahflmid, ahflm);
  IO_put_var_double (restart->nc_file_id, nc_wlm1mid, wlm1m);
  IO_put_var_double (restart->nc_file_id, nc_tsm1mid, tsm1m);
  IO_put_var_double (restart->nc_file_id, nc_tdm1mid, tdm1m);
  IO_put_var_double (restart->nc_file_id, nc_wsm1mid, wsm1m);
  IO_put_var_double (restart->nc_file_id, nc_wdm1mid, wdm1m);
  IO_put_var_double (restart->nc_file_id, nc_snm1mid, snm1m);
  IO_put_var_double (restart->nc_file_id, nc_emtermid, emterm);
  IO_put_var_double (restart->nc_file_id, nc_trsolmid, trsolm);
  IO_put_var_double (restart->nc_file_id, nc_runoffmid, runoffm);
  IO_put_var_double (restart->nc_file_id, nc_srad0umid, srad0um);
  IO_put_var_double (restart->nc_file_id, nc_sradsumid, sradsum);
  IO_put_var_double (restart->nc_file_id, nc_tradsumid, tradsum);
  IO_put_var_double (restart->nc_file_id, nc_albedomid, albedom);
  IO_put_var_double (restart->nc_file_id, nc_tsurfmid, tsurfm);
  IO_put_var_double (restart->nc_file_id, nc_tsnmid, tsnm);
  IO_put_var_double (restart->nc_file_id, nc_td3mid, td3m);
  IO_put_var_double (restart->nc_file_id, nc_td4mid, td4m);
  IO_put_var_double (restart->nc_file_id, nc_td5mid, td5m);
  IO_put_var_double (restart->nc_file_id, nc_tsnm1mid, tsnm1m);
  IO_put_var_double (restart->nc_file_id, nc_td3m1mid, td3m1m);
  IO_put_var_double (restart->nc_file_id, nc_td4m1mid, td4m1m);
  IO_put_var_double (restart->nc_file_id, nc_td5m1mid, td5m1m);
  IO_put_var_double (restart->nc_file_id, nc_tdclm1mid, tdclm1m);
  IO_put_var_double (restart->nc_file_id, nc_seaicemid, seaicem);
  IO_put_var_double (restart->nc_file_id, nc_sicedmid, sicedm);
  IO_put_var_double (restart->nc_file_id, nc_ustar3mid, ustar3m);
  IO_put_var_double (restart->nc_file_id, nc_teffmid, teffm);
  IO_put_var_double (restart->nc_file_id, nc_glacmid, glacm);
  IO_put_var_double (restart->nc_file_id, nc_aclcmid, aclcm);
  IO_put_var_double (restart->nc_file_id, nc_aclcacmid, aclcacm);
  IO_put_var_double (restart->nc_file_id, nc_snmelmid, snmelm);
  IO_put_var_double (restart->nc_file_id, nc_runtocmid, runtocm);
  IO_put_var_double (restart->nc_file_id, nc_tslinmid, tslinm);
  IO_put_var_double (restart->nc_file_id, nc_dsnacmid, dsnacm);
  IO_put_var_double (restart->nc_file_id, nc_t2maxmid, t2maxm);
  IO_put_var_double (restart->nc_file_id, nc_t2minmid, t2minm);
  IO_put_var_double (restart->nc_file_id, nc_tsmaxmid, tsmaxm);
  IO_put_var_double (restart->nc_file_id, nc_tsminmid, tsminm);
  IO_put_var_double (restart->nc_file_id, nc_wimaxmid, wimaxm);
  IO_put_var_double (restart->nc_file_id, nc_topmaxmid, topmaxm);
  IO_put_var_double (restart->nc_file_id, nc_aclcvmid, aclcvm);
  IO_put_var_double (restart->nc_file_id, nc_qvimid, qvim);
  IO_put_var_double (restart->nc_file_id, nc_alwcvimid, alwcvim);
  IO_put_var_double (restart->nc_file_id, nc_runlndmid, runlndm);
  IO_put_var_double (restart->nc_file_id, nc_rgcgnmid, rgcgnm);
  IO_put_var_double (restart->nc_file_id, nc_sodifmid, sodifm);
  IO_put_var_double (restart->nc_file_id, nc_srafsmid, srafsm);
  IO_put_var_double (restart->nc_file_id, nc_trafsmid, trafsm);
  IO_put_var_double (restart->nc_file_id, nc_sraf0mid, sraf0m);
  IO_put_var_double (restart->nc_file_id, nc_traf0mid, traf0m);
  IO_put_var_double (restart->nc_file_id, nc_emtefmid, emtefm);
  IO_put_var_double (restart->nc_file_id, nc_trsofmid, trsofm);
  IO_put_var_double (restart->nc_file_id, nc_tclfsmid, tclfsm);
  IO_put_var_double (restart->nc_file_id, nc_sclfsmid, sclfsm);
  IO_put_var_double (restart->nc_file_id, nc_tclf0mid, tclf0m);
  IO_put_var_double (restart->nc_file_id, nc_sclf0mid, sclf0m);
  IO_put_var_double (restart->nc_file_id, nc_tkemid, tkem);
  IO_put_var_double (restart->nc_file_id, nc_tkem1mid, tkem1m);
  IO_put_var_double (restart->nc_file_id, nc_auxil1mid, auxil1m);
  IO_put_var_double (restart->nc_file_id, nc_auxil2mid, auxil2m);

  IO_close (restart->nc_file_id);

  return;
}
