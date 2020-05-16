#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/stat.h>

#include "prep.h"
#include "io_nc.h"


double *
read_sp (FILE * file, int ilen, HEADER * info, size_t * rb)
{
  static double *f;
  UINT8 *in, *out;

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));
  *rb = fread (in, sizeof (UINT8), ilen, file);
  out = (UINT8 *) f;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (UINT8 *));
    }
#else
  memcpy (f, in, ilen * sizeof (UINT8 *));
#endif

  free (in);

  return (f);
}

double *
read_gl (FILE * file, int ilen, HEADER * info, DIMENSIONS * model,
	 size_t * rb)
{
  static double *f;
  int n, npos, ilenpp;
  size_t lrb;
  UINT8 *in, *out;

  *rb = 0;
  ilenpp = ilen / model->ngl;

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));
  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb = fread (in + npos * ilenpp, sizeof (UINT8), ilenpp, file);
      *rb += lrb;
    }
  out = (UINT8 *) f;

#ifndef CRAY
  if (!strncmp (info->def, "CRAY", 4))
    {
      cray2ieee (in, out, &ilen);
    }
  else
    {
      memcpy (f, in, ilen * sizeof (UINT8 *));
    }
#else
  memcpy (f, in, ilen * sizeof (UINT8 *));
#endif

  free (in);

  return (f);
}

double *
read_g3 (int no, long data_offset, FILE * file, int ilen, HEADER * info,
	 DIMENSIONS * model, size_t * rb)
{
  static double *f;
  int n, npos, ilenpp;
  long noffset;
  size_t lrb;
  UINT8 *in, *out;

  *rb = 0;
  ilenpp = ilen / model->ngl;

  /* now  go to begin of data section */

  rewind (file);
  fseek (file, data_offset, SEEK_CUR);

  /* select begin of field */

  noffset = (long) ((no - 1) * ilenpp * sizeof (UINT8));
  fseek (file, noffset, SEEK_CUR);

  noffset = (long) (14 * ilenpp * sizeof (UINT8));

  in = (UINT8 *) malloc (ilen * sizeof (UINT8));
  f = (double *) malloc (ilen * sizeof (double));
  for (n = 0; n < model->ngl; n++)
    {
      npos =
	((n % 2) * (model->ngl - (n - 1) / 2 - 1) + ((n + 1) % 2) * n / 2);
      lrb = fread (in + npos * ilenpp, sizeof (UINT8), ilenpp, file);
      if (n < model->ngl - 1)
	fseek (file, noffset, SEEK_CUR);
      *rb += lrb;
    }
  out = (UINT8 *) f;

  if (no == 9)
    {
      memcpy (f, in, ilen * sizeof (UINT8));
    }
  else
    {
#ifndef CRAY
      if (!strncmp (info->def, "CRAY", 4))
	{
	  cray2ieee (in, out, &ilen);
	}
      else
	{
	  memcpy (f, in, ilen * sizeof (UINT8));
	}
#else
      memcpy (f, in, ilen * sizeof (UINT8));
#endif
    }

  free (in);

  return (f);
}


/*=========================================================================== */

void
write_anc (netCDF_file * initial, HEADER * info, DIMENSIONS * d,
	   GRID * g, LABELS * run, TIMES * t, char *ofile,
	   double *svo, double *sd, double *stp, double *q, int fileformat)
{
  int nlatid, nlonid, nlevid, nlevp1id, nspid, nvclevid, n2id, nlp2id;
  int nc_glid, nc_lonid, nc_vctaid, nc_vctbid;
  int nc_svoid, nc_sdid, nc_stpid;
  int nc_qid, ilev, ilat, ioff, ioff2;
  int IO_file_id;
  time_t current_time;
  struct passwd *user;
  int dlen2, nlon, nlp2, nlat, nlev;
  double *dtmp;
  extern int Debug;
  extern char *Commandline;
  extern Table Table_specini[];
  /* extern int MaxCode_dim;  not used yet */
  extern Table Table_dim[];
  extern char institute[];

  /* setup of global attributes for netCDF file (creation information) */

  /* strcpy (ncheader.nc_conventions, "COARDS");  the spectral files are NOT COARDS conform ! */
  strcpy (initial->nc_model, "ECHAM");
  strcpy (initial->nc_creation_program, Commandline);
  strcpy (initial->nc_institution, institute);
  current_time = time (NULL);
  strcpy (initial->nc_creation_date, ctime (&current_time));
  initial->nc_creation_date[strlen (initial->nc_creation_date) - 1] = '\0';
  user = getpwuid (getuid ());
  strcpy (initial->nc_creation_user, user->pw_gecos);

  initial->nc_access_mode = NC_WRITE;
  IO_create (ofile, initial->nc_access_mode, &(initial->nc_file_id));

  IO_file_id = initial->nc_file_id;

  strcpy (initial->nc_binary_source, info->def);
  strcpy (initial->nc_file_type, "Initial file spectral");


  if (Debug)
    fprintf (stderr, " writing ECHAM netCDF Fileformat:  %d \n", fileformat);

  if (fileformat == 1)
    {
      IO_def_dim (IO_file_id, "ngl", d->ngl, &nlatid);
      IO_def_dim (IO_file_id, "nlon", d->nlon, &nlonid);
      IO_def_dim (IO_file_id, "nlp2", d->nlp2, &nlp2id);
    }
  else if (fileformat == 2)
    {
      IO_def_dim (IO_file_id, Table_dim[1].name, d->ngl, &nlatid);
      IO_def_dim (IO_file_id, Table_dim[0].name, d->nlon, &nlonid);
    }
  IO_def_dim (IO_file_id, "nlev", d->nlev, &nlevid);
  IO_def_dim (IO_file_id, "nlevp1", d->nlevp1, &nlevp1id);
  IO_def_dim (IO_file_id, "nsp", d->nsp, &nspid);
  IO_def_dim (IO_file_id, "nvclev", g->nvclev, &nvclevid);
  IO_def_dim (IO_file_id, "n2", 2, &n2id);

  /* IO_put_att_text (IO_file_id, NC_GLOBAL, "Conventions",
     strlen (ncheader.nc_conventions), initial.nc_conventions); */
  IO_put_att_text (IO_file_id, NC_GLOBAL, "model",
		   strlen (initial->nc_model), initial->nc_model);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "file_type",
		   strlen (initial->nc_file_type), initial->nc_file_type);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "source_type",
		   strlen (initial->nc_binary_source),
		   initial->nc_binary_source);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "title", 29,
		   "No information - initial file");
  IO_put_att_text (IO_file_id, NC_GLOBAL, "institution",
		   strlen (initial->nc_institution), initial->nc_institution);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "user",
		   strlen (initial->nc_creation_user),
		   initial->nc_creation_user);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "history",
		   strlen (initial->nc_creation_program),
		   initial->nc_creation_program);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "created",
		   strlen (initial->nc_creation_date),
		   initial->nc_creation_date);

  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_1",
		   strlen (run->label[0]), run->label[0]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_2",
		   strlen (run->label[1]), run->label[1]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_3",
		   strlen (run->label[2]), run->label[2]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_4",
		   strlen (run->label[3]), run->label[3]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_5",
		   strlen (run->label[4]), run->label[4]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_6",
		   strlen (run->label[5]), run->label[5]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_7",
		   strlen (run->label[6]), run->label[6]);
  IO_put_att_text (IO_file_id, NC_GLOBAL, "label_8",
		   strlen (run->label[7]), run->label[7]);

  IO_put_att_int (IO_file_id, NC_GLOBAL, "fdate", NC_INT, 1, &(t->fdate));
  IO_put_att_int (IO_file_id, NC_GLOBAL, "ftime", NC_INT, 1, &(t->ftime));
  IO_put_att_int (IO_file_id, NC_GLOBAL, "vdate", NC_INT, 1, &(t->vdate));
  IO_put_att_int (IO_file_id, NC_GLOBAL, "vtime", NC_INT, 1, &(t->vtime));

  IO_put_att_int (IO_file_id, NC_GLOBAL, "spherical_truncation_n",
		  NC_INT, 1, &(d->nn));
  IO_put_att_int (IO_file_id, NC_GLOBAL, "spherical_truncation_m",
		  NC_INT, 1, &(d->nm));
  IO_put_att_int (IO_file_id, NC_GLOBAL, "spherical_truncation_k",
		  NC_INT, 1, &(d->nk));

  /*
   * Coordinates, but not following the netCDF convention,
   * do it an other way ?
   *
   */

  initial->nc_dims[0] = nlatid;
  IO_def_var (IO_file_id, Table_dim[1].name, NC_DOUBLE, 1, initial->nc_dims,
	      &nc_glid);
  IO_put_att_text (IO_file_id, nc_glid, "units", strlen (Table_dim[1].unit),
		   Table_dim[1].unit);
  IO_put_att_text (IO_file_id, nc_glid, "long_name",
		   strlen (Table_dim[1].longname), Table_dim[1].longname);

  initial->nc_dims[0] = nlonid;
  IO_def_var (IO_file_id, Table_dim[0].name, NC_DOUBLE, 1, initial->nc_dims,
	      &nc_lonid);
  IO_put_att_text (IO_file_id, nc_lonid, "units", strlen (Table_dim[0].unit),
		   Table_dim[0].unit);
  IO_put_att_text (IO_file_id, nc_lonid, "long_name",
		   strlen (Table_dim[0].longname), Table_dim[0].longname);

  initial->nc_dims[0] = nvclevid;
  IO_def_var (IO_file_id, "vct_a", NC_DOUBLE, 1, initial->nc_dims,
	      &nc_vctaid);
  IO_put_att_text (IO_file_id, nc_vctaid, "units", 1, "");
  IO_put_att_text (IO_file_id, nc_vctaid, "long_name", 35,
		   "vertical-coordinate parameter set A");

  initial->nc_dims[0] = nvclevid;
  IO_def_var (IO_file_id, "vct_b", NC_DOUBLE, 1, initial->nc_dims,
	      &nc_vctbid);
  IO_put_att_text (IO_file_id, nc_vctbid, "units", 1, "");
  IO_put_att_text (IO_file_id, nc_vctbid, "long_name", 35,
		   "vertical-coordinate parameter set B");

  /* global spectral fields */

  initial->nc_dims[0] = nspid;
  initial->nc_dims[1] = n2id;
  initial->nc_dims[2] = nlevid;


  IO_def_var (IO_file_id, Table_specini[0].name, NC_DOUBLE, 3,
	      initial->nc_dims, &nc_svoid);
  IO_put_att_text (IO_file_id, nc_svoid, "long_name",
		   strlen (Table_specini[0].longname),
		   Table_specini[0].longname);
  IO_put_att_text (IO_file_id, nc_svoid, "units",
		   strlen (Table_specini[0].unit), Table_specini[0].unit);

  IO_def_var (IO_file_id, Table_specini[1].name, NC_DOUBLE, 3,
	      initial->nc_dims, &nc_sdid);
  IO_put_att_text (IO_file_id, nc_sdid, "long_name",
		   strlen (Table_specini[1].longname),
		   Table_specini[1].longname);
  IO_put_att_text (IO_file_id, nc_sdid, "units",
		   strlen (Table_specini[1].unit), Table_specini[1].unit);

  initial->nc_dims[2] = nlevp1id;

  IO_def_var (IO_file_id, Table_specini[2].name, NC_DOUBLE, 3,
	      initial->nc_dims, &nc_stpid);
  IO_put_att_text (IO_file_id, nc_stpid, "long_name",
		   strlen (Table_specini[2].longname),
		   Table_specini[2].longname);
  IO_put_att_text (IO_file_id, nc_stpid, "units",
		   strlen (Table_specini[2].unit), Table_specini[2].unit);

  /* gl1-buffer q field */

  initial->nc_dims[0] = nlatid;
  initial->nc_dims[1] = nlevid;

  if (fileformat == 1)
    {
      initial->nc_dims[2] = nlp2id;
    }
  else if (fileformat == 2)
    {
      initial->nc_dims[2] = nlonid;
    }

  IO_def_var (IO_file_id, Table_specini[3].name, NC_DOUBLE, 3,
	      initial->nc_dims, &nc_qid);
  IO_put_att_text (IO_file_id, nc_qid, "long_name",
		   strlen (Table_specini[3].longname),
		   Table_specini[3].longname);
  IO_put_att_text (IO_file_id, nc_qid, "units",
		   strlen (Table_specini[3].unit), Table_specini[3].unit);
  IO_put_att_text (IO_file_id, nc_qid, "axis", 3, "YZX");

  IO_enddef (IO_file_id);

  IO_put_var_double (IO_file_id, nc_glid, g->lat);
  IO_put_var_double (IO_file_id, nc_lonid, g->lon);
  IO_put_var_double (IO_file_id, nc_vctaid, g->vct_a);
  IO_put_var_double (IO_file_id, nc_vctbid, g->vct_b);

  IO_put_var_double (IO_file_id, nc_svoid, svo);
  IO_put_var_double (IO_file_id, nc_sdid, sd);
  IO_put_var_double (IO_file_id, nc_stpid, stp);


  /* add zeroes for old fileformat 1 , leave it as is for fileformat 2 */
  if (fileformat == 1)
    {
      nlev = d->nlev;
      nlat = d->ngl;
      nlon = d->nlon;
      nlp2 = d->nlp2;
      dlen2 = nlev * nlat * nlp2;
      dtmp = (double *) malloc (dlen2 * sizeof (double));

      /* ==============    add zeroes to Q at nlp2  ============== */

      for (ilev = 0; ilev < nlev; ilev++)
	{
	  for (ilat = 0; ilat < nlat; ilat++)
	    {
	      ioff = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlon);
	      ioff2 = get_offset (0, ilat, ilev, 0, 1, nlat, nlev, nlp2);
	      memcpy (&dtmp[ioff2], &q[ioff], nlon * sizeof (double));
	    }
	}

      IO_put_var_double (IO_file_id, nc_qid, dtmp);
    }
  else if (fileformat == 2)
    {
      IO_put_var_double (IO_file_id, nc_qid, q);
    }
  IO_close (IO_file_id);

  return;
}

/* =============================================================================   */

void
read_snc (GRID * g, TIMES * t, char *ifile, SINITIAL * vars)
{
  int nlatid, nlonid, nlp2id, IO_file_id, dlen, d2len;
  int status, nc_varid, ivar, icode, nlp2, nlon, nlat, i;
  int oldformat;
  double *d2tmp;
  int varcode;
  char *varname;
  extern Table Table_surfini[MAXCODE_PTsurfini];
  extern Table Table_dim[];

  netCDF_file ncheader;		/* global accessible netCDF file information */


  ncheader.nc_access_mode = NC_NOWRITE;
  IO_open (ifile, ncheader.nc_access_mode, &(ncheader.nc_file_id));

  IO_file_id = ncheader.nc_file_id;

  status = nc_inq_dimid (IO_file_id, "nlp2", &nlp2id);
  if (status == NC_NOERR)
    {
      oldformat = 1;
      IO_inq_dimid (IO_file_id, "nlon", &nlonid);
      IO_inq_dimid (IO_file_id, "ngl", &nlatid);
      IO_inq_dimlen (IO_file_id, nlp2id, &nlp2);
    }
  else
    {
      oldformat = 0;
      IO_inq_dimid (IO_file_id, Table_dim[0].name, &nlonid);
      IO_inq_dimid (IO_file_id, Table_dim[1].name, &nlatid);
    }

  IO_inq_dimlen (IO_file_id, nlonid, &nlon);
  IO_inq_dimlen (IO_file_id, nlatid, &nlat);

  g->nlon = nlon;
  g->nlat = nlat;

  d2tmp = NULL;
  if (oldformat)
    {
      d2len = nlat * nlp2;
      d2tmp = (double *) malloc (d2len * sizeof (double));
    }

  dlen = nlat * nlon;

  IO_get_att_int (IO_file_id, NC_GLOBAL, "fdate", &t->fdate);
  IO_get_att_int (IO_file_id, NC_GLOBAL, "ftime", &t->ftime);

  icode = 0;
  for (ivar = 0; ivar < MAXCODE_PTsurfini; ivar++)
    {
      varname = Table_surfini[ivar].name;
      varcode = Table_surfini[ivar].code;
      status = nc_inq_varid (IO_file_id, varname, &nc_varid);
      if (status == NC_NOERR)
	{
	  vars->code[icode] = varcode;
	  vars->var[icode] = (double *) malloc (dlen * sizeof (double));
	  if (oldformat)
	    {
	      IO_get_var_double (IO_file_id, nc_varid, d2tmp);
	      for (i = 0; i < nlat; i++)
		{
		  memcpy (&vars->var[icode][nlon * i], &d2tmp[nlp2 * i],
			  nlon * sizeof (double));
		}
	    }
	  else
	    {
	      IO_get_var_double (IO_file_id, nc_varid, vars->var[icode]);
	    }
	  icode++;
	}
    }
  vars->ncodes = icode;

  IO_close (IO_file_id);

  if (oldformat)
    free (d2tmp);

  return;
}

/* ===========================================================================   */

void
write_snc (GRID * g, char *tfile, char *ofile, SINITIAL * vars,
	   int fileformat, DIMENSIONS * d, LABELS * run, TIMES * t)
{
  int nc_varid, nlatid, nlonid, nlp2id, nlevid, nlevp1id, nspid, nvclevid,
    n2id;
  int status, nlp2, nlon, nlat, nc_glid, nc_lonid, nc_vctaid, nc_vctbid;
  time_t current_time;
  struct passwd *user;
  int nvars, icode;
  int I_file_id=0, O_file_id, i, dlen, d2len, ivar, irec;
  netCDF_file template, initial;
  int varcode;
  char *varname, *varlongname, *varunit;
  int oldformat = 1, templ = 1, labels = 1;
  int tdate, ttime;
  double *d2tmp=0;
  extern char *Commandline;
  extern int Verbose;
  extern Table Table_surfini[MAXCODE_PTsurfini];
  extern Table Table_dim[];
  extern char institute[];
  extern int Debug;
  extern int Verbose;

  if (tfile == NULL)
    {
      if (atoi (run->label[0]) == 9999)
	{
	  if (Debug)
	    fprintf (stderr, "will not print labels \n");
	  labels = 0;
	  strcpy (initial.nc_binary_source, "IEEE");
	}
      else if (*run->label[0] == ' ')
	{
	  labels = 1;
	  strcpy (initial.nc_binary_source, "IEEE");
	}
      else if (*run->label[0] == 'O')
	{
	  strcpy (initial.nc_binary_source, "IEEE");
	}
      else
	{
	  strcpy (initial.nc_binary_source, "CRAY");
	}

      templ = 0;
      /* this meeans we are writing from inisurf2snc   */
      if (Debug)
	{
	  fprintf (stderr, "==>  no template-file will be used  ! \n ");
	  fprintf (stderr, " File-format of output-file is: %d\n",
		   fileformat);
	}
      nlon = d->nlon;
      nlat = d->ngl;
      nlp2 = d->nlp2;
    }
  else
    strcpy (initial.nc_binary_source, "IEEE");


  /* setup of global attributes for netCDF file (creation information) */

  strcpy (initial.nc_creation_program, Commandline);
  current_time = time (NULL);
  strcpy (initial.nc_creation_date, ctime (&current_time));
  initial.nc_creation_date[strlen (initial.nc_creation_date) - 1] = '\0';
  user = getpwuid (getuid ());
  strcpy (initial.nc_creation_user, user->pw_gecos);


  initial.nc_access_mode = NC_WRITE;
  IO_create (ofile, initial.nc_access_mode, &(initial.nc_file_id));
  strcpy (initial.nc_file_type, "Initial file surface");

  O_file_id = initial.nc_file_id;

  /* ++++++++++++++++++++++   begin Template    ++++++++++++++++++++++++++++++++++ */
  if (templ)
    {
      template.nc_access_mode = NC_NOWRITE;
      IO_open (tfile, template.nc_access_mode, &(template.nc_file_id));
      I_file_id = template.nc_file_id;

      status = nc_inq_dimid (I_file_id, "nlp2", &nlp2id);
      if (status == NC_NOERR)
	{
	  oldformat = 1;
	  IO_inq_dimid (I_file_id, "nlon", &nlonid);
	  IO_inq_dimid (I_file_id, "ngl", &nlatid);
	  IO_inq_dimlen (I_file_id, nlp2id, &nlp2);

	  IO_inq_dimid (I_file_id, "nlev", &nlevid);
	  IO_inq_dimid (I_file_id, "nlevp1", &nlevp1id);
	  IO_inq_dimid (I_file_id, "nsp", &nspid);
	  IO_inq_dimid (I_file_id, "nvclev", &nvclevid);
	  IO_inq_dimid (I_file_id, "n2", &n2id);

	  IO_inq_dimlen (I_file_id, nlevid, (size_t *) & d->nlev);
	  IO_inq_dimlen (I_file_id, nlevp1id, (size_t *) & d->nlevp1);
	  IO_inq_dimlen (I_file_id, nspid, (size_t *) & d->nsp);
	  IO_inq_dimlen (I_file_id, nvclevid, (size_t *) & (g->nvclev));
	  IO_inq_dimlen (I_file_id, n2id, (size_t *) & d->n2);

	  IO_inq_varid (I_file_id, "vct_a", &nc_vctaid);
	  IO_inq_varid (I_file_id, "vct_b", &nc_vctbid);
	  g->vct_a = (double *) malloc (g->nvclev * sizeof (double));
	  g->vct_b = (double *) malloc (g->nvclev * sizeof (double));
	  IO_get_var_double (I_file_id, nc_vctaid, g->vct_a);
	  IO_get_var_double (I_file_id, nc_vctbid, g->vct_b);
	}
      else
	{
	  oldformat = 0;
	  IO_inq_dimid (I_file_id, Table_dim[0].name, &nlonid);
	  IO_inq_dimid (I_file_id, Table_dim[1].name, &nlatid);
	}

      if (!oldformat)
	fileformat = 2;

      if (!fileformat)
	{
	  if (oldformat)
	    fileformat = 1;
	  else
	    fileformat = 2;
	}

      if (Debug)
	{
	  if (oldformat)
	    fprintf (stderr, " Template-file has old file-format! \n");
	  else
	    fprintf (stderr, " Template-file has new file-format! \n");

	  fprintf (stderr, " File-format of output-file is: %d\n",
		   fileformat);
	}


      IO_inq_dimlen (I_file_id, nlonid, &nlon);
      IO_inq_dimlen (I_file_id, nlatid, &nlat);

      d->nlon = nlon;
      d->ngl = nlat;


      if (vars->ncodes)
	{
	  /* check if LAT is the same as in SRV-file */
	  if (d->ngl != g->nlat)
	    {
	      fprintf (stderr,
		       "number of latitudes of Template-file and SERVICE-file differ: !\n"
		       "  ====>   %i : %i \a \n ", d->ngl, g->nlat);
	      exit (1);
	    }
	  /* check if LON is the same as in SRV-file */
	  if (d->nlon != g->nlon)
	    {
	      fprintf (stderr,
		       "number of longitudes of Template-file and SERVICE -file differ: ! \n"
		       "   ====>   %i : %i \a \n ", d->nlon, g->nlon);
	      exit (1);
	    }
	}

      IO_inq_varid (I_file_id, Table_dim[1].name, &nc_glid);
      IO_inq_varid (I_file_id, Table_dim[0].name, &nc_lonid);

      g->lat = (double *) malloc (d->ngl * sizeof (double));
      g->lon = (double *) malloc (d->nlon * sizeof (double));

      IO_get_var_double (I_file_id, nc_glid, g->lat);
      IO_get_var_double (I_file_id, nc_lonid, g->lon);

      /*   now read each  variables that we didn't got from the srv-file    */
      d2tmp = NULL;
      if (oldformat)
	{
	  d2len = nlat * nlp2;
	  d2tmp = (double *) malloc (d2len * sizeof (double));
	}

      dlen = nlat * nlon;

      irec = 1;
      IO_get_att_int (I_file_id, NC_GLOBAL, "fdate", &tdate);
      IO_get_att_int (I_file_id, NC_GLOBAL, "ftime", &ttime);

      if (Verbose)
	{
	  fprintf (stderr, "\n");
	  fprintf (stderr, " Used data from Template file\n");
	  fprintf (stderr, " ============================\n");
	  fprintf (stderr, " Rec :     Date   Time Code  Level  Lon  Lat :"
		   "     Minimum        Mean     Maximum\n");
	}
      for (ivar = 0; ivar < MAXCODE_PTsurfini; ivar++)
	{
	  varname = Table_surfini[ivar].name;
	  varcode = Table_surfini[ivar].code;
	  for (i = 0; i < vars->ncodes; i++)
	    {
	      if (vars->code[i] == varcode)
		break;
	    }
	  if (i == vars->ncodes)
	    {
	      status = nc_inq_varid (I_file_id, varname, &nc_varid);
	      if (status == NC_NOERR)
		{
		  nvars = vars->ncodes;
		  vars->code[nvars] = varcode;
		  vars->var[nvars] =
		    (double *) malloc (dlen * sizeof (double));
		  if (oldformat)
		    {
		      IO_get_var_double (I_file_id, nc_varid, d2tmp);
		      for (i = 0; i < nlat; i++)
			{
			  memcpy (&vars->var[nvars][nlon * i],
				  &d2tmp[nlp2 * i], nlon * sizeof (double));
			}
		    }
		  else
		    {
		      IO_get_var_double (I_file_id, nc_varid,
					 vars->var[nvars]);
		    }
		  if (Verbose)
		    {
		      fprintf (stderr, "%4d : %8d %6d %4d%7d %4d %4d :", irec,
			       tdate, ttime / 100, varcode, 0, nlon, nlat);
		      fprintf (stderr, "%12.6g%12.6g%12.6g\n",
			       min_darray (vars->var[nvars], dlen),
			       mean_darray (vars->var[nvars], dlen),
			       max_darray (vars->var[nvars], dlen));
		    }
		  irec++;
		  vars->ncodes++;
		}
	    }
	}

    }
  /* ++++++++++++++++++++++   end Template    ++++++++++++++++++++++++++++++++++ */
  if (fileformat == 1)
    {
      IO_def_dim (O_file_id, "ngl", d->ngl, &nlatid);
      IO_def_dim (O_file_id, "nlon", d->nlon, &nlonid);
      IO_def_dim (O_file_id, "nlp2", nlp2, &nlp2id);
      IO_def_dim (O_file_id, "nlev", d->nlev, &nlevid);
      IO_def_dim (O_file_id, "nlevp1", d->nlevp1, &nlevp1id);
      IO_def_dim (O_file_id, "nsp", d->nsp, &nspid);
      IO_def_dim (O_file_id, "nvclev", g->nvclev, &nvclevid);
      IO_def_dim (O_file_id, "n2", 2, &n2id);
    }
  else if (fileformat == 2)
    {
      IO_def_dim (O_file_id, Table_dim[1].name, d->ngl, &nlatid);
      IO_def_dim (O_file_id, Table_dim[0].name, d->nlon, &nlonid);
    }

  /* IO_put_att_text (O_file_id, NC_GLOBAL, "Conventions", 6, "COARDS"); */
  if (fileformat == 2)
    {
      strcpy (initial.nc_conventions, "COARDS");
      strcpy (initial.nc_model, "ECHAM");

      IO_put_att_text (O_file_id, NC_GLOBAL, "Conventions",
		       strlen (initial.nc_conventions),
		       initial.nc_conventions);
      IO_put_att_text (O_file_id, NC_GLOBAL, "model",
		       strlen (initial.nc_model), initial.nc_model);

    }

  strcpy (initial.nc_institution, institute);
  IO_put_att_text (O_file_id, NC_GLOBAL, "institution",
		   strlen (initial.nc_institution), initial.nc_institution);
  IO_put_att_text (O_file_id, NC_GLOBAL, "file_type",
		   strlen (initial.nc_file_type), initial.nc_file_type);
  IO_put_att_text (O_file_id, NC_GLOBAL, "source_type",
		   strlen (initial.nc_binary_source),
		   initial.nc_binary_source);
  IO_put_att_text (O_file_id, NC_GLOBAL, "user",
		   strlen (initial.nc_creation_user),
		   initial.nc_creation_user);
  IO_put_att_text (O_file_id, NC_GLOBAL, "history",
		   strlen (initial.nc_creation_program),
		   initial.nc_creation_program);
  IO_put_att_text (O_file_id, NC_GLOBAL, "created",
		   strlen (initial.nc_creation_date),
		   initial.nc_creation_date);

  /* ++++++++++++++++++++++   begin Template  att  ++++++++++++++++++++++++++++++++++ */
  if (templ)
    {
      IO_copy_att (I_file_id, NC_GLOBAL, "title", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_1", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_2", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_3", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_4", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_5", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_6", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_7", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "label_8", O_file_id, NC_GLOBAL);

      IO_copy_att (I_file_id, NC_GLOBAL, "fdate", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "ftime", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "vdate", O_file_id, NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "vtime", O_file_id, NC_GLOBAL);

      IO_copy_att (I_file_id, NC_GLOBAL, "spherical_truncation_n", O_file_id,
		   NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "spherical_truncation_m", O_file_id,
		   NC_GLOBAL);
      IO_copy_att (I_file_id, NC_GLOBAL, "spherical_truncation_k", O_file_id,
		   NC_GLOBAL);

      IO_close (I_file_id);

    }
  else
    {

      IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "title", 29,
		       "No information - initial file");

      if (labels)
	{
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_1",
			   strlen (run->label[0]), run->label[0]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_2",
			   strlen (run->label[1]), run->label[1]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_3",
			   strlen (run->label[2]), run->label[2]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_4",
			   strlen (run->label[3]), run->label[3]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_5",
			   strlen (run->label[4]), run->label[4]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_6",
			   strlen (run->label[5]), run->label[5]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_7",
			   strlen (run->label[6]), run->label[6]);
	  IO_put_att_text (initial.nc_file_id, NC_GLOBAL, "label_8",
			   strlen (run->label[7]), run->label[7]);

	}
      IO_put_att_int (initial.nc_file_id, NC_GLOBAL, "fdate", NC_INT, 1,
		      &(t->fdate));
      IO_put_att_int (initial.nc_file_id, NC_GLOBAL, "ftime", NC_INT, 1,
		      &(t->ftime));
      IO_put_att_int (initial.nc_file_id, NC_GLOBAL, "vdate", NC_INT, 1,
		      &(t->vdate));
      IO_put_att_int (initial.nc_file_id, NC_GLOBAL, "vtime", NC_INT, 1,
		      &(t->vtime));
      if (labels)
	{
	  IO_put_att_int (initial.nc_file_id, NC_GLOBAL,
			  "spherical_truncation_n", NC_INT, 1, &(d->nn));
	  IO_put_att_int (initial.nc_file_id, NC_GLOBAL,
			  "spherical_truncation_m", NC_INT, 1, &(d->nm));
	  IO_put_att_int (initial.nc_file_id, NC_GLOBAL,
			  "spherical_truncation_k", NC_INT, 1, &(d->nk));
	}
    }
  /* ++++++++++++++++++++++   end Template att  ++++++++++++++++++++++++++++++++++ */


  initial.nc_dims[0] = nlatid;
  IO_def_var (O_file_id, Table_dim[1].name, NC_DOUBLE, 1, initial.nc_dims,
	      &nc_glid);
  IO_put_att_text (O_file_id, nc_glid, "units", strlen (Table_dim[1].unit),
		   Table_dim[1].unit);
  IO_put_att_text (O_file_id, nc_glid, "long_name",
		   strlen (Table_dim[1].longname), Table_dim[1].longname);

  initial.nc_dims[0] = nlonid;
  IO_def_var (O_file_id, Table_dim[0].name, NC_DOUBLE, 1, initial.nc_dims,
	      &nc_lonid);
  IO_put_att_text (O_file_id, nc_lonid, "units", strlen (Table_dim[0].unit),
		   Table_dim[0].unit);
  IO_put_att_text (O_file_id, nc_lonid, "long_name",
		   strlen (Table_dim[0].longname), Table_dim[0].longname);

  if (fileformat == 1)
    {
      initial.nc_dims[0] = nvclevid;
      IO_def_var (O_file_id, "vct_a", NC_DOUBLE, 1, initial.nc_dims,
		  &nc_vctaid);
      IO_put_att_text (O_file_id, nc_vctaid, "units", 1, "");
      IO_put_att_text (O_file_id, nc_vctaid, "long_name", 35,
		       "vertical-coordinate parameter set A");

      initial.nc_dims[0] = nvclevid;
      IO_def_var (O_file_id, "vct_b", NC_DOUBLE, 1, initial.nc_dims,
		  &nc_vctbid);
      IO_put_att_text (O_file_id, nc_vctbid, "units", 1, "");
      IO_put_att_text (O_file_id, nc_vctbid, "long_name", 35,
		       "vertical-coordinate parameter set B");
    }

  /* So long no more real data fields are defined ... */

  initial.nc_dims[0] = nlatid;
  if (fileformat == 1)
    initial.nc_dims[1] = nlp2id;
  else if (fileformat == 2)
    initial.nc_dims[1] = nlonid;


  for (ivar = 0; ivar < MAXCODE_PTsurfini; ivar++)
    {
      varcode = Table_surfini[ivar].code;
      varname = Table_surfini[ivar].name;
      varlongname = Table_surfini[ivar].longname;
      varunit = Table_surfini[ivar].unit;

      for (icode = 0; icode < vars->ncodes; icode++)
	{
	  if (vars->code[icode] == varcode)
	    break;

	}

      if (icode < vars->ncodes)
	{
	  IO_def_var (O_file_id, varname, NC_DOUBLE, 2, initial.nc_dims,
		      &vars->varid[icode]);
	  IO_put_att_text (O_file_id, vars->varid[icode], "long_name",
			   strlen (varlongname), varlongname);
	  IO_put_att_text (O_file_id, vars->varid[icode], "units",
			   strlen (varunit), varunit);
	}
    }


  IO_enddef (O_file_id);

  IO_put_var_double (O_file_id, nc_glid, g->lat);
  IO_put_var_double (O_file_id, nc_lonid, g->lon);
  if (fileformat == 1)
    {
      IO_put_var_double (O_file_id, nc_vctaid, g->vct_a);
      IO_put_var_double (O_file_id, nc_vctbid, g->vct_b);
    }

  if (!templ && oldformat)
    {
      d2len = nlat * nlp2;
      d2tmp = (double *) malloc (d2len * sizeof (double));
    }
  for (icode = 0; icode < vars->ncodes; icode++)
    {
      if (fileformat == 1)
	{
	  for (i = 0; i < nlat; i++)
	    {
	      memcpy (&d2tmp[nlp2 * i], &vars->var[icode][nlon * i],
		      nlon * sizeof (double));
	    }
	  IO_put_var_double (O_file_id, vars->varid[icode], d2tmp);
	}
      else if (fileformat == 2)
	{
	  IO_put_var_double (O_file_id, vars->varid[icode], vars->var[icode]);
	}
    }

  IO_close (O_file_id);

  if (oldformat)
    free (d2tmp);

  return;
}
