#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>
#include <sys/stat.h>


#include "prep.h"


void
ncheader (char *argv[], netCDF_file * restart, HEADER * info, DIMENSIONS * d,
	  GRID * g, LABELS * run, TIMES * t)
{
  time_t current_time;
  struct passwd *user;

  /* setup of global attributes for netCDF file (creation information) */

  strcpy (restart->nc_creation_program, argv[0]);
  current_time = time (NULL);
  strcpy (restart->nc_creation_date, ctime (&current_time));
  restart->nc_creation_date[strlen (restart->nc_creation_date) - 1] = '\0';
  user = getpwuid (getuid ());
  strcpy (restart->nc_creation_user, user->pw_gecos);


  strcpy (restart->nc_binary_source, info->def);
  strcpy (restart->nc_file_type, "Restart file");

  /* IO_put_att_text (restart->nc_file_id, NC_GLOBAL, "Conventions", 6, "COARDS"); */
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "file_type",
		   strlen (restart->nc_file_type), restart->nc_file_type);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "source_type",
		   strlen (restart->nc_file_type), restart->nc_binary_source);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "user",
		   strlen (restart->nc_creation_user),
		   restart->nc_creation_user);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "history",
		   strlen (restart->nc_creation_program),
		   restart->nc_creation_program);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "created",
		   strlen (restart->nc_creation_date),
		   restart->nc_creation_date);


  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_1",
		   strlen (run->label[0]), run->label[0]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_2",
		   strlen (run->label[1]), run->label[1]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_3",
		   strlen (run->label[2]), run->label[2]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_4",
		   strlen (run->label[3]), run->label[3]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_5",
		   strlen (run->label[4]), run->label[4]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_6",
		   strlen (run->label[5]), run->label[5]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_7",
		   strlen (run->label[6]), run->label[6]);
  nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "label_8",
		   strlen (run->label[7]), run->label[7]);

  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "fdate", NC_INT, 1,
		  &(t->fdate));
  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "ftime", NC_INT, 1,
		  &(t->ftime));
  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "vdate", NC_INT, 1,
		  &(t->vdate));
  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "vtime", NC_INT, 1,
		  &(t->vtime));

  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "spherical_truncation_n",
		  NC_INT, 1, &(d->nn));
  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "spherical_truncation_m",
		  NC_INT, 1, &(d->nm));
  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "spherical_truncation_k",
		  NC_INT, 1, &(d->nk));

  nc_put_att_int (restart->nc_file_id, NC_GLOBAL, "nstep", NC_INT, 1,
		  &(d->nstep));

  nc_put_att_double (restart->nc_file_id, NC_GLOBAL, "timestep", NC_DOUBLE, 1,
		     &(d->dt));

  if (d->nhtrac == 0)
    {
      nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "tracer_definition",
		       19, "No tracer available");
    }
  else
    {
      nc_put_att_text (restart->nc_file_id, NC_GLOBAL, "tracer_definition",
		       16, "Tracer available");
    }

}
