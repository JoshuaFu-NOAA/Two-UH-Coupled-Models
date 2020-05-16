/*
 * This file contains a set of routines to provide operating system support
 * for f90 on Unix machines. Unfortunatelly there is no implemented set of
 * such functions available. To restrict the f90 source code to f90 standard
 * only this mostly Posix C compliant implementations are developed. One
 * basic convention is to name all routines with a leading util_ to make
 * this visible to the developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 *
 * $Id: util_system.c,v 1.7 1999/09/17 11:18:47 m214089 Exp $
 *
 */

#include "util_fortran.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <sys/time.h>
#include <sys/resource.h>

#ifndef HAVE_VALLOC
#define valloc malloc
#endif

#if defined(FORTRANCAPS)
#define util_cgetenv_          UTIL_CGETENV
#define util_igetenv_          UTIL_IGETENV
#define util_abort_            UTIL_ABORT
#define util_cputime_          UTIL_CPUTIME
#define util_walltime_         UTIL_WALLTIME
#define util_system_           UTIL_SYSTEM
#define util_address_bit_size_ UTIL_ADDRESS_BIT_SIZE
#define util_indef_            UTIL_INDEF
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define util_cgetenv_          util_cgetenv__
#define util_igetenv_          util_igetenv__
#define util_abort_            util_abort__
#define util_cputime_          util_cputime__
#define util_walltime_         util_walltime__
#define util_system_           util_system__
#define util_address_bit_size_ util_address_bit_size__
#define util_ihpstat_          util_ihpstat__
#elif !defined(FORTRANUNDERSCORE)
#define util_cgetenv_          util_cgetenv
#define util_igetenv_          util_igetenv
#define util_abort_            util_abort
#define util_cputime_          util_cputime
#define util_walltime_         util_walltime
#define util_system_           util_system
#define util_address_bit_size_ util_address_bit_size
#define util_ihpstat_          util_ihpstat
#endif

#define ALLOC(n, size)  valloc((n) * (size))
#define FREE(x) { free((x)); (x) = NULL; }

#define COPY_TO(to,from_ftn,len) \
{ strncpy((to),_fcdtocp((from_ftn)),(len)); (to)[(len)] = '\0'; }

#define COPY_FROM(to_ftn,from,len) \
{ strncpy(_fcdtocp(to_ftn),(from),(len)); }

#define STRLEN(x) strlen(_fcdtocp(x))

/*
 * util_cgetenv - this routine gives access to character environment
 * variables. envname is the environment variable name, def a possible 
 * default string, output the retrieved name, and actual_len the
 * length of the output string.
 *
 */

FORTRAN_CALL
void util_cgetenv_(_fcd envname, _fcd def, _fcd output,
         	   INT *actual_len
#ifndef CRAY
		   , INT envname_len, INT def_len, INT output_len
#endif
)
{
#ifdef CRAY
  INT envname_len = _fcdlen(envname);
  INT def_len     = _fcdlen(def);
  INT output_len  = _fcdlen(output);
#endif
  char *p_envname = (char *) ALLOC((size_t) envname_len+1,sizeof(*p_envname));
  char *p;
  char *p_def = NULL;

  COPY_TO(p_envname, envname, (size_t) envname_len);
  p = getenv(p_envname);
  FREE(p_envname);

  if (!p) {
    p = p_def = (char *) ALLOC((size_t) def_len+1,sizeof(*p_def));
    COPY_TO(p_def, def, (size_t) def_len);
  }
  
  COPY_FROM(output, p, (size_t) output_len);

  *actual_len = (strlen(p) < output_len) ? STRLEN(output) : output_len;
  if (*actual_len < output_len) {
    memset(&_fcdtocp(output)[*actual_len],' ',(size_t) (output_len - *actual_len));
  }

  if (p_def) FREE(p_def);
}

/*
 * util_igetenv - this routine gives access to integer environment
 * variables. envname is the environment variable name, def a possible 
 * default value, and output the retrieved value.
 *
 */

FORTRAN_CALL
void util_igetenv_(_fcd envname, INT *def, INT *output
#ifndef CRAY
		   , INT envname_len
#endif
)
{
#ifdef CRAY
  INT envname_len = _fcdlen(envname);
#endif
  char *p_envname = (char *) ALLOC((size_t) envname_len+1,sizeof(*p_envname));
  char *p;

  COPY_TO(p_envname, envname, (size_t) envname_len);
  p = getenv(p_envname);

  *output = p ? atoi(p) : *def;

  FREE(p_envname);
}

/*
 * util_abort - the hard way out ...
 *
 */

FORTRAN_CALL
void util_abort_()
{
  /*  abort(); */
  exit(1);
}

/* Portable CPU-timer (User + Sys); also WALL CLOCK-timer */

#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/param.h>

extern clock_t times (struct tms *buffer);
#define clock_ticks ( (double) sysconf(_SC_CLK_TCK) )

FORTRAN_CALL
INT util_cputime_(REAL *user_time, REAL *system_time)
{
  struct tms tbuf;
  if (times(&tbuf) == -1) return ((INT) (-1)); 

  *user_time   = ((REAL) tbuf.tms_utime) / clock_ticks; 
  *system_time = ((REAL) tbuf.tms_stime) / clock_ticks;

  return (0);
}

#ifndef __uxp__
FORTRAN_CALL
REAL util_walltime_()
{
  static REAL time_init = 0.;
  REAL time_in_secs;
  struct timeval tbuf;
  if (gettimeofday(&tbuf,NULL) == -1) perror("UTIL_WALLTIME");

  if (time_init == 0.) time_init =
    (REAL) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0);

  time_in_secs =
  (REAL) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0) - time_init;

  return (time_in_secs);
}
#else
FORTRAN_CALL
REAL util_walltime_()
{
  REAL time_in_secs;
  double w;
  static double wallref = 0;
  extern FORTRAN_CALL gettod_(double *);
  if (wallref == 0) gettod_(&wallref);
  gettod_(&w);
  time_in_secs = (w - wallref) * 0.000001;
  return time_in_secs;
}
#endif

#ifdef __uxp__

#include <sys/types.h>
#include <sys/param.h>
#include <sys/signal.h>
#include <sys/fault.h>
#include <sys/syscall.h>
#include <sys/procfs.h>
#include <sys/proc.h>
#include <fcntl.h>

static int fujitsu_getrusage(int who, struct rusage *rusage)
{
  int rc = -1;

  if (rusage) rusage->ru_maxrss = 0;

  if (who == RUSAGE_SELF && rusage) {
    static int maxrss =  0;
    static int oldpid = -1;
    static char procfile[20] = "";
    static char *pf = NULL;
    static proc_t proc;
    int pid = getpid();
    static int fildes = -1;
    unsigned int size;

    if (oldpid != pid) {
      oldpid = pid;
      maxrss = 0;
      pf = NULL;
    }

    if (!pf) {
      sprintf(procfile,"/proc/%d",pid);
      pf = procfile;
      fildes = open(procfile, O_RDONLY);
    }

    if (fildes == -1) return rc;

    if (ioctl(fildes, PIOCGETPR, &proc) == -1) {
      perror("ioctl@fujitsu_getrusage(PIOCGETPR)");
      return rc;
    }

    size  = proc.p_brksize + proc.p_stksize;
    if (size > maxrss) maxrss = size*8;
    rusage->ru_maxrss = maxrss;

    rc = 0;
  }
  return rc;
}
#endif 

FORTRAN_CALL
INT util_system_ (_fcd s
#ifndef CRAY
		   , INT s_len
#endif
)
{
#ifdef CRAY
  INT s_len = _fcdlen(s);
#endif

  char *buff = (char *) ALLOC((size_t) s_len+1, sizeof(*buff));
  int rv;

  COPY_TO(buff, s, (size_t) s_len);
  rv = system(buff);
  FREE(buff);

  return rv;
}

/* returns the size of a (FILE *) in bit */

INT util_address_bit_size_ (void)
{
  return ( (INT) 8*sizeof (FILE *));
}

#ifndef CRAY

FORTRAN_CALL
INT util_ihpstat_(INT *option)
{
  /* 
   * *option =  1 -> Current heap length
   * *option = 12 -> Amount by which the heap can grow (without paging) 
   */ 
  INT ret_value = 0;
  INT pagesize;
 
#ifdef __uxp__
  pagesize = 1;
#else
  pagesize = (INT) sysconf(_SC_PAGESIZE);
#endif
  if (*option == 1) {
    struct rusage rusage;
#ifdef __uxp__
    fujitsu_getrusage(0, &rusage);
#else
    getrusage(RUSAGE_SELF, &rusage);
#endif
#ifdef sun
    fprintf (stderr, "This option is not supported on Solaris 2.5.1.\n");
#else    
    ret_value = (rusage.ru_maxrss * pagesize + 7); /* in byte */
#endif
  } else if (*option == 12) {
#ifdef sun
    ret_value = pagesize*sysconf(_SC_AVPHYS_PAGES);
#else
    fprintf (stderr, "This option is not supported.\n");
#endif
  }
  return ret_value;
}
 
#endif /* not def CRAY */
 



