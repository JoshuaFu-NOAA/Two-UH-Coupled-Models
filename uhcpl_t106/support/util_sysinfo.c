/*
 * This file contains a set of routines to provide operating system information
 * for f90 on Unix machines. Unfortunatelly there is no implemented set of
 * such functions available. To restrict the f90 source code to f90 standard
 * only this mostly POSIX C compliant implementations are developed. One
 * basic convention is to name all routines with a leading util_ to make
 * this visible to the developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 *
 * $Id: util_sysinfo.c,v 1.9 2000/02/03 10:21:39 m214003 Exp $
 *
 */

#include "util_fortran.h"

#include <sys/param.h>
#include <sys/utsname.h>
#include <netdb.h>
 
#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif HAVE_SYS_UNISTD_H
#include <sys/unistd.h>
#endif
 
#include <pwd.h>
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#ifdef __uxp__
#include <sys/systeminfo.h>
#endif

#ifndef HAVE_VALLOC
#define valloc malloc
#endif

#if defined(FORTRANCAPS)
#define util_os_system_ UTIL_OS_SYSTEM
#define util_user_name_ UTIL_USER_NAME
#define util_node_name_ UTIL_NODE_NAME
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define util_os_system_ util_os_system__
#define util_user_name_ util_user_name__
#define util_node_name_ util_node_name__
#elif !defined(FORTRANUNDERSCORE)
#define util_os_system_ util_os_system
#define util_user_name_ util_user_name
#define util_node_name_ util_node_name
#endif


#define ALLOC(n, size)  valloc((n) * (size))
#define FREE(x) { free((x)); (x) = NULL; }

#define COPY_TO(to,from_ftn,len) \
{ strncpy((to),_fcdtocp((from_ftn)),(len)); (to)[(len)] = '\0'; }

#define COPY_FROM(to_ftn,from,len) \
{ strncpy(_fcdtocp(to_ftn),(from),(len)); }

#define STRLEN(x) strlen(_fcdtocp(x))

/* funcion implemetations */ 

FORTRAN_CALL
void util_os_system_  (_fcd name, INT *actual_len
#ifndef CRAY
		   , INT name_len
#endif
)
{
#ifdef CRAY
  INT name_len = _fcdlen(name);
#endif

  char *p_name = (char *) ALLOC((size_t)name_len+1,sizeof(*p_name));

  struct utsname utname;
 
  uname( &utname );

  /* Cray doesn't handle operating system information proper */

#ifdef CRAY
  strcpy (p_name, "UNICOS");
#else
  strcpy (p_name, utname.sysname);
#endif
  strcat (p_name, " ");
  strcat (p_name, utname.release);
  strcat (p_name, " on ");
  strcat (p_name, utname.machine);

  COPY_FROM(name, p_name, (size_t)name_len);

  *actual_len = (strlen(p_name) < name_len) ? STRLEN(name) : name_len;
  if (*actual_len < name_len) {
    memset(&_fcdtocp(name)[*actual_len],' ',(size_t)(name_len - *actual_len));
  }
  
  FREE (p_name);
}

FORTRAN_CALL  
void util_user_name_ (_fcd name, INT *actual_len
#ifndef CRAY
		   , INT name_len
#endif
) 
{
#ifdef CRAY
  INT name_len   = _fcdlen(name);
#endif

  char *p_name = (char *) ALLOC((size_t)name_len+1,sizeof(*p_name));

  struct passwd *current;
 
  current = getpwuid(getuid());
  if (current == NULL) {
    strcpy (p_name, "unknown user name");
  } else {
    if (strlen(current->pw_name) == 0) {
      strcpy (p_name, "unknown user name");
    } else {
      strcpy (p_name, current->pw_gecos);
      strcat (p_name, " (");
      strcat (p_name, current->pw_name);
      strcat (p_name, ")");
    }
  }    

  COPY_FROM(name, p_name, (size_t)name_len);

  *actual_len = (strlen(p_name) < name_len) ? STRLEN(name) : name_len;
  if (*actual_len < name_len) {
    memset(&_fcdtocp(name)[*actual_len],' ',(size_t)(name_len - *actual_len));
  }
  
  FREE (p_name);
}

FORTRAN_CALL 
void util_node_name_ ( _fcd name, INT  *actual_len
#ifndef CRAY
                       ,INT name_len
#endif
)
{
#ifdef CRAY
  INT name_len = _fcdlen(name);
#endif

  char *p_name = (char *) ALLOC((size_t)name_len+1,sizeof(*p_name));

  char *hostname;

  if ((hostname = getenv ("HOST")) == NULL) {
    strcpy (p_name, "unknown");
  } else {
    strcpy (p_name, hostname);
  }

  COPY_FROM(name, p_name, (size_t)name_len);

  *actual_len = (strlen(p_name) < name_len) ? STRLEN(name) : name_len;
  if (*actual_len < name_len) {
    memset(&_fcdtocp(name)[*actual_len],' ',(size_t)(name_len - *actual_len));
  }
  
  FREE (p_name);
}
