/*
 * This files contains the binary C-IO interface to read and write unblocked
 * files originally provided through the EMOS library from ECMWF. 
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:   30.5.1998
 *
 * $Id: util_pbio.c,v 1.4 1998/12/11 13:23:58 m214003 Exp $
 *
 */

#ifndef EMOS

#include "util_fortran.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef HAVE_VALLOC
#define valloc malloc
#endif

#if defined(FORTRANCAPS)
#define pbopen_  PBOPEN  
#define pbread_  PBREAD  
#define pbwrite_ PBWRITE 
#define pbseek_  PBSEEK  
#define pbclose_ PBCLOSE 
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define pbopen_  pbopen__  
#define pbread_  pbread__  
#define pbwrite_ pbwrite__ 
#define pbseek_  pbseek__  
#define pbclose_ pbclose__ 
#elif !defined(FORTRANUNDERSCORE)
#define pbopen_  pbopen  
#define pbread_  pbread  
#define pbwrite_ pbwrite 
#define pbseek_  pbseek  
#define pbclose_ pbclose 
#endif


#define ALLOC(n, size)  valloc((n) * (size))
#define FREE(x) { free((x)); (x) = NULL; }
 
#define COPY_TO(to,from_ftn,len) \
{ strncpy((to),_fcdtocp((from_ftn)),(len)); (to)[(len)] = '\0'; }
 
#define COPY_FROM(to_ftn,from,len) \
{ strncpy(_fcdtocp(to_ftn),(from),(len)); }
 
#define STRLEN(x) strlen(_fcdtocp(x))

/*
 * Purpose:	Opens file, return UNIX FILE pointer.
 *
 * Function returns: long iret:  -1 = Could not open file.
 *                               -2 = Invalid file name.
 *                               -3 = Invalid open mode specified
 *                                0 = OK.
 */

FORTRAN_CALL
void pbopen_ (FILE **unit, _fcd name, _fcd mode, INT *iret
#ifndef CRAY
              , INT name_len, INT mode_len
#endif
)
{
#ifdef CRAY
  INT name_len = _fcdlen(name);
  INT mode_len = _fcdlen(mode);
#endif
  char *p_name = (char *) ALLOC((size_t) name_len+1,sizeof(*p_name));
  char *p_mode = (char *) ALLOC((size_t) mode_len+1,sizeof(*p_mode));
  char *p;
  char flags[4];

  COPY_TO(p_name, name, (size_t) name_len);
  COPY_TO(p_mode, mode, (size_t) mode_len);

  strcpy (flags, "");

  *unit = NULL;
  *iret = 0;

/* strip trailing blanks */

  p  = p_name + name_len - 1 ;
  while(*p == ' ') {
    *p = 0;
    p--;
  }
		
/* build open flags from "modes" */

  p = p_mode;
  while (*p && (strlen(flags) < 3 )) {
    switch (*p) {
    case 'a':
    case 'A':
      strcat (flags, "a");
      break;
    case 'c':
    case 'C':
    case 'w':
    case 'W':
      strcat (flags, "w");
      break;
    case 'r':
    case 'R':
      strcat (flags, "r");
      break;
    default:
      *iret = -3;
      return;
    }
    p++;
  }

/* if read/write change flags */

  if (!strcmp (flags, "wr") || !strcmp (flags, "rw")) {
    strcpy (flags, "r+w" );
  }
  *unit = fopen (p_name, flags);
  
  if (*unit == NULL) {
    perror (p_name);
    *iret = -1;
  }
  
  FREE (p_name);
  FREE (p_mode);

  return;
}

/*
 *
 * Purpose:	Seeks to specified location in file.
 *
 * Function returns:	int status : 	-2 = error in handling file,
 *					-1 = end-of-file
 *                              otherwise,  = byte offset from start of file.
 *
 *	whence	= 0, from start of file
 *		= 1, from current position
 *		= 2, from end of file.	
 */

FORTRAN_CALL
void pbseek_ (FILE **unit, INT *offset, INT *whence, INT *iret)
{
  long my_offset = *offset;

/* must use negative offset if working from end-of-file	*/

  if ( *whence == 2) my_offset = - labs(my_offset);

  *iret = fseek(*unit, my_offset, (int) *whence);

  if(*iret != 0) {
    if ( ! feof(*unit) ) {
      *iret = -2;		/* error in file-handling */
      perror("pbseek");
    } else {
      *iret = -1;		/* end-of-file	*/
    }
    clearerr(*unit);
    return;
  }

  *iret = (INT) ftell(*unit);		/* byte offset from start of file */

  return;
}

/*
 *
 * Purpose:      Reads a block of bytes from a file.
 *
 * Function returns:     int status :   -2 = error in reading file,
 *                                       -1 = end-of-file,
 *                               otherwise, = number of bytes read.
 */

FORTRAN_CALL
void pbread_ (FILE ** unit, char *buffer, INT *nbytes, INT *iret)
{
 
  if ( (*iret = fread(buffer, 1, (size_t) *nbytes, *unit) ) != *nbytes) {
    /*      Read problem */
    if ( ! feof(*unit) ) {
      *iret = -2;             /*  error in file-handling  */
      perror("pbread");
      clearerr(*unit);
      return;
    } else {
      *iret = -1;             /*  end-of-file */
      clearerr(*unit);
    }
  }

  return;
}

/*
 * Purpose:	Writes a block of bytes to a file.
 *
 * Function returns:	int status : -1 = Could not write to file.
 *                                    >=0 = Number of bytes written.
 */

FORTRAN_CALL
void pbwrite_ (FILE **unit, char *buffer, INT *nbytes, INT *iret)
{
  if ( (*iret = fwrite(buffer, 1, (size_t) *nbytes, *unit) ) != *nbytes) {	
    /* Problem with write */   
    perror("pbwrite");
    *iret = -1;
  }

  return;
}

/*
 *
 * Purpose:	Closes file.
 *
 * Function returns:	int status : non-0 = error in handling file.
 *                                     	  0 = OK.
 */

FORTRAN_CALL
void pbclose_ (FILE **unit, INT *iret)
{
  *iret = fclose(*unit);
  
  if(*iret != 0) perror("pbclose");

  return;
}

#else
#define UNUSED
#endif
