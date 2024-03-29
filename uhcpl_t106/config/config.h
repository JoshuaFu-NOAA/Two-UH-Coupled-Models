/* config/config.h.  Generated automatically by configure.  */
#ifndef CONFIG_H
#define CONFIG_H

/* config.h.in.  Generated automatically from configure.in by autoheader.  */

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define if you have the ANSI C header files.  */
#define STDC_HEADERS 1

/* Define if you can safely include both <sys/time.h> and <time.h>.  */
#define TIME_WITH_SYS_TIME 1

/* Define if your processor stores words with the most significant
   byte first (like Motorola and SPARC, unlike Intel and VAX).  */
/* #undef WORDS_BIGENDIAN */

/* The number of bytes in a char.  */
#define SIZEOF_CHAR 1

/* The number of bytes in a double.  */
#define SIZEOF_DOUBLE 8

/* The number of bytes in a float.  */
#define SIZEOF_FLOAT 4

/* The number of bytes in a int.  */
#define SIZEOF_INT 4

/* The number of bytes in a long.  */
#define SIZEOF_LONG 8

/* The number of bytes in a long double.  */
#define SIZEOF_LONG_DOUBLE 16

/* The number of bytes in a long long.  */
#define SIZEOF_LONG_LONG 8

/* The number of bytes in a short.  */
#define SIZEOF_SHORT 2

/* The number of bytes in a int *.  */
#define SIZEOF_INT_P 8

/* Define if you have the getrusage function.  */
#define HAVE_GETRUSAGE 1

/* Define if you have the gettimeofday function.  */
#define HAVE_GETTIMEOFDAY 1

/* Define if you have the sysconf function.  */
#define HAVE_SYSCONF 1

/* Define if you have the uname function.  */
#define HAVE_UNAME 1

/* Define if you have the valloc function.  */
#define HAVE_VALLOC 1

/* Define if you have the <fcntl.h> header file.  */
#define HAVE_FCNTL_H 1

/* Define if you have the <fortran.h> header file.  */
/* #undef HAVE_FORTRAN_H */

/* Define if you have the <limits.h> header file.  */
#define HAVE_LIMITS_H 1

/* Define if you have the <malloc.h> header file.  */
#define HAVE_MALLOC_H 1

/* Define if you have the <netdb.h> header file.  */
#define HAVE_NETDB_H 1

/* Define if you have the <pwd.h> header file.  */
#define HAVE_PWD_H 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define if you have the <sys/param.h> header file.  */
/* #undef HAVE_PARAM_TIME_H */

/* Define if you have the <sys/time.h> header file.  */
#define HAVE_SYS_TIME_H 1

/* Define if you have the <sys/unistd.h> header file.  */
#define HAVE_SYS_UNISTD_H 1

/* Define if you have the <sys/utsname.h> header file.  */
#define HAVE_SYS_UTSNAME_H 1

/* Define the fortran library name mangling type */

/* #undef FORTRANDOUBLEUNDERSCORE */
#define FORTRANUNDERSCORE 1
/* #undef FORTRANCAPS */
/* #undef FORTRANNOUNDERSCORE */

/* Define native Fortran INTEGER and REAL byte size */

#define FORT_INTEGER_LEN 4
#define FORT_REAL_LEN 8

/*
 * Define corresponding Fortran C datatyp
 *  0 not supported
 *  1 char
 *  4 short
 *  6 int
 *  8 long      ! anyway, int and long should have the same size in byte
 * 10 float
 * 11 double
 * 13 long long
 */

#define FORT_INT 6
#define FORT_INT1 1
#define FORT_INT2 4
#define FORT_INT4 6
#define FORT_INT8 13
#define FORT_INT16 0

#define FORT_REAL 11
#define FORT_REAL4 10
#define FORT_REAL8 11
#define FORT_REAL16 0

#endif
