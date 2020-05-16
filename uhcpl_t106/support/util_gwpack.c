/*
 * This file contains a set of routines to provide 8 byte integer support
 * for unpacking gravity drag information fields in a fast way.           
 * To restrict the f90 source code to f90 standard only this implementations 
 * are developed. This C routines interface to f90. A basic convention is to 
 * name all routines with a leading util_ to make this visible to the 
 * developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 * 
 * $Id: util_gwpack.c,v 1.4 1999/10/04 09:26:11 m214003 Exp $
 *
 */

#include "util_fortran.h"

#include <string.h>

#if defined(FORTRANCAPS)
#define util_gwpack_    UTIL_GWPACK
#define util_gwunpk_    UTIL_GWUNPK
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define util_gwpack_    util_gwpack__
#define util_gwunpk_    util_gwunpk__
#elif !defined(FORTRANUNDERSCORE)
#define util_gwpack_    util_gwpack
#define util_gwunpk_    util_gwunpk
#endif

/*
 * util_gwpack - packs four input fields into one output field
 * each input word is rounded to the nearest integer
 * and truncated on the right by 7 bits.
 * the bits are packed into a subfield of size 16 bits
 * therefore, all input numbers must be less than 8388480
 *
 * a, b, c, d are input arrays of length klen words
 * p is the packed output array of length klen words
 */

FORTRAN_CALL
void util_gwpack_ (REAL *a, REAL *b, REAL *c, REAL *d, REAL *p, INT *klen)
{
  ULONG ia, ib, ic, id, ip;

  const ULONG mask57 = 0xffffffffffffff80;

  const int itrunc = 7;
  const double rounder = 64.0; /* 2**(itrunc-1) */

  int jl;

  for (jl = 0; jl < *klen; jl++) {

    ia = (ULONG) (*(a+jl) + rounder);
    ib = (ULONG) (*(b+jl) + rounder);
    ic = (ULONG) (*(c+jl) + rounder);
    id = (ULONG) (*(d+jl) + rounder);

    ip =      ((ia & mask57) << (48-itrunc));
    ip = ip | ((ib & mask57) << (32-itrunc));
    ip = ip | ((ic & mask57) << (16-itrunc));
    ip = ip | (id >> itrunc);

    memcpy(p+jl, &ip, sizeof(ULONG));
  }

  return;
}

/*
 * util_gwunpk - unpacks one input array into four output arrays
 * each packed field is held in 16 bits
 * it is shifted by 7 bits and floated
 *
 * p is the packed input array of length klen words
 * a, b, c, d are gwunpked output arrays of length klen words
 */

FORTRAN_CALL
void util_gwunpk_ (REAL *a, REAL *b, REAL *c, REAL *d, REAL *p, INT *klen)
{
  ULONG ia, ib, ic, id, ip;
  ULONG imask;

  const ULONG mask16 = 0xffff000000000000;

  const int itrunc = 7;

  int jl;

  imask = mask16 >> (48-itrunc);

  for (jl = 0; jl < *klen; jl++) {

    memcpy(&ip, p+jl, sizeof(ULONG));

    ia= ip >> (48-itrunc) & imask;
    ib= ip >> (32-itrunc) & imask;
    ic= ip >> (16-itrunc) & imask;
    id= ip << itrunc & imask;

    *(a+jl) = (REAL) ia;
    *(b+jl) = (REAL) ib;
    *(c+jl) = (REAL) ic;
    *(d+jl) = (REAL) id;
  }

  return;
}

