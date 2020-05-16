/*
 * This file contains a set of routines to provide CRAY to IEEE binary and 
 * vice versa conversion for f90 on Unix machines. Unfortunatelly there is 
 * no implemented set of such functions in general available. To restrict
 * the f90 source code to f90 standard only this implementations are developed.
 * One basic convention is to name all routines with a leading util_ to make
 * this visible to the developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 *
 * $Id: util_convert.c,v 1.5 1999/10/04 09:26:10 m214003 Exp $
 *
 */

#include "util_fortran.h"

#include <stdio.h>
#include <string.h>

#if defined(FORTRANCAPS)
#define util_cray2ieee_ UTIL_CRAY2IEEE
#define util_ieee2cray_ UTIL_IEEE2CRAY
#define util_i4toi8_    UTIL_I4TOI8
#define util_i8toi4_    UTIL_I8TOI4
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define util_cray2ieee_ util_cray2ieee__
#define util_ieee2cray_ util_ieee2cray__
#define util_i4toi8_    util_i4toi8__
#define util_i8toi4_    util_i8toi4__
#elif !defined(FORTRANUNDERSCORE)
#define util_cray2ieee_ util_cray2ieee
#define util_ieee2cray_ util_ieee2cray
#define util_i4toi8_    util_i4toi8
#define util_i8toi4_    util_i8toi4
#endif

/*
 *  Convert Floating point, Cray to IEEE 64-bit
 *
 *
 *	input : cf	Cray Floating point numbers
 *		nf	Number of elements in cf
 *	output: ieeef	IEEE Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */

FORTRAN_CALL
void util_cray2ieee_ (ULONG *crayf, ULONG *ieeef, INT *nf)
{

  const ULONG cray_expbias = 0x4000;   /* Cray exponent bias */
  const ULONG ieee_expbias = 0x3ff;    /* IEEE exponent bias */

  const ULONG mask1 = 0x8000000000000000;
  const ULONG mask2 = 0x7fff;
  const ULONG mask3 = 0xfffffffffffff;
  const ULONG mask4 = 0x7ff;
  const ULONG indef = 0xfff8000000000000;   /* IEEE indefinite, 64 bit form */

  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < *nf; i++) {
    *(ieeef+i) = (*(crayf+i) & mask1)
       |
        ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) << 52)
       |
    	((*(crayf+i) << 5) & mask3);

  }
  /*
   * handle 0.0, underflow and overflow :
   *
   * if the exponent field goes negative then the Cray number was
   * either 0.0 or too small to represent in IEEE, in either case
   * set the IEEE result to all 0's which is the IEEE representation for 0.0.
   *
   * if the exponent field is too large for the IEEE field, set
   * the result to the indefinite value.
   */

  /* Cray's internal conversion seems to skip this underflow part and instead
     of indef sets 0 .... */

  for (i = 0; i < *nf; i++) {
    if ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) < 0) {
      *(ieeef+i) = 0;
    } else if ((((*(crayf+i) >> 48) & mask2)-cray_expbias+ieee_expbias-1) 
	       > mask4) {
/*    *(ieeef+i) = indef;  */ 
      *(ieeef+i) = 0;
    }
  }
  

  return;
}

/*
 *  Convert Floating point, IEEE 64-bit, to Cray Floating point
 *
 *	input : ieeef	IEEE Floating point numbers (double precision)
 *		nf	Number of elements in ieeef
 *	output: crayf	Cray Floating point numbers
 *
 *  Format :
 *            sign  exponent  mantissa	unused
 *  IEEE :     1      11        52
 *  Cray :     1      15        48
 */
FORTRAN_CALL
void util_ieee2cray_ (ULONG *ieeef, ULONG *crayf, INT *nf)
{
  const ULONG cray_expbias = 0x4000;      /* Cray exponent bias */
  const ULONG ieee_expbias = 0x3ff;       /* IEEE exponent bias */
  const ULONG implied = 0x10000000000000; /* implied bit in IEEE mantissa */

  const ULONG mask1 = 0x8000000000000000;
  const ULONG mask3 = 0xfffffffffffff;
  const ULONG mask4 = 0x7ff;

  int i;

  /* Set sign bit, exponent and mantissa in one vector loop */

  for (i = 0; i < *nf; i++) {
    if (*(ieeef+i) == 0) {
      *(crayf+i) = 0;
    } else {
      *(crayf+i) = (*(ieeef+i) & mask1)
	 |
 	  ((((*(ieeef+i) >> 52) & mask4)-ieee_expbias+cray_expbias+1) << 48)
 	 |
 	  (((*(ieeef+i) & mask3)+implied) >> 5);
    }
  }

  return;
}

FORTRAN_CALL
void util_i8toi4_ (LONG *int64, INT *int32, INT *ni)
{
#if SIZEOF_INT == 8
  memcpy (int32, int64, (size_t) *ni*8);
#else
  int i;

  for (i = 0; i < *ni; i++) {
#ifdef DEBUG
    if ( (*(int64+i) > INT_MAX) ||  (*(int64+i) < INT_MIN)) {
      fprintf (stderr, 
	       "Integer out of range. Converted number contains nonsens\n");
    }
#endif
    /*
    *(int32+i) = (INT) *(int64+i);
    */
    *(int32+i) = (int) *(int64+i);
  }
#endif
  return;
}

FORTRAN_CALL
void util_i4toi8_ (INT *int32, ULONG *int64, INT *ni)
{
#if SIZEOF_INT == 8
  memcpy (int64, int32, (size_t) *ni*8);
#else
  const ULONG mask1 = 0xffffffff00000000;
  const ULONG mask2 = 0x0;
  int i;

  for (i = 0; i < *ni; i++) {
    if ( *(int32+i) < 0 ) {
      *(int64+i) =  *(int32+i) | mask1;
    } else {
      *(int64+i) =  *(int32+i) | mask2;
    }
  }
#endif
  return;
}

