#include "util_fortran.h"

#include <string.h>

#if defined(FORTRANCAPS)
#define util_reshape_ UTIL_RESHAPE
#elif defined(FORTRANDOUBLEUNDERSCORE)
#define util_reshape_ util_reshape__
#elif !defined(FORTRANUNDERSCORE)
#define util_reshape_ util_reshape
#endif

FORTRAN_CALL
void util_reshape_(REAL *dest, REAL *source, INT *amount)
{
  size_t iamount;

  iamount = (size_t) *amount * sizeof(REAL);

  memcpy (dest, source, iamount);

  return;
}
