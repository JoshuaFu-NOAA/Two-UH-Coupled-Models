MODULE mo_kind

  ! $Id: mo_kind.f90,v 1.8 1999/09/27 12:29:01 m214089 Exp $

  IMPLICIT NONE

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
#ifdef SX
#ifdef CRAY
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
#else
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(15)
#endif
#else
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(18)
#endif

#ifdef CP4
#define CP_LEN 4
#else
#define CP_LEN 8
#endif

#if (CP_LEN == 4)
  INTEGER, PARAMETER :: cp = i4
#elif (CP_LEN == 8)
  INTEGER, PARAMETER :: cp = i8
#else
  #error "This memory layout is not supported!"
#endif


END MODULE mo_kind
