!+ computes the shortwave transmission functions
!+ $Id: swtt.f90,v 1.5 1999/12/10 10:11:02 m214003 Exp $

SUBROUTINE swtt(kdlon,knu,ka,pu,ptr)

  ! Description:
  !
  ! This routine computes the transmission functions for all the
  ! absorbers (H2O, uniformly mixed gases, and O3) in the two spectral
  ! intervals.
  !
  ! Method:
  !
  ! *swtt* is called from *sw*.
  !
  ! Explicit arguments:
  ! knu    :                     ; index of the spectral interval
  ! ka     :                     ; index of the absorber
  ! pu     : (kdlon)             ; absorber amount
  ! ==== outputs ===
  ! ptr    : (kdlon)             ; transmission function
  !
  ! Transmission function are computed using pade approximants
  ! and horner's algorithm.
  !
  ! Reference:
  ! See radiation's part of the model's documentation and
  ! ECMWF research department documentation of the "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, December 1988, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_shortwave

  IMPLICIT NONE

  !  Scalar arguments 
  INTEGER :: ka, kdlon, knu

  !  Array arguments 
  REAL :: ptr(kdlon), pu(kdlon)

  !  Local scalars: 
  REAL :: zr1, zr2
  INTEGER :: jl


  !  Executable statements 

!-- 1. Horner's algorithm to compute transmission function

  DO jl = 1, kdlon
    zr1 = apad(knu,ka,1) + pu(jl)*(apad(knu,ka,2) + pu(jl)*(apad(knu,ka,3) &
                         + pu(jl)*(apad(knu,ka,4) + pu(jl)*(apad(knu,ka,5) &
                         + pu(jl)*(apad(knu,ka,6) + pu(jl)*(apad(knu,ka,7)))))))

    zr2 = bpad(knu,ka,1) + pu(jl)*(bpad(knu,ka,2) + pu(jl)*(bpad(knu,ka,3) &
                         + pu(jl)*(bpad(knu,ka,4) + pu(jl)*(bpad(knu,ka,5) &
                         + pu(jl)*(bpad(knu,ka,6) + pu(jl)*(bpad(knu,ka,7)))))))

!-- 2. Add the background transmission

    ptr(jl) = (zr1/zr2)*(1.-d(knu,ka)) + d(knu,ka)
  END DO

  RETURN
END SUBROUTINE swtt
