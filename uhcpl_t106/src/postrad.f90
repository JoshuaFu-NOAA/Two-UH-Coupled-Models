!+ reallocate space for the matrix for the helmoltz equation and recompute
!  it after a radiation step. it also complets and prints diagnostics
!  for radiation.
!+ $Id: postrad.f90,v 1.5 1998/10/28 12:32:02 m214003 Exp $

SUBROUTINE postrad

  ! Description:
  !
  ! Reallocate space for the matrix for the helmoltz equation and recompute
  ! it after a radiation step. It also complets and prints diagnostics
  ! for radiation.
  !
  ! Method:
  !
  ! *postrad* is called from *nnsc1*.
  !
  ! Externals:
  ! *helmo*   compute matrix needed to invert *helmoltz equation.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, June 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tmp_buffer
  USE mo_parameters
  USE mo_control
  USE mo_doctor
  USE mo_rad_switches
  USE mo_start_dataset

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zalbpr, zdift1, zdift2, zdift3, zdift4, zdift5
  INTEGER :: jk
  LOGICAL :: lodia

  !  Local arrays: 
  REAL, ALLOCATABLE :: zdiat(:)

  !  Intrinsic functions 
  INTRINSIC MOD


  !  Executable statements 

  IF ( .NOT. lrad) RETURN
  IF (MOD(nstep,nradfr)/=0) RETURN

!-- 1. Complete diagnostics

  IF (nradpfr/=0) THEN
    lodia = MOD(nstep,nradpfr) == 0
  ELSE
    lodia = .FALSE.
  END IF
  IF (lodia) THEN
    ALLOCATE (zdiat(nlevp1))
    DO jk = 1, nlevp1
      zdiat(jk) = diag(jk,3) + diag(jk,4)
    END DO
    zdift1 = dia1(1) - dia1(2)
    zdift2 = dia1(4) - dia1(5)
    zdift3 = dia1(6) - dia1(7)
    zdift4 = zdift1 - dia1(3)
    zdift5 = zdift2 - zdift3
    IF (dia1(1)>0.01) THEN
      zalbpr = dia1(2)/dia1(1)*100.
    ELSE
      zalbpr = 0.
    END IF
    WRITE (nout, '(/,a,/,a,/)') &
&        ' Global Radiation horizontally averaged:', &
&        '   (Fluxes - T=Top, B=Bottom, D=Down, U=Up, S=Solar, L=Thermal, N=Net)'
    WRITE (nout, '(a)') ' Temperature  [K]'
    WRITE (nout, '(10f6.1)') diag(1:nlev,1)
    WRITE (nout, '(a)') ' Cloud Cover     '
    WRITE (nout, '(10f6.1)') diag(1:nlev,2)
    WRITE (nout, '(a)') ' Shortwave heating rates [K/day]  '
    WRITE (nout, '(10f6.1)') diag(1:nlev,3)
    WRITE (nout, '(a)') ' Longwave heating rates [K/day]  '
    WRITE (nout, '(10f6.1)') diag(1:nlev,4)
    WRITE (nout, '(a)') ' Net heating rates [K/day]  ' 
    WRITE (nout, '(10f6.1)') zdiat(1:nlev)

    WRITE (nout, '(a)') ' Fluxes:'
    WRITE (nout, '(3(a,f6.1))') &
&        '   TSD = ', dia1(1), ' TSU = ', dia1(2), ' TLU = ', dia1(3)
    WRITE (nout, '(4(a,f6.1))') &
&        '   BSD = ', dia1(4),  ' BSU = ', dia1(5), &
&        ' BLU = ', dia1(6), ' BLD = ', dia1(7)
    WRITE (nout, '(5(a,f6.1))') '   TS  = ', zdift1, &
&      ' BS  = ',zdift2 ,' BL  = ', zdift3, ' TND = ', zdift4,' BND = ', zdift5
    WRITE (nout, '(a,f5.1,a)') '   Planetary Albedo = ', zalbpr, ' %'
    WRITE (nout, '(2/)')

    DEALLOCATE (zdiat)

  END IF

  RETURN
END SUBROUTINE postrad
