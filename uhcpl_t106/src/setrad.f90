!+ preset and modify constants in radiation modules.
!+ $Id: setrad.f90,v 1.8 1999/08/17 12:58:06 m214003 Exp $

SUBROUTINE setrad

  ! Description:
  !
  ! Preset and modify constants in radiation modules
  !
  ! Method:
  !
  ! This subroutine preset and modify constants in modules
  ! used in the parameterisation of radiative processes.
  ! it also allocates,if requested,space for some diagnostics.
  !
  ! *setrad* is called from *inimod*.
  !
  ! Externals:
  ! *radclc*    compute critical relative humidity at each
  !             vertical level for radiation.
  ! *aerdis*    compute aerosol distributions.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_tmp_buffer
  USE mo_parameters
  USE mo_mpi
  USE mo_doctor
  USE mo_control
  USE mo_hyb
  USE mo_param_switches
  USE mo_rad_switches
  USE mo_rad2
  USE mo_aerosols
  USE mo_longwave

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: iradfr

  !  External subroutines 
  EXTERNAL aerdis, suaerx, suradi, susw

  !  Intrinsic functions 
  INTRINSIC MOD

  INCLUDE 'radctl.inc'


  !  Executable statements 

  !-- 1. Preset constants in *mo_rad_switches*

  lrad = .TRUE.
  lsolc = .TRUE.
  ldiur = .TRUE.
  nradfr = -2
  co2fac = 1.
  nradpfr = 126
  nradpla = ngl/3
  nmonth = 0
  laer = .TRUE.
  lgadsrh = .FALSE.
  lcfc = .FALSE.

  ! Preset aerosol definition array in yomaer

  DO ja = 1, 12
     ndfaer(ja) = 0
  END DO

  !-- 2. Read namelist

  IF (p_parallel) THEN
     IF (p_parallel_io) THEN
        READ (nin,radctl)
     ENDIF
     CALL p_bcast (lrad, p_io)
     CALL p_bcast (ldiur, p_io)
     CALL p_bcast (lsolc, p_io)
     CALL p_bcast (laer, p_io)
     CALL p_bcast (lcfc, p_io)
     CALL p_bcast (lgadsrh, p_io)
     CALL p_bcast (co2fac, p_io)
     CALL p_bcast (nmonth, p_io)
     CALL p_bcast (nradfr, p_io)
     CALL p_bcast (nradpfr, p_io)
     CALL p_bcast (nradpla, p_io)
     CALL p_bcast (ndfaer, p_io)
  ELSE
     READ (nin,radctl)
  ENDIF

  !-- 3. Modify constants in *mo_rad_switches*

  IF ( .NOT. lphys) THEN
     lrad = .FALSE.
     ldiur = .FALSE.
  END IF

  IF ( .NOT. lrad) THEN
     ldiur = .FALSE.
  END IF
  iradfr = nradfr
  IF (nradfr<0) nradfr = -nradfr*3600./dtime + 0.1
  nradpfr = nradpfr*nradfr
  IF (MOD(nradpla,2)==0 .AND. nradpla/=0) nradpla = nradpla + 1

  ! Set more radiation parameters

  CALL suradi
  CALL suaerx
  CALL susw

  !-- 4. Allocate space for diagnostics

  IF (lrad .AND. nradpfr/=0) THEN
     ALLOCATE (diag(nlevp1,4))
  END IF

  !-- 5. Compute aerosol distribution

  CALL aerdis(cetah,cvdaes,cvdael,cvdaeu,cvdaed,nlevp1,ctrbga,cvobga,cstbga, &
       &      caeops,caeopl,caeopu,caeopd,ctrpt,caeadk,caeadm)

  !-- 6. Write namelist

#ifdef DEBUG
  WRITE (nerr,radctl)
#endif

  !-- 7. Print message

  IF (p_pe == p_io) THEN
     IF (lrad) THEN
        IF (lsolc) THEN
           WRITE (nout,'(3(a,/))') &
                ' Diagnostic of solar cloud forcing is switched on,',        &
                '   to save 3-5% cpu time set lsolc = .false. in namelist ', &
                '   radctl if this diagnostic is not necessary.'
        ELSE
           WRITE (nout, '(a)') &
                ' Diagnostic of solar cloud forcing is switched off.' 
        END IF
     END IF
  END IF

END SUBROUTINE setrad
