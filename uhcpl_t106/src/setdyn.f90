!+ preset and modify constants in dynamics,initialisation and general
!  purposes modules.
!+ $Id: setdyn.f90,v 1.18 1999/11/10 08:21:44 m214003 Exp $

SUBROUTINE setdyn

  ! Description:
  !
  ! Preset and modify constants in dynamics,initialisation
  ! and general purposes modules.
  !
  ! Method:
  !
  ! *setdyn* is called from *initialise*.
  !
  ! Externals.:
  !
  ! *inicon*   preset constants in *mo_constants*.
  ! *inictl*   preset constants in *mo_control*.
  ! *inifft*   preset constants in *mo_fft*.
  ! *inigau*   preset constants in *mo_gaussgrid*.
  ! *inihyb*   preset constants in *mo_hyb*.
  ! *inhysi*   modify constants in *mo_hyb*.
  ! *helmo*    compute matrix for the *helmoltz equation.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, December 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, November 1998, nudging
  ! L. Kornblueh, MPI, June 1999, parallel version (MPI based)
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_exception
  USE mo_tmp_buffer
  USE mo_parameters
  USE mo_mpi
  USE mo_doctor
  USE mo_control
  USE mo_truncation,       ONLY: ntrk, ntrm, ntrn, scpar
  USE mo_semi_impl
  USE mo_hdiff
  USE mo_hyb
  USE mo_param_switches
  USE mo_forecast_switches
  USE mo_constants
  USE mo_start_dataset
  USE mo_midatm,            ONLY: damhih, spdrag, enspodi, nlvspd1, nlvspd2
  USE mo_fft,               ONLY: inifft
  USE mo_gaussgrid,         ONLY: inigau ! module subroutine

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zaa, zp0, zt, zt0
  INTEGER :: idt, jk, jl

  !  Local arrays: 
  REAL :: zp0a(1), zpf(nlev), zph(nlevp1)

  !  External subroutines 

  EXTERNAL helmo, inhysi, pres, presf, sudif

  INCLUDE 'dynctl.inc'

  !  Executable statements 

  !-- 1. Preset constants in modules

  !-- 1.1 Preset constants in *mo_constants*

  !! call now in INICTL
  !!  CALL inicon

  !-- 1.2 Preset constants in *mo_forecast_switches*

  ndiadfr = -120
  ndiavfr = 0
  lsimdt = .TRUE.
  lsimzq = .TRUE.
  lvtmpc1 = .TRUE.
  lvtmpc2 = .TRUE.
  lumax = .FALSE.
  lzondia = .FALSE.
  ldiahdf = .FALSE.

  !-- 1.3 Preset constants in *mo_fft*

  CALL inifft(nlon)

  !-- 1.4 Preset constants in *mo_gaussgrid*

  CALL inigau

  !-- 1.5 Preset constants in *mo_hdiff*

  IF (nn==21) THEN
     IF (lmidatm) THEN
        dampth = 192.
        damhih = 1000.
        enstdif = 1.
        nlvstd1 = 1
        nlvstd2 = 1
        spdrag  = 0.0
        enspodi = 1.
        nlvspd1 = 1
        nlvspd2 = 1
     ELSE
        dampth = 6.
        enstdif = 0.1
        nlvstd1 = 1
        nlvstd2 = 1
     ENDIF
  ELSE IF (nn==30) THEN
     dampth = 12.
     IF (lmidatm) THEN
        damhih = 100.
        enstdif = 1.
        nlvstd1 = 1
        nlvstd2 = 1
        spdrag  = 0.0
        enspodi = 1.
        nlvspd1 = 1
        nlvspd2 = 1
     ELSE
        enstdif = 0.4
        nlvstd1 = 2
        nlvstd2 = 2
     ENDIF
  ELSE IF (nn==106) THEN
     dampth = 3.
     IF (lmidatm) THEN
        damhih = 100.
        enstdif = 1.
        nlvstd1 = 1
        nlvstd2 = 1
        spdrag  = 0.0
        enspodi = 1.
        nlvspd1 = 1
        nlvspd2 = 1
     ELSE
        enstdif = 1.25
        nlvstd1 = 2
        nlvstd2 = 2
     ENDIF
  ELSE
     dampth = 9.
     IF (lmidatm) THEN
        damhih = 100.
        enstdif = 1.
        nlvstd1 = 1
        nlvstd2 = 1
        spdrag  = 0.0
        enspodi = 1.
        nlvspd1 = 1
        nlvspd2 = 1
     ELSE
        enstdif = 0.6
        nlvstd1 = 2
        nlvstd2 = 2
     ENDIF
  END IF

  IF (nlev>=19 .OR. lgwdrag) THEN
     ldrag = .FALSE.
  ELSE
     ldrag = .TRUE.
  END IF
  cdrag = 14.*dayl

  !-- 1.5.1 Preset constants in *mo_diff*

  CALL sudif

  !-- 1.6 Preset constants in *mo_truncation*

  CALL scpar (nm, nn, nk)

  IF (nlev/=16 .AND. nlev/=19) THEN
     CALL message('setdyn','')
     CALL message('setdyn','')
     CALL message('setdyn','')
     CALL message('setdyn','check setting of stratospheric diffusion')
     CALL message('setdyn','')
     CALL message('setdyn','')
     CALL message('setdyn','')
  END IF

  IF (nm==nn .AND. nm==nk) THEN
     IF (nm>=106 .AND. nlev>=19) THEN
        ntrn(1) = nm - 24
        ntrn(2) = nm - 22
        ntrn(3) = nm - 20
        ntrn(4) = nm - 18
        ntrn(5) = nm - 16
        ntrn(6) = nm - 13
        ntrn(7) = nm - 10
        ntrn(8) = nm - 6
        ntrn(9) = nm - 3
        ntrn(10) = nm - 1
     ELSE IF (nm>=106 .AND. nlev>=16) THEN
        ntrn(1) = nm - 22
        ntrn(2) = nm - 18
        ntrn(3) = nm - 14
        ntrn(4) = nm - 10
        ntrn(5) = nm - 6
        ntrn(6) = nm - 3
        ntrn(7) = nm - 1
     ELSE IF (nm==63 .AND. nlev>=19) THEN
        ntrn(1) = 56
        ntrn(2) = 56
        ntrn(3) = 57
        ntrn(4) = 57
        ntrn(5) = 58
        ntrn(6) = 58
        ntrn(7) = 59
        ntrn(8) = 60
        ntrn(9) = 61
        ntrn(10) = 62
     ELSE IF (nm==63 .AND. nlev>=16) THEN
        ntrn(1) = 56
        ntrn(2) = 57
        ntrn(3) = 58
        ntrn(4) = 59
        ntrn(5) = 60
        ntrn(6) = 61
        ntrn(7) = 62
     END IF

  END IF

  !-- 1.7 Preset constants in *mo_semi_impl*

  apr = 80000.
  tr = 300.
  IF (nn==106) THEN
     vcrit = 68.
  ELSE
     vcrit = 85.
  END IF
  hdamp = 1.
  vcheck = 200.

  !-- 1.8 Preset constants in *mo_hyb*

  CALL inihyb

  !-- 1.9 Compute icao-based correction profile for
  !       temperature diffusion

  zp0 = 101320.
  zp0a(1) = zp0
  zt0 = 288.
  zaa = 1./5.256

  CALL pres(zph,1,zp0a,1)
  CALL presf(zpf,1,zph,1)

  DO jk = 1, nlev
     zt = zt0*((zpf(jk)/zp0)**zaa)
     IF (zt>216.5) THEN
        diftcor(jk) = 0.5*(vct(nvclev+jk)+vct(nvclev+jk+1))*zaa*zt*zp0/zpf(jk)
     ELSE
        diftcor(jk) = 0.
     END IF
  END DO

  !-- 2. Read namelists

  IF (p_parallel) THEN
     IF (p_parallel_io) THEN
        READ (nin,dynctl)
     ENDIF
     CALL p_bcast (ndiadfr, p_io)
     CALL p_bcast (ndiavfr, p_io)
     CALL p_bcast (ntrn, p_io)
     CALL p_bcast (nlvstd1, p_io)
     CALL p_bcast (nlvstd2, p_io)
     CALL p_bcast (lumax, p_io)
     CALL p_bcast (lzondia, p_io)
     CALL p_bcast (ldrag, p_io)
     CALL p_bcast (ldiahdf, p_io)
     CALL p_bcast (vcrit, p_io)
     CALL p_bcast (hdamp, p_io)
     CALL p_bcast (enstdif, p_io)
     CALL p_bcast (apsurf, p_io)
     CALL p_bcast (vcheck, p_io)
     CALL p_bcast (eps, p_io)
     CALL p_bcast (dampth, p_io)
     CALL p_bcast (cdrag, p_io)
     CALL p_bcast (damhih, p_io)
     CALL p_bcast (spdrag, p_io)
     CALL p_bcast (enspodi, p_io)
     CALL p_bcast (nlvspd1, p_io)
     CALL p_bcast (nlvspd2, p_io)
  ELSE
     READ (nin,dynctl)
  ENDIF

  !-- 3. Modify and derive some constants

  !-- 3.1 Modify constants in *mo_constants*

  IF ( .NOT. lvtmpc1) vtmpc1 = 0.
  IF ( .NOT. lvtmpc2) vtmpc2 = 0.

  !-- 3.2 Modify constants in *mo_truncation*

  IF (nm==nn .AND. nm==nk) THEN

     DO jl = 1, nlev
        ntrm(jl) = ntrn(jl)
        ntrk(jl) = ntrn(jl)
     END DO

  ELSE

     DO jl = 1, nlev
        ntrn(jl) = nn
     END DO

  END IF

  !-- 3.3 Modify constants in *mo_forecast_switches*

  IF (ndiadfr<0) ndiadfr = -ndiadfr*3600./dtime + 0.1
  IF (ndiavfr<0) ndiavfr = -ndiavfr*3600./dtime + 0.1

  !-- 3.4 Modify constants in *mo_semi_impl*

  IF (lmidatm) THEN
     ! vcrit is defined as the courant velocity: a/(nn*dtime)
     ! multiplied by the truncation nn. it is used in hdiff.
     vcrit = a/dtime
  ELSE
     ! The following statement is based on the empirical
     ! constatation that a 20 minutes time step was always safe
     ! for t63 if the wind does not exceed 85(vcrit)m/s.
     vcrit = vcrit*1200.*63./dtime
  ENDIF

  IF (p_pe == p_io) THEN
     WRITE (nout, '(/,a,e12.3,/,a,e12.3,/)') &
          ' Damping factor for strong stratospheric damping      (hdamp) = ', &
          hdamp, &
          ' Critical velocity for horizontal diffusion enhancing (vcrit) = ', &
          vcrit
  END IF


  IF (lsimdt) THEN
     betadt = 0.75
  ELSE
     betadt = 0.
  END IF

  IF (lsimzq) THEN
     betazq = 1.
  ELSE
     betazq = 0.
  END IF

  !-- 3.5 Modify constants in *mo_hyb*

  CALL inhysi

  !-- 3.6 Compute matrix for the *helmoltz equation

  ALLOCATE (cn(nlev,nlev,nkp1))
  idt = 1
  IF (nstep>nstart) idt = 2
  CALL helmo(idt)

  !-- 3.7 Modify constants in dynctl

  IF (cdrag<0) cdrag = -cdrag*dayl

  !-- 4. Write namelist values

#ifdef DEBUG
  WRITE (nerr,dynctl)
#endif

END SUBROUTINE setdyn
