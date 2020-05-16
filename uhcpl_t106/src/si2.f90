!+ 2nd part of the semi-implicit scheme (done in fourier space).
!+ $Id: si2.f90,v 1.14 2000/06/05 12:35:25 m214003 Exp $
!OCL NOALIAS

SUBROUTINE si2

  ! Description:
  !
  ! 2nd part of the semi-implicit scheme (done in fourier space).
  !
  ! Method:
  !
  ! This subroutine computes the contribution in
  ! *Fourier space to the semi-implicit scheme (mainly for
  ! the vorticity and humidity equations).
  !
  ! *si2* is called from *fcc1*.
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, January 1982, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! I. Kirchner, MPI, December 1998, vorticity semi-implicit diagnosticsterm modified
  ! I. Kirchner, MPI, May 2000, revision of tendency diagnostics
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: lptime, ltdiag, twodt
  USE mo_gaussgrid,     ONLY: racst, twomu
  USE mo_truncation,    ONLY: am
  USE mo_semi_impl,     ONLY: betazq
  USE mo_constants,     ONLY: a
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_diag_tendency, ONLY: pdiga, pdigb, pdigs
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_buffer_fft,    ONLY: fdm1, fvol, fvom, fu0, fdu0, fvom1

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: z1, z2, z3, zbmi, zbmr, zcmi, zcmr, zdl, zdt, zdtma2, zdtmda2, &
&      zdu0, zrcst, zu0, zvoli, zvolr, zvomi, zvomr &
    , zdigsr, zdigsi, zdigsmr, zdigsmi, zdigslr, zdigsli

  INTEGER :: irow, jlev, jm, jrow, krow, jglat
  INTEGER :: nflev, nmp1, ngl, nflat

  !  Executable statements 

  ! loop bounds on this PE

  nflev = dc% nflev  ! local (Fourier space) : number of levels 
  nflat = dc% nflat  ! local (Fourier space) : number of latitudes  
  nmp1  = dc% nm+1   ! global                : number of coefficients m
  ngl   = dc% nlat   ! global                : number of latitudes

!-- 1. Locate and allocate storage

!-- 1.1 Fourier components

  ! latitude loop indices

  DO jrow = 1, nflat
! jrow  = nrow(2)                        ! local  continuous north to south
  jglat = dc% glat(jrow)                 ! global continuous north -> south
  irow  = MIN(2*jglat-1,2*(ngl+1-jglat)) ! global ping pong index
  krow  = (irow+1)/2                     ! half latitude index

!  ALLOCATE (dl(nlp2,nflev))

!-- 2. Skip over *si2* during initialisation iterations
!      or compute temporary quantities.

!-- 2.1 Compute temporary quantities

  zdt = twodt*.5
  IF (nstep==nstart) zdt = zdt*.5
  z1 = betazq*zdt*racst(krow)
  z2 = z1*a
  z3 = twomu(irow)*racst(krow)
  zrcst = racst(krow)*a

!-- 3. Semi implicit modifications

!-- 3.1 Initiate scans

  DO jlev = 1, nflev
    zu0  = z1*fu0(jlev,jrow)
    zdu0 = z2*(fdu0(jlev,jrow)+z3*fu0(jlev,jrow))

!DIR$ IVDEP
!OCL NOVREC
    DO jm = 1, nmp1
      zdtma2 = am(jm)*zu0
      zdtmda2 = am(jm)*zdu0

      zdl = am(jm)*zrcst

      zbmr = 1./(1.+zdtma2**2)
      zbmi = -zdtma2*zbmr

      zcmr = -2.*zdtmda2*    zdtma2    *zbmr**2
      zcmi =    -zdtmda2*(1.-zdtma2**2)*zbmr**2

!-- 3.2 Divergence equation

      fdm1(2*jm-1,jlev,jrow) =fdm1(2*jm-1,jlev,jrow)+zdl*fvom(2*jm,  jlev,jrow)
      fdm1(2*jm,  jlev,jrow) =fdm1(2*jm,  jlev,jrow)-zdl*fvom(2*jm-1,jlev,jrow)

!-- 3.3 Vorticity equation
      zvolr = zbmr*(fvom1(2*jm-1,jlev,jrow) - am(jm)*fvol(2*jm,  jlev,jrow))  &
&           - zbmi*(fvom1(2*jm,  jlev,jrow) + am(jm)*fvol(2*jm-1,jlev,jrow))  &
&           - zcmr* fvom (2*jm-1,jlev,jrow) &
&           + zcmi* fvom (2*jm,  jlev,jrow)

      zvoli = zbmi* (fvom1(2*jm-1,jlev,jrow) - am(jm)*fvol(2*jm,  jlev,jrow)) &
&           + zbmr* (fvom1(2*jm,  jlev,jrow) + am(jm)*fvol(2*jm-1,jlev,jrow)) &
&           - zcmi* fvom (2*jm-1,jlev,jrow) &
&           - zcmr* fvom (2*jm,  jlev,jrow)

      zvomr = zbmr*fvom(2*jm-1,jlev,jrow) - zbmi*fvom(2*jm,jlev,jrow)
      zvomi = zbmi*fvom(2*jm-1,jlev,jrow) + zbmr*fvom(2*jm,jlev,jrow)

      IF (ltdiag) THEN
         ! store L-Term before adjustment
         zdigslr = fvol(2*jm-1,jlev,jrow)
         zdigsli = fvol(2*jm  ,jlev,jrow)

         ! store M-Term before adjustment
         zdigsmr = fvom(2*jm-1,jlev,jrow)
         zdigsmi = fvom(2*jm  ,jlev,jrow)
      ENDIF

      fvol(2*jm-1,jlev,jrow) = zvolr
      fvol(2*jm,  jlev,jrow) = zvoli
      fvom(2*jm-1,jlev,jrow) = zvomr
      fvom(2*jm,  jlev,jrow) = zvomi

      IF (ltdiag) THEN
         ! compose explicit part, accounted in SI1
         zdigsr = - am(jm)*pdigs(2*jm  ,jlev,1,jrow)
         zdigsi =   am(jm)*pdigs(2*jm-1,jlev,1,jrow)

         ! adjustment of semi-implicit part of vorticity, L-Term
         !-------------------------------------
         ! I. Kirchner
         ! NOTE: there is a small hole in the accumulation of
         !       all right hand side terms, it can be
         !       related to a missing part of semi-implicit ter
         !       .. only for vorticity
         !       subject of investigation !!
         !-------------------------------------
         ! add implicit part

         pdigs(2*jm-1,jlev,1,jrow) = zdigsr + zdtma2*zvoli
         pdigs(2*jm,  jlev,1,jrow) = zdigsi - zdtma2*zvolr

         pdigs(2*jm-1,jlev,4,jrow) =   zdtma2*zvomi
         pdigs(2*jm,  jlev,4,jrow) = - zdtma2*zvomr
      ENDIF

    END DO

    IF (ltdiag .AND. lptime) THEN
      ! zonal derivatives with factor m/(1-mu**2)
      DO jm = 1,nmp1
        zdl = am(jm)*zrcst
        pdigb(2*jm-1,jlev,:,jrow) = - zdl * pdiga(2*jm,  jlev,1:10,jrow)
        pdigb(2*jm,  jlev,:,jrow) =   zdl * pdiga(2*jm-1,jlev,1:10,jrow)
      ENDDO
    END IF

  END DO
END DO

END SUBROUTINE si2
