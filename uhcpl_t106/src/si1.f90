!OCL NOALIAS

SUBROUTINE si1

  ! Description:
  !
  ! 1st part of the semi-implicit scheme (done in grid point space). 
  !
  ! Method:
  !
  ! This subroutine computes the contribution in
  ! grid points to the semi-implicit scheme.
  !
  ! *si1* is called from *gpc*.
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
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_sc1,           ONLY: alps, alpse, d, qe, rh, t, te, u0, vo, vol, &
                              vom, xe, xte
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xm1, xtm1
  USE mo_tmp_buffer,    ONLY: dm
  USE mo_gaussgrid,     ONLY: racst
  USE mo_control,       ONLY: lptime, ltdiag, nlev, nrow, twodt
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_constants,     ONLY: a
  USE mo_semi_impl,     ONLY: betazq
  USE mo_tracer,        ONLY: ntrac
  USE mo_diag_tendency, ONLY: pdiga, pdsga, pdigs, pdsgs
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: z1, z1u0, z2, z3, zdt, ztwodt
  INTEGER :: irow, jlev, jlon, jrow, krow, jt
  INTEGER :: nglon, nglpx
  LOGICAL :: loperm

  !  Local arrays: 
  REAL, TARGET  :: zd(dc% nglpx, nlev)
  REAL          :: zp(dc% nglpx), zt(dc% nglpx, nlev)
  REAL, POINTER :: zr(:,:)

  !  External subroutines 

  EXTERNAL conteq, pgrad


  !  Executable statements 

  !-- 1. Locate and allocate storage

  irow  = nrow(1)    ! global ping pong index
  krow  = (irow+1)/2 ! global half latitude index
  jrow  = nrow(2)    ! local  continuous index
  nglon = dc% nglon  ! local  number of latitudes
  nglpx = dc% nglpx  ! local  number of latitudes allocated

  ALLOCATE (dm(nglpx,nlev))

  !   zp(:)   = 0.0
  !   zt(:,:) = 0.0
  !   zd(:,:) = 0.0

  !-- 1.2 Equivalence arrays *zd* and *zr*

  zr => zd

  !-- 2. Skip over *si1* during initialisation iterations
  !      or prepare some temporary constants.

  !-- 2.1 Compute temporary constants

  zdt = twodt*.5
  IF (nstep==nstart) zdt = zdt*.5
  ztwodt = 2.*zdt
  z1 = betazq*zdt*racst(krow)
  z2 = ztwodt/a
  z3 = ztwodt*racst(krow)

  IF (ltdiag .AND. lptime) THEN
     ! correction of pdiga only at accumulation time steps
     ! rescaling with 
     !   pdiga(*,1...10,*) with 2*dt/a
     !   pdiga(*,11..23,*) with 2*dt
     pdiga(:,:, 1:10,jrow) = pdiga(:,:, 1:10,jrow)*z2
     pdiga(:,:,11:20,jrow) = pdiga(:,:,11:20,jrow)*ztwodt
     pdiga(:,:,24:25,jrow) = pdiga(:,:,24:25,jrow)*ztwodt
     pdsga(:,      1,jrow) = pdsga(:,      1,jrow)*ztwodt
  ENDIF

  !-- 3. Semi implicit computations

  !-- 3.1 Vorticity and humidity equations
  !       and term dmu of divergence equation.
  DO jlev = 1, nlev

     z1u0 = z1*u0(jlev)

     ! get explicit part of time schema, (2*dt) included
     ! explicit part of vorticity
     IF (ltdiag) pdigs(1:nglon,jlev,1,jrow) = &
          -z1u0*(vom1(1:nglon,jlev,jrow)-2.*vo(1:nglon,jlev))

!DIR$ IVDEP
     DO jlon = 1, nglon

        dm(jlon,jlev)  =  z2  * vol (jlon,jlev)

        vol(jlon,jlev) =   z3 * vol (jlon,jlev) &
                          -z1u0*(vom1(jlon,jlev,jrow) - 2.*vo(jlon,jlev))
        vom(jlon,jlev) = -z2  * vom (jlon,jlev)

        qm1(jlon,jlev,jrow) = qm1(jlon,jlev,jrow) + ztwodt*qe(jlon,jlev)
        xm1(jlon,jlev,jrow) = xm1(jlon,jlev,jrow) + ztwodt*xe(jlon,jlev)

     END DO
  END DO

  DO jt = 1, ntrac
     DO jlev = 1, nlev
        DO jlon = 1, nglon
           xtm1(jlon,jlev,jt,jrow) = xtm1(jlon,jlev,jt,jrow) + ztwodt*xte(jlon,jlev,jt)
        END DO
     END DO
  END DO

  !-- 3.2 Compute implicit contribution of divergence
  !       to temperature and surface equations.

  DO jlev = 1, nlev
     DO jlon = 1, nglon
        zd(jlon,jlev) = .5*dm1(jlon,jlev,jrow) - d(jlon,jlev)
     END DO
  END DO
  loperm = .FALSE.
  CALL conteq(zt,nglpx,zp,zd,nglpx,nglon,loperm)

  IF (ltdiag) THEN
     pdigs(:,:,3,jrow) = 2.*zt(:,:)   ! explicit part of temperature
     pdsgs(:    ,jrow) = 2.*zp(:)     ! explicit part of pressure
  ENDIF

  !-- 3.3 Update *zt* and *zp* to compute the contribution
  !       of temperature and surface pressure to the
  !       divergence equation.

  DO jlev = 1, nlev
     DO jlon = 1, nglon
        zt(jlon,jlev) = zt(jlon,jlev) + tm1(jlon,jlev,jrow) - t(jlon,jlev) + &
                        zdt*te(jlon,jlev)
     END DO
  END DO
  DO jlon = 1, nglon
     zp(jlon) = zp(jlon) + alpsm1(jlon,jrow) - alps(jlon) + zdt*alpse(jlon)
  END DO
  CALL pgrad(zr,zt,nglpx,zp,nglon)

  ! explicit part of divergence
  IF (ltdiag) pdigs(:,:,2,jrow) = -zr(:,:)*ztwodt

  !-- 3.4 Complete computation of the terms to be
  !       passed to *fftd*.

  DO jlev = 1, nlev
     DO jlon = 1, nglon
        tm1(jlon,jlev,jrow) = 2.*zt(jlon,jlev) - tm1(jlon,jlev,jrow) + 2.*t(jlon,jlev)
        rh(jlon,jlev) = -ztwodt*(rh(jlon,jlev)+zr(jlon,jlev))
     END DO
  END DO

  DO jlon = 1, nglon
     alpsm1(jlon,jrow) = 2.*zp(jlon) - alpsm1(jlon,jrow) + 2.*alps(jlon)
  END DO

  RETURN
END SUBROUTINE si1
