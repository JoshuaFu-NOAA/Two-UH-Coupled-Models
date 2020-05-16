!+ fills 3dim arrays for the semi lagrangian scheme
!+ $Id: sltfill.f90,v 1.19 2000/03/17 09:13:01 m214003 Exp $

!OCL NOALIAS

SUBROUTINE sltfill(ub,vb,fb,pdel)

  ! Description:
  !
  ! Fills 3dim arrays for the semi lagrangian scheme
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, February 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_decomposition, ONLY: dc => local_decomposition
  USE mo_control,       ONLY: ngl
  USE mo_gaussgrid,     ONLY: sqcst
  USE mo_tracer,        ONLY: ntrac, lslt
  USE mo_grid,          ONLY: plond, plev, plat, pcnst, platd, plon, plato2, &
                              js, nxpt, plevp1
  USE mo_scan_buffer,   ONLY: u_scb, v_scb
  USE mo_memory_g1a,    ONLY: qm1, xm1, xtm1, alpsm1

  IMPLICIT NONE

  !  Array arguments 
  REAL :: pdel(plond,plev,plat)             ! pressure differences
  REAL :: fb(plond,plev,pcnst,platd,2)      ! constituent fields
  REAL :: ub(plond,plev,platd,2)            ! u-component of wind
  REAL :: vb(plond,plev,platd,2)            ! v-component of wind

  !  Local scalars: 
  REAL :: zlimit, zrcst
  INTEGER :: i, ita, j, jk, jl, jt, k, l, jglat, igprow, ih

  !  Local arrays: 
  REAL :: zphm1(plon,plevp1,plat), zpsm1(plon,plat)

  !  External subroutines 
  EXTERNAL pres

  !  Intrinsic functions 
  INTRINSIC ABS


  !  Executable statements 

  ! Set scalar fields to zero if the absolute value
  ! is smaller than a hardware dependent value
  ! to prevent slt scheme from underflow.

#if (defined CRAY) && (! defined _CRAYMPP)
  zlimit = 1.E-2200
#else
  zlimit = 1.E-200
#endif

  WHERE (ABS(qm1(:,:,:)) < zlimit) 
     qm1(:,:,:) = 0.0
  END WHERE

  WHERE (ABS(xm1(:,:,:)) < zlimit)
     xm1(:,:,:) = 0.0
  END WHERE

  WHERE (ABS(xtm1(:,:,:,:)) < zlimit) 
     xtm1(:,:,:,:) = 0.0
  END WHERE

  zpsm1(:,:) = EXP(alpsm1(1:plon,:))

  ! Compute half level pressure at t-dt for pdel.

  DO i = 1, plat
     CALL pres(zphm1(:,:,i), plon, zpsm1(:,i), plon)
  END DO

  ! Fill up 3d arrays for slt
  ! and divide winds by cos(latitude)

  DO i = 1, plat
     IF (MOD(i,2)==1) THEN
        j = plat-i/2
     ELSE
        j = i/2
     END IF

     IF (i <= plato2) THEN
        ih = 1
        k = i + js - 1
     ELSE
        ih = 2
        k = i + js - 1 - plato2
     END IF

     !     l = MOD((i+1),2)*(plat-(i-1)/2)+MOD(i,2)*(i+1)/2
     l=plat-i+1

     ! the latitudes in sqcst are interleaved
     jglat  = dc%glat(l)   ! global index (continuous from north to south)
     igprow = MIN(2*jglat-1,2*(ngl+1-jglat))    ! global ping pong index
     zrcst  = 1.0/sqcst(igprow)                 ! 1./cos(latitude)

     ! arrays for SLT must have the format south -> north
     DO jk = 1, plev
        DO jl = 1, plon
           pdel(jl+nxpt,jk,i)    = zphm1(jl,jk+1,l)-zphm1(jl,jk,l)
           ub(jl+nxpt,jk,k,ih)   = u_scb(jl,jk,l)*zrcst
           vb(jl+nxpt,jk,k,ih)   = v_scb(jl,jk,l)*zrcst
           fb(jl+nxpt,jk,1,k,ih) = qm1(jl,jk,l)
           fb(jl+nxpt,jk,2,k,ih) = xm1(jl,jk,l)
        END DO
     END DO

     ita = 0
     DO jt = 1, ntrac
        IF (lslt(jt)) THEN
           ita = ita + 1
           DO jk = 1, plev
              DO jl = 1, plon
                 fb(jl+nxpt,jk,2+ita,k,ih) = xtm1(jl,jk,jt,l)
              END DO
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE sltfill
