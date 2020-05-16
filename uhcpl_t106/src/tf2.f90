!OCL NOALIAS

SUBROUTINE tf2

  ! Description:
  !
  ! This subroutine completes the second part of the
  ! time filtering:xm1=xm1+eps*x.
  !
  ! Method:
  !
  ! *tf2* is called from *scans1l*.
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, May 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! I. Kirchner, MPI, October 1998, tendency diagnostics
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,       ONLY: eps, ltdiag, nlev, nrow
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_tracer,        ONLY: ntrac
  USE mo_scan_buffer,   ONLY: alps_scb, d_scb, t_scb, u_scb, v_scb, vo_scb
  USE mo_memory_gl,     ONLY: q, x, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xm1, xtm1
  USE mo_memory_g2a,    ONLY: um1, vm1
  USE mo_diag_tendency, ONLY: ptfh1, ptfh2, pdiga, pdsga, ldinit
  USE mo_decomposition, ONLY: dc => local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL    :: zeps
  INTEGER :: jk, jl, jt, jrow
  INTEGER :: nglon, nglpx, nglat


  !  Executable statements 

  ! local array bounds

  nglon = dc%nglon ! number of longitudes
  nglpx = dc%nglpx ! number of longitudes allocated
  nglat = dc%nglat ! number of latitudes

  ! local latitude index, continuous north -> south

  !  jrow = nrow(2)
  DO jrow = 1, nglat

     !-- 1. Second part of time filtering

     zeps = eps
     IF (nstep-nstart<=1) zeps = 0.

     DO jk = 1, nlev
!DIR$ IVDEP
        DO jl = 1, nglpx
           vom1(jl,jk,jrow) = vom1(jl,jk,jrow) + zeps*vo_scb(jl,jk,jrow)
           dm1 (jl,jk,jrow) = dm1 (jl,jk,jrow) + zeps* d_scb(jl,jk,jrow)
           qm1 (jl,jk,jrow) = qm1 (jl,jk,jrow) + zeps* q    (jl,jk,jrow)
           xm1 (jl,jk,jrow) = xm1 (jl,jk,jrow) + zeps* x    (jl,jk,jrow)
           tm1 (jl,jk,jrow) = tm1 (jl,jk,jrow) + zeps* t_scb(jl,jk,jrow)
           um1 (jl,jk,jrow) = um1 (jl,jk,jrow) + zeps* u_scb(jl,jk,jrow)
           vm1 (jl,jk,jrow) = vm1 (jl,jk,jrow) + zeps* v_scb(jl,jk,jrow)
        END DO
     END DO

     DO jt = 1, ntrac
        DO jk = 1, nlev
!DIR$ IVDEP
           DO jl = 1, nglon
              xtm1(jl,jk,jt,jrow) = xtm1(jl,jk,jt,jrow) + zeps*xt(jl,jk,jt,jrow)
           END DO
        END DO
     END DO

     DO jl = 1, nglon
        alpsm1(jl,jrow) = alpsm1(jl,jrow) + zeps*alps_scb(jl,jrow)
     END DO

     IF (ltdiag) THEN
        ! time filter correction for (t-1)
        IF (ldinit) THEN
           ptfh1(:,:,1,jrow) = 0.0
           ptfh1(:,:,2,jrow) = 0.0
           ptfh1(:,:,3,jrow) = 0.0
           ptfh2(:    ,jrow) = 0.0
        ELSE
           ptfh1(1:nglon,:,1,jrow) = vom1  (1:nglon,:,jrow) - ptfh1(1:nglon,:,1,jrow)
           ptfh1(1:nglon,:,2,jrow) = dm1   (1:nglon,:,jrow) - ptfh1(1:nglon,:,2,jrow)
           ptfh1(1:nglon,:,3,jrow) = tm1   (1:nglon,:,jrow) - ptfh1(1:nglon,:,3,jrow)
           ptfh2(1:nglon    ,jrow) = alpsm1(1:nglon,jrow)   - ptfh2(1:nglon    ,jrow)
        ENDIF
        ! accumulate time filter part
        pdiga(1:nglon,:,21,jrow) = pdiga(1:nglon,:,21,jrow) + ptfh1(1:nglon,:,1,jrow)
        pdiga(1:nglon,:,22,jrow) = pdiga(1:nglon,:,22,jrow) + ptfh1(1:nglon,:,2,jrow)
        pdiga(1:nglon,:,23,jrow) = pdiga(1:nglon,:,23,jrow) + ptfh1(1:nglon,:,3,jrow)
        pdsga(1:nglon,   2,jrow) = pdsga(1:nglon,   2,jrow) + ptfh2(1:nglon    ,jrow)
     ENDIF

  END DO

END SUBROUTINE tf2
