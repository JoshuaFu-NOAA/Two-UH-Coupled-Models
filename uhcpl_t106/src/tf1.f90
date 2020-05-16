!OCL NOALIAS

SUBROUTINE tf1

  ! Description:
  !
  ! It performs the first part of time filtering for next time step:
  ! xf=x+eps*(xm1-2*x).
  !
  ! Method:
  !
  ! *tf1* is called from *scan1sl* (subscan 2).
  !
  ! Reference:
  ! 1-appendix *b1:organisation of the spectral model.
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

  USE mo_control,       ONLY: eps, ltdiag, nrow
  USE mo_start_dataset, ONLY: nstart, nstep
  USE mo_tracer,        ONLY: ntrac
  USE mo_sc1,           ONLY: alps, d, t, u, v, vo
  USE mo_memory_gl,     ONLY: q, x, xt
  USE mo_memory_g1a,    ONLY: alpsm1, dm1, qm1, tm1, vom1, xm1, xtm1
  USE mo_memory_g1b,    ONLY: alpsf, df, qf, tf, vof, xf, xtf
  USE mo_memory_g2a,    ONLY: um1, vm1
  USE mo_memory_g2b,    ONLY: uf, vf
  USE mo_diag_tendency, ONLY: ptfh1, ptfh2
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL    :: z1m2eps, zeps
  INTEGER :: jk, jl, jt, jrow
  INTEGER :: nlev, nglon, nglpx


  !  Executable statements 

  nlev  = dc%nlev   ! global: number of levels
  nglon = dc%nglon  ! local:  number of longitudes
  nglpx = dc%nglpx  ! local:  number of longitudes allocated
  jrow  = nrow(2)   ! local:  latitude index, continuous

  !-- 1. First part of time filtering

  IF (nstep==nstart) THEN
     zeps = 0.
     z1m2eps = 1.
  ELSE
     zeps = eps
     z1m2eps = 1. - 2.*eps
  END IF

  DO jk = 1, nlev
!DIR$ IVDEP
     DO jl = 1, nglpx
        vof(jl,jk,jrow) = z1m2eps*vo(jl,jk) + zeps*vom1(jl,jk,jrow)
        df (jl,jk,jrow) = z1m2eps*d(jl,jk) + zeps*dm1(jl,jk,jrow)
        qf (jl,jk,jrow) = z1m2eps*q(jl,jk,jrow) + zeps*qm1(jl,jk,jrow)
        xf (jl,jk,jrow) = z1m2eps*x(jl,jk,jrow) + zeps*xm1(jl,jk,jrow)
        tf (jl,jk,jrow) = z1m2eps*t(jl,jk) + zeps*tm1(jl,jk,jrow)
        uf (jl,jk,jrow) = z1m2eps*u(jl,jk) + zeps*um1(jl,jk,jrow)
        vf (jl,jk,jrow) = z1m2eps*v(jl,jk) + zeps*vm1(jl,jk,jrow)
     END DO
  END DO

  DO jt = 1, ntrac
     DO jk = 1, nlev
!DIR$ IVDEP
        DO jl = 1, nglon
           xtf(jl,jk,jt,jrow) = z1m2eps*xt(jl,jk,jt,jrow) + zeps*xtm1(jl,jk,jt,jrow)
        END DO
     END DO
  END DO

!DIR$ IVDEP
  DO jl = 1, nglon
     alpsf(jl,jrow) = z1m2eps*alps(jl) + zeps*alpsm1(jl,jrow)
  END DO

  IF (ltdiag) THEN          ! store at (t) used at (t+1) in TF2
     ptfh1(1:nglon,:,1,jrow) = vo  (1:nglon,:)
     ptfh1(1:nglon,:,2,jrow) = d   (1:nglon,:)
     ptfh1(1:nglon,:,3,jrow) = t   (1:nglon,:)
     ptfh2(1:nglon    ,jrow) = alps(1:nglon)
  ENDIF

  RETURN
END SUBROUTINE tf1
