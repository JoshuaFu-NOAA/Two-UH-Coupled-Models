!+ mass fixing of constituents and air mass, tendencies
!+ $Id: slt2.f90,v 1.15 2000/03/17 09:13:01 m214003 Exp $

#ifndef SLDIAG
SUBROUTINE slt2(fb,qfcst,alpha,psm1cor,pscor,etamid,kftype)
#else
SUBROUTINE slt2(fb,qfcst,alpha,psm1cor,pscor,etamid,kftype,hqfcst)
#endif

  ! Description:
  !
  ! Mass fixing of constituents and air mass, tendencies
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, March 1994, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! T. Diehl, DKRZ, July 1999, parallel version
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters,    ONLY: jps
  USE mo_control,       ONLY: twodt
  USE mo_start_dataset, ONLY: nstep, nstart
  USE mo_tracer,        ONLY: ntrac, lslt
  USE mo_grid,          ONLY: pcnst, plev, plond, platd, plat, plato2, js, &
                              plon, nxpt
  USE mo_scan_buffer,   ONLY: alps_scb, qe_scb, xe_scb, xte_scb
  USE mo_memory_g1a,    ONLY: alpsm1
  USE mo_decomposition, ONLY: dc => local_decomposition

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: pscor, psm1cor

  !  Array arguments 
  REAL :: alpha(pcnst), etamid(plev), fb(plond,plev,pcnst,platd,2), &
          qfcst(plond,plev,pcnst,platd,2)
  INTEGER :: kftype(pcnst)

  !  Local scalars: 
  REAL :: zlpsc, zlpsm1c, ztmst
  INTEGER :: ilat, jlat, ita, jcen, jk, jl, jt, ih

  !  Local arrays: 
#ifdef SLDIAG
  REAL :: qf3m(plon,plev,pcnst), vqfm(plon,plev,pcnst)
  REAL :: hqfcst(plond,plev,pcnst,plat)
#endif

  !  External subroutines 
  EXTERNAL fixer


  !  Executable statements 

  IF (nstep/=nstart) THEN
     ztmst = twodt
  ELSE
     ztmst = 0.5*twodt
  END IF

  DO jlat  = 1, dc% nglat  ! local  index north -> south
     ilat = plat + 1 - jlat  ! continuous local index for hqfcst
     IF (ilat > plato2) THEN
        ih = 2
        jcen = js - 1 + ilat - plato2
     ELSE
        ih = 1
        jcen = js - 1 + ilat
     END IF

     !-- 1. Fix mass of semi lagrangian transported scalars

#ifdef SLDIAG
     CALL fixer(ztmst,alpha,fb(1,1,1,jcen,ih),qfcst(1,1,1,jcen,ih), &
                etamid,kftype,hqfcst(1,1,1,ilat),vqfm,qf3m)
#else
     CALL fixer(ztmst,alpha,fb(1,1,1,jcen,ih),qfcst(1,1,1,jcen,ih), &
                etamid,kftype)
#endif

     !-- 2. Compute tendencies of semi lagrangian transport

     DO jk = 1, plev
!DIR$ IVDEP
!OCL NOVREC
        DO jl = 1, plon
           qe_scb(jl,jk,jlat)                                     &
                       = qe_scb(jl,jk,jlat) + (qfcst(jl+nxpt,jk,1,jcen,ih)-     &
                                                  fb(jl+nxpt,jk,1,jcen,ih))/ztmst
           xe_scb(jl,jk,jlat)                                     &
                       = xe_scb(jl,jk,jlat) + (qfcst(jl+nxpt,jk,2,jcen,ih)-     &
                                                  fb(jl+nxpt,jk,2,jcen,ih))/ztmst
        END DO
     END DO

     ita = 0
     DO jt = 1, ntrac
        IF (lslt(jt)) THEN
           ita = ita + 1
           DO jk = 1, plev
!DIR$ IVDEP
!OCL NOVREC
              DO jl = 1, plon
                 xte_scb(jl,jk,jt,jlat) =                         &
                            xte_scb(jl,jk,jt,jlat) +              &
                            (qfcst(jl+nxpt,jk,ita+jps,jcen,ih) -  &
                             fb(jl+nxpt,jk,ita+jps,jcen,ih))/ztmst
              END DO
           END DO
        END IF
     END DO

     !-- 3. Fix mass of air

     zlpsm1c = LOG(psm1cor)
     zlpsc = LOG(pscor)

     DO jl = 1, plon
        alpsm1  (jl,jlat) = alpsm1  (jl,jlat) + zlpsm1c
        alps_scb(jl,jlat) = alps_scb(jl,jlat) + zlpsc
     END DO
  END DO

END SUBROUTINE slt2
