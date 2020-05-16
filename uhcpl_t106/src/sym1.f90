!+ computes symmetric and antisymmetric parts of fourier components.
!+ $Id: sym1.f90,v 1.16 2000/02/10 12:59:03 m214003 Exp $
!OCL NOALIAS

SUBROUTINE sym1

  ! Description:
  !
  ! This subroutine computes symmetric and antisymmetric
  ! parts of *fourier components.
  !
  ! Method:
  !
  ! *sym1* is called from *fcc1*.
  !
  ! *sym1* computes symetric and antisymmetric parts of
  ! *Fourier components in 2 steps:The contribution of northern
  ! is added when *nrow* is odd and of southern hemisphere when
  ! *nrow* is even. For multi-tasking, real contributions for
  ! each of two rows are computed for even tasks, the imm-
  ! aginary contributions being computed in odd tasks.
  !
  ! Results:
  ! The contributions are returned as symmetric and antisymmetric
  ! *Fourier coefficients.
  !
  ! Reference:
  ! 1-appendix *b1:Organisation of the spectral model.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, March 1982, original source
  ! J. K. Gibson, ECMWF, August 1983, Adapted for multi-tasking
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_memory_f,      ONLY: fadl, fadm, far, fatp1, faul, fazl, fazm, &
                              fsdl, fsdm, fsr, fstp1, fsul, fszl, fszm
  USE mo_buffer_fft,    ONLY: ldl, ldm, ltm1, lrh, lvol, lvom ,lul, lalpsm1


  IMPLICIT NONE

  !  Local scalars: 

  INTEGER :: nlev, nlevp1, nmp1, nlnm0, ngl, nhgl
  INTEGER :: jm, jl, jg, jgr


  !  Executable statements 

  ! local array bounds on this PE

  nlev   = dc% nflev   ! number of levels
  nlevp1 = dc% nflevp1 ! number of levels including highest level
  nmp1   = dc% nlm     ! number of spectral coefficients m
  nlnm0  = dc% nlnm0   ! number of spectral coefficients n for m=0
  ngl    = dc% nlat    ! number of latitudes
  nhgl   = ngl / 2

  DO jg = 1, nhgl
    jgr = ngl + 1 - jg
    DO jl = 1, nlev
!DIR$ IVDEP
      DO jm = 1, nmp1

        fszl (jl,1,jm,jg) = (lvol(2*jm-1,jl,jg) + lvol(2*jm-1,jl,jgr))*.5
        fszl (jl,2,jm,jg) = (lvol(2*jm  ,jl,jg) + lvol(2*jm  ,jl,jgr))*.5
        fazm (jl,1,jm,jg) = (lvom(2*jm-1,jl,jg) - lvom(2*jm-1,jl,jgr))*.5
        fazm (jl,2,jm,jg) = (lvom(2*jm  ,jl,jg) - lvom(2*jm  ,jl,jgr))*.5
        fsdl (jl,1,jm,jg) = (ldl (2*jm-1,jl,jg) + ldl (2*jm-1,jl,jgr))*.5
        fsdl (jl,2,jm,jg) = (ldl (2*jm  ,jl,jg) + ldl (2*jm  ,jl,jgr))*.5
        fadm (jl,1,jm,jg) = (ldm (2*jm-1,jl,jg) - ldm (2*jm-1,jl,jgr))*.5
        fadm (jl,2,jm,jg) = (ldm (2*jm  ,jl,jg) - ldm (2*jm  ,jl,jgr))*.5
        fsr  (jl,1,jm,jg) = (lrh (2*jm-1,jl,jg) + lrh (2*jm-1,jl,jgr))*.5
        fsr  (jl,2,jm,jg) = (lrh (2*jm  ,jl,jg) + lrh (2*jm  ,jl,jgr))*.5
        fstp1(jl,1,jm,jg) = (ltm1(2*jm-1,jl,jg) + ltm1(2*jm-1,jl,jgr))*.5
        fstp1(jl,2,jm,jg) = (ltm1(2*jm  ,jl,jg) + ltm1(2*jm  ,jl,jgr))*.5
        fazl (jl,1,jm,jg) = (lvol(2*jm-1,jl,jg) - lvol(2*jm-1,jl,jgr))*.5
        fazl (jl,2,jm,jg) = (lvol(2*jm  ,jl,jg) - lvol(2*jm  ,jl,jgr))*.5
        fszm (jl,1,jm,jg) = (lvom(2*jm-1,jl,jg) + lvom(2*jm-1,jl,jgr))*.5
        fszm (jl,2,jm,jg) = (lvom(2*jm  ,jl,jg) + lvom(2*jm  ,jl,jgr))*.5
        fadl (jl,1,jm,jg) = (ldl (2*jm-1,jl,jg) - ldl (2*jm-1,jl,jgr))*.5
        fadl (jl,2,jm,jg) = (ldl (2*jm  ,jl,jg) - ldl (2*jm  ,jl,jgr))*.5
        fsdm (jl,1,jm,jg) = (ldm (2*jm-1,jl,jg) + ldm (2*jm-1,jl,jgr))*.5
        fsdm (jl,2,jm,jg) = (ldm (2*jm  ,jl,jg) + ldm (2*jm  ,jl,jgr))*.5
        far  (jl,1,jm,jg) = (lrh (2*jm-1,jl,jg) - lrh (2*jm-1,jl,jgr))*.5
        far  (jl,2,jm,jg) = (lrh (2*jm  ,jl,jg) - lrh (2*jm  ,jl,jgr))*.5
        fatp1(jl,1,jm,jg) = (ltm1(2*jm-1,jl,jg) - ltm1(2*jm-1,jl,jgr))*.5
        fatp1(jl,2,jm,jg) = (ltm1(2*jm  ,jl,jg) - ltm1(2*jm  ,jl,jgr))*.5

      END DO
    END DO
  END DO

  IF (nlevp1>nlev) THEN
    DO jm = 1, nmp1

      fstp1(nlevp1,1,jm,:)=(lalpsm1(2*jm-1,:nhgl)+lalpsm1(2*jm-1,ngl:nhgl+1:-1))*.5
      fstp1(nlevp1,2,jm,:)=(lalpsm1(2*jm  ,:nhgl)+lalpsm1(2*jm  ,ngl:nhgl+1:-1))*.5
      fatp1(nlevp1,1,jm,:)=(lalpsm1(2*jm-1,:nhgl)-lalpsm1(2*jm-1,ngl:nhgl+1:-1))*.5
      fatp1(nlevp1,2,jm,:)=(lalpsm1(2*jm  ,:nhgl)-lalpsm1(2*jm  ,ngl:nhgl+1:-1))*.5

    END DO
  END IF

  IF (nlnm0>0) fsul(:,:) = (lul(:,1:nhgl)+lul(:,ngl:nhgl+1:-1))*.5
  IF (nlnm0>0) faul(:,:) = (lul(:,1:nhgl)-lul(:,ngl:nhgl+1:-1))*.5

END SUBROUTINE sym1
