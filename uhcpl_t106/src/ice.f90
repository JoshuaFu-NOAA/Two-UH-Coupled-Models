!+ computes seaice cover and depth for uncoupled runs

SUBROUTINE atmice

  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, July 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_g3a,    ONLY: auxil2m, slmm, tsm, tsm1m
  USE mo_memory_g3b,    ONLY: seaice, siced
  USE mo_tmp_buffer,    ONLY: ahfli
  USE mo_control,       ONLY: nrow, lamip2
  USE mo_physc2,        ONLY: ctfreez
  USE mo_decomposition, ONLY: dc=>local_decomposition

  IMPLICIT NONE

  !  Local scalars: 
  REAL    :: zalpha, zdice
  INTEGER :: irow, jrow, jl
  INTEGER :: nglon, nglpx
  LOGICAL :: lo, lo2

  !  Local arrays: 
  REAL, POINTER :: tskin1m(:)


  !  Executable statements 

  ! computational constant
  ! diffusion coefficient zalpha

  zalpha = 2.0

!-- 1. Locate and position space

  nglon = dc% nglon ! local number of longitudes
  nglpx = dc% nglpx ! local number of longitudes allocated
  irow  = nrow(1)   ! global ping pong latitude index
  jrow  = nrow(2)   ! local continuous latitude index

  tskin1m => auxil2m(:,jrow)

  ALLOCATE (ahfli(nglpx))

!-- 2. Compute seaice

  lo2 = ((irow/2)*2 == irow) ! true in southern hemisphere

  DO jl = 1, nglon
    lo = (tsm1m(jl,jrow) <= ctfreez) .AND. (slmm(jl,jrow) < 0.5)
    IF (lo) THEN
      seaice(jl,jrow)  = 1.
      IF (lo2) THEN
        siced(jl,jrow) = 1.
      ELSE
        IF (lamip2) THEN
          siced(jl,jrow) = 1.5
        ELSE
          siced(jl,jrow) = 2.
        END IF
      END IF
    ELSE
      seaice(jl,jrow) = 0.
      siced(jl,jrow)  = 0.
    END IF
  END DO

!-- 3. Set skin temperature

  zdice = 0.10
  DO jl = 1, nglon
    ahfli(jl) = 0.
    lo = (siced(jl,jrow) > zdice)
    IF (lo) THEN
      tsm1m(jl,jrow) = tskin1m(jl)
      tsm(jl,jrow)   = tsm1m(jl,jrow)
      ahfli(jl)      = zalpha*(tsm(jl,jrow)-ctfreez)/siced(jl,jrow)
    END IF
  END DO

  RETURN
END SUBROUTINE atmice
