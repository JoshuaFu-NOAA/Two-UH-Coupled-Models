!+ resets accumulated g4x-fileds
!+ $Id: setzerox.f90,v 1.7 1999/08/17 09:13:03 m214003 Exp $

SUBROUTINE setzerox

  ! Description:
  !
  ! All accumulated g4x-fields are reset to default values
  ! this routine is called from the incore version of the model.
  !
  ! Method:
  !
  ! Locate every latitude row in the model's wokspace
  ! then set all values to zero
  !
  ! *setzerox* is called from *scan1sl*.
  !
  ! Results:
  ! The accumulated fields are reset
  !
  ! Authors:
  !
  ! U. Schlese, DKRZ, May 1993, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_parameters
  USE mo_start_dataset
  USE mo_control
  USE mo_io_tables
  USE mo_post

  IMPLICIT NONE

  !  Local scalars: 
!  INTEGER :: i1dim, i1xdim, ilen, irow, jl, jx


  !  Executable statements 

!-- 1. Process extra g4-fields

!!$  IF (lg4x) THEN
!!$
!!$    ! Call locate(ig4x,'g4x01   ',igptype,-1)
!!$
!!$    ! I1XDIM=IG4X
!!$    DO jx = 1, ng4xp
!!$      ! I1DIM=I1XDIM-NPTRB(NG4BBK)+NPTGLO(NG1A)+NDISP4
!!$      ilen = nlon*ng4xl(jx)
!!$      IF (l4xaccu(jx)) THEN
!!$        DO irow = 1, ngl
!!$          ! Z1DIM(JL)=0.
!!$          DO jl = 1, ilen
!!$          END DO
!!$!          i1dim = i1dim + hdntab(ng1a) %length*jmbyte
!!$        END DO
!!$      END IF
!!$!      i1xdim = i1xdim + ilen*jmbyte
!!$    END DO
!!$
!!$  END IF

  RETURN
END SUBROUTINE setzerox
