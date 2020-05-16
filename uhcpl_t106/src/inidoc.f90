!+ preset constants in mo_doctor.
!+ $Id: inidoc.f90,v 1.9 1999/09/20 09:10:20 m214089 Exp $

SUBROUTINE inidoc

  ! Description:
  !
  ! Preset constants in mo_doctor.
  !
  ! Method:
  !
  ! *inidoc* is called from *initialise*.
  !
  ! Authors:
  !
  ! J. K. Gibson, ECMWF, April 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_doctor
  USE mo_mpi, ONLY: p_pe, p_io, p_bcast

  IMPLICIT NONE

  !  Local scalars: 
  INTEGER :: nlena, nlenb, nlenc
  REAL :: zvers
  CHARACTER (3)   :: yovers
  CHARACTER (256) :: os_name, user_name, host_name
  CHARACTER (8)   :: ydate
  CHARACTER (10)  :: ytime

  !  External subroutines
  EXTERNAL util_os_system, util_user_name, util_node_name

  !  Executable statements 

  yovers = '4.0'
  READ (yovers,'(f8.1)') zvers
  nvers = zvers*10.0 + 10
  IF (p_pe == p_io) THEN
     CALL util_os_system (os_name, nlena)
     CALL util_user_name (user_name, nlenb)
     CALL util_node_name (host_name, nlenc)
  END IF
  CALL p_bcast (os_name, p_io)
  CALL p_Bcast (nlena, p_io)
  CALL p_bcast (user_name, p_io)
  CALL p_Bcast (nlenb, p_io)
  CALL p_bcast (host_name, p_io)
  CALL p_Bcast (nlenc, p_io)

  ylabel1(:) = ' '
  ylabel2(:) = ' '
  ylabel3(:) = ' '
  ylabel4(:) = ' '
  ylabel5(:) = ' '
  ylabel6(:) = ' '
  ylabel7(:) = ' '
  ylabel8(:) = ' '

  CALL date_and_time(ydate, ytime)

  ylabel1(1:) = ' Atmospheric model version ' // yovers
  ylabel2(1:) = ' Library 05/24/96'
  ylabel3(1:) = ' Spectral and semi Lagrangian dynamics'
  ylabel4(1:) = ' Modified ECMWF physics'
  ylabel5(1:) = ' Modified ECMWF radiation'
  ylabel6(1:) = ' Date - ' // ydate(1:8) // ' Time - ' // ytime(1:6)
  ylabel7(2:) = user_name(1:nlenb) // ' on ' // host_name(1:nlenc) 
  ylabel8(2:) = os_name(1:nlena)

END SUBROUTINE inidoc


