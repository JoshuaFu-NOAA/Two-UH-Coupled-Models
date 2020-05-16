!*     NAMELIST *NMICTL* to control NMI
!
! Authors:
!   I. Kirchner, MPI, August 1998
!
! description of parameter see in *MO_NMI*
!
NAMELIST /nmictl/ &
&     ntpre, &  ! INTEGER number of timesteps skiped before accumulation of tendencies
&     ntdia, &  ! INTEGER number of accumulation time steps for diabatic tendencies
&     ntiter, & ! INTEGER number of iteration time steps
&     pcut, &   ! INTEGER cut off period for fast gravity modes in hours
&     pcutd     ! INTEGER cut off period for filtering of diabatic tendencies
!
! version 10-oct-1999
