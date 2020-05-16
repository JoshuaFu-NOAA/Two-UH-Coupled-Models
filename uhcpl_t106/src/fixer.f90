!+ computes  1) time tendency due to vertical advection
!            2) time tendency due to mass adjustment of the constituents
!+ $Id: fixer.f90,v 1.5 1999/08/27 17:00:21 m214030 Exp $

#ifndef SLDIAG
SUBROUTINE fixer(ztodt,alpha,qin,qout,etamid,kftype)
#else
SUBROUTINE fixer(ztodt,alpha,qin,qout,etamid,kftype,hqfm,vqfm,qf3m)
#endif

  ! Description:
  !
  ! Computes  1) time tendency due to vertical advection.
  !           2) time tendency due to mass adjustment of the constituents.
  !
  ! Method:
  !
  ! Also, the constituent fields are updated based upon the fixer (mass
  ! adjustment) tendency.
  !
  ! qfix=alpha*F*q2*|q2 - q1|**Beta.
  ! Two options are available:
  !  1. kftype=1 : F=1.  and Beta=1.5
  !  2. kftype=2 : F=eta and Beta=1.
  !
  ! .On input
  ! lat     latitude index (values run from 1 to nlat; from southern-most
  !         to northern-most latitude).
  ! ztodt   length of time step (in seconds)
  ! alpha   array of fixer coefficients (one for each constituent) used
  !         in computing the mass adjustment tendency
  ! qin     constituent fields from the previous time step
  ! hqfm    horizontal time tendency  (optional)
  ! qout    set of FORECASTED constituent fields
  ! etamid  vert. coord. at full levels
  ! kftype  type of mass fixer
  !
  ! .On return.
  ! vqfm    vertical time tendency  (optional)
  ! qfm3    fixer time tendency     (optional)
  ! qout    set of FORECASTED AND FIXED constituent fields
  !
  ! .Local parameters: required.
  ! plev    Number of levels in global grid.
  ! plon    Number of longitudes in global grid.
  ! plond   number of longitudes in extended grid.
  ! pcnst   Number of constituents

! USE mo_slt
  USE mo_grid, ONLY: pcnst, plev, plond, istart, istop

  IMPLICIT NONE

  !  Scalar arguments 
  REAL :: ztodt

  !  Array arguments 
  REAL :: alpha(pcnst), etamid(plev), &
&      qin(plond,plev,pcnst), qout(plond,plev,pcnst) 
#ifdef SLDIAG
  REAL :: qf3m(plon,plev,pcnst), hqfm(plond,plev,pcnst), vqfm(plon,plev,pcnst)
#endif


  INTEGER :: kftype(pcnst)

  !  Local scalars: 
  INTEGER :: i, jc, k

  !  Intrinsic functions 
  INTRINSIC ABS, MAX, SQRT


  !  Executable statements 
#ifdef SLDIAG
  DO jc = 1, pcnst
     DO k = 1, plev
        DO i = istart, istop
           vqfm(i-nxpt,k,jc) = (qout(i,k,jc)-qin(i,k,jc))/ztodt - &
&                              hqfm(i-nxpt,k,jc)
        END DO
     END DO
  END DO
#endif
  DO jc = 1, pcnst

    IF (kftype(jc)==1) THEN

      DO k = 1, plev
        DO i = istart, istop
#ifdef SLDIAG
          qf3m(i-nxpt,k,jc) = alpha(jc)*qout(i,k,jc)*(SQRT(ABS(qout(i,k,jc)- &
&                             qin(i,k,jc))))**3/ztodt
          qout(i,k,jc) = qout(i,k,jc) + ztodt*qf3m(i-nxpt,k,jc)
#else
          qout(i,k,jc) = qout(i,k,jc) + alpha(jc)*qout(i,k,jc)*(SQRT(ABS(qout &
&                        (i,k,jc)-qin(i,k,jc))))**3
#endif
        END DO
      END DO

    ELSE IF (kftype(jc)==2) THEN

      DO k = 1, plev
        DO i = istart, istop
#ifdef SLDIAG
          qf3m(i-nxpt,k,jc) = alpha(jc)*qout(i,k,jc)*etamid(k)* &
&              ABS(qout(i,k,jc)-qin(i,k,jc))/ztodt
          qout(i,k,jc) = qout(i,k,jc) + ztodt*qf3m(i-nxpt,k,jc)
#else
          qout(i,k,jc) = qout(i,k,jc) + alpha(jc)*qout(i,k,jc)*etamid(k)*ABS( &
&              qout(i,k,jc)-qin(i,k,jc))
#endif
        END DO
      END DO
    END IF

    ! Fill up negative values

    DO k = 1, plev
      DO i = istart, istop
        qout(i,k,jc) = MAX(qout(i,k,jc),0.)
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE fixer
