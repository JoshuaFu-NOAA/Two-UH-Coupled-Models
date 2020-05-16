!+ initialize module radint controlling radint

SUBROUTINE suradi

  ! Description:
  !
  ! Initialize radint, the module that controls the
  ! radiation interface
  !
  ! Method:
  !
  ! Reference:
  ! Ecmwf research department documentation of the
  ! "in core model"
  !
  ! Authors:
  !
  ! J.-J. Morcrette, ECMWF, December 1988, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_control,     ONLY: lamip2
  USE mo_radint
  USE mo_rad_switches

  IMPLICIT NONE

  !  Local scalars: 
  REAL :: zairmwg, zch4mwg, zco2mwg, zn2omwg


  !  Executable statements 

!-- 1. Set default values

  zsnowal = 0.01

  calbsea = 0.07
  cemiss  = 0.996
  calbsno = 0.80

  !  Concentration of various trace gases (ipcc/sacc values for 1990)
  !    CO2     CH4       N2O     CFC11    CFC12    HCFC113  CFC114
  !  353ppmv 1.72ppmv  310ppbv  280pptv  484pptv   60pptv   15pptv
  !  CFC115   HCFC22   HCFC123  HCFC124   HFC125   HFC134A HCFC141B
  !    5pptv  122pptv     ?        ?        ?         ?       ?
  ! HCFC142B  HFC143A  HFC152A   CCL4     H3CCL3
  !    ?        ?         ?     146pptv     ?

  zairmwg = 28.970
  zco2mwg = 44.011
  zch4mwg = 16.043
  zn2omwg = 44.013

  IF (lamip2) THEN
    !
    !     CO2       CH4       N2O   
    !  348 ppmv  1.65 ppmv  306 ppbv
    !
    ccardi   = 348.E-06*zco2mwg/zairmwg*co2fac
    zmetha   = 1.65E-06*zch4mwg/zairmwg
    znitox   = 306.E-09*zn2omwg/zairmwg
  ELSE
    ccardi   = 353.E-06*zco2mwg/zairmwg*co2fac
    zmetha   = 1.72E-06*zch4mwg/zairmwg
    znitox   = 310.E-09*zn2omwg/zairmwg
  END IF

  zcfc(1)  = 280.E-12
  zcfc(2)  = 484.E-12
  zcfc(3)  =  60.E-12
  zcfc(4)  =  15.E-12
  zcfc(5)  =   5.E-12
  zcfc(6)  = 122.E-12
  zcfc(7)  =   0.
  zcfc(8)  =   0.
  zcfc(9)  =   0.
  zcfc(10) =   0.
  zcfc(11) =   0.
  zcfc(12) =   0.
  zcfc(13) =   0.
  zcfc(14) =   0.
  zcfc(15) = 146.E-12
  zcfc(16) =   0.

  zepsec = 1.E-12
  zepclc = 1.E-12
  zeph2o = 1.E-12
  zepalb = 1.E-12

  RETURN
END SUBROUTINE suradi
