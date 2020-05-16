MODULE mo_aerosols

  IMPLICIT NONE

  !  ------------------------------------------
  ! 
  !  radiative characteristics of the aerosols
  ! 
  !  ------------------------------------------

  !    ----------------------
  ! 1. shortwave coefficients
  !    ----------------------

  !  normalized optical thickness at 0.55 micron
  REAL :: taua(2,5) = RESHAPE ( (/ .730719, .730719,    &
                                   .912819, .912819,    &
                                   .725059, .725059,    &
                                   .745405, .745405,    &
                                   .682188, .682188 /), &
                                   (/2,5/))
  !  single scattering albedo
  REAL :: piza(2,5) = RESHAPE ( (/ .872212, .872212,    &
                                   .982545, .982545,    &
                                   .623143, .623143,    &
                                   .944887, .944887,    &
                                   .997975, .997975 /), &
                                   (/2,5/))
  !  assymetry factor
  REAL :: cga(2,5) = RESHAPE ( (/  .647596, .647596,    &
                                   .739002, .739002,    &
                                   .580845, .580845,    &
                                   .662657, .662657,    &
                                   .624246, .624246 /), &
                                   (/2,5/))  

  !    ---------------------
  ! 2. longwave coefficients
  !    ---------------------

  !  absorption coefficients
  !  opical parameters for gads components, 8 rh classes
  REAL :: caer(5,5) = RESHAPE ( (/ &
       .038520, .037196, .040532, .054934, .038520,    &
       .12613 , .18313 , .10357 , .064106, .126130,    &
       .012579, .013649, .018652, .025181, .012579,    &
       .011890, .016142, .021105, .028908, .011890,    &
       .013792, .026810, .052203, .066338, .013792 /), &
       (/5,5/))

  REAL :: tauan(8,2,12)  !  extinction cross section in cm2/particle
  REAL :: pizan(8,2,12)  !  single scattering albedo
  REAL :: cgan(8,2,12)   !  assymetry factor
                         !    longwave:
  REAL :: caern(8,5,12)  !  extinction cross section in cm2/particle
  REAL :: fcvaer(12)     !  conversion factor in part/g

  INTEGER :: ndfaer(12)  !  contains the number of aerosol component, corresponding
                         !    with the gads number. element k=(1,...,12) corresponds
                         !    with the aerosol mixing ratio written in paer(...,5+k)
  INTEGER :: newaer      !  total number of additional aerosol components

  INTEGER :: ja, in

END MODULE mo_aerosols
