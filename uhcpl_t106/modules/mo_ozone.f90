MODULE mo_ozone

  ! readozo: read ozone distribution  (called by control) 
  !----------------------------------------------------------
  ! *** Version for the Fortuin ozone from Holland **** 
  ! --------------------------------------------------------    

  ! Note: mo_midatm.f90 contains the code to read ozone data

  ! jpozlev :  number of vertical levels of ozone distribution
  ! nuoz    :  unit number associated with ozone file

  INTEGER, PARAMETER :: jpozlev=19, nuoz=41

END MODULE mo_ozone
