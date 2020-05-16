MODULE mo_timeint

  IMPLICIT NONE

  REAL    :: wgt1, wgt2   ! weighting factors for interpolation in time
  INTEGER :: nmw1, nmw2   ! month indices

CONTAINS

  SUBROUTINE timeintold

    ! *TIMEINTOLD* - Calculates weighting factores for CLSST, CLALAI and CLAVGRAT
    ! TIMEINT use nstep+1 to calculate zday
    !

    USE mo_parameters
    USE mo_constants
    USE mo_control
    USE mo_start_dataset
    USE mo_year

    REAL :: zday, zdayl, zdfrac, zmonthl, zmonthlp1, zmonthlm1, &
         &         zmohlf, zmohlfp1, zmohlfm1, zmomid, zsec
    INTEGER :: id, iday, idayl, im, imlen, imlenp1, imlenm1, imm1, imp1, isec, iy
    INTEGER :: kmm1, kmp1

    LOGICAL :: lday

    INTRINSIC NINT


    ! Set calendar related parameters

    zdayl = dayl
    idayl = NINT(dayl)

    ! Determine century day

    zday   = ncbase+(ntbase+dtime*(nstep))/dayl
    iday   = INT(zday+1.E-12)
    zdfrac = zday-iday

    ! Determine date

    CALL cd2dat(iday,id,im,iy)

    ! Determine length of months and position within current month

    kmp1 = im + 1
    kmm1 = im - 1

    imp1 = im + 1
    imm1 = im - 1

    IF (imp1 > 12) imp1 = 1
    IF (imm1 <  1) imm1 = 12

    IF (ly365) THEN
      imlen   = im2day(im,iy)
      imlenp1 = im2day(imp1,iy)
      imlenm1 = im2day(imm1,iy)
    ELSE
      imlen   = 30
      imlenp1 = 30
      imlenm1 = 30
    END IF

    zmonthl   = imlen*zdayl
    zmonthlp1 = imlenp1*zdayl
    zmonthlm1 = imlenm1*zdayl

    zmohlf   = zmonthl*0.5
    zmohlfp1 = zmonthlp1*0.5
    zmohlfm1 = zmonthlm1*0.5
    zmomid   = zmohlf

    zsec = ((id-1)+zdfrac)*zdayl
    isec = NINT(zsec)

    ! Set logical switch for global statistics output in postatd

    lday = .FALSE.
    IF (nstep /= nstart) THEN
       lday = (MOD(isec,idayl) == 0)
    END IF

    ! Weighting factors for first/second half of month

    nmw1 = im
    IF (zsec <= zmomid) THEN
       wgt1 = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
       wgt2 = 1.-wgt1
       nmw2 = kmm1
    ELSE
       wgt2 = (zsec-zmohlf)/(zmohlf+zmohlfp1)
       wgt1 = 1.-wgt2
       nmw2 = kmp1
    END IF

    RETURN
  END SUBROUTINE timeintold

  SUBROUTINE timeint

    ! *TIMEINT* - Calculates weighting factores for CLSST2 and OZONE
    !
    ! U. SCHLESE  -  SEP-98

    USE mo_parameters
    USE mo_constants
    USE mo_control
    USE mo_start_dataset
    USE mo_year

    REAL :: zday, zdayl, zdfrac, zmonthl, zmonthlp1, zmonthlm1, &
         &         zmohlf, zmohlfp1, zmohlfm1, zmomid, zsec
    INTEGER :: id, iday, idayl, im, imlen, imlenp1, imlenm1, imm1, imp1, isec, iy
    INTEGER :: kmm1, kmp1

    LOGICAL :: lday

    INTRINSIC NINT


    ! Set calendar related parameters

    zdayl = dayl
    idayl = NINT(dayl)

    ! Determine century day

    zday   = ncbase+(ntbase+dtime*(nstep+1))/dayl
    iday   = INT(zday+1.E-12)
    zdfrac = zday-iday

    ! Determine date

    CALL cd2dat(iday,id,im,iy)

    ! Determine length of months and position within current month

    kmp1 = im + 1
    kmm1 = im - 1

    imp1 = im + 1
    imm1 = im - 1

    IF (imp1 > 12) imp1 = 1
    IF (imm1 <  1) imm1 = 12

    IF (ly365) THEN
      imlen   = im2day(im,iy)
      imlenp1 = im2day(imp1,iy)
      imlenm1 = im2day(imm1,iy)
    ELSE
      imlen   = 30
      imlenp1 = 30
      imlenm1 = 30
    END IF

    zmonthl   = imlen*zdayl
    zmonthlp1 = imlenp1*zdayl
    zmonthlm1 = imlenm1*zdayl

    zmohlf   = zmonthl*0.5
    zmohlfp1 = zmonthlp1*0.5
    zmohlfm1 = zmonthlm1*0.5
    zmomid   = zmohlf

    zsec = ((id-1)+zdfrac)*zdayl
    isec = NINT(zsec)

    ! Set logical switch for global statistics output in postatd

    lday = .FALSE.
    IF (nstep /= nstart) THEN
       lday = (MOD(isec,idayl) == 0)
    END IF

    ! Weighting factors for first/second half of month

    nmw1 = im
    IF (zsec <= zmomid) THEN
       wgt1 = (zmohlfm1+zsec)/(zmohlfm1+zmohlf)
       wgt2 = 1.-wgt1
       nmw2 = kmm1
    ELSE
       wgt2 = (zsec-zmohlf)/(zmohlf+zmohlfp1)
       wgt1 = 1.-wgt2
       nmw2 = kmp1
    END IF

    RETURN
  END SUBROUTINE timeint

END MODULE mo_timeint
