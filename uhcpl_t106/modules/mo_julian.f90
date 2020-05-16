MODULE mo_julian

  IMPLICIT NONE

  ! Routines to handle Julian Dates and Times
  !
  ! References:
  !
  !   Montenbruck, Oliver, "Practical Ephmeris Calculations", Ch. 2, pp 33-40.
  !   The 1992 Astronomical Almanac, page B4.
  !
  ! The Julian Date is defined as the the number of days which have elapsed
  ! since the 1st of January of the year 4713 BC 12:00 Universal Time.
  !
  ! Up to 4th October 1582 AD the Julian Calendar was in force with leap
  ! years every four years, but from then on the Gregorian Calendar carried 
  ! on from 15th October 1582 with the leap years defined by the rule:
  !
  !  "Leap year is every year whose yearly number is divisible by four, but
  !   not by a hundred, or is divisible by four hundred."
  !
  ! At midday on 4th October 1582, 2,299,160.5 Julian days had elapsed.
  ! The Gregorian Calendar then carried on at this point from 15th October
  ! 1582, so its beginning occured on the Julian date 2,299,160.5.
  ! 
  ! Note: the astronomical year -4712 corresponds to the year 4713 BC, the
  !       year 0 to the year 1 BC; thereafter the astronomical year match
  !       the year AD, e.g. year 1 = year 1 AD.
  !
  !       This routines work for the years -5877402 BC until 5868098 AD. This 
  !       dates are neveetheless coverable by the current GRIB edition 1. Which
  !       can cover dates between 1 AD and 25599 AD.  
  !
  !       The "Modified Julian Date" is the Julian Date - 2400000.5 and so 
  !       has a zero point of 17th November 1858 AD 00:00 Universal Time.
  ! 
  !       Due to the small area coverable by the GRIB output there is no need
  !       to use a "Modified Julian Date".
  !
  ! The Julian day number is stored as two doubles to guarantee sufficient
  ! precision on all machines. The complete value is (day + fraction)
  ! although this addition will sometimes lose precision. Note that day
  ! might not be an integer value and fraction might be greater than one.

  ! This is the clean version !!!!!!!

  !INTEGER, PRIVATE :: days_in_month(0:12) = (/ &
!&      0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! This is the version compatible with the old code ...

  INTEGER, PRIVATE :: days_in_month(0:12) = (/ &
&      365, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)



  TYPE julian_date
     REAL :: day
     REAL :: fraction
  END TYPE julian_date

CONTAINS

  FUNCTION SetYMD (ky, km, kd, zfraction) RESULT (julian_day)

    TYPE (julian_date) :: julian_day
    
    REAL, OPTIONAL :: zfraction

    INTEGER, INTENT(IN) :: ky, km, kd

    INTEGER :: ib, iy, im

    REAL :: zf

    IF (.NOT. PRESENT(zfraction)) THEN
       julian_day%fraction = 0.0
    ELSE 
       julian_day%fraction = zfraction
    ENDIF

    IF (km <= 2) THEN
        iy = ky-1
        im = km+12
     ELSE 
        iy = ky
        im = km
    ENDIF

    IF (ky > 1582 .OR. (ky == 1582 .AND. km > 10 &
&        .OR. (km == 10 .AND. kd >= 15))) THEN

       ! 15th October 1582 AD or later

       ib = INT(iy/400)-INT(iy/100)
    ELSE 

       ! 4th October 1582 AD or earlier

       ib = -2
    ENDIF

    julian_day%day = FLOOR(365.25*iy)+INT(30.6001*(im+1))+ib+1720996.5+kd

    zf = julian_day%day-AINT(julian_day%day)+julian_day%fraction
    julian_day%day = AINT(julian_day%day)
    if (zf >= 1.0) THEN
       julian_day%fraction = zf-AINT(zf)
       julian_day%day = julian_day%day+AINT(zf)
    ELSE
       julian_day%fraction = zf-AINT(zf)
    ENDIF 

  END FUNCTION SetYMD

  SUBROUTINE YMD (julian_day, ky, km, kd, zfraction)
    
    TYPE (julian_date), INTENT(IN) :: julian_day
    INTEGER, INTENT(OUT) :: km, kd, ky
    REAL, INTENT(OUT) :: zfraction
    
    REAL :: za, zb, zc, zd, ze, zf

    za = FLOOR(julian_day%day+julian_day%fraction+0.5)

    IF (za < 2299161.0) THEN
       zc = za+1524.0
    ELSE 
       zb = FLOOR((za-1867216.25)/36524.25)
       zc = za+zb-FLOOR(zb/4)+1525
    ENDIF

    zd = FLOOR((zc-122.1)/365.25)
    ze = FLOOR(365.25*zd)
    zf = FLOOR((zc-ze)/30.6001)

    kd = INT(zc-ze- FLOOR(30.6001*zf))
    km = INT(zf-1-12*FLOOR((zf+0.0001)/14))
    ky = INT(zd-4715-((7+km)/10))
    zfraction = (julian_day%day+0.5-za)+julian_day%fraction;

  END SUBROUTINE YMD

  SUBROUTINE YD (julian_day, ky, kd, zfraction)
    
    TYPE (julian_date), INTENT(IN) :: julian_day
    INTEGER, INTENT(OUT) :: kd, ky
    REAL, INTENT(OUT) :: zfraction
    
    INTEGER :: i, im
    REAL :: za, zb, zc, zd, ze, zf

    za = FLOOR(julian_day%day+julian_day%fraction+0.5)

    IF (za < 2299161.0) THEN
       zc = za+1524.0
    ELSE 
       zb = FLOOR((za-1867216.25)/36524.25)
       zc = za+zb-FLOOR(zb/4)+1525
    ENDIF

    zd = FLOOR((zc-122.1)/365.25)
    ze = FLOOR(365.25*zd)
    zf = FLOOR((zc-ze)/30.6001)

    kd = INT(zc-ze- FLOOR(30.6001*zf))
    im = INT(zf-1-12*FLOOR((zf+0.0001)/14))
    ky = INT(zd-4715-((7+im)/10))
    zfraction = (julian_day%day+0.5-za)+julian_day%fraction

    DO i = 1, im
       IF (im == 2 .AND. (MOD(ky,4) == 0 .AND. MOD(ky,100) /= 0) &
&           .OR. MOD(ky,400) == 0) THEN
            kd = kd+29
            CYCLE
        ENDIF
        kd = kd + days_in_month(i-1)
    ENDDO

  END SUBROUTINE YD
        
END MODULE mo_julian

