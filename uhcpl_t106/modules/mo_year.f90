MODULE mo_year

  USE mo_julian
  USE mo_start_dataset, ONLY: ly365

  ! Different treatment for years with 360 or 365/366 days

  IMPLICIT NONE

  INTEGER, PARAMETER :: idays(12)=(/0,31,59,90,120,151,181,212,243,273,304,334/)
  INTEGER, PARAMETER :: jdays(12)=(/0,31,60,91,121,152,182,213,244,274,305,335/)
  INTEGER, PARAMETER :: kdays(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

  ! cd2dat    - convert century day/julian day to standard date
  ! cd2dy     - convert century day/julian day to day of year and year
  ! ic2ymd    - returns an integer yymmdd/yyyymmdd given the 
  !             century day/julian day
  ! iymd2c    - returns the century day/julian day given an 
  !             integer yymmdd/yyyymmdd
  ! idat2c    - returns forecast day starting at 1st of Jan of first 
  !           - forecast-year
  ! isec2hms  - convert seconds of day into HHMMSS format
  ! ihms2sec  - convert HHMMSS format in seconds of day

CONTAINS

  FUNCTION idat2c(kd,km,ky)

    ! Description:
    !
    ! This function returns forecast day
    ! starting at 1st of jan  of first forecast-year
    ! kd        day
    ! km        month
    ! ky        year
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, March 1989
    ! L. Kornblueh, MPI, May 1998
    ! U. Schulzweida, MPI, May 1998

    !  Function return value 
    INTEGER :: idat2c

    !  Scalar arguments 
    INTEGER :: kd, km, ky

    !  Local scalars: 
    INTEGER :: iy, iyday, ily

    !  Intrinsic functions 
    INTRINSIC MOD


    ! External statements

    IF (ly365) THEN  ! 365 days per year
       iy=MOD(ky,100)
       iyday=idays(km)+kd         
       IF (MOD(iy,4) == 0) iyday=jdays(km)+kd
       ily=(iy-1)/4
       idat2c=iy*365+ily+iyday
    ELSE             ! 360 days per year
       iy = MOD(ky,100)
       iyday = (km-1)*30 + kd
       idat2c = (iy-1)*360 + iyday
    ENDIF

  END FUNCTION idat2c

  SUBROUTINE cd2dat(kcd,kd,km,ky)

    ! Description:
    !
    ! This subroutine converts
    !    kcd       either century day or julian day
    ! back to
    !    kd        day
    !    km        month
    !    ky        year
    !
    ! Authors:
    !
    ! U. Schulzweida, MPI, July 1998
    ! L. Kornblueh, MPI, December 1998
    !

    !  Scalar arguments 
    INTEGER, INTENT(IN)  :: kcd
    INTEGER, INTENT(OUT) :: kd, km, ky

    !  Intrinsic functions 
    INTRINSIC MOD
    ! Local scalars
    TYPE (julian_date) :: julian_day 
    REAL :: zfraction
    INTEGER :: idays


    !  Executable statements 

    IF (ly365) THEN

       ! to modify for noon shifting of julian date
       
       julian_day%day      = REAL(kcd-0.5)
       julian_day%fraction = 0.0         
       
       CALL YMD (julian_day, ky, km, kd, zfraction) 

    ELSE

       ky = (kcd-1)/360 + 1
       idays = MOD(kcd-1,360) + 1
       km = (idays-1)/30 + 1
       kd = MOD(idays-1,30) + 1

    ENDIF

  END SUBROUTINE cd2dat

  SUBROUTINE cd2dy(kcd,kd,ky)

    ! Description:
    !
    ! This subroutine converts
    !    kcd       either century day or julian day in day of year
    ! back to
    !    kd        day
    !    ky        year
    !
    ! Authors:
    !
    ! U. Schulzweida, MPI, July 1998
    ! L. Kornblueh, MPI, December 1998
    !

    !  Scalar arguments 
    INTEGER, INTENT(IN)  :: kcd
    INTEGER, INTENT(OUT) :: kd, ky

    !  Intrinsic functions 
    INTRINSIC MOD

    ! Local scalars
    TYPE (julian_date) :: julian_day 
    REAL :: zfraction


    !  Executable statements 

    IF (ly365) THEN

       ! to modify for noon shifting of julian date
       
       julian_day%day      = REAL(kcd-0.5)
       julian_day%fraction = 0.0         
       
       CALL YD (julian_day, ky, kd, zfraction) 

    ELSE

       ky = (kcd-1)/360 + 1
       kd = MOD(kcd-1,360) + 1

    ENDIF

  END SUBROUTINE cd2dy

  FUNCTION ic2ymd(kcd)

    ! Description:
    !
    ! This function returns an integer yymmdd given the century day
    ! where yy is the year
    !       mm is the month
    !       dd is the day
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, March 1989
    ! L. Kornblueh, MPI, May 1998
    ! U. Schulzweida, MPI, May 1998
    ! H.-S. Bauer, MPI, July 1998

    !  Function Return Value 
    INTEGER :: ic2ymd

    !  Scalar arguments 
    INTEGER :: kcd

    !  Local scalars: 
    INTEGER :: id, im, iy


    !  Executable Statements 

    CALL cd2dat(kcd,id,im,iy)

    ic2ymd = id + (im+iy*100)*100

  END FUNCTION ic2ymd

  FUNCTION iymd2c(kymd)

    ! Description:
    ! This function returns the century day given an integer yymmdd/yyyymmdd
    ! where yy is the year
    !       mm is the month
    !       dd is the day

    !  Function return value 
    INTEGER :: iymd2c

    !  Scalar arguments 
    INTEGER :: kymd

    !  Local types:
    TYPE (julian_date) :: julian_day


    !  Intrinsic functions 
    INTRINSIC MOD


    !  Executable Statements 
    
    IF (ly365) THEN
       julian_day = SetYMD(kymd/10000,ABS(MOD(kymd/100,100)),ABS(MOD(kymd,100)))
       iymd2c = INT(julian_day%day+julian_day%fraction+0.5)
    ELSE
       iymd2c = idat2c(ABS(MOD(kymd,100)),ABS(MOD(kymd/100,100)),kymd/10000)
    END IF

  END FUNCTION iymd2c

  FUNCTION im2day(km,ky)

    ! Description:
    !
    ! Returns number of days of a given month
    !
    ! km        month
    ! ky        year (19xx or xx)
    !
    ! Authors:
    !
    ! U. Schlese, DKRZ, Jan 1996
    ! L. Kornblueh, MPI, May 1998
    ! U. Schulzweida, MPI, May 1998
    !

    ! Function return value
    INTEGER :: im2day

    ! Scalar arguments
    INTEGER :: km,ky

    ! Intrinsic functions
    INTRINSIC MOD


    ! Executable statements

    IF (ly365) THEN
      IF (km /= 2) THEN
         im2day = kdays(km)
      ELSE
         !
         ! check for a leap year
         !
         IF(MOD(ky,400) == 0) THEN
            im2day=29
         ELSE IF (MOD(ky,100) == 0) THEN
            im2day=28
         ELSE IF (MOD(ky,4) /= 0) THEN   
            im2day=28
         ELSE
            im2day=29
         END IF
      END IF
    ELSE
      im2day = 30
    END IF

  END FUNCTION im2day

  FUNCTION isec2hms(isec) 
 
    ! convert seconds of day into HHMMSS format
  
    INTEGER :: isec, isec2hms, ihh, imm, iss
 
    ! Intrinsic functions
    INTRINSIC MOD
 
    ihh = isec/3600
    imm = MOD(isec,3600)/60
    iss = MOD(isec,60)
 
    isec2hms = ihh*10000 + imm*100 + iss
    
  END FUNCTION isec2hms

  FUNCTION ihms2sec(hms) 

    ! convert HHMMSS format into seconds of day

    INTEGER :: hms, ihms2sec, ihh, imm, iss

    ! Intrinsic functions
    INTRINSIC MOD

    ihh = hms/10000
    imm = MOD(hms,10000)/100
    iss = MOD(hms,100)

    ihms2sec = ihh*3600 + imm*60 + iss

  END FUNCTION ihms2sec

END MODULE mo_year
