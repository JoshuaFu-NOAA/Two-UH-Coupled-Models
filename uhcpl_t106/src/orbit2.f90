SUBROUTINE ORBIT2 (PCLOCK,PVETIM,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
  !
  !**** *ORBIT2* - COMPUTES SOLAR ORBITAL PARAMETERS.
  !
  !     S.J.LORENZ  UNIVERSITY OF BREMEN                     31-JUL-96.
  !     LAST CHANGE  S.LORENZ      *UNIHB*                   03-JUL-98.
  !       "          U.SCHLESE      DKRZ                     25-SEP-98
  !
  !     PURPOSE.
  !     --------
  !
  !          THIS ROUTINE COMPUTES THE SOLAR CONSTANT NORMALISED BY ITS
  !     ANNUAL MEAN AND THE DECLINATION OF THE SUN AS A FUNCTION OF TIME
  !     FROM VERNAL EQUINOX IN RADIANS.
  !         DIFFERENT PARAMETERS OF THE EARTH'S ORBIT FOR PALEOCLIMATIC
  !     CALCULATIONS CAN BE USED.
  !
  !**   INTERFACE.
  !     ----------
  !
  !     *ORBIT* IS CALLED FROM *PRERAD*.
  !
  !     THERE ARE SEVEN DUMMY ARGUMENTS:
  !
  !     INPUT:
  !     ------
  !     *PCLOCK* TIME OF THE DAY
  !     
  !     *PVETIM* TIME OF THE YEAR FROM THE VERNAL EQUINOX(!) IN RADIANS
  !   
  !     OUTPUT: 
  !     ------  
  !     *PDISSE* RATIO OF THE SOLAR CONSTANT TO ITS ANNUAL MEAN
  !     
  !     *PZEN1*, *PZEN2* AND *PZEN3* ARE ZENITHAL PARAMETERS.
  !
  !     *PRAE*   RATIO OF THE HEIGHT OF THE EQUIVALENT ATMOSPHERE
  !              TO THE RADIUS OF THE EARTH.
  !     
  !
  !     METHOD.
  !     -------
  !
  !          STAIGHTFORWARD. INTERMEDIATE VARIABLES ARE THE ECCENTRIC
  !     ANOMALY AND THE TRUE LONGITUDE ON THE EARTH'S ORBIT.
  !
  !     EXTERNALS.
  !     ----------
  !
  !          NONE.
  !
  !     REFERENCE.
  !     ----------
  !
  !          MONIN,A.S. AN INTRODUCTION TO THE THEORY OF CLIMATE, D.REIDEL
  !     PUBLISHING COMPANY, DORDRECHT, 1986 (PP 10-12).
  !
  !
  !*    DATA STATEMENTS.
  !     ---- -----------
  !
  !     *ZECC*     ECCENTRICITY OF THE EARTH'S ORBIT
  !     *ZOBLD*    OBLIQUITY IN DEGREES
  !     *ZLONP*    LONGITUDE OF PERIHELION MEASURED FROM VERNAL EQUINOX
  !     *ZRAE*     RATIO OF THE HEIGHT OF THE EQUIVALENT ATMOSPHERE
  !                TO THE RADIUS OF THE EARTH
  !
  !
  USE mo_parameters
  USE mo_constants
  !
  !
  DATA ZECC    / 0.016715 /
  DATA ZOBLD   / 23.441 /
  DATA ZLONP   / 282.7 /
  DATA ZRAE    / 0.1277E-02/
  !
  !
  !     ------------------------------------------------------------------
  !
  !*         1.     PRELIMINARY SETTINGS.
  !                 ----------- ---------
  !
100 CONTINUE
  !
  !
  ZOBLR=ZOBLD*API/180.
  ZLONPR=ZLONP*API/180.
  ZSQECC=SQRT((1+ZECC)/(1-ZECC))
  !
  ZEPS=1.E-6
  !
  !     CALCULATION OF ECCENTRIC ANOMALY OF VERNAL EQUINOX
  !
  ZEVE=2.*ATAN(TAN(0.5*ZLONPR)/ZSQECC)
  !
  !     CALCULATION OF TIME ANGLE IN RADIANS OF LONGITUDE OF PERIHELION
  !       FROM VERNAL EQUINOX
  !
  ZTLONPR= ZEVE - ZECC*SIN(ZEVE)
  !
  ZCLOCK=PCLOCK
  ZVETIM=PVETIM
  !
  !
  !     ------------------------------------------------------------------
  !
  !*         2.     COMPUTATIONS.
  !                 -------------
  !
200 CONTINUE
  !
  !     CALCULATE ECCENTRIC ANOMALY:
  !
  !     USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION FOR
  !     EQUATION OF ECCENTRIC ANOMALY *ZENEW*
  !
  ZZTIME=ZVETIM-ZTLONPR
  !
  ZEOLD=ZZTIME/(1.-ZECC)
  ZENEW=ZZTIME
  ITER=0
  !
250 CONTINUE
  ZZEPS=ZEOLD-ZENEW
  IF (ITER.GE.10) GO TO 270
  IF (ABS(ZZEPS).LT.ZEPS) GO TO 280
  ITER=ITER+1
  ZEOLD=ZENEW
  ZCOSE=COS(ZENEW)
  ZENEW=(ZZTIME+ZECC*(SIN(ZENEW)-ZENEW*ZCOSE))/(1.-ZECC*ZCOSE)
  GO TO 250
  !
270 PRINT*,' SUBROUTINE *ORBIT* - ECCENTRIC ANOMALY NOT FOUND!'
  PRINT*,' ERROR IN   *ORBIT* -- STOP'
  !     CALL ERREXIT
  !
280 CONTINUE
  ZDISSE=(1./(1.-ZECC*COS(ZENEW)))**2
  !
  !     CALCULATION OF THE DECLINATION.
  ! 
  ZTGEAN=TAN(ZENEW*0.5)
  !
  !     *ZNU*: TRUE ANOMALY 
  !            (ACTUAL ANGLE OF EARTH'S POSITION FROM PERIHELION)
  !     *ZLAMBDA*: TRUE LONGITUDE OF THE EARTH 
  !                (ACTUAL ANGLE FROM VERNAL EQUINOX)
  !
  ZNU=2.*ATAN(ZSQECC*ZTGEAN)
  ZLAMBDA=ZNU+ZLONPR
  ZSINDE=SIN(ZOBLR)*SIN(ZLAMBDA)
  ZDECLI=ASIN(ZSINDE)
  !
  ZZEN1=SIN(ZDECLI)
  ZZEN2=COS(ZDECLI)*COS(ZCLOCK)
  ZZEN3=COS(ZDECLI)*SIN(ZCLOCK)

  !
  !     ------------------------------------------------------------------
  !
  !*         3.     RETURN.
  !                 -------
  !
300 CONTINUE
  !
  PDISSE=ZDISSE
  PZEN1=ZZEN1
  PZEN2=ZZEN2
  PZEN3=ZZEN3
  PRAE=ZRAE
  !
  !
  RETURN
END SUBROUTINE ORBIT2