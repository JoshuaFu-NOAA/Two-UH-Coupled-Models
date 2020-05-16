MODULE mo_soil_impl

  IMPLICIT NONE

  !      -----------------------------------------------------------------
  !
  ! *    Parameters used in implicit soil temperature scheme
  !
  !      -----------------------------------------------------------------

  INTEGER, PARAMETER :: ngrndmx=5         ! number of soil layers

  REAL, PARAMETER :: snowcri=1.5          ! only sublimation above this value
  REAL, PARAMETER :: sneige=snowcri/1000. ! lower limit of snow amount

  !      ---------------------------------------------------------------
  !
  ! *    *common* *comsoil* - quantities needed for the implicit soil scheme
  !
  !      ---------------------------------------------------------------

  REAL :: sn_capa
  REAL :: sn_cond
  REAL :: sn_dens

  REAL :: dz1(ngrndmx)
  REAL :: dz2(ngrndmx)
  REAL ::  zz(ngrndmx)

  REAL, ALLOCATABLE :: cgrnd(:,:,:)
  REAL, ALLOCATABLE :: dgrnd(:,:,:)
  REAL, ALLOCATABLE :: pcapcal(:,:)
  REAL, ALLOCATABLE :: pfluxgrd(:,:)
  REAL, ALLOCATABLE :: rnatur(:,:)
  REAL, ALLOCATABLE :: ptn(:,:,:)

CONTAINS

  SUBROUTINE inisoil2ech(jt, ngrid, ptimestep, ptsol, td4, td5, &
                         tdp, tdcl, so_tdif, so_capa, snow)

    !   Author:  Frederic Hourdin     30/01/92
    !   -------
    !            Adapted to the LMD-GCM by Jan Polcher  26/02/92
    !            Adapted to ECHAM-GCM,  Jan-Peter Schulz   April 1996
    !
    !            J.-P. Schulz   MPI - October 1997 :
    !               Routine used for implementation of full-implicit
    !               coupling of land surface and atomsphere in the
    !               ECHAM4 GCM.
    !
    !   inputs:
    !   -------
    !   jt                       Index of the line being treated
    !   ngrid                    Number of gridpoints per line
    !   ptimestep                Time-step (s)
    !   ptsol(ngrid)             Initial temperature at soil surface
    !   td4(ngrid)                  "          "     of second ECHAM layer
    !   td5(ngrid)                  "          "     of third  ECHAM layer
    !   tdp(ngrid)                  "          "     of fourth ECHAM layer
    !   tdcl(ngrid)                 "          "     of fifth  ECHAM layer
    !   so_tdif(ngrid)           Soil temperature diffusivity [m**2/s] (FAO)
    !   so_capa(ngrid)           Soil vol. heat capacity    [J/m**3/K] (FAO)
    !   snow(ngrid)              Snow depth (liquid water equivalent) [m]
    !
    !   outputs:
    !   --------
    !   cgrnd(nlon,ngl,ngrndmx) coefficient of the soil temperature scheme
    !   dgrnd(nlon,ngl,ngrndmx) coefficient of the soil temperature scheme
    !   pcapcal(nlon,ngl)       specific heat (full grid) (W*m-2*s-1*K-1)
    !   pfluxgrd(nlon,ngl)      surface diffusive flux from ground
    !                            (full grid) (W*m-2)
    !

    USE mo_parameters
    USE mo_control
    USE mo_gaussgrid
    USE mo_truncation
    USE mo_semi_impl
    USE mo_hdiff
    USE mo_hyb
    USE mo_diagnostics
    USE mo_rad_switches
    USE mo_param_switches
    USE mo_forecast_switches
    USE mo_physc1
    USE mo_physc2
    USE mo_constants
    USE mo_start_dataset

    IMPLICIT NONE

    INTEGER :: jt, ngrid
    REAL :: ptimestep
    REAL :: ptsol(ngrid)
    REAL :: td4(ngrid), td5(ngrid)
    REAL :: tdp(ngrid), tdcl(ngrid)
    REAL :: so_tdif(ngrid), so_capa(ngrid)
    REAL :: snow(ngrid)

    !  local:

    INTEGER :: ig, jk
    REAL :: so_cond(ngrid)
    REAL :: pkappa(nlon,ngrndmx)
    REAL :: pcapa(nlon,ngrndmx)
    REAL :: pt(nlon,ngrndmx)
    REAL :: zdz2(nlon,ngrndmx)
    REAL :: zdz1(nlon,ngrndmx)
    REAL :: snow_h, zx1, zx2, zsncri, z1
    REAL :: lambda


    !  Set a few variables to zero.

    DO ig = 1,ngrid
      pcapcal(ig,jt)  = 0.0
      pfluxgrd(ig,jt) = 0.0
    END DO

    cgrnd(:,jt,:) = 0.
    dgrnd(:,jt,:) = 0.

    ! Initialize pt

    DO ig = 1,ngrid

      !  Land surface or glacier

      IF (rnatur(ig,jt) > 0.0) THEN
        pt(ig,1) = ptsol(ig)
        pt(ig,2) = td4(ig)
        pt(ig,3) = td5(ig)
        pt(ig,4) = tdp(ig)
        pt(ig,5) = tdcl(ig)
      ELSE

        !  Sea surface temperature (SST) or sea ice

        DO jk=1,ngrndmx
          pt(ig,jk) = 0.0
        END DO
      END IF

    END DO


    !*    1.  Specifying the depths of the temperature levels.
    !         ------------------------------------------------

    !*    1.1 Some constants used in the temperature scheme.

    sn_cond = 0.22              ! Snow thermal conductivity [J/s/m/K]
    sn_dens = 300.0             ! Snow density              [kg/m**3]
    sn_capa = 634500.0          ! Snow vol. heat capacity   [J/m**3/K]
    zsncri  = csncri            ! Critical snow height  [m water equ.]

    DO jk = 1, ngrndmx
      dz2(jk) = cdel(jk)     ! Thickness of soil layers
      zz(jk)  = cmid(jk)     ! Depth of mids of soil layers
    ENDDO

    !*    1.2 Computing some useful constants.

    DO jk = 1,ngrndmx-1
      dz1(jk) = 1./(zz(jk+1)-zz(jk))
    END DO
    lambda = zz(1)*dz1(1)

    !*    1.3 Compute of the soil thermal conductivity [J/s/m/K] from
    !*        the soil temperature diffusivity [m**2/s].

    DO ig = 1,ngrid
      so_cond(ig) = so_capa(ig)*so_tdif(ig)
    END DO

    !*    1.4 Pre-set thermal conductivity at all levels.

    DO jk = 1,ngrndmx
      DO ig = 1,ngrid
        pkappa(ig,jk) = so_cond(ig)
        pcapa(ig,jk)  = so_capa(ig)
      END DO
    END DO

    !   Computation of the Cgrd and Dgrd coefficient for the next step:

    DO ig = 1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN

        IF (rnatur(ig,jt)<2.0 .AND. rnatur(ig,jt)>0.0) THEN
          snow_h = snow(ig)*1000. / sn_dens
        ELSE
          snow_h = 10.0
        ENDIF

        !*       Traitement special pour le premiere couche

        IF ( snow_h > zz(2) ) THEN
          pcapa(ig,1)  = sn_capa
          pkappa(ig,1) = sn_cond
        ELSE IF (snow_h > 0.0 .AND. snow_h <= zz(2)) THEN 
          zx1 = snow_h / zz(2)
          zx2 = ( zz(2) - snow_h) / zz(2)
          pcapa(ig,1)  = zx1*sn_capa + zx2*so_capa(ig)
          pkappa(ig,1) = 1.0/(zx1/sn_cond + zx2/so_cond(ig)) 
        ELSE
          pcapa(ig,1)  = so_capa(ig)
          pkappa(ig,1) = so_cond(ig)
        ENDIF

        DO jk = 2, ngrndmx - 2
          IF ( snow_h > zz(jk+1) ) THEN
            pcapa(ig,jk)  = sn_capa
            pkappa(ig,jk) = sn_cond
          ELSE IF ( snow_h > zz(jk) .AND. snow_h <= zz(jk+1) ) THEN
            zx1 = (snow_h - zz(jk)) * dz1(jk)
            zx2 = ( zz(jk+1) - snow_h) * dz1(jk)
            pcapa(ig,jk)  = zx1*sn_capa+zx2*so_capa(ig)
            pkappa(ig,jk) = 1.0/(zx1/sn_cond + zx2/so_cond(ig))
          ELSE
            pcapa(ig,jk)  = so_capa(ig)
            pkappa(ig,jk) = so_cond(ig)
          ENDIF
        ENDDO

      END IF
    ENDDO

    DO jk=1,ngrndmx
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          zdz2(ig,jk) = pcapa(ig,jk) * dz2(jk)/ptimestep 
        END IF
      ENDDO
    ENDDO

    DO jk=1,ngrndmx-1
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          zdz1(ig,jk) = dz1(jk) * pkappa(ig,jk)
        END IF
      ENDDO
    ENDDO

    DO ig=1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN
        z1 = zdz2(ig,ngrndmx) + zdz1(ig,ngrndmx-1)
        cgrnd(ig,jt,ngrndmx-1) = zdz2(ig,ngrndmx)*pt(ig,ngrndmx)/z1 
        dgrnd(ig,jt,ngrndmx-1) = zdz1(ig,ngrndmx-1)/z1 
      END IF
    ENDDO

    DO jk=ngrndmx-1,2,-1
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          z1 = 1./(zdz2(ig,jk) + zdz1(ig,jk-1) + &
               zdz1(ig,jk) * (1.-dgrnd(ig,jt,jk)))
          cgrnd(ig,jt,jk-1) =  &
               (pt(ig,jk)*zdz2(ig,jk)+zdz1(ig,jk)*cgrnd(ig,jt,jk))*z1 
          dgrnd(ig,jt,jk-1) = zdz1(ig,jk-1)*z1
        END IF
      ENDDO
    ENDDO

    !   computation of the surface diffusive flux from ground and
    !   calorific capacity of the ground:

    DO ig=1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN
        pfluxgrd(ig,jt) = zdz1(ig,1) * (cgrnd(ig,jt,1) &
                        + (dgrnd(ig,jt,1)-1.) * pt(ig,1))
        pcapcal(ig,jt)  = (zdz2(ig,1) * ptimestep + &
                          ptimestep * (1.-dgrnd(ig,jt,1)) * zdz1(ig,1))

        IF (snow(ig) >= zsncri) THEN
          z1 = lambda * (1.-dgrnd(ig,jt,1)) + 1.
          pcapcal(ig,jt)  = pcapcal(ig,jt) / z1
          pfluxgrd(ig,jt) = pfluxgrd(ig,jt) + pcapcal(ig,jt) * &
               (pt(ig,1) * z1 - lambda * cgrnd(ig,jt,1) - ptsol(ig)) / ptimestep
        END IF

      END IF
    ENDDO

    RETURN
  END SUBROUTINE inisoil2ech

  SUBROUTINE soil2ech(jt, ngrid, ptimestep, ptsol, so_tdif, so_capa, snow)

    !   Author:  Frederic Hourdin     30/01/92
    !   -------
    !            Adapted to the LMD-GCM by Jan Polcher  26/02/92
    !            Adapted to the ECHAM-GCM by Jan-Peter Schulz, MPI  03/02/96
    !
    !            J.-P. Schulz   MPI - October 1997 :
    !               Routine used for implementation of full-implicit
    !               coupling of land surface and atomsphere in the
    !               ECHAM4 GCM.
    !
    !   Object:  computation of : the ground temperature evolution
    !   ------                    the ground specific heat "Capcal"
    !                             the surface diffusive flux from ground "F0"
    !
    !
    !   Method:  implicit time integration
    !   -------
    !   Consecutives ground temperatures are related by:
    !           T(k+1) = C(k) + D(k)*T(k)  (1)
    !   the coefficients C and D are computed at the t-dt time-step.
    !   Routine structure:
    !   1)new temperatures are computed  using (1)
    !   2)C and D coefficients are computed from the new temperature
    !     profile for the t+dt time-step
    !   3)the coefficients A and B are computed where the diffusive
    !     fluxes at the t+dt time-step is given by
    !            Fdiff = A + B Ts(t+dt)
    !     or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
    !            with F0 = A + B (Ts(t))
    !                 Capcal = B*dt
    !           
    !   Interface:
    !   ----------
    !
    !   Arguments:
    !   ----------
    !   Inputs:
    !
    !   jt                       Index of the line being treated
    !   ptimestep                Time-step (s)
    !   ptsol(ngrid)             Initial temperature at soil surface
    !   so_tdif(ngrid)           Soil temperature diffusivity [m**2/s] (FAO)
    !   so_capa(ngrid)           Soil vol. heat capacity    [J/m**3/K] (FAO)
    !
    !   Outputs:
    !   --------
    !   ptn(nlon,ngl,ngrndmx)   ground temperatures of ngrndmx layers
    !   cgrnd(nlon,ngl,ngrndmx) coefficient of the soil temperature scheme
    !   dgrnd(nlon,ngl,ngrndmx) coefficient of the soil temperature scheme
    !   pcapcal(nlon,ngl)       specific heat (full grid) (W*m-2*s-1*K-1)
    !   pfluxgrd(nlon,ngl)      surface diffusive flux from ground
    !                            (full grid) (W*m-2)
    !   

    USE mo_parameters
    USE mo_control
    USE mo_gaussgrid
    USE mo_truncation
    USE mo_semi_impl
    USE mo_hdiff
    USE mo_hyb
    USE mo_diagnostics
    USE mo_rad_switches
    USE mo_param_switches
    USE mo_forecast_switches
    USE mo_physc1
    USE mo_physc2
    USE mo_constants
    USE mo_start_dataset

    INTEGER :: jt, ngrid
    REAL :: ptimestep
    REAL :: ptsol(ngrid), so_tdif(ngrid), so_capa(ngrid), snow(ngrid)

    !  local arrays

    INTEGER :: ig, jk
    REAL :: so_cond(ngrid)
    REAL :: pkappa(nlon,ngrndmx)
    REAL :: pcapa(nlon,ngrndmx)
    REAL :: zdz2(nlon,ngrndmx)
    REAL :: zdz1(nlon,ngrndmx)
    REAL :: snow_h, zx1, zx2, zsncri, z1
    REAL :: lambda


    sn_cond = 0.22              ! Snow thermal conductivity [J/s/m/K]
    sn_dens = 300.0             ! Snow density              [kg/m**3]
    sn_capa = 634500.0          ! Snow vol. heat capacity   [J/m**3/K]
    zsncri  = csncri          ! Critical snow height  [m water equ.]

    DO jk = 1, ngrndmx
      dz2(jk) = cdel(jk)     ! Thickness of soil layers
      zz(jk)  = cmid(jk)     ! Depth of mids of soil layers
    ENDDO

    !*    1.2 Computing some useful constants.

    DO jk = 1,ngrndmx-1
      dz1(jk) = 1./(zz(jk+1)-zz(jk))
    END DO
    lambda = zz(1)*dz1(1)

    !*    1.3 Compute of the soil thermal conductivity [J/s/m/K] from
    !*        the soil temperature diffusivity [m**2/s].

    DO ig = 1,ngrid
      so_cond(ig) = so_capa(ig)*so_tdif(ig)
    ENDDO

    !*    1.4 Pre-set thermal conductivity at all levels.

    DO jk = 1,ngrndmx
      DO ig = 1,ngrid
        pkappa(ig,jk) = so_cond(ig)
        pcapa (ig,jk) = so_capa(ig)
      ENDDO
    ENDDO

    !   A few variables are set to zero

    DO jk = 1,ngrndmx
      DO ig = 1,ngrid
        ptn(ig,jt,jk) = 0.0
      END DO
    END DO

    DO ig = 1,ngrid
      pcapcal (ig,jt)  = 0.0
      pfluxgrd(ig,jt)  = 0.0
    ENDDO

    !   Computation of the ground temperatures using the Cgrd and Dgrd
    !   coefficients computed at the previous time-step:

    !    surface temperature

    DO ig = 1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN
        IF (snow(ig) >= zsncri) THEN
          ptn(ig,jt,1) = (lambda*cgrnd(ig,jt,1) + &
                         ptsol(ig))/(lambda*(1.-dgrnd(ig,jt,1))+1.)
        ELSE
          ptn(ig,jt ,1) = ptsol(ig)
        END IF
      END IF
    ENDDO

    !   other temperatures

    DO jk = 1,ngrndmx-1
      DO ig = 1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          ptn(ig,jt,jk+1)=cgrnd(ig,jt,jk)+dgrnd(ig,jt,jk)*ptn(ig,jt,jk)
        END IF
      ENDDO
    ENDDO

    !   Computation of the Cgrd and Dgrd coefficients for the next step:

    DO ig = 1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN

        IF (rnatur(ig,jt) < 2.0 .AND. rnatur(ig,jt) > 0.0) THEN
          snow_h = snow(ig)*1000. / sn_dens
        ELSE
          snow_h = 10.0
        ENDIF

        !*       Traitement special pour le premiere couche

        IF ( snow_h > zz(2) ) THEN
          pcapa(ig,1) = sn_capa
          pkappa(ig,1) = sn_cond
        ELSE IF (snow_h > 0.0 .AND. snow_h <= zz(2)) THEN
          zx1 = snow_h / zz(2)
          zx2 = ( zz(2) - snow_h) / zz(2)
          pcapa(ig,1) = zx1 * sn_capa + zx2 * so_capa(ig)
          pkappa(ig,1) = 1.0/(zx1/sn_cond + zx2/so_cond(ig))
        ELSE
          pcapa(ig,1)  = so_capa(ig)
          pkappa(ig,1) = so_cond(ig)
        ENDIF

        DO jk = 2, ngrndmx - 2
          IF ( snow_h > zz(jk+1) ) THEN
            pcapa(ig,jk)  = sn_capa
            pkappa(ig,jk) = sn_cond
          ELSE IF ( snow_h > zz(jk) .AND. snow_h <= zz(jk+1)) THEN
            zx1 = (snow_h - zz(jk)) * dz1(jk)
            zx2 = ( zz(jk+1) - snow_h) * dz1(jk)
            pcapa(ig,jk)  = zx1*sn_capa + zx2*so_capa(ig)
            pkappa(ig,jk) = 1.0/(zx1/sn_cond + zx2/so_cond(ig))
          ELSE
            pcapa(ig,jk)  = so_capa(ig)
            pkappa(ig,jk) = so_cond(ig)
          ENDIF
        ENDDO

      END IF
    ENDDO

    DO jk=1,ngrndmx
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          zdz2(ig,jk) = pcapa(ig,jk) * dz2(jk)/ptimestep
        END IF
      ENDDO
    ENDDO

    DO jk=1,ngrndmx-1
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          zdz1(ig,jk) = dz1(jk) * pkappa(ig,jk)
        END IF
      ENDDO
    ENDDO

    DO ig=1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN
        z1 = zdz2(ig,ngrndmx)+zdz1(ig,ngrndmx-1)
        cgrnd(ig,jt,ngrndmx-1) = zdz2(ig,ngrndmx)*ptn(ig,jt,ngrndmx)/z1
        dgrnd(ig,jt,ngrndmx-1) = zdz1(ig,ngrndmx-1)/z1
      END IF
    ENDDO

    DO jk=ngrndmx-1,2,-1
      DO ig=1,ngrid
        IF (rnatur(ig,jt) > 0.0) THEN
          z1 = 1./(zdz2(ig,jk) + zdz1(ig,jk-1) + &
               zdz1(ig,jk) * (1.-dgrnd(ig,jt,jk)))
          cgrnd(ig,jt,jk-1) = (ptn(ig,jt,jk) * zdz2(ig,jk) + &
               zdz1(ig,jk)*cgrnd(ig,jt,jk))*z1
          dgrnd(ig,jt,jk-1) = zdz1(ig,jk-1)*z1
        END IF
      ENDDO
    ENDDO

    !   computation of the surface diffusive flux from ground and
    !   calorific capacity of the ground:

    DO ig=1,ngrid
      IF (rnatur(ig,jt) > 0.0) THEN
        pfluxgrd(ig,jt) = zdz1(ig,1) * (cgrnd(ig,jt,1) &
                        + (dgrnd(ig,jt,1)-1.) * ptn(ig,jt,1))
        pcapcal(ig,jt)  = (zdz2(ig,1) * ptimestep + &
                           ptimestep * (1.-dgrnd(ig,jt,1)) * zdz1(ig,1))

        IF (snow(ig) >= zsncri) THEN
          z1 = lambda * (1.-dgrnd(ig,jt,1)) + 1.
          pcapcal(ig,jt)  = pcapcal(ig,jt) / z1
          pfluxgrd(ig,jt) = pfluxgrd(ig,jt) +                     &
                            pcapcal(ig,jt)*(ptn(ig,jt,1) * z1 -   &
                            lambda * cgrnd(ig,jt,1) - ptsol(ig)) / ptimestep
        END IF
      END IF
    ENDDO

    RETURN
  END SUBROUTINE soil2ech

  SUBROUTINE evapotr4(nlon, nlp2, laland,                   &
                      pcvs, pcvw, pvgrat, zeveff, zrelhum,  &
                      zqnatm, zqs, cdrag, qsol, qair,       &
                      zqsnfl, zqlfl, zqvfl, zqgfl)

    ! calculates sublimation, interception loss,
    ! transpiration and evaporation from bare soil.

    !  Authors:
    !
    !  J.-P. Schulz                       *MPI*   9/10/97
    !  based on SECHIBA subroutines by
    !  N. Ducoudre and J. Polcher         *LMD*   7/25/91

    !     Purpose.
    !
    !       This subroutine gets the atmospheric conditions 
    !  calculated by the PBL at that timestep and diagnoses 
    !  the different components of evaporation of the
    !  land surface model.
    !
    !         NLON    -  Number of points per line.(I)
    !         NLP2    -  NLON + 2.(I)
    !         LALAND  -  Land sea mask.(I)
    !         PCVS    -  Fractional snow cover.(I)
    !         PCVW    -  Fractional skin reservoir content.(I)
    !         PVGRAT  -  Fractional vegetation cover.(I)
    !         ZEVEFF  -  Evaporation efficiency (transpiration).(I)
    !         ZRELHUM -  Relative surface air humidity
    !                    (bare soil evaporation).(I)
    !         ZQNATM  -  Old atm. spec. humidity for old temp. PTM1.(I)
    !         ZQS     -  Surface sat. spec. humidity for old temp. PTM1.(I)
    !         CDRAG   -  Transfer coefficient.(I)
    !         QSOL    -  New saturation specific humidity.(I)
    !         QAIR    -  New air specific humidity.(I)
    !         ZQSNFL  -  Sublimation.(O)           [mm/s] \
    !         ZQLFL   -  Interception loss.(O)     [mm/s]  \ as they were used
    !         ZQVFL   -  Transpiration.(O)         [mm/s]  / in the implicit
    !         ZQGFL   -  Bare soil evaporation.(O) [mm/s] /  scheme of surftemp

    !  Method.
    !
    !       This subroutine calculates the different evaporations with
    !  the resistances which were calculated in the subroutine VDIFF.
    !  But this time the computations are done with the new values for
    !  temperature and specific humidity as they are calculated by the
    !  PBL (in subroutine VDIFF) and the surface routine SURFTEMP.

    USE mo_parameters

    IMPLICIT NONE

    ! -- INPUT
    INTEGER :: nlon, nlp2
    LOGICAL :: laland(nlp2)
    REAL :: pcvs(nlp2), pcvw(nlp2), pvgrat(nlp2),                              &
            zeveff(nlon), zrelhum(nlon), zqnatm(nlon), zqs(nlon), cdrag(nlon), &
            qsol(nlon), qair(nlon)

    ! -- OUTPUT
    REAL :: zqsnfl(nlon), zqlfl(nlon), zqvfl(nlon), zqgfl(nlon)

    ! -- Local scalars:
    INTEGER :: jl
    LOGICAL :: lo
    !
    ! CALCULATION OF THE MOISTURE FLUXES FROM THE 4 DIFFERENT RESERVOIRS
    !
    ! ----------------------------------------------------------------------
    ! RESERVOIR :   SNOW      SKIN RESERVOIR      VEGETATION     BARE SOIL
    ! ----------------------------------------------------------------------
    ! FLUX      :   ZQSNFL    ZQLFL               ZQVFL          ZQGFL
    ! ----------------------------------------------------------------------
    !
    !*    1.   Main loop.

    DO jl = 1,nlon
      IF (laland(jl)) THEN

        ! 2.   Sublimation on non ocean points.

        zqsnfl(jl) = cdrag(jl)*pcvs(jl)*(qair(jl)-qsol(jl))

        ! 3.   Interception loss.

        zqlfl(jl) = cdrag(jl)*(1.-pcvs(jl))*pcvw(jl) &
                  *(qair(jl)-qsol(jl))

        ! 4.   Transpiration.

        zqvfl(jl)  = cdrag(jl)*(1.-pcvs(jl))*(1.-pcvw(jl)) &
                   *pvgrat(jl)*zeveff(jl)*(qair(jl)-qsol(jl))

        ! 5.   Bare soil evaporation.

        IF ((zrelhum(jl) <= zqnatm(jl)/zqs(jl)) .AND. &
             (zqnatm(jl) <= zqs(jl))) THEN
          zqgfl(jl) = 0.
        ELSE
          lo=zqnatm(jl) > zqs(jl)
          zrelhum(jl) = MERGE(1.,zrelhum(jl),lo)
          zrelhum(jl) = MERGE(zrelhum(jl),1.,laland(jl))

          zqgfl(jl)   = cdrag(jl)*(1.-pcvs(jl))*(1.-pcvw(jl)) &
                      *(1.-pvgrat(jl))*(qair(jl)-zrelhum(jl)*qsol(jl))
        END IF
      ELSE
        zqsnfl(jl) = 0.
        zqlfl(jl)  = 0.
        zqvfl(jl)  = 0.
        zqgfl(jl)  = 0.
      END IF
    END DO

    RETURN
  END SUBROUTINE evapotr4

  SUBROUTINE surftemp(klon, klp2, laland, pdt, pemi, pboltz, pcp, pc16,     &
                      platev, platsu, pfscoe, pescoe, pfqcoe, peqcoe,       &
                      psold, pqsold, pdqsold, pnetrad, pgrdfl, pcfh, pcair, &
                      pcsat, pfracsu, pgrdcap, psnew, pqsnew) 

    !  Authors:
    !
    !  J. POLCHER  *LMD*  and  J.-P. SCHULZ  *MPI*,  MAY 1995
    !
    !  Computes the energy balance at the surface with an implicit scheme 
    !  that is connected to the Richtmyer and Morton algorithm of the PBL.
    !
    !  Modifications
    !  J.-P. Schulz  *MPI*,  October 1997:
    !     Modify according to latent heat flux formulation in vdiff
    !     using zcair and zcsat coefficients.
    !  J.-P. Schulz  *MPI*,  August 1998:
    !     Modify according to sensible heat flux formulation in vdiff.


    !  Input
    !  -----
    !  klon     : Length of arrays to be used
    !  pdt      : timestep in seconds
    !  pemi     : surface emissivity
    !  pboltz   : Stefan-Boltzman constant
    !  pcp      : Specific heat of air
    !  pc16     : cpd*vtmpc2=cpd*(delta-1) for sens. heat flux (cf. vdiff)
    !  platev   : Latent heat of evaporation
    !  platsu   : Latent heat of sublimation 
    !  pfscoe, pescoe : Coefficients of the Richtmyer and Morton scheme 
    !                   for dry static energy
    !  pfqcoe, peqcoe : as above but for specific humidity
    !  psold   : Old surface dry static energy (Ts * Cp)
    !  pqsold  : Saturated  specific humidity for old temperature
    !  pdqsold : Derivative of saturated  specific humidity at the 
    !            old temperature
    !  pnetrad : Net radiation at the surface (upward longwave is 
    !            included but for the old surface temperature) 
    !  pgrdfl  : Ground heat flux
    !  pcfh    : Diffusion coefficient for static energy and moisture
    !  pcair   : Coefficient in latent heat flux formula (see VDIFF)
    !  pcsat   : Coefficient in latent heat flux formula (see VDIFF)
    !  pfracsu : Fraction of surface for sublimation
    !  pgrdcap : Surface heat capacity 

    !  Output
    !  ------
    !  psnew   : New surface static energy
    !  pqsnew  : New saturated surface air moisture

    IMPLICIT NONE

    INTEGER :: klon, klp2
    REAL :: pdt, pemi, pboltz, pcp(klon), pc16, platev, platsu
    REAL :: pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL :: psold(klon), pqsold(klon), pdqsold(klon)
    REAL :: pnetrad(klon), pgrdfl(klon)
    REAL :: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL :: pgrdcap(klon)
    REAL :: psnew(klon), pqsnew(klon)
    LOGICAL :: laland(klp2)

    ! Local variables
    INTEGER :: jl
    REAL :: zcolin, zcohfl, zcoind, zicp, zca, zcs


    DO jl = 1,klon
      IF (laland(jl)) THEN
        zicp = 1./pcp(jl)
        zca  = platsu*pfracsu(jl)+platev*(pcair(jl)-pfracsu(jl))
        zcs  = platsu*pfracsu(jl)+platev*(pcsat(jl)-pfracsu(jl))

        zcolin = pgrdcap(jl)*zicp +                                &
                 pdt*(zicp*4.*pemi*pboltz*((zicp*psold(jl))**3.) - &
                 pcfh(jl)*(zca*peqcoe(jl)-zcs-zicp*pc16*psold(jl)* &
                (pcair(jl)*peqcoe(jl)-pcsat(jl)))*zicp*pdqsold(jl))

        zcohfl = -pdt*pcfh(jl)*(pescoe(jl)-1.)

        zcoind = pdt*(pnetrad(jl)+pcfh(jl)*pfscoe(jl)+                     &
                 pcfh(jl)*((zca*peqcoe(jl)-zcs)*pqsold(jl)+zca*pfqcoe(jl)- &
                 zicp*pc16*psold(jl)*                                      &
                 ((pcair(jl)*peqcoe(jl)-pcsat(jl))*pqsold(jl)+             &
                 pcair(jl)*pfqcoe(jl)))+pgrdfl(jl))

        psnew(jl)  = (zcolin*psold(jl)+zcoind)/(zcolin+zcohfl)
        pqsnew(jl) = pqsold(jl)+zicp*pdqsold(jl)*(psnew(jl)-psold(jl))
      ELSE
        ! Over sea points we have to use the imposed values of the sea
        ! surface temperature to compute the new atmospheric values and
        ! the surface fluxes.

        psnew(jl)  = psold(jl)
        pqsnew(jl) = pqsold(jl)
      END IF
    END DO

    RETURN
  END SUBROUTINE surftemp

END MODULE mo_soil_impl
