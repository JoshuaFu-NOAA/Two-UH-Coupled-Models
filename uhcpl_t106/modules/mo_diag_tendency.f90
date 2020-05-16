MODULE mo_diag_tendency

  ! Authors:

  ! I. Kirchner, MPI, October 1998 -
  ! I. Kirchner, MPI, May-2000, patch v3.22p1

  !-------------------------------------------------------------------
  ! overview about the separation of terms in diagnostic equations
  !
  ! ***************** Interface to ECHAM *************************
  !
  ! *ECHAM*       *DIAGNOSTICS*
  !
  ! CONTROL -+
  !          +--> DIAG_Init(1)                general initialization
  !                  |
  !                  +---> DIAG_Rerun('R')    read rerun data
  !
  ! STEPON --+
  !          +--> DIAG_SumAll                 count total tendency
  !          |
  !          +--> DIAG_Write                  store data in model output
  !          |       |
  !          |       +---> DIAG_Check
  !          |       |
  !          |       +---> DiagWriteSP
  !          |
  !          +--> DIAG_Init(-1)               clear diagnostic memory
  !                  |
  !                  +---> DIAG_Rerun('W')    store rerun data
  !
  ! SCAN1SL -+
  !          +--> DIAG_fftd
  !          |
  !          SI2
  !          |
  !          +--> DIAG_SpecTrans
  !                  |
  !                  +--> DIAG_sym1
  !                  |
  !                  +--> DIAG_ltd
  !
  !---------------------------------------------------------------------
  ! count terms in additional arrays
  !   *spectral*    *Gaussian*
  !
  !                 DYN            adiabatic terms
  !                  |
  !                 TF2            time filter
  !                  |
  !                 TF1            time filter
  !                  |
  !                 GPC
  !                  +--> PHYSC    diabatic terms
  !                  |      |
  !                  |      +---> M_RADHEAT long/short wave
  !                  |
  !                  +--> SI1      semi-implicit terms
  !                  |
  !                 DIAG_fftd
  !                 |
  !               SI2                semi-implicit terms
  !               |
  !             DIAG_SpecTrans
  !                     |
  !                   DIAG_sym1
  !                     |
  !                   DIAG_ltd
  ! SCCD                           semi-implicit terms
  !   |
  ! SCCTP                          semi-implicit terms
  !   |
  ! HDIFF                          horizontal diffusion
  !   |
  ! DIAG_SumAll
  !   |
  ! DIAG_Write
  !
  !**************************************************************

  USE mo_doctor,        ONLY: nout
  USE mo_exception,     ONLY: finish, message
  USE mo_start_dataset, ONLY: lres, nstep
  USE mo_year,          ONLY: ic2ymd

  USE mo_truncation,    ONLY: nmp, nnp
  USE mo_control, ONLY: &
&       nlev, nlevp1, nsp, n2sp, nlp2, ngl, nmp1, nhgl, &
&       lptime, twodt, nrow, nlon, nptime, nvclev, vct, dtime

  USE mo_fft,           ONLY: trig, nfax, fft991cy
  USE mo_legendre,      ONLY: pnmd, anmd, rnmd, legmod

  USE mo_post,          ONLY: nunitsp, nunitdf
  USE mo_grib
  USE mo_memory_sp,     ONLY: sd, svo, stp
  USE mo_time_control

  IMPLICIT NONE

  !  LTDIAG      additional switch in &RUNCTL
  !              .true. ==> addional diagnostics
  !
  PUBLIC  :: DIAG_Init      ! allocate/deallocate memory, initialize diagnostic arrays
  PUBLIC  :: DIAG_SpecTrans ! second part of spectral transform
  PUBLIC  :: DIAG_fftd      ! Fourier transformation of diagnostics arrays
  PRIVATE :: DIAG_sym1      ! composition of symetric/antisymetric part
  PRIVATE :: DIAG_ltd       ! legendre transform of diagnostic arrays
  PUBLIC  :: DIAG_Write     ! store diagnostic fields
  PRIVATE :: DIAG_Rerun     ! handle rerun data
  PRIVATE :: DIAG_Check     ! global check of accounting error

  ! special diagnostic fields

  REAL, ALLOCATABLE :: &    ! ARRAYS in spectral space
&       pdvor (:,:,:,:), &  !   terms for vorticity equation
&       pddiv (:,:,:,:), &  !   terms for divergence equation
&       pdtem (:,:,:,:), &  !   terms for temperature equation
&       pdprs (:,:,:), &    !   terms for surface pressure equation
&       pdprl (:,:,:)       !   surface pressure term vertical integration

  REAL, ALLOCATABLE :: &    ! old value for accumulation of tendency
&       pdovor (:,:,:,:), pdodiv (:,:,:,:), pdotep (:,:,:,:)
  LOGICAL :: lpdo1, lpdo2   ! marker for correct storage
  INTEGER :: ipdoaccu       ! counter for total tendency part
  INTEGER :: ipdaccu        ! general counter for  other parts

  LOGICAL :: ldinit         ! state of initialization

  INTEGER, PARAMETER :: NO_PTFH = 3
  REAL, ALLOCATABLE ::   &  ! ARRAYS in grid point space
&       ptfh1 (:,:,:,:), &  !   time integration parts
&       ptfh2 (:,:)

  INTEGER, PARAMETER :: NO_PDIGA=25  ! no. of 3-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDSGA=2   ! no. of 2-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDIGB=10  ! no. of 3-dim fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS=4   ! no. of 2-dim mixed grid point and fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS_FRAC=3   ! no. of 2-dim grid point space fields fraction

  REAL, ALLOCATABLE :: &    ! ARRAYS for counting
&       pdiga (:,:,:,:), pdigaa(:,:,:,:), pdigas(:,:,:,:), &
&       pdsga (:,:,:),   pdsgaa(:,:,:),   pdsgas(:,:,:),   &
&       pdigb (:,:,:,:), pdigba(:,:,:,:), pdigbs(:,:,:,:), &
&       pdigs (:,:,:,:), pdigsa(:,:,:,:), pdigss(:,:,:,:), &
&       pdsgs (:,:),     pdsgsa(:,:),     pdsgss(:,:)


!  INTEGER,PARAMETER :: ldunit=90 ! unit for rerun tendency data
  INTEGER,PARAMETER :: ldunit=33 ! unit for rerun tendency data

  ! insert the new codes here, if GRIB1 works fine into the postprecessing
  !
  ! codes for diagnostic terms
  ! vorticity
  INTEGER, PARAMETER :: ndvor=9 ! number of terms in equation
  INTEGER, PARAMETER :: ncvor(ndvor)=(/41,42,43,44,45,46,47,48,49/)
  !integer, parameter :: ncvor(ndvor)=(/155,156,157,158,159,160,161,162/)
  !     1     041/155 dynamic tendencies horizontal part
  !     2     042/156 vertical advection
  !     3     043/157 vertical dissipation due to VDIFF 
  !     4     044/158 vertical dissipation due to GWDRAG
  !     5     045/159 vertical dissipation due to CUCALL
  !     6     046/160 timefilter
  !     7     047/161 semi-implicit time integration part
  !     8     048/162 horizontal diffusion
  !     9     049/163 total tendency, from step to step only (t+1)-(t-1)
  !
  ! divergence
  INTEGER, PARAMETER :: nddiv=9 ! number of terms in equation
  INTEGER, PARAMETER :: ncdiv(nddiv)=(/61,62,63,64,65,66,67,68,69/)
  !integer, parameter :: ncdiv(nddiv)=(/170,171,172,173,174,175,176,177,178/)
  !     1     061/170 dynamic tendencies horizontal part
  !     2     062/171 vertical advection
  !     3     063/172 vertical dissipation due to VDIFF 
  !     4     064/173 vertical dissipation due to GWDRAG
  !     5     065/174 vertical dissipation due to CUCALL
  !     6     066/175 timefilter
  !     7     067/176 semi-implicit time integration part
  !     8     068/177 horizontal diffusion
  !     9     069/178 total tendency, from step to step only (t+1)-(t-1)
  !
  ! temperature
  INTEGER, PARAMETER :: ndtem=14 ! number of terms in equation
  INTEGER, PARAMETER :: nctem(ndtem)=(/81,82,83,84,85,86,87,88,89,90,91,92,93,94/)
  !integer, parameter :: nctem(ndtem)=(/185,186,187,188,189,190,191,192,193,194,195,196,197,198/)
  !     1     081/185 dynamic tendencies horizontal part
  !     2     082/186 vertical advection
  !     3     083/187 energy conversion
  !     4     084/188 radiation processes
  !     5     085/189 vertical dissipation due to VDIFF 
  !     6     086/190 vertical dissipation due to GWDRAG
  !     7     087/191 cond, prec, eva due to CUCALL
  !     8     088/192 cond, prec, eva due to CUCOND
  !     9     089/193 timefilter
  !     10    090/194 semi-implicit time integration part
  !     11    091/195 horizontal diffusion
  !     12    092/196 longwave radiation
  !     13    093/197 shortwave radiation
  !     14    094/198 total tendency, from step to step only (t+1)-(t-1)
  !
  ! surface pressure
  INTEGER, PARAMETER :: ndprs=4 ! number of terms in equation
  INTEGER, PARAMETER :: ncprs(ndprs)=(/101,102,103,104/)
  !integer, parameter :: ncprs(ndprs)=(/201,202,203,204/)
  !     1     101/201 divergence term
  !     2     102/202 timefilter
  !     3     103/203 semi-implicit time integration part
  !     4     104/204 total tendency, from step to step only (t+1)-(t-1)
  !
  INTEGER, PARAMETER :: ncprl = 100
  !integer, parameter :: ncprl = 200
  !     1     100/200 horizontal divergence in the layers

  INTEGER, PARAMETER :: IDIAG_INIT = 1
  INTEGER, PARAMETER :: IDIAG_INI_GBUF = 40
  INTEGER, PARAMETER :: IDIAG_INI_PREV = 100
  INTEGER, PARAMETER :: IDIAG_RER_WR = -1

  CHARACTER(256) :: diag_mess

CONTAINS


!======================================================================
  SUBROUTINE DIAG_Init(itype)

    INTEGER, INTENT(in) :: itype       ! function separator

    SELECT CASE(itype)
    CASE(IDIAG_INIT)             ! startup initialization

       ldinit   = .TRUE.

       ipdoaccu = 0
       ipdaccu  = 0

       CALL message('DIAG_Init', &
&            ' Tendency Diagnostics E4AMIP2-V3.22-p2 23-May-2000 (kirchner@dkrz.de)')
       ! -----------------------------------------------------------
       ! allocate spectral space arrays
       !
       !   diagnostic arrays all accumulated
       !
       !     g .... accumulated in grid point space
       !     gs ... parts from gridpoint space, accumulated in spectral space
       !     s .... accumulated in spectral space
       !
       !     spectral terms
       !
       !   PDVOR  vorticity equation
       ALLOCATE (pdvor (nlev,2,nsp,ndvor)); pdvor(:,:,:,:) = 0.0
       ALLOCATE (pdovor(nlev,2,nsp,2)); pdovor(:,:,:,:) = 0.0
       !
       ! g    1  horizontal advection, pressure gradient, coriolisterm (DYN)
       ! g    2  vertical advection (DYN)
       ! g    3  vertical diffusion due to impuls (VDIFF)
       ! g    4  gravity wave drag (GWDRAG)
       ! g    5  moisture mass flux (CUCALL)
       ! g    6  timefilter
       ! gs   7  semi-implicit part of time integration
       ! s    8  horizontal diffusion
       ! s    9  total step to step
       !
       !   PDDIV  divergence equation
       ALLOCATE (pddiv (nlev,2,nsp,nddiv)); pddiv(:,:,:,:) = 0.0
       ALLOCATE (pdodiv(nlev,2,nsp,2)); pdodiv(:,:,:,:) = 0.0
       !
       ! g    1  horizontal advection, coriolisterm, pressure gradient, G-term
       ! g    2  vertical advection
       ! g    3  vertical diffusion due to imuls (VDIFF)
       ! g    4  gravity wave drag (GWDRAG)
       ! g    5  moisture mass flux (CUCALL)
       ! g    6  timefilter
       ! gs   7  semi-implicit part of time integration
       ! s    8  horizontal diffusion
       ! s    9  total step to step
       !
       !   PDTEM  temperature
       ALLOCATE (pdtem (nlev,2,nsp,ndtem)); pdtem(:,:,:,:) = 0.0
       ALLOCATE (pdotep(nlevp1,2,nsp,2)); pdotep(:,:,:,:) = 0.0
       !
       ! g    1  horizontal advection
       ! g    2  vertical advection
       ! g    3  energy conversion
       ! g    4  radiation (RADHEAT)
       ! g    5  vertical diffusion due to turbulence (VDIFF)
       ! g    6  gravity wave drag (GWDRAG)
       ! g    7  convection (CUCALL)
       ! g    8  large scale cloud processes (COND)
       ! g    9  timefilter
       ! gs   10 semi-implicit part of time integration
       ! s    11 horizontal diffusion
       ! g    12 longwave radiation
       ! g    13 shortwave radiation
       ! s    14 total step to step
       !
       !   PDPRL  divergence for each layer
       ALLOCATE (pdprl(nlev,2,nsp)); pdprl(:,:,:) = 0.0
       ! g       convergence in each layer
       !
       !   PDPRS  pressure equation
       ALLOCATE (pdprs     (2,nsp,ndprs)); pdprs(:,:,:) = 0.0
       !
       ! g    1  vertical integrated convergence
       ! g    2  timefilter
       ! s    3  semi-implicit part of time integration
       ! s    4  total step to step
       !
       ! -----------------------------------------------------------
       ! allocate gridpoint space arrays
       !
       !   PTFH1   memory of timefilter
       ! updated in TF1
       ALLOCATE (ptfh1 (nlp2,nlev,NO_PTFH,ngl)); ptfh1(:,:,:,:) = 0.0
       !           1 ... vorticity
       !           2 ... divergence
       !           3 ... temperature
       !   PTFH2   pressure
       ! updated in TF1
       ALLOCATE (ptfh2 (nlp2,       ngl)); ptfh2(:,:) = 0.0
       !
       ! -----------------------------------------------------------
       ! workspace for accumulation of terms in grid point space
       !
       !     Index
       !     ............. vorticity and divergenc equation
       !     VOM = Fu, VOL = Fv
       !
       !     Fu (VOM parts in GPC)
       ! DYN    1   coriolisterm and pressure tendency of Fu
       ! DYN    2   vertical advection of Fu
       ! PHYSC  3   diffusion due to impuls Pu VDIFF
       ! PHYSC  4   diffusion due to gravity wave drag Pu GWDRAG
       ! PHYSC  5   diffusion due to mass flux Pu CUCALL
       !
       !     Fv (VOL parts in GPC)
       ! DYN    6   coriolisterm and pressure tendency of Fv
       ! DYN    7   vertical advection of Fv
       ! PHYSC  8   diffusion due to impuls Pv VDIFF
       ! PHYSC  9   diffusion due to gravity wave drag Pv GWDRAG
       ! PHYSC  10  diffusion due to mass flux Pv CUCALL
       !
       ! DYN    11  potential and kinetic energy term G
       !
       !     ............. temperature tendency equation
       ! DYN    12  horizontal advection term
       ! DYN    13  vertical advection term
       ! DYN    14  energy conversion
       ! PHYSC  15  RADHEAT radiation tendency
       ! PHYSC  16  VDIFF turbulence
       ! PHYSC  17  GWDRAG gravity wave drag
       ! PHYSC  18  CUCALL convective condensation
       ! PHYSC  19  COND large+scale condensation
       !
       !     ............. pressure tendency equation
       ! DYN    20 level dependend divergence part, 
       !
       ! TF2    21  timefilter of vorticity
       ! TF2    22  timefilter of divergence
       ! TF2    23  timefilter of temperature
       !
       ! RADHEAT 24 longwave radiation
       ! RADHEAT 25 shortwave radiation
       !
       ALLOCATE (pdiga    (nlp2,nlev,NO_PDIGA,ngl)); pdiga(:,:,:,:) = 0.0
       ALLOCATE (pdigaa (nlev,2,nmp1,NO_PDIGA))
       ALLOCATE (pdigas (nlev,2,nmp1,NO_PDIGA))
       !
       ! one level arrays of surface pressure
       !      1     last level total integral
       !      2     timefilter for pressure
       !
       ALLOCATE (pdsga    (nlp2,NO_PDSGA,ngl)); pdsga(:,:,:) = 0.0
       ALLOCATE (pdsgaa (2,nmp1,NO_PDSGA))
       ALLOCATE (pdsgas (2,nmp1,NO_PDSGA))
       !
       !     array for d/dlambda derivatives calculated in fourierspace
       !     corresponding to pdiga(*,1...10,*)
       !
       ALLOCATE (pdigb    (nlp2,nlev,NO_PDIGB,ngl))
       ALLOCATE (pdigba (nlev,2,nmp1,NO_PDIGB))
       ALLOCATE (pdigbs (nlev,2,nmp1,NO_PDIGB))
       !
       !     local memory buffer for accumulation of semi-implicit parts
       !     fraction solved in grid point space
       !
       !    1 vorticity and implicit and explicit part (L)
       !    2 divergence
       !    3 temperatur and pressure
       !    4 implicit part of vorticity (M) used in Fourierspace
       ALLOCATE (pdigs    (nlp2,nlev,NO_PDIGS,ngl))
       ALLOCATE (pdigsa (nlev,2,nmp1,NO_PDIGS))
       ALLOCATE (pdigss (nlev,2,nmp1,NO_PDIGS))
       !
       !     semi-implicit part of log surface pressure
       ALLOCATE (pdsgs    (nlp2,ngl))
       ALLOCATE (pdsgsa (2,nmp1))
       ALLOCATE (pdsgss (2,nmp1))
       !
       lpdo1 = .FALSE.
       lpdo2 = .FALSE.
       !
       ! special initialization for restart
       IF (lres) CALL DIAG_Rerun('R')     ! read old rerun data

    CASE(IDIAG_RER_WR)        ! get memory free

       CALL DIAG_Rerun('W')     ! store rerun data

       DEALLOCATE (pdvor)
       DEALLOCATE (pddiv)
       DEALLOCATE (pdtem)
       DEALLOCATE (pdprs)
       DEALLOCATE (pdprl)

       DEALLOCATE (pdovor)
       DEALLOCATE (pdodiv)
       DEALLOCATE (pdotep)

       DEALLOCATE (ptfh1)
       DEALLOCATE (ptfh2)

       DEALLOCATE (pdiga)
       DEALLOCATE (pdigaa)
       DEALLOCATE (pdigas)

       DEALLOCATE (pdsga)
       DEALLOCATE (pdsgaa)
       DEALLOCATE (pdsgas)

       DEALLOCATE (pdigb)
       DEALLOCATE (pdigba)
       DEALLOCATE (pdigbs)

       DEALLOCATE (pdigs)
       DEALLOCATE (pdigsa)
       DEALLOCATE (pdigss)

       DEALLOCATE (pdsgs)
       DEALLOCATE (pdsgsa)
       DEALLOCATE (pdsgss)

    CASE (IDIAG_INI_GBUF)
       ! reset local accumulation buffer
       pdigb (:,:,:,:) = 0.0
       pdigs (:,:,:,:) = 0.0
       pdsgs (:,:)     = 0.0
       
    CASE (IDIAG_INI_PREV)

       ! calculate tendency
       IF (lpdo1) THEN
          ! memory is filled, the tendency can be calculated as [X(t+1) - X(t-1)]
          pdvor(:,:,:,ndvor) = pdvor(:,:,:,ndvor)+(svo(:,:,:)      - pdovor(:,:,:,1))
          pddiv(:,:,:,nddiv) = pddiv(:,:,:,nddiv)+(sd (:,:,:)      - pdodiv(:,:,:,1))
          pdtem(:,:,:,ndtem) = pdtem(:,:,:,ndtem)+(stp(1:nlev,:,:) - pdotep(1:nlev,:,:,1))
          pdprs(  :,:,ndprs) = pdprs(  :,:,ndprs)+(stp(nlevp1,:,:) - pdotep(nlevp1,:,:,1))
          ipdoaccu = ipdoaccu + 1
       END IF
       IF (lpdo2) THEN
          ! move (t) into (t-1) memory
          pdovor(:,:,:,1) = pdovor(:,:,:,2)
          pdodiv(:,:,:,1) = pdodiv(:,:,:,2)
          pdotep(:,:,:,1) = pdotep(:,:,:,2)
          lpdo1 = .TRUE.
       END IF
       ! store (t+1) for next integration loop
       ! it is the unfiltered valu ein the spectral space
       ! in TF1 the same field is available at grid points
       pdovor(:,:,:,2) = svo(:,:,:)
       pdodiv(:,:,:,2) = sd (:,:,:)
       pdotep(:,:,:,2) = stp(:,:,:)
       lpdo2 = .TRUE.

       ipdaccu = ipdaccu + 1

    CASE default
      WRITE (diag_mess,*) 'type not implemented, ITYPE= ',itype
      CALL finish('DIAG_Init',diag_mess)

    END SELECT
    
  END SUBROUTINE DIAG_Init



!======================================================================
  SUBROUTINE DIAG_SpecTrans    ! -->>> insert in SCAN1SL after CALL SI2
    ! second part of spectral transform

    INTEGER :: jlat

    DO jlat  = 1, nhgl                         ! global  index north -> south

       CALL DIAG_sym1(jlat)
       CALL legmod(jlat)
       CALL DIAG_ltd

    END DO

  END SUBROUTINE DIAG_SpecTrans



!======================================================================
  SUBROUTINE DIAG_fftd    ! -->>> insert in SCAN1SL after CALL FFTD
    ! transform the diagnostic arrays into fourierspace

    INTEGER :: irow                 ! actual latitude index
    INTEGER :: inc, isign

    REAL :: zinp(nlp2*nlev*NO_PDIGA) ! work space for transformation
    REAL :: zwrk(nlp2*nlev*NO_PDIGA) ! work space for transformation

    INTEGER :: iamount

#if (defined CRAY) || (defined SX)
    EXTERNAL util_reshape
#else
    INTRINSIC RESHAPE
#endif

    inc   = 1
    isign = -1
    DO irow = 1,ngl

       ! transform 3-dim fields semi-implicit parts counted in grid space
       iamount = nlp2*nlev*NO_PDIGS_FRAC
#if (defined CRAY) || (defined SX)
       CALL util_reshape(zinp,pdigs(1,1,1,irow),iamount)
#else
       zinp(:) = RESHAPE(     pdigs(:,:,1:NO_PDIGS_FRAC,irow),(/iamount/))
#endif
       CALL fft991cy(zinp,zwrk,trig,nfax,inc,nlp2,nlon,nlev*NO_PDIGS_FRAC,isign)
       iamount = nlp2*nlev*NO_PDIGS_FRAC
#if (defined CRAY) || (defined SX)
       CALL util_reshape(pdigs(1,1,1,irow),zinp,iamount)
#else
       pdigs(:,:,1:NO_PDIGS_FRAC,irow) = RESHAPE(zinp(:),(/nlp2,nlev,NO_PDIGS_FRAC/))
#endif

       ! transform 2-dim fields
       iamount = nlp2
#if (defined CRAY) || (defined SX)
       CALL util_reshape(zinp,pdsgs(1,irow),iamount)
#else
       zinp(:) = RESHAPE(     pdsgs(:,irow),(/iamount/))
#endif
       CALL fft991cy(zinp,zwrk,trig,nfax,inc,nlp2,nlon,1,isign)
       iamount = nlp2
#if (defined CRAY) || (defined SX)
       CALL util_reshape(pdsgs(1    ,irow),zinp,iamount)
#else
       pdsgs(:       ,irow) = RESHAPE(zinp(:),(/nlp2/))
#endif

       IF (lptime) THEN       ! transform other terms during output time step

          ! 3-dim fields
          iamount = nlp2*nlev*NO_PDIGA
#if (defined CRAY) || (defined SX)
          CALL util_reshape(zinp,pdiga(1,1,1,irow),iamount)
#else
          zinp(:) = RESHAPE(pdiga(:,:,:,irow),(/iamount/))
#endif
          CALL fft991cy(zinp,zwrk,trig,nfax,inc,nlp2,nlon,nlev*NO_PDIGA,isign)
          iamount = nlp2*nlev*NO_PDIGA
#if (defined CRAY) || (defined SX)
          CALL util_reshape(pdiga(1,1,1,irow),zinp,iamount)
#else
          pdiga(:,:,:,irow) = RESHAPE(zinp(:),(/nlp2,nlev,NO_PDIGA/))
#endif
          ! 2-dim fields
          iamount = nlp2*NO_PDSGA
#if (defined CRAY) || (defined SX)
          CALL util_reshape(zinp,pdsga(1,1  ,irow),iamount)
#else
          zinp(:) = RESHAPE(pdsga(:,:  ,irow),(/iamount/))
#endif
          CALL fft991cy(zinp,zwrk,trig,nfax,inc,nlp2,nlon,NO_PDSGA,isign)
          iamount = nlp2*NO_PDSGA
#if (defined CRAY) || (defined SX)
          CALL util_reshape(pdsga(1,1  ,irow),zinp,iamount)
#else
          pdsga(:,:  ,irow) = RESHAPE(zinp(:),(/nlp2,NO_PDSGA/))
#endif
      ENDIF

    ENDDO

  END SUBROUTINE DIAG_fftd


!======================================================================
  SUBROUTINE DIAG_sym1(ihrow)    ! -->>> insert in SCAN1SL after CALL SYM1
    ! separate the diagnostics arrays

    INTEGER :: jl, jm, ihrow, irow_n, irow_s

    !     even and odd components
    irow_n = ihrow
    irow_s = ngl + 1 - ihrow

    DO jm = 1,nmp1
       DO jl = 1,nlev
          pdigss(jl,:,jm,:) = 0.5*(pdigs(2*jm-1:2*jm,jl,:,irow_n) + pdigs(2*jm-1:2*jm,jl,:,irow_s))
          pdigsa(jl,:,jm,:) = 0.5*(pdigs(2*jm-1:2*jm,jl,:,irow_n) - pdigs(2*jm-1:2*jm,jl,:,irow_s))
       ENDDO
       pdsgss(:,jm) = 0.5*(pdsgs(2*jm-1:2*jm,irow_n) + pdsgs(2*jm-1:2*jm,irow_s))
       pdsgsa(:,jm) = 0.5*(pdsgs(2*jm-1:2*jm,irow_n) - pdsgs(2*jm-1:2*jm,irow_s))
    ENDDO
    IF (lptime) THEN
       DO jm = 1,nmp1
          DO jl = 1,nlev
             pdigas(jl,:,jm,:) = 0.5*(pdiga(2*jm-1:2*jm,jl,:,irow_n) + pdiga(2*jm-1:2*jm,jl,:,irow_s))
             pdigaa(jl,:,jm,:) = 0.5*(pdiga(2*jm-1:2*jm,jl,:,irow_n) - pdiga(2*jm-1:2*jm,jl,:,irow_s))
             pdigbs(jl,:,jm,:) = 0.5*(pdigb(2*jm-1:2*jm,jl,:,irow_n) + pdigb(2*jm-1:2*jm,jl,:,irow_s))
             pdigba(jl,:,jm,:) = 0.5*(pdigb(2*jm-1:2*jm,jl,:,irow_n) - pdigb(2*jm-1:2*jm,jl,:,irow_s))
          ENDDO
          pdsgas(:,jm,:) = 0.5*(pdsga(2*jm-1:2*jm,:,irow_n) + pdsga(2*jm-1:2*jm,:,irow_s))
          pdsgaa(:,jm,:) = 0.5*(pdsga(2*jm-1:2*jm,:,irow_n) - pdsga(2*jm-1:2*jm,:,irow_s))
       ENDDO
    ENDIF

  END SUBROUTINE DIAG_sym1


!======================================================================
  SUBROUTINE DIAG_ltd    ! -->>> insert in SCAN1SL after CALL LTD
    ! perform legendre transform for diagnostic arrays

    INTEGER :: jm, ims, ins, is, jn, jh, iu

    DO jh = 1,2 ! 1: north, 2:south
       iu = 2-jh
       DO jm = 1,nmp1
          ims = nmp(jm)-iu
          ins = nnp(jm)+iu
          DO jn=2,ins,2
             is = ims+jn
             IF (jh == 1) THEN   !     calculations for northern hemisphere
                !
                ! semi-implicit parts of vorticity
                pdvor(:,:,is, 7) = pdvor(:,:,is, 7) + pdigss(:,:,jm,1)*pnmd(is) &
&                                                   - pdigsa(:,:,jm,4)*anmd(is)

                ! explicit part of divergence, temperature and pressure
                pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigss(:,:,jm,2)*rnmd(is)
                pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigss(:,:,jm,3)*pnmd(is)
                pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgss  (:,jm  )*pnmd(is)
                !
                IF (lptime) THEN
                   ! dynamic and physical tendencies
                   ! vorticity
                   pdvor(:,:,is,1:5) = pdvor (:,:,is,1:5) + pdigbs(:,:,jm,6:10)*pnmd(is) &
&                                                         + pdigaa(:,:,jm,1:5 )*anmd(is)
                   pdvor(:,:,is,  6) = pdvor (:,:,is,  6) + pdigas(:,:,jm,21  )*pnmd(is)
                   !
                   ! divergence
                   pddiv(:,:,is,1:5) = pddiv (:,:,is,1:5) + pdigbs(:,:,jm,1:5 )*pnmd(is) &
&                                                         - pdigaa(:,:,jm,6:10)*anmd(is)
                   pddiv(:,:,is,  6) = pddiv (:,:,is,  6) + pdigas(:,:,jm,22  )*pnmd(is)
                   ! laplacian operation with G-term
                   pddiv(:,:,is,  1) = pddiv (:,:,is,  1) + pdigas(:,:,jm,11  )*rnmd(is)
                   !
                   ! temperature
                   pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigas(:,:,jm,12:19)*pnmd(is)
                   pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigas(:,:,jm,23   )*pnmd(is)
                   pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigas(:,:,jm,24:25)*pnmd(is)
                   !
                   ! pressure
                   pdprl(:,:,is)     = pdprl(:,:,is)      + pdigas(:,:,jm,20)*pnmd(is)
                   pdprs  (:,is,1:2) = pdprs  (:,is,1:2)  + pdsgas  (:,jm, :)*pnmd(is)
                ENDIF
                !
             ELSE                ! calculations for southern hemisphere
                !
                ! semi-implicit parts of vorticity
                pdvor(:,:,is, 7) = pdvor (:,:,is, 7) + pdigsa(:,:,jm,1)*pnmd(is) &
&                                                    - pdigss(:,:,jm,4)*anmd(is)
                ! explicit part divergence, temperature and pressure
                pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigsa(:,:,jm,2)*rnmd(is)
                pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigsa(:,:,jm,3)*pnmd(is)
                pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgsa  (:,jm  )*pnmd(is)
                !
                IF (lptime) THEN
                   ! dynamic and physical tendencies
                   ! vorticity
                   pdvor(:,:,is,1:5) = pdvor(:,:,is,1:5) + pdigba(:,:,jm,6:10)*pnmd(is) &
&                                                        + pdigas(:,:,jm,1:5 )*anmd(is)
                   pdvor(:,:,is,  6) = pdvor(:,:,is,  6) + pdigaa(:,:,jm,21  )*pnmd(is)
                   !
                   ! divergence
                   pddiv(:,:,is,1:5) = pddiv(:,:,is,1:5) + pdigba(:,:,jm,1:5 )*pnmd(is) &
&                                                        - pdigas(:,:,jm,6:10)*anmd(is)
                   pddiv(:,:,is,  6) = pddiv(:,:,is,  6) + pdigaa(:,:,jm,22  )*pnmd(is)
                   pddiv(:,:,is,  1) = pddiv(:,:,is,  1) + pdigaa(:,:,jm,11  )*rnmd(is)
                   !
                   ! temperature
                   pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigaa(:,:,jm,12:19)*pnmd(is)
                   pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigaa(:,:,jm,23   )*pnmd(is)
                   pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigaa(:,:,jm,24:25)*pnmd(is)
                   !
                   ! pressure
                   pdprl(:,:,is)     = pdprl(:,:,is)     + pdigaa(:,:,jm,20)*pnmd(is)
                   pdprs  (:,is,1:2) = pdprs  (:,is,1:2) + pdsgaa  (:,jm, :)*pnmd(is)
                ENDIF
                !
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE DIAG_ltd

!======================================================================
  SUBROUTINE DIAG_Write
    REAL    :: tfact1, tfact2, tfact
    INTEGER :: iterm, jlev, ierr, iret
    REAL    :: worksp(2,nsp)

    EXTERNAL codegb5, pbwrite

    ! proportionality factor
    ! the data accumulated over different intervalls
    ! all tendencies are stored as VALUE/sec
    IF (ipdaccu > 0) THEN
       tfact1 = 0.5/(REAL(ipdaccu)*dtime)     ! all tendency terms
    ELSE
       tfact1 = 0.0
    END IF
    IF (ipdoaccu > 0) THEN
       tfact2 = 0.5/(REAL(ipdoaccu)*dtime) ! the total tendency
       ! perform an internal check
       CALL DIAG_Check
    ELSE
       tfact2 = 0.0
    END IF
    !tfact = 0.5/time_inc_sec(nptime,TIME_INC_DAYS)

    ! Finish definition GRIB block 1

!    ksec1(1) = nudging_table
    ksec1(1) = local_table

    CALL set_output_time

    ksec1(10) = year
    ksec1(11) = month
    ksec1(12) = day
    ksec1(13) = hour
    ksec1(14) = minute
    ksec1(21) = century

    ksec4(1) = n2sp

    WRITE (diag_mess,*) 'century year mon day hour min ', &
&         century,year,month,day,hour,minute
    CALL message('DIAG_Write',diag_mess)

    level_type = 109
    ksec1(7) = level_type


    DO iterm = 1,ndvor    ! 1. vorticity equation
       tfact = tfact1
       IF (iterm==ndvor) tfact=tfact2
       ksec1(6) = ncvor(iterm)
       DO jlev = 1 , nlev
          level_p2 = jlev
          ksec1(8) = level_p2
          worksp(:,:) = pdvor(jlev,:,:,iterm)*tfact
          CALL DiagWriteSP
       ENDDO
    ENDDO
    pdvor(:,:,:,:) = 0.0

    DO iterm = 1,nddiv    ! 2. divergence equation
       tfact = tfact1
       IF (iterm==nddiv) tfact=tfact2
       ksec1(6) = ncdiv(iterm)
       DO jlev = 1 , nlev
          level_p2 = jlev
          ksec1(8) = level_p2
          worksp(:,:) = pddiv(jlev,:,:,iterm)*tfact
          CALL DiagWriteSP
       ENDDO
    ENDDO
    pddiv(:,:,:,:) = 0.0

    DO iterm = 1,ndtem    ! 3. temperature equation
       tfact = tfact1
       IF (iterm==ndtem) tfact=tfact2
       ksec1(6) = nctem(iterm)
       DO jlev = 1 , nlev
          level_p2 = jlev
          ksec1(8) = level_p2
          worksp(:,:) = pdtem(jlev,:,:,iterm)*tfact
          CALL DiagWriteSP
       ENDDO
    ENDDO
    pdtem(:,:,:,:) = 0.0


    ksec1(6) = ncprl    ! 4. pressure equation
    DO jlev = 1 , nlev
       level_p2 = jlev
       ksec1(8) = level_p2
       worksp(:,:) = pdprl(jlev,:,:)*tfact1
       CALL DiagWriteSP
    ENDDO
    pdprl(:,:,:) = 0.0

    level_type = 1    ! surface data
    ksec1(7) = level_type
    level_p2   = 0
    ksec1(8)   = level_p2
    DO iterm = 1,ndprs
       tfact = tfact1
       IF (iterm==ndprs) tfact=tfact2
       ksec1(6) = ncprs(iterm)
       worksp(:,:) = pdprs(:,:,iterm)*tfact
       CALL DiagWriteSP
    ENDDO
    pdprs(:,:,:) = 0.0

    ! reset accumulated grid point arrays
    pdiga(:,:,:,:) = 0.0
    pdsga(:,:,:)   = 0.0

    ! reset accumulation counter
    ipdoaccu = 0
    ipdaccu  = 0

  CONTAINS

    SUBROUTINE DiagWriteSP
#ifdef EMOS
      CALL gribex (ksec0, ksec1, ksec2_sp, psec2, ksec3, psec3, &
&                   ksec4, worksp, n2sp, kgrib, kleng, kword, 'C', ierr)
#else
      CALL codegb5(worksp,n2sp,16,nbit,ksec1,ksec2_sp,vct,2*nvclev,&
&               kgrib,klengo,kword,0,ierr)
#endif
      CALL pbwrite(nunitdf,kgrib,kword*iobyte,iret)
      IF (iret /= kword*iobyte) &
&           CALL finish('DIAG_Write','I/O error on output - disk full?')

    END SUBROUTINE DiagWriteSP

  END SUBROUTINE DIAG_Write


!======================================================================
  SUBROUTINE DIAG_Rerun(ctyp)

    INTEGER   :: ilevels, ilons2, ilats, insp, jl, flag
    CHARACTER :: ctyp*1

    REAL, ALLOCATABLE ::  rvct(:)

    SELECT CASE(ctyp)
    CASE('r','R')      ! read rerun data from file
       REWIND(ldunit,err=100)

       READ(ldunit,err=100,END=100) ilevels, ilons2, ilats, insp

       IF (ilevels /= nlev) THEN
          CALL message('DIAG_Rerun','Warning: number of levels inconsistent')
       ELSE IF (ilons2 /= nlp2) THEN
          CALL message('DIAG_Rerun','Warning: number of longitudes inconsistent')
       ELSE IF (ilats /= ngl) THEN
          CALL message('DIAG_Rerun','Warning: number of latitudes inconsistent')
       ELSE IF (insp /= nsp) THEN
          CALL message('DIAG_Rerun','Warning: number of spectral coefficients inconsistent')
       ELSE
          ! get output accumulation numbers
          READ(ldunit) ipdoaccu, ipdaccu
          ! get coordiante system from file
          ALLOCATE (rvct(2*nlev+2))
          READ(ldunit,err=100) rvct
          flag = 0
          DO jl = 1,2*nlev+2
             IF (ABS(rvct(jl)-vct(jl)) > 0.) flag = flag + 1
          ENDDO
          DEALLOCATE (rvct)
          IF (flag /= 0) THEN
             CALL message('DIAG_Rerun','Warning: inconsistent vertical model structure')
          ELSE
             READ(ldunit) pdvor
             READ(ldunit) pddiv
             READ(ldunit) pdtem
             READ(ldunit) pdprl
             READ(ldunit) pdprs
             READ(ldunit) ptfh1
             READ(ldunit) ptfh2
             READ(ldunit) pdovor
             READ(ldunit) pdodiv
             READ(ldunit) pdotep
             lpdo1 = .TRUE.
             lpdo2 = .TRUE.
             ldinit = .FALSE.    ! no additional initialization for timefilter
             CALL message('DIAG_Rerun',' ... read restart data')
          ENDIF
       ENDIF

    CASE('w','W')       ! write restart data to file
       REWIND(ldunit)
       WRITE(ldunit) nlev,nlp2,ngl,nsp
       WRITE(ldunit) ipdoaccu, ipdaccu
       WRITE(ldunit) vct(1:2*nlev+2)
       WRITE(ldunit) pdvor
       WRITE(ldunit) pddiv
       WRITE(ldunit) pdtem
       WRITE(ldunit) pdprl
       WRITE(ldunit) pdprs
       WRITE(ldunit) ptfh1
       WRITE(ldunit) ptfh2
       WRITE(ldunit) pdovor
       WRITE(ldunit) pdodiv
       WRITE(ldunit) pdotep

       CALL message('DIAG_Rerun',' ... wrote restart data')

    END SELECT

100 CONTINUE

  END SUBROUTINE DIAG_Rerun


!======================================================================
  SUBROUTINE DIAG_Check
    REAL :: sumt, sumv, sumd, sump, val1, val2,&
&         sumt1, sumt2, sumd1, sumd2, sumv1, sumv2, sump1, sump2, &
&         relv, reld, relt, relp
    INTEGER :: jl, jk, js, jv
    ! compare new minus old with tendency

    sumv = 0.0; sumv1 = 0.0; sumv2 = 0.0; relv = 0.0
    sumd = 0.0; sumd1 = 0.0; sumd2 = 0.0; reld = 0.0
    sumt = 0.0; sumt1 = 0.0; sumt2 = 0.0; relt = 0.0
    sump = 0.0; sump1 = 0.0; sump2 = 0.0; relp = 0.0

    DO js=1,nsp
       DO jk=1,2
          DO jl=1,nlev
             IF (js > 1) THEN
                val1  = pdvor(jl,jk,js,ndvor)
                val2 = SUM(pdvor(jl,jk,js,1:ndvor-1))
                sumv  = sumv  + val1*val2
                sumv1 = sumv1 + val1*val1
                sumv2 = sumv2 + val2*val2

                val1  = pddiv(jl,jk,js,nddiv)
                val2 = SUM(pddiv(jl,jk,js,1:nddiv-1))
                sumd  = sumd  + val1*val2
                sumd1 = sumd1 + val1*val1
                sumd2 = sumd2 + val2*val2
             END IF

             val1  = pdtem(jl,jk,js,ndtem)
             val2 = SUM(pdtem(jl,jk,js,1:ndtem-3))
             sumt  = sumt  + val1*val2
             sumt1 = sumt1 + val1*val1
             sumt2 = sumt2 + val2*val2

          END DO

          IF (js > 1) THEN
             val1 = pdprs(jk,js,ndprs)
             val2 = SUM(pdprs(jk,js,1:ndprs-1))
             sump  = sump  + val1*val2
             sump1 = sump1 + val1*val1
             sump2 = sump2 + val2*val2
          END IF

       END DO

    END DO

    IF (sumv1 > 0.0 .AND. sumv2 > 0.0) relv = sumv/SQRT(sumv1*sumv2)
    IF (sumd1 > 0.0 .AND. sumd2 > 0.0) reld = sumd/SQRT(sumd1*sumd2)
    IF (sumt1 > 0.0 .AND. sumt2 > 0.0) relt = sumt/SQRT(sumt1*sumt2)
    IF (sump1 > 0.0 .AND. sump2 > 0.0) relp = sump/SQRT(sump1*sump2)

    WRITE(diag_mess,'(4(a,f8.5))') 'VOR ',relv,' DIV ',reld,' TEM ',relt,' PRS ',relp
    CALL message('DIAG_Check',diag_mess)

  END SUBROUTINE DIAG_Check

END MODULE mo_diag_tendency
