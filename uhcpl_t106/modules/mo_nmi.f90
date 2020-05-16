#ifdef CRAY
#define dgemm sgemm
#define dgeev sgeev
#define dgesv sgesv
#define dsyev ssyev
#endif

MODULE mo_nmi

  ! Authors:

  ! W. Wergen, ECMWF, July 1982 - December 1982
  ! W. Heckley, ECMWF, July 1986 - November 1988
  ! C. Temperton, ECMWF, January 1989 - May 1991
  ! L. Kornblueh, MPI, June 1998
  ! I. Kirchner, MPI, August 1998 - January 1999
  ! I. Kirchner, MPI, October 1999, March 2000

  ! ******************** Interface to ECHAM **********************************
  !
  ! *ECHAM*       *NMI*
  !
  ! CONTROL ----> NMI_Init
  !                 +---> NMI_Normal_Modes
  !                         +---> NMI_Horizontal_Modes
  !
  ! STEPON -----> NMI_Make(5)
  ! Nudging ----> NMI_Make(3)
  !                 +---> NMI_LoadData
  !                 |       +---> NMI_Diagnose
  !                 +---> NMI_Filter

  USE mo_doctor,        ONLY: nout, nin
  USE mo_start_dataset, ONLY: nstep, lres, nstart
  USE mo_constants,     ONLY: a, omega, api, g, rd
  USE mo_exception,     ONLY: finish, message
  USE mo_rad_switches,  ONLY: ldiur

  USE mo_nudging_buffer,ONLY: sdobs3, svobs3, stobs3, sdobs, svobs, stobs,&
&       sdfast_a,    svfast_a,    stfast_a, &
&       sdfast_b,    svfast_b,    stfast_b, &
&       sdfast_accu, svfast_accu, stfast_accu, &
&       ifast_accu, lfill_a, lfill_b

  USE mo_memory_sp,     ONLY: svo, sd, stp

  USE mo_semi_impl,     ONLY: apr, betadt, tr
  USE mo_hyb,           ONLY: bb, delpr, rdelpr, nlevm1
  ! rdtr /= rd*tr, see INHYSI for details

  USE mo_truncation,    ONLY: nnp, nmp
  USE mo_control,       ONLY: &
&         nm, nmp1, nn, nnp1, nk, nlev, nlevp1, nsp, n2sp, &
&         lnudge, lnmi, nlp2, ngl, vct, nvclev, dtime, twodt, &
&         lptime

  IMPLICIT NONE

  PUBLIC :: NMI_Init
  PUBLIC :: NMI_Make

  PRIVATE :: NMI_Normal_Modes
  PRIVATE :: NMI_Horizontal_Modes
  PRIVATE :: NMI_Filter
  PRIVATE :: unitmx
  PRIVATE :: df, bf, cf, ef, hf, of, pf, qf 
  PRIVATE :: NMI_LoadData
  PRIVATE :: NMI_Diagnose

  ! development steps
  ! 1.1     revision
  ! 1.0     final version with explcit filtering only
  ! 0.10    last checks global average
  ! 0.9     new methode for separation, Lagrangian multiplier
  ! 0.8     check of separation of H in T and lnPS
  ! 0.4     filtering in MNI_Filter
  ! 0.3     calculation of modes checked
  !         no tides, no initialisation, no diabtic tendencies
  ! 3.0     revision march 2000

  ! define the filter characteristics
  ! normal filter
  REAL                 :: pcut
  INTEGER, ALLOCATABLE :: iplim(:,:,:) ! cut off limit for fast modes
  INTEGER              :: nvm          ! number of vertical modes included
  ! diabatic tendencies
  REAL                 :: pcutd
  INTEGER, ALLOCATABLE :: iplimd(:,:,:)! cut off limit for fast modes
  INTEGER              :: nvmd         ! number of vertical modes included

  ! accumulation of physical tendencies and iteration
  INTEGER  :: ntpre, ntdia, ntiter

  LOGICAL  :: lpsrs = .FALSE.  ! switch for pressure restore.

  ! internal used parameters and arrays

  ! vertical modes

  REAL, ALLOCATABLE :: zdb(:), &    ! delpr/b
&                      phib(:), &   ! phibar, equivalent depths
&                      vm_ev(:,:), &! m (vertical modes) eigenvectors
&                      vm_evi(:,:)  ! inverse(m)

  ! horizontal modes

  INTEGER  :: ing, &                     ! number of gravity modes
&             inr                        ! number of rossby modes
  REAL, ALLOCATABLE    :: zhmod(:)       ! eigenvectors and eigenvalues
  INTEGER, ALLOCATABLE :: ihmod(:,:,:,:) ! pointer array for horizontal modes
                                         !  (1,.,.,.) block starting pointer
                                         !  (2,.,.,.) number of Gravity modes (east+west)
                                         !  (3,.,.,.) number of Rossby modes
                                         !  (.,X,.,.) 1 .. symmetric, 2 .. antisymmetric part
                                         !  (.,.,X,.) zonal wave number
                                         !  (.,.,.,X) vertical mode

  ! work space for NMI calculation

  REAL, ALLOCATABLE, DIMENSION(:) :: &
&       svo_nmi, sd_nmi, stp_nmi, &  ! buffer A absolut
&       svo_dt,  sd_dt,  stp_dt      ! buffer B tendencies

  CHARACTER (256) :: nmi_mess

  ! options for NMI_make
  INTEGER, PARAMETER :: NMI_MAKE_TEST=0
  INTEGER, PARAMETER :: NMI_MAKE_FOD =1   ! Filter Observed Data
  INTEGER, PARAMETER :: NMI_MAKE_FMMO=3   ! Filter Model Minus Observed
  INTEGER, PARAMETER :: NMI_MAKE_NMI =5   ! Normal Mode Initialization

  ! options for data buffer handling
  INTEGER, PARAMETER :: NMI_LOAD_CLEAR   =0  ! clear buffer
  INTEGER, PARAMETER :: NMI_LOAD_ECH_BUFA=1  ! echam data into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_ECH_BUFB=2  ! echam data into buffer B
  INTEGER, PARAMETER :: NMI_LOAD_ORD_BUFA=3  ! observed read buffer into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_ORD_BUFB=4  ! observed read buffer into buffer B
  INTEGER, PARAMETER :: NMI_LOAD_OIN_BUFA=11 ! interpolated data into buffer A
  INTEGER, PARAMETER :: NMI_LOAD_OIN_BUFB=5  ! interpolated data into buffer B

  INTEGER, PARAMETER :: NMI_COPY_AB   =10  ! exchange buffer
  INTEGER, PARAMETER :: NMI_COPY_BA   =15
  INTEGER, PARAMETER :: NMI_CALC_TEND =20  ! calculate tendencies
  INTEGER, PARAMETER :: NMI_CALC_DIFF =25  ! difference A minus B
  INTEGER, PARAMETER :: NMI_CALC_ADD  =30  ! add both buffer
  INTEGER, PARAMETER :: NMI_CALC_DIAG =99  ! calculates diagnostics

  INTEGER, parameter :: NMI_STORE_FAST=50  ! store the fast mode part

  ! filter options, type of mode separation
  INTEGER, PARAMETER :: NMI_SEL_ALL     =0   ! use all modes without global average
  INTEGER, PARAMETER :: NMI_SEL_ALL_W0  =10  ! use all modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_GM      =1   ! use GW modes without global average
  INTEGER, PARAMETER :: NMI_SEL_GM_W0   =8   ! use GW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_RM      =2   ! use RW modes without global average
  INTEGER, PARAMETER :: NMI_SEL_RM_IMPL =5   ! subtract GW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_RM_W0   =9   ! use RW modes, result includes global average
  INTEGER, PARAMETER :: NMI_SEL_FILT    =6   ! GW modes without global average inside window
  INTEGER, PARAMETER :: NMI_SEL_FILT_W0 =3   ! GW modes partly, result includes global average, inside window
  INTEGER, PARAMETER :: NMI_SEL_FILTD   =7   ! GW modes without global average, inside second window
  INTEGER, PARAMETER :: NMI_SEL_FILTD_W0=4   ! GW modes partly, result includes global average
CONTAINS

!======================================================================
  SUBROUTINE NMI_Init
    ! initialization of NMI for ECHAM
    INCLUDE 'nmictl.inc'    ! define parameters for external manipulation of NMI

    INTEGER :: imem, ii, imo

    WRITE (nmi_mess,* ) 'NMI Version Rev. 3.0 17-March-2000 (kirchner@dkrz.de)'
    CALL message('NMI_Init',nmi_mess)

    ! iteration control
    ntpre  = 2
    ntdia  = 8
    ntiter = 2

    ! filter time limits [hours]
    pcut   = 12.0
    pcutd  =  5.0

    READ (nin,nmictl)       ! read namelist *NMICTL*

    ! general filter characteristic
    nvm     = nlev      ! maximal number of vertical modes normal filter
    nvmd    = nlev      ! maximal number of vertical modes diabatic filter

    ! check switches for initialization process
    ntpre = MIN(MAX(ntpre,2),1)
    ntdia = MIN(MAX(ntdia,8),1)
    IF (nk >= 106) THEN
       ! For resolutions greater or equal to t106 make the default
       ! number of timesteps used in the accumulation of physical
       ! tendencies for the initialisation equal to 2 hours as a
       ! function of the timestep
       ntdia = NINT(REAL(7200)/dtime+0.499)+2
    ENDIF
    ntiter = MIN(ntiter,2)

    pcut   = MAX(pcut, 0.0)
    pcutd  = MAX(pcutd,0.0)

    ! allocate space for vertical modes
    ALLOCATE (zdb (nlev));      zdb(:)   = 0.0
    ALLOCATE (phib(nlev));      phib(:)  = 0.0
    ALLOCATE (vm_ev (nlev,nlev)); vm_ev (:,:)  = 0.0
    ALLOCATE (vm_evi(nlev,nlev)); vm_evi(:,:) = 0.0

    ! allocate space for horizontal mode index table
    ALLOCATE (ihmod(3,2,nmp1,nlev)); ihmod(:,:,:,:) = 0

    ! estimate the memory for horizontal modes
    imem = 0
    DO ii = 1,nmp1
       imo  = 3*nnp(ii)
       imem = imem + (1+imo)*imo
    ENDDO
    ALLOCATE (zhmod(nlev*imem)); zhmod(:) = 0.0

    ! allocate memory for cut off limit of NMI filter
    ALLOCATE (iplim(2,nmp1,nlev)); iplim(:,:,:) = 0

    ! filter characteristic of diabatic tendencies
    ALLOCATE (iplimd(2,nmp1,nlev)); iplimd(:,:,:) = 0

    ! allocate spectral buffer A
    ALLOCATE (svo_nmi(nlev  *n2sp)); svo_nmi(:) = 0.0
    ALLOCATE (sd_nmi (nlev  *n2sp)); sd_nmi(:) = 0.0
    ALLOCATE (stp_nmi(nlevp1*n2sp)); stp_nmi(:) = 0.0

    ! allocate spectral buffer B
    ALLOCATE (svo_dt(nlev  *n2sp)); svo_dt(:) = 0.0
    ALLOCATE (sd_dt (nlev  *n2sp)); sd_dt(:) = 0.0
    ALLOCATE (stp_dt(nlevp1*n2sp)); stp_dt(:) = 0.0

    CALL message('','Optional parameters:')
    WRITE(nmi_mess,*) '   cut off period (hours)           PCUT  = ',pcut
    CALL message('',nmi_mess)
    IF (.NOT. lnudge) THEN
       WRITE(nmi_mess,*) '   cut off period (hours) diabatic  PCUTD = ',pcutd
       CALL message('',nmi_mess)
       WRITE(nmi_mess,*) '   start up period (steps)         NTPRE  = ',ntpre
       CALL message('',nmi_mess)
       WRITE(nmi_mess,*) '   accumulation period (steps)     NTDIA  = ',ntdia
       CALL message('',nmi_mess)
       WRITE(nmi_mess,*) '   iteration period (steps)        NTITER = ',ntiter
       CALL message('',nmi_mess)
       CALL message('','NMI in INITIALIZATION mode')
    ELSE
       CALL message('','NMI in FILTER mode')
    ENDIF

    CALL NMI_Normal_Modes    ! Initialisation is required - establish modes

    CALL message('','NMI initialized successful')

  END SUBROUTINE NMI_Init

!======================================================================

  SUBROUTINE NMI_Make(ityp)
    ! Control the NMI procedure

    INTEGER   :: ityp

    SELECT CASE (ityp)
    CASE(NMI_MAKE_TEST) ! test loop
       CALL finish('TEST BENCH ENDE','Run aborted.')

    CASE (NMI_MAKE_FOD)            ! *** filter ERA data
       CALL NMI_LoadData(NMI_LOAD_ORD_BUFA)    ! observed --> buf A
       CALL NMI_Filter(.TRUE.,NMI_SEL_FILT_W0) ! filter ERA
       CALL NMI_LoadData(-NMI_LOAD_ORD_BUFA)   ! buf A --> observed

!ik    CASE (2)                       ! *** compose nudging data
!ik       CALL NMI_LoadData(NMI_LOAD_ECH_BUFA)    ! model --> buf A
!ik       CALL NMI_Filter(.TRUE.,NMI_SEL_FILT)    ! calculate fast GW mode part
!ik       CALL NMI_LoadData(NMI_LOAD_OIN_BUFB)    ! interpolated --> buf B
!ik       CALL NMI_LoadData(NMI_CALC_ADD)         ! buf A + buf B --> buf B
!ik       CALL NMI_LoadData(-NMI_LOAD_OIN_BUFB)   ! buf B --> interpolated

    CASE (NMI_MAKE_FMMO)           ! *** filter anomaly (model-obs)
       CALL NMI_LoadData(NMI_LOAD_ECH_BUFA)    ! model --> buf A
       CALL NMI_LoadData(NMI_LOAD_OIN_BUFB)    ! interpolated --> buf B
       CALL NMI_LoadData(NMI_CALC_DIFF)        ! buf A - buf B --> buf A
       CALL NMI_Filter(.TRUE.,NMI_SEL_FILT)    ! calculate fast part of model-era
!ik       if (lptime) call NMI_LoadData(NMI_STORE_FAST)  ! store fast part
       call NMI_LoadData(NMI_STORE_FAST)       ! store fast part
       CALL NMI_LoadData(NMI_CALC_ADD)         ! buf A + buf B --> buf B
       CALL NMI_LoadData(-NMI_LOAD_OIN_BUFB)   ! buf B --> interpolated

    CASE (NMI_MAKE_NMI)            ! *** initialsation loop

       WRITE(nmi_mess,*) 'not implemented now'
       CALL finish('',nmi_mess)

       IF (nstep == nstart + ntpre) THEN     ! Start diabatic accumulation
          WRITE(nmi_mess,*) 'start diabatic accumulation at step ',nstep
          ! perform anomalies
          ! (n.b. the first two timesteps are not used in the accumulation
          ! of the physical tendencies).
          CALL NMI_LoadData(NMI_LOAD_ECH_BUFB)     ! load data in second NMI buffer

       ELSE IF (nstep < nstart + ntpre + ntdia) THEN     ! Initialize iteration
          WRITE(nmi_mess,*) 'accumulate forcing at step ',nstep
          ! interation mit tagesgang aktiv, auswirkungen auf initialisierung des modells??
          ! ** accumulate filtered forcing over NDIFFST timesteps
          ! store diabatic forcing as spectral arrays
          ! accumulate terms

       ELSE IF (nstep >= nstart + ntpre + ntdia) THEN    ! Perform iteration steps
          ldiur = .FALSE.                                ! no diurnal cycle
          ! rueckspeichern der startwerte welche notwendig ???
          IF (nstep == nstart + ntpre + ntdia) THEN          ! Initialize iteration
             WRITE(nmi_mess,*) 'filter accumulated forcing at step ',nstep
             ! accumulate diabatic forcing in grid point space
             ! problem with filtering, before filter spectral transform neccessary
             ! ** at the end of NDIFFST timesteps filter the diabatic forcing
             !<<< load DIABATIC forcing in NMI buffer
             CALL NMI_Filter(.TRUE.,NMI_SEL_FILTD_W0) ! diabatic filter
             !<<< ** average diabatic part
             !<<< reload diabatic terms for iteration
             CALL NMI_LoadData(-NMI_LOAD_ECH_BUFB)    ! restore data for iteration start
             ! start iteration with unfiltered data
             ! now only adiabatic part calculated, add diabatic part in gridpoint space
          ELSE
             WRITE(nmi_mess,*) 'interation at step ',nstep
             ! ** add adiabatic and averaged diabatic forcing
             CALL NMI_LoadData(NMI_LOAD_ECH_BUFA)     ! load new spectral arrays
             CALL NMI_LoadData(NMI_CALC_TEND)         ! calculate tendencies
             CALL NMI_Filter(.FALSE.,NMI_SEL_FILT_W0) ! filter tendencies and perform iteration
             CALL NMI_LoadData(-NMI_LOAD_ECH_BUFA)    ! reset spectral arrays
             CALL NMI_LoadData(NMI_COPY_BA)           ! load filtered data in second array for next iteration
          ENDIF
          CALL message('NMI_Make',nmi_mess)
          IF (nstep == nstart + ntpre + ntdia + ntiter ) THEN
             WRITE(nmi_mess,*) 'stop iteration at step ',nstep
             CALL message('NMI_Make',nmi_mess)
             ! at the end of the iteration initialized arrays OK
             ldiur = .TRUE.              ! switch diurnal cycle on
             lnmi  = .FALSE.             ! reset the NMI switch
             nstep = nstart + ntpre      ! restart the forecast
             WRITE(nmi_mess,'(a,i10)') 'restart forecast at step ',nstep
             CALL message('NMI_Make',nmi_mess)
             ! clear NMI memory
             DEALLOCATE (zdb)
             DEALLOCATE (phib)
             DEALLOCATE (vm_ev)
             DEALLOCATE (vm_evi)
             DEALLOCATE (ihmod)
             DEALLOCATE (zhmod)
             DEALLOCATE (iplim)
             DEALLOCATE (iplimd)
             DEALLOCATE (svo_nmi)
             DEALLOCATE (sd_nmi)
             DEALLOCATE (stp_nmi)
             DEALLOCATE (svo_dt)
             DEALLOCATE (sd_dt)
             DEALLOCATE (stp_dt)
             WRITE(nmi_mess,*) 'end of INITIAL NMI'
         ENDIF
      ENDIF
      CALL message('NMI_Make',nmi_mess)

    END SELECT

  END SUBROUTINE NMI_Make

!======================================================================
  SUBROUTINE NMI_Normal_Modes
    ! calculate normal modes

    REAL, ALLOCATABLE :: &
&           zwrk1(:), zwrk2(:,:), zwrk3(:,:), zwrk4(:,:) &   ! workspace
&         , zvr(:,:), zvi(:,:)                   ! real/imaginary part of eigenvectors
    REAL    :: zevima, zevimi, zsum, zscale, zdtim2 &
&         , renorm1                    ! time unit factor

    INTEGER,ALLOCATABLE :: ipivot(:)            ! permutation index
    INTEGER :: &
&           ifail &                  ! failure code for external procedures
&         , imar(1), imai(1) &
&         , ipos &                   ! pointer in horizontal mode array
&         , ings &       ! no. of modes symmetric part
&         , jk, js, jm, jj, jk1, jk2 &
&         , iposr, iposge, iposgw    ! pointer for different wave types

    CALL message('NMI_Normal_Modes','Vertical mode setup:')

    ! allocate memory for vertical mode calculation
    ALLOCATE (ipivot(nlev))
    ALLOCATE (zvr(nlev,nlev))
    ALLOCATE (zvi(nlev,nlev))
    ALLOCATE (zwrk1(4*nlev))
    ALLOCATE (zwrk2(nlev,nlev)) 
    ALLOCATE (zwrk3(nlev,nlev))
    ALLOCATE (zwrk4(nlev,nlev))

    ! Find, sort and print eigenvalues and vectors of bb
    IF (betadt == 0.0) THEN
       zdtim2 = (dtime*0.5)**2         ! neccessary for full explicit scheme
    ELSE
       zdtim2 = (dtime*0.5*betadt)**2  ! use this factor to get meaningfull data range
    ENDIF

    ! the vertival structure matrix is a function of the Basic State
    ! calculate eigenvalues and eigenvectors
    ! rescale structur matrix (gravitation-matrix)
    zwrk2(1:nlev,1:nlev) = TRANSPOSE(bb(1:nlev,1:nlev)) / zdtim2

    ! all other possibilities are not correct
    ! bb^T or bb --> eigenvalues larger then expected
    ! bb or bb/zdtim2 --> eigenvektor others then given in paper

    ifail = 0
    CALL dgeev (&    ! LAPACK eigenvalue solver A*X = lamdba*X
&         'n',&          ! calculate no left eigenvectors
&         'v',&          ! calculate right eigenvectors
&         nlev,&         ! order of matrix A
&         zwrk2,&        ! inp: structure matrix A, out: overwritten
&         nlev,&         ! leading dimension of A
&         zdb,&          ! real part of eigenvalues
&         phib,&         ! imaginary part of eigenvalues
&         zvr,&          ! contains nothing
&         nlev,&         ! leading dimension of left eigenvectors
&         zvi,&          ! right eigenvectors in columns
&         nlev,&         ! leading dimension of right eigenvectors
&         zwrk1,&        ! work space
&         4*nlev,&       ! dimension of workspace
&         ifail)         ! error return code
    IF (ifail /= 0) THEN
       CALL message('','Calculation of eigenvectors/values failed.')
       WRITE (nmi_mess,'(a,i4)') ' Failure in s/dgeev, ifail = ', ifail
       CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    imai   = MAXLOC(ABS(phib(1:nlev)))    ! find eigenvalue with largest imaginary part
    zevima = phib(imai(1))
    imai   = MINLOC(ABS(phib(1:nlev)))    ! find eigenvalue with smallest imaginary part
    zevimi = phib(imai(1))
    IF ((ABS(zevima) > 0.0).OR.(ABS(zevimi) > 0.0)) &
&         CALL finish ('NMI_Normal_Modes',&
&         'Found complex eigenvalues. Failure in vertical structure matrix')

    ! sort eigenvalues and eigenvectors
    jk = 1
    DO
       IF (ABS(phib(jk)) > 0.0) THEN     ! found complex and complex conjugate eigenvalue
          vm_ev(:,jk)    =  zvi(:,jk)
          vm_evi(:,jk)   =  zvi(:,jk+1)
          vm_ev(:,jk+1)  =  zvi(:,jk)
          vm_evi(:,jk+1) = -zvi(:,jk+1)
          jk = jk+2
       ELSE                          ! real eigenvalue
          vm_ev(:,jk)  = zvi(:,jk)
          vm_evi(:,jk) = 0.0
          jk = jk+1
       ENDIF
       IF (jk > nlev) EXIT
    ENDDO

    DEALLOCATE (zvr)
    DEALLOCATE (zvi)
    DEALLOCATE (zwrk1)

    DO jk=1,nlev
       imar         = MAXLOC(zdb(:))      ! find next largest real part of all eigenvalues
       phib(jk)     = zdb(imar(1))
       zdb(imar(1)) = 0.0
       vm_evi(:,jk) = vm_ev(:,imar(1))       ! store real part of eigenvector
    ENDDO

    ! equivalent depths in PHIB
    ! corresponding eigenvectors in VM_EVI

    WRITE (nmi_mess,'(a,f6.1,a,f8.0,a)') &
&         '  Vertical modes for tr = ', tr, ' k;   pr = ', apr,' pa'
    CALL message('',nmi_mess)
    WRITE (nmi_mess,'(a)') &
&         '   level       a(k)        b(k)     eigenvalues    phase velocity   EV/g[gpm]'
    CALL message('',nmi_mess)
    DO jk = 1,nlev
       WRITE (nmi_mess,'(3x,i4,2x,f13.6,f12.8,f13.2,f14.2,f13.2)') &
&         jk,vct(jk),vct(nvclev+jk), phib(jk), SQRT(phib(jk)), phib(jk)/g
       CALL message('',nmi_mess)
    END DO

    ! eigenvectors normalized
    ! convert vertical modes to "universal" scaling

    zevima = SUM(delpr(:))                ! check the normalization
    IF (ABS(zevima-apr) > 1.e-5) THEN
       CALL message('','Universal scaling failure.')
       WRITE (nmi_mess,'(2(a,e20.10))') &
&            '  surface pressure ... expected = ',apr,' calculated = ',zevima
       CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    ! for rescaling and orthogonalization of vertical modes see paper (198?) by C.Temperton
    !
    !    WRITE (nout,'(a)') ' **NMI_Normal_Modes** Rescaling factor for: '
    !    DO  jk=1,nlev
    !       zsum      = SUM(vm_evi(:,jk)*vm_evi(:,jk)*delpr(:))
    !       zscale    = SQRT(apr/zsum)
    !       !       vm_ev(:,jk)  = zscale*vm_evi(:,jk)
    !       vm_ev(:,jk)  = vm_evi(:,jk)
    !       zrscale   = 1.0/zscale
    !       WRITE (nout,'(a,i3,a,f14.10,a,f14.10,a)') &
    !&            '   vertical mode', jk, ' = ', zscale,&
    !&            '   ( ~ for coefficients = ', zrscale, ')'
    !    ENDDO
    !    WRITE (nout,'(a)') '  do not use the universal rescaling now !!!'
    !
    ! scaled eigenvectors in ZM

    ! Compute h**(-1) and delpr*b**(-1)
    ! calculate inverse of ZM

    zwrk2(:,:) = vm_ev(:,:)
    CALL unitmx(zwrk3,nlev)
    ifail = 0
    CALL dgesv (&       ! BLAS linear equation solver
&         nlev,&            ! number of linear equations NEQ [int]
&         nlev,&            ! number of right hand sides NRHS [int]
&         zwrk2,&           ! inp: coefficients matrix A(LDA,NEQ), out: l-u-factorization [real]
&         nlev,&            ! leading dimension LDA [int]
&         ipivot,&          ! out: permutation matrix IPIV(NEQ) [int]
&         zwrk3,&           ! inp: right hand side X(LDX,RHS), out: solution matrix [real]
&         nlev,&            ! leading dimension LDX [int]
&         ifail)            ! error code [int]
    IF (ifail /= 0) THEN
       WRITE (nmi_mess,'(a,i4)') 'Inverting VM_EV, failure in s/dgesv, ifail = ', ifail
       CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF
    vm_evi(:,:) = zwrk3(:,:)    ! inverse eigenvector matrix

    ! calculate inverse of BB
    zwrk2(1:nlev,1:nlev) = TRANSPOSE(bb(1:nlev,1:nlev)) / zdtim2
    CALL unitmx(zwrk3,nlev)
    ifail = 0
    CALL dgesv (nlev,nlev,zwrk2,nlev,ipivot,zwrk3,nlev,ifail)
    IF (ifail /= 0) THEN
       WRITE (nmi_mess,'(a,i4)') 'Inverting BB, failure in s/dgesv, ifail = ', ifail
       CALL finish ('NMI_Normal_Modes',nmi_mess)
    ENDIF

    DO  jk=1,nlev                ! final calculation of ZDB (S^T*BB^-1)
       zdb(jk) = DOT_PRODUCT(delpr(:),zwrk3(:,jk))/apr
    ENDDO

    ! delpr/bb/apr stored in ZDB
    DEALLOCATE (ipivot)
    DEALLOCATE (zwrk2)
    DEALLOCATE (zwrk3)
    DEALLOCATE (zwrk4)

    CALL message('NMI_Normal_Modes','Calculate horizontal modes:')
    ipos = 0
    DO  jm = 0, nm                            ! calculation for all zonal wave numbers
       CALL NMI_Horizontal_Modes (jm,ipos)
                                              ! Print frequencies of horizontal modes
       renorm1 = api/omega/3600.0             ! get period in hours
       ! Set the filter limits for the NMI
       DO jk = 1,nlev
          DO js = 1,2                      ! loop over symmetric/antisymmetric part
             ings = ihmod(2,js,jm+1,jk)/2  ! get number of one set of gravity modes
             DO jj = 1,ings                ! loop over one set of gravity modes (westerly part)
                iposgw = ihmod(1,js,jm+1,jk) -1 +ings +jj  ! position of eigenvalue
                ! normal filter ... NVM and IPLIM
                IF (renorm1/zhmod(iposgw) < pcut) THEN
                   ! found period smaller then the cut-off period
                   nvm = jk
                   iplim(js,jm+1,jk) = iplim(js,jm+1,jk) + 1
                END IF
                !  diabatic tendency filter ... NVMD and IPLIMD
                IF (renorm1/zhmod(iposgw) < pcutd) THEN
                   ! found period smaller then the cut-off period
                   nvmd = jk
                   iplimd(js,jm+1,jk) = iplimd(js,jm+1,jk) + 1
                END IF
             END DO
          END DO
       END DO

       WRITE (nmi_mess,'(a,i5,a)') &
&            'No of horizontal modes for wavenumber ', jm, ' in cut-off window:'
       CALL message('',nmi_mess)
       jk1 = 1
       DO
          jk2 = jk1 + 8
          jk2 = MIN(jk2,nvm)
          WRITE(nmi_mess,'(20(a,i3.3,a,i5,1x))') &
               ('VM_(',js,'): ',iplim(1,jm+1,js)+iplim(2,jm+1,js),js=jk1,jk2)
          CALL message('',nmi_mess)
          IF (jk2 == nvm) EXIT
          jk1 = jk2+1
       END DO
       IF (.NOT.lnudge) THEN
          jk1 = 1
          DO
             jk2 = jk1 + 8
             jk2 = MIN(jk2,nvm)
             WRITE(nmi_mess,'(20(a,i3.3,a,i5,1x))') &
                  ('VMD(',js,'): ',iplim(1,jm+1,js)+iplim(2,jm+1,js),js=jk1,jk2)
             CALL message('',nmi_mess)
             IF (jk2 == nvm) EXIT
             jk1 = jk2+1
          END DO
       END IF

    ENDDO

  END SUBROUTINE NMI_Normal_Modes
!======================================================================
  SUBROUTINE NMI_Horizontal_Modes(km, ipos)
    ! calculate horizontal modes for one zonal wave number

    INTEGER, INTENT(IN) :: km        ! actual zonal wave numbers for which
                                     ! modes are to be computed
    INTEGER, INTENT(INOUT) :: ipos   ! record position in ihmod

    REAL, ALLOCATABLE :: &
&         za(:,:), &         ! horizontal struktur matrix
&         zr(:), &           ! eigenvalues of horizontal structur matrix
&         work(:)

    INTEGER :: &
&         ipk, &             ! maximal number of meridional index
&         kmo, &             ! maximal field dimension
&         jk, &              ! vertical loop index
&         jsym, &            ! symmetry loop
&         jr, il, &
&         imo, &             ! number of modes
&         jms, &
&         jn, &
&         ifail, inge, j, jj, jl, ingt

    REAL    :: zed        ! local equivalent depths
    INTEGER :: imo2

    ! allocate local memory
    ipk = MIN(nn+km,nk)
    kmo = 3*nnp1 + 6
    ALLOCATE (za  (kmo,kmo))
    ALLOCATE (zr  (kmo))
    ALLOCATE (work(3*kmo-1))

    ! start loop on equivalent depths

    DO jk = 1, nlev       ! all vertical modes are included
       zed = phib(jk)     ! get the equivalent depths

       DO jsym = 0, 1     ! separate loop for symmetric and antisymmetric case
          za(:,:) = 0.0
          zr(:)   = 0.0
          work(:) = 0.0

          ! Define matrix A (za)
          ! fill in upper triangel + diagonal of A
          ! IMO ... number of modes can be solved

          imo = 2*INT((ipk-km+jsym)/2)+INT((ipk-km+1-jsym)/2)
          IF (km == 0) THEN
             jms = 2
             il  = 2 + jsym
          ELSE
             imo = imo + 2 - jsym
             jms = km
             il  = 1
          ENDIF

          imo2 = 0
          IF (jsym == 0) THEN                ! symmetric flow
             DO jn = jms, ipk, 2             ! loop for velocity potential and mass
                za(il,il)   =  cf(km,jn)
                za(il,il+1) =  ef(jn,zed)
                za(il,il+2) = -bf(km,jn+1)
                il   = il   + 3
                imo2 = imo2 + 2
             ENDDO
             il = 3                          ! loop for streamfct.
             IF (km == 0) il = 1
             DO jn = km+1, ipk, 2
                za(il,il)   =  cf(km,jn)
                za(il,il+1) = -bf(km,jn+1)
                il   = il   + 3
                imo2 = imo2 + 1
             ENDDO
          ELSE                               ! antisymmetric flow
             DO jn = jms, ipk, 2             ! loop for streamfct.
                za(il,il)   =  cf(km,jn)
                za(il,il+1) = -bf(km,jn+1)
                il   = il   + 3
                imo2 = imo2 + 1
             ENDDO
             il = 2                          ! loop for velocity pot. and mass
             IF (km == 0) il = 1
             DO jn = km+1, ipk, 2
                za(il,il)   =  cf(km,jn)
                za(il,il+1) =  ef(jn,zed)
                za(il,il+2) = -bf(km,jn+1)
                il   = il   + 3
                imo2 = imo2 + 2
             ENDDO
          ENDIF
          ! setup of matrix A following W.Wergen (OK 02-feb-99)
          IF (imo/=imo2) THEN
             WRITE(nout,*) ' ERROR: count modes ',imo2,' expected ',imo
             CALL finish ('NMI_Horizontal_Modes', 'Run aborted.')
          ENDIF

          DO jr = 1, imo                   ! fill in lower triangal
             DO jn = jr+1, imo
                za(jn,jr) = za(jr,jn)
             ENDDO
          ENDDO

          ! Solve eigenproblem
          ifail = 0
          CALL dsyev (&  ! LAPACK get eigenvalues and eigenvectors of real symmetric matrix
&               'V',&        ! get eigenvalues and eigenvectors
&               'L',&        ! lower triangle is in  A stored
&               imo,&        ! order of matrix, N
&               za,&         ! inp: symmetric matrix A(LDA,N), out: orthonormal eigenvectors
&               kmo,&        ! leading dimension of matrix A LDA
&               zr,&         ! out: eigenvalues in ascending order
&               work,&       ! workspace
&               3*kmo-1,&    ! dimension of workspace
&               ifail)       ! return code
          IF (ifail /= 0) THEN
             WRITE (nmi_mess,'(a,i2,a,g13.5,a,i3,a,i2)') &
&                  ' Error in s/dsyev - ifail = ', ifail, &
&                  ' for phi = ', zed, ' m = ', km, ' and isym = ', jsym
             CALL finish ('NMI_Horizontal_Modes',nmi_mess)
          ENDIF

          ! Store sorted frequencies and modes in buffer ZHMOD
          !
          ! zr == free frequencies of normal modes
          ! order is (gravity east, rossby, gravity west)
          ! in ZR (GE[max..0] R[max..0] GW[0..max])
          ! following the frequency [meridionalindex]
          !
          ! reorder in ZHMOD
          ! (GE[0..max] GW[max..0] R[0..max]
          ! slow......fast ........ slow

          inge = 0
          DO  j = 1,imo      ! count eastward moving gravity waves
             IF (zr(j) > -1e-10) EXIT
             inge = inge+1
          ENDDO
          ing = 2*inge       ! all gravity modes
          inr = imo-ing      ! residue is number of rossby modes
          ! check the number of modes
          IF (jsym == 0) THEN
             ingt = (nk-km+2)/2
             IF (km == 0) ingt = ingt-1
          ELSE
             ingt = (nk-km+1)/2
          ENDIF
          IF (inge /= ingt) THEN
             WRITE (nmi_mess,'(a,/,a,3i3,a,i4,a,i4)') &
&                  '   Possible problem in counting eastward gravity waves', &
&                  '   - for jm, jk, jsym (', km, jk, jsym, ') theory: ', &
&                  inge, ' reality: ', ingt
             CALL finish ('NMI_Horizontal_Modes',nmi_mess)
          ENDIF

          ! set pointer array for horizontal modes
          ! ipos is the last index of the previous array
          ! ihmod(1,...) is the position of the first place in the actual array

          ihmod(1,jsym+1,km+1,jk) = ipos + 1    ! start position
          ihmod(2,jsym+1,km+1,jk) = ing         ! no. of gravity modes
          ihmod(3,jsym+1,km+1,jk) = inr         ! no. of rossby modes
          DO jj=1,inge                          ! store eigenvalues ZR
             zhmod(ipos     +jj) = zr(inge+1-jj)  ! GW east
             zhmod(ipos+inge+jj) = zr(imo +1-jj)  ! GW west
          END DO
          ipos = ipos + ing

          IF (km==0) THEN
             ! reset frequencies of rossby modes for wave number zero
             ! these Rossby modes stationary

             zr(inge+1:inge+inr) = 0.0

             ! the eigenvectors are non zero 
             ! the mixed rossby-gravity mode and the kelvin mode are
             ! gravity modes here
          END IF

          DO jj=1,inr
             zhmod(ipos+jj) = zr(inge+inr+1-jj)   ! Rossby modes
          END DO
          ipos = ipos + inr
          DO jj=1,inge                            ! store eigenvectors ZA
             DO jl = 1,imo
                zhmod(ipos+jl) = za(jl,inge+1-jj)
             END DO
             ipos = ipos + imo
          END DO
          DO jj=1,inge
             DO jl = 1,imo
                zhmod(ipos+jl) = za(jl,imo+1 -jj)
             END DO
             ipos = ipos + imo
          END DO
          DO jj=1,inr
             DO jl=1,imo
                zhmod(ipos+jl) = za(jl,inge+inr+1-jj)
             END DO
             ipos = ipos + imo
          END DO

          ! IPOS is now the last memory position in ZHMOD used
       ENDDO
    ENDDO

    DEALLOCATE (work)
    DEALLOCATE (zr)
    DEALLOCATE (za)

    RETURN
  END SUBROUTINE NMI_Horizontal_Modes

!======================================================================
  SUBROUTINE NMI_Filter(lproj,iiproj)
    ! perform the NMI filter or the non-linear iteration

    LOGICAL, INTENT(in) :: lproj
                            ! switch for projection or initialization
                            ! .true.    only projection of absolut values and filtering
                            ! .false.   solve the tendecy equation for next iteration
    INTEGER, INTENT(in) :: iiproj              ! type of projection

    INTEGER :: iproj
                             !   0= all modes
                             !   1= gravity modes
                             !   2= rossby modes
                             !   3= normal filter window (see pcut)
                             !   4= diabatic filter window (see pcutd)
                             !   5= all minus gravity part

    INTEGER :: invm          ! number of modes actual used

    REAL, ALLOCATABLE :: &
&         zzdb(:), &
&         zph(:), &                  ! half level pressure
&         zrlnpr(:), &               ! log (delta p)
&         zralphr(:), &              ! alpha
&         zralprr(:), &              ! 1/alpha
&         zbs(:), &    ! spectral coefficients T, lnPS
&         zworks(:), & ! spectral coefficients VO, D, H
&         zbv(:), &    ! parts of different waves in vertical space
&         zworkv(:), & ! subset of coefficients in vertical space
&         zworkm(:), & ! projection in horizontal mode space
&         zbb(:), &    ! separation matrix
&         zlmm(:,:), zlmy(:) ! workspace for Lagragian multiplier method

    INTEGER,ALLOCATABLE :: ipivot(:)            ! permutation index

    REAL    :: &
&         zedra, &   ! sqrt(equivalent depth) / (radius of earth)
&         zs, zo     ! helpvariables explicit solution

    INTEGER :: &
&         jpfs3d, jpfsm, jpbnm, &   ! array dimensions
&         isc1, isc2, isc3, &       ! index for scaling, rescaling
&         zdtim, &         ! scaling of equivalent depths
&         ijump, &       ! record length of actual vertical mode data set
&         isp, &         ! local pointer in spectral array (nlev)
&         ismp, &        ! local pointer in spectral array (nlevp1)
&         jm, &          ! index of wavenumber
&         jn, &          ! zonal index number
&         ipk, &         ! maximal zonal index for local wave number
&         innp2, &       ! actual number of one set of spectral coefficients real+imaginary
&         ivmo, &        ! one set of spectral coefficients (VO+D+H)
&         innp2nm, &     ! one set of components for all zonal index
&         ilt, &         ! actual record length in VO and D
&         iltp, &        ! actual record length in TP
&         ist, &         ! pointer in T and lnPS
&         ish, &         ! pointer in H
&         jk, &          ! loop index of level and vertical modes
&         izr, izi, idr, idi, ihr, ihi, &     ! pointer in vertical space arrays
&         ipos, &        ! pointer in horizontal mode array
&         imo, &         ! total number of horizontal modes
&         ingrx, &       ! number of modes used
&         iposx, &       ! local position in mode array
&         iposv, &       ! local position in eigenvector array
&         jsym, &        ! index for symmetry loop
&         idx1, idx2, idx3, idx4, idx5, idx6 ! symmetry loop steering 

    REAL    :: oro

    ! initialization and checks

    IF (lproj) THEN
       iproj = iiproj
       lpsrs = .FALSE.
    ELSE
       iproj = NMI_SEL_FILT_W0     ! only one filter for initialization allowed
    ENDIF

    SELECT CASE (iproj)
    CASE (NMI_SEL_FILT_W0,NMI_SEL_FILT)  ; invm = nvm  ! normal filter
    CASE (NMI_SEL_FILTD_W0,NMI_SEL_FILTD); invm = nvmd ! diabatic filter
    CASE default                         ; invm = nlev ! all vertical modes included
    END SELECT

    IF (betadt == 0.0) THEN
       zdtim = dtime*0.5         ! neccessary for full explicit scheme
    ELSE
       zdtim = dtime*0.5*betadt  ! use this factor to get meaningfull data range
    ENDIF

    ! vertikal matrix components for basic state

    ALLOCATE (zph     (nlevp1))
    ALLOCATE (zrlnpr  (nlev  ))
    ALLOCATE (zralphr (nlev  ))
    ALLOCATE (zralprr (nlev  ))

    ALLOCATE (zlmm(2*nlev+1,2*nlev+1))
    ALLOCATE (zlmy(2*nlev+1))
    ALLOCATE (ipivot(2*nlev+1))

    ! rdtr is modified in model setup  (rdtr .ne. rd*tr) here

    ! ===========================================================================

    jpfs3d = 2 * nlev   * nnp1  ! max length of spectral array for D and VO
    jpfsm  = 2 * nlevp1 * nnp1  ! max length of spectral array for T and lnPS
    jpbnm  = 3 * nnp1           ! max number of horizontal modes

    ! workspace in spectral space
    ALLOCATE ( zbs    (2 * jpfsm) )        ! ZBS(TP)
    ALLOCATE ( zworks (3 * jpfs3d) )       ! ZWORKS(V+D+H)
    ALLOCATE ( zbb    (2* nnp1) )          ! ZBB(SPECTCOEF)

    ! workspace in vertical mode space
    ALLOCATE ( zbv (invm * 2 * jpbnm) )    ! ZBV(invm,3*(REAL+IMAG)
    ALLOCATE ( zworkv     (2 * jpbnm) )    ! ZWORKV(3*(REAL+IMAG))

    ! workspace in horizontal mode space
    ALLOCATE ( zworkm (2*3*nnp1) )         ! ZWORKM(imo,2)

    ALLOCATE (zzdb (nlev))
    ! ===========================================================================

    ! setup basic state, only one column
    CALL pres(zph,1,apr,1) ! pressure levels for basic state
    zrlnpr(:) = LOG(apr)
    CALL auxhyb(&
&         delpr,&      ! OUT: pressure difference
&         rdelpr,&     ! OUT: reziproce of DELPR
&         zrlnpr,&     ! OUT: logarithm of pressure differences
&         zralphr,&    ! OUT: ALPHAs for integration
&         zph,&        ! INP: half-level pressure
&         1,&          ! first dimension of arrays
&         1)           ! number of points
    zralprr(:) = 1./zralphr(:)

    ! initialize pointer in spectral arrays

    ijump = 2*invm

    nmi_waveno_loop : DO jm = 0, nm    !>>>>>>>>>>>>>>>>>>> wavenumber loop START

       ! define local dimensions and pointer in work space
       isp     = nmp(jm+1)*2*nlev   + 1
       ismp    = nmp(jm+1)*2*nlevp1 + 1

       ipk     = MIN(nn+jm,nk)
       innp2   = nnp(jm+1) * 2   ! spectral components of given wavenumber (R+I)
       ilt     = innp2 * nlev    ! offset in VO and D spectral field
       iltp    = innp2 * nlevp1  ! offset in T/lnPS spectral field
       ivmo    = innp2 * 3       ! normal mode components (VO+D+H)
       innp2nm = innp2 * invm    ! vertical space components

       zworks(    1:  ilt) = svo_nmi(isp :isp +ilt -1)   ! Fill in VO and D
       zworks(ilt+1:2*ilt) = sd_nmi (isp :isp +ilt -1)
       zbs(       1:iltp)  = stp_nmi(ismp:ismp+iltp-1)   ! Fill in T and lnPS

       ! Calculate H with T' and lnPS'
       !
       ! compute h for lowest level
       DO jn = 1, innp2
          ist =         jn*nlevp1                ! pointer  lnPS
          ish = 2*ilt + jn*nlev                  ! pointer  H

          ! orography correction not necessary ... test only
          oro = 0.0
          zworks(ish) = &                                    ! ==> H(nlev)
&               zralphr(nlev) * zbs(ist-1) &                 ! g * T(nlev)
&               + rd*tr       * zbs(ist)   &                 ! R*Tr*lnPS
&               + oro    ! orography correction
       END DO

       ! compute H for remaining levels
       DO jn = 1,innp2
          ist =         (jn-1)*nlevp1            ! pointer  T
          ish = 2*ilt + (jn-1)*nlev              ! pointer  H
          DO jk = nlevm1,1,-1
             zworks(ish+jk) = &                                  ! ==> H(level)
&                  zworks(ish+jk+1) &                            ! H(level+1)
&                  + zbs(ist+jk)  *              zralphr(jk) &   ! T(level)
&                  + zbs(ist+jk+1)*(zrlnpr(jk+1)-zralphr(jk+1))  ! T(level+1)
          END DO
       END DO

       IF (jm == 0) THEN

          ! conserve wave 0 first meridional index (= global average unchanged)
          ! this should be the basic state
          zworks(      1:      2*nlev)   = 0.0
          zworks(  ilt+1:  ilt+2*nlev)   = 0.0
          zworks(2*ilt+1:2*ilt+2*nlev)   = 0.0
       ENDIF

       ! ZWORKS complete defined (VO+D+H) structure is (NLEV,2,3*NNP)
       !
       ! ========================================================================
       ! Expansion into vertical mode space (page 12 W.Wergen)
       !
       ! M == VM_EV and invers(M) == vm_evi
       ! use vm_evi for projection, organizes as VM_EVI(nlev,invm+)
       !
       ! transpose[VM_EVI(nlev,invm)]*ZWORKS(nlev,ivmo) = ZBV(invm,ivmo)

       CALL dgemm &
&            ('N', 'N', invm, ivmo, nlev, 1.0, vm_evi, nlev, zworks, nlev, 0.0, zbv,invm)

       ! ZBV(invm,ivmo) now defined
       !>>>>>>>>>>>>>>>>> vertical mode loop START
       ! Treat horizontal dependence for first invm vertical modes

       vmode_loop : DO jk = 1, invm

          zedra  = a/SQRT(phib(jk))

          ! phib is calculated using rescaled BB
          !zedra  = a/(SQRT(phib(jk))*zdtim)
          ! rescaling not necessary
          ! set pointer for accumulation of wavenumber fraction

          izr    = jk                      ! ZETA real part
          izi    = jk             + invm   !      imaginary part 
          idr    = jk +   innp2nm          ! DIV  real part 
          idi    = jk +   innp2nm + invm   !      imaginary part
          ihr    = jk + 2*innp2nm          ! H    real part
          ihi    = jk + 2*innp2nm + invm   !      imaginary part

          ! Scale fields in ZBV (will be multiplied by i later)
          ! the global average is skipped

          IF (jm == 0) THEN   ! set up indicies for scaling loop
             isc1 = 1
             isc2 = ijump
          ELSE
             isc1 = jm
             isc2 = 0
          ENDIF

          isc3 = isc2
          DO jn = isc1,ipk
             zbv(izr+isc3) = zbv(izr+isc3) * qf(jn)
             zbv(izi+isc3) = zbv(izi+isc3) * qf(jn)
             zbv(idr+isc3) = zbv(idr+isc3) * qf(jn)
             zbv(idi+isc3) = zbv(idi+isc3) * qf(jn)
             zbv(ihr+isc3) = zbv(ihr+isc3) * zedra
             zbv(ihi+isc3) = zbv(ihi+isc3) * zedra
             isc3 = isc3 + ijump
          ENDDO

          ! Separate the type of the solution for nonlinear terms
          ! no separation, implicit excluded
          ! in projection mode all times explicit

          IF ( lproj .OR. &      ! separation during initialization
&               ((iplim(1,jm+1,jk)+iplim(2,jm+1,jk)) < nnp(jm+1)) ) THEN

             ! case A: projection only, then always explicit method
             ! case B: filtering
             !         only if gravity waves of all meridional index
             !         inside the filter window the implicit method is applicable ?
             !
             ! Explicit solution
             !
             ! explicit nmi for this wavenumber & vertical mode            
             ! the following code is passed through twice,
             ! symmetric flow (jsym = 0), antisymmetric flow (jsym = 1)

             sym_loop : DO jsym = 0,1

                ! Get eigenvalues and eigenvectors of horizontal modes

                ipos = ihmod(1,jsym+1,jm+1,jk)
                ing  = ihmod(2,jsym+1,jm+1,jk)
                inr  = ihmod(3,jsym+1,jm+1,jk)
                imo  = ing + inr

                !======== load ZWORKV (REAL+IMAG) == (IMO,2)
                ! real part      == (VOr,Dr,Hr)*(No of spectral coefficients)
                ! imaginary part == (VOi,Di,Hi)*(No of spectral coefficients)

                zworkv(1:2*imo) = 0.
                   
                ! get number of horizontal normal modes used
                ! ZHMOD(ipos+) contains gravity and rossby eigenvectors ==> matrix A
                ! gravity part = ZHMOD(1:imo,    1:ing)
                ! rossby part  = ZHMOD(1:imo,ing+1:ing+inr)
                ! for projection ZHMOD should be transposed using xGEMM moduls

                SELECT CASE (iproj)
                CASE (NMI_SEL_ALL, NMI_SEL_ALL_W0)
                   ingrx = imo                      ! number of modes used
                   iposx = ipos+imo                 ! start pointer for eigenvectors
                   iposv = ipos                     ! start pointer for eigenvalues

                CASE (NMI_SEL_GM, NMI_SEL_RM_IMPL, NMI_SEL_GM_W0)                   
                   ingrx = ing
                   iposx = ipos+imo
                   iposv = ipos

                CASE (NMI_SEL_RM, NMI_SEL_RM_W0)
                   ingrx = inr
                   iposx = ipos+imo+imo*ing
                   iposv = ipos+ing

                CASE (NMI_SEL_FILT, NMI_SEL_FILT_W0)
                   ingrx = 2*iplim(jsym+1,jm+1,jk)
                   iposx = ipos+imo+imo*(0.5*ing-iplim(jsym+1,jm+1,jk))
                   iposv = ipos+(0.5*ing-iplim(jsym+1,jm+1,jk))

                CASE (NMI_SEL_FILTD, NMI_SEL_FILTD_W0)
                   ingrx = 2*iplimd(jsym+1,jm+1,jk)
                   iposx = ipos+imo+imo*(0.5*ing-iplimd(jsym+1,jm+1,jk))
                   iposv = ipos+(0.5*ing-iplimd(jsym+1,jm+1,jk))

                CASE default
                   WRITE (nmi_mess,'(a,i5)') 'not implemented IPROJ = ',iproj
                   CALL finish ('NMI_Filter',nmi_mess)
                END SELECT

                IF (ingrx /= 0) THEN   ! perform filtering explicit

                   ! Compose vector x
                   CALL setidx ( jsym, ijump, jm,  idx1, idx2, idx3, idx4, idx5, idx6 )
                   DO jn = idx5,ipk,2                       ! load ZETA
                      zworkv(idx2)     = -zbv(izr+idx1)
                      zworkv(idx2+imo) = -zbv(izi+idx1)
                      idx1 = idx1 + 2*ijump
                      idx2 = idx2 + 3
                   ENDDO
                   DO jn = idx6,ipk,2                       ! load DIV and H
                      zworkv(idx4)       = -zbv(idi+idx3)
                      zworkv(idx4+imo)   =  zbv(idr+idx3)
                      zworkv(idx4+1)     =  zbv(ihr+idx3)
                      zworkv(idx4+1+imo) =  zbv(ihi+idx3)
                      idx3 = idx3 + 2*ijump
                      idx4 = idx4 + 3
                   ENDDO

                   ! Project into horizontal normal mode space
                   CALL dgemm &
&                        ('T','N',ingrx,2,imo,1.0,zhmod(iposx),imo,zworkv,imo,0.0,zworkm,ingrx)

                   ! ZWORKV(imo,2) ==>> ZWORKM(ingrx,2)

                   IF (.NOT.lproj) THEN  ! explicit solution for initialisation
                      DO jn = 1,ingrx
                         zs               =  zworkm(ingrx+jn)
                         zo               =  2.*omega*zhmod(iposv+jn-1)
                         zworkm(ingrx+jn) =  zworkm(jn) / zo
                         zworkm(jn)       =  - zs / zo
                      ENDDO
                      ! ZWORKM contains now the nonlinear forcing increment
                   ENDIF

                   ! Modify normal mode coefficients in normal mode space
                   ! Exclude tidal signal from initialisation ... removed
                   ! Project back into vertical space
                   !
                   !  ZWORKM(ingrx,2) ==>> ZWORKV(imo,2)

                   CALL dgemm &
&                        ('N','N',imo,2,ingrx,1.0,zhmod(iposx),imo,zworkm,ingrx,0.0,zworkv,imo)

                ENDIF

                ! Decompose vector x
                CALL setidx ( jsym, ijump, jm, idx1, idx2, idx3, idx4, idx5, idx6 )

                ! store explicit part in ZBV
                DO jn = idx5,ipk,2
                   zbv(izr+idx1) = -zworkv(idx2)
                   zbv(izi+idx1) = -zworkv(idx2+imo)
                   idx1 = idx1 + 2*ijump
                   idx2 = idx2 + 3
                ENDDO
                DO jn = idx6,ipk,2
                   zbv(idi+idx3) = -zworkv(idx4)
                   zbv(idr+idx3) =  zworkv(idx4+imo)
                   zbv(ihr+idx3) =  zworkv(idx4+1)
                   zbv(ihi+idx3) =  zworkv(idx4+1+imo)
                   idx3 = idx3 + 2*ijump
                   idx4 = idx4 + 3
                ENDDO

             END DO sym_loop

             ! end of explicit calculation for given zonal wavenumber
             !==========================
             ! no implicit calculation included
             ! Implicit nmi for this wavenumber & vertical mode
          ENDIF

          ! now ZBV new calculated
          isc3 = isc2                  ! Rescaling
          DO jn = isc1,ipk
             zbv(izr+isc3) = zbv(izr+isc3) / qf(jn)
             zbv(izi+isc3) = zbv(izi+isc3) / qf(jn)
             zbv(idr+isc3) = zbv(idr+isc3) / qf(jn)
             zbv(idi+isc3) = zbv(idi+isc3) / qf(jn)
             zbv(ihr+isc3) = zbv(ihr+isc3) / zedra
             zbv(ihi+isc3) = zbv(ihi+isc3) / zedra
             isc3 = isc3 + ijump
          ENDDO

       END DO vmode_loop

       ! Inverse vertical transform
       CALL dgemm &
&            ('N', 'N', nlev, ivmo, invm, 1.0, vm_ev, nlev, zbv, invm, 0.0, zworks,nlev)

       ! correct global average calculation... removed
       ! alternative separation methode ... removed
       ! calculate changes of surface pressure --> eq. 61 in W.Wergen
       ! ZDB(1,nlev)*ZWORKS(nlev,innp2) = ZBB(1,innp2)
       ! ZBB(1,innp2) ... changes of log(surface pressure)
       ! correct H with orography ... removed
       !
       ! Split H into lnPS and T
       zzdb(:)=zdb(:)
       CALL dgemm &
&            ('T', 'N', 1,innp2,nlev, 1.0, zzdb, nlev, zworks(2*ilt+1), nlev, 0.0, zbb,1)

       ! possible modifications
       ! 1. no pressure changes
       ! 2. calculate changes of lnPs from changes of H

       DO jn = 1,innp2
          !zbs(jn*nlevp1) = 0.0
          zbs(jn*nlevp1) = zbb(jn)
       ENDDO

       IF (jm==0) THEN
          ! no changes in global averaged pressure
          zbs(  nlevp1) = 0.0
          zbs(2*nlevp1) = 0.0
       END IF

       ! Compute changes in T

       DO jn = 1,innp2       !     lowest level
          ish = 2*ilt + jn*nlev          ! pointer  H
          ist =         jn*nlevp1        ! pointer  lnPS
          zbs(ist-1) = zralprr(nlev) * &  ! ==> T(nlev)
&               ( zworks(ish) - &         ! H(nlev)
&               rd*tr*zbs(ist) )          ! R*Tr*lnPS
       ENDDO
       DO jn = 1,innp2       !     remaining levels
          ish = (jn-1)*nlev   + 2*ilt    ! pointer  H
          ist = (jn-1)*nlevp1            ! pointer  T
          DO  jk = nlevm1,1,-1
             zbs(ist+jk) =  zralprr(jk) * ( &                            ! ==> T(level)
&                    zworks(ish+jk) &                                    ! H(level)
&                  - zworks(ish+jk+1) &                                  ! H(level+1)
&                  - zbs(ist+jk+1) * ( zrlnpr(jk+1) - zralphr(jk+1) ) &  ! T(level+1)
&                  )
          ENDDO
       ENDDO

       ! Compute new values

       IF (lpsrs) THEN       !        zero surface pressure changes of all components
          DO jn=1,innp2
             zbs(jn*nlevp1) = 0.
          ENDDO
       ENDIF

       SELECT CASE(iproj)
       CASE(NMI_SEL_FILT, NMI_SEL_FILTD,&
            NMI_SEL_GM,   NMI_SEL_RM,   NMI_SEL_ALL,&
            NMI_SEL_GM_W0,NMI_SEL_RM_W0,NMI_SEL_ALL_W0)        ! store filtered fields
          svo_nmi(isp :isp +ilt -1) = zworks(    1:  ilt)
          sd_nmi (isp :isp +ilt -1) = zworks(ilt+1:2*ilt)
          SELECT CASE(iproj)
          CASE(NMI_SEL_GM_W0, NMI_SEL_RM_W0, NMI_SEL_ALL_W0)
             IF (jm==0) THEN
                ! global average T and lnPs unchanged
                stp_nmi(ismp+2*nlevp1:ismp+iltp-1) = zbs(1+2*nlevp1:iltp)
             ELSE
                stp_nmi(ismp:ismp+iltp-1) = zbs(1:iltp)
             END IF
          CASE default
             stp_nmi(ismp:ismp+iltp-1) = zbs(1:iltp)
          END SELECT

       CASE(NMI_SEL_FILT_W0, NMI_SEL_FILTD_W0, NMI_SEL_RM_IMPL) ! subtract the calculated fraction
          ! VO and D
          svo_nmi(isp :isp +ilt -1) = svo_nmi(isp :isp +ilt -1) - zworks(    1:ilt)
          sd_nmi (isp :isp +ilt -1) = sd_nmi (isp :isp +ilt -1) - zworks(ilt+1:2*ilt)
          ! ln(ps) and T
          stp_nmi(ismp:ismp+iltp-1) = stp_nmi(ismp:ismp+iltp-1) - zbs   (    1:iltp)

       END SELECT

    END DO nmi_waveno_loop

    ! clear work space
    DEALLOCATE (zph)
    DEALLOCATE (zrlnpr)
    DEALLOCATE (zralphr)
    DEALLOCATE (zralprr)
    DEALLOCATE (zbs)
    DEALLOCATE (zworks)
    DEALLOCATE (zbb)
    DEALLOCATE (zbv)
    DEALLOCATE (zworkv)
    DEALLOCATE (zworkm)
    DEALLOCATE (zlmm)
    DEALLOCATE (zlmy)
    DEALLOCATE (ipivot)
    DEALLOCATE (zzdb)

  CONTAINS

    SUBROUTINE setidx ( jsym, jump, m, i1, i2, i3, i4, i5, i6 )
      ! set flow dependent index table
      INTEGER, INTENT(in)  :: jsym, jump, m
      INTEGER, INTENT(out) :: i1, i2, i3, i4, i5, i6

      i1 = (1-jsym)*jump      ! i1 ok
      i2 = 3 - 2*jsym         ! i2 ok
      i3 = jsym*jump          ! i3 ok
      i4 = 1 + jsym           ! i4 ok
      i5 = m + (1-jsym)       ! i5 ok
      i6 = m + jsym           ! i6 ok
      IF (m == 0) THEN
         i4 = 2 - jsym ! ok
         IF (jsym == 0) THEN  ! symmetric flow
            i2 = 1            ! i2 ok
            i3 = 2 * jump     ! i3 ok
            i6 = 2            ! i6 ok
         ELSE                 ! antisymmetric flow
            i1 = 2 * jump     ! i1 ok
            i2 = 3            ! i2 ok
            i5 = 2            ! i5 ok
         ENDIF
      ENDIF
    END SUBROUTINE setidx

  END SUBROUTINE NMI_Filter

!======================================================================
  ! set a unit matrix
  SUBROUTINE unitmx(array,dim)
    REAL                :: array(:,:)   ! out: unit matrix
    INTEGER, INTENT(in) :: dim          ! order of matrix
    INTEGER             :: j
    array(1:dim,1:dim) = 0.0
    DO j=1,dim
       array(j,j) = 1.0
    ENDDO
  END SUBROUTINE unitmx
!======================================================================
  FUNCTION df (m, n) RESULT (d)
    INTEGER, INTENT(IN) :: m, n
    REAL :: d
    d = SQRT(REAL((n**2-m**2))/(4*n**2-1))
  END FUNCTION df
!======================================================================
  FUNCTION bf (m, n) RESULT (b)
    INTEGER, INTENT(IN) :: m, n
    REAL :: b
    b = df(m,n)*SQRT(REAL(n-1)*(n+1)/n**2)
  END FUNCTION bf
!======================================================================
  FUNCTION cf (m, n) RESULT (c)
    INTEGER, INTENT(IN) :: m, n
    REAL :: c
    c = REAL(m)/(n*(n+1))
  END FUNCTION cf
!======================================================================
  FUNCTION ef (n, ped) RESULT (e)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN)    :: ped
    REAL :: e
    e = SQRT(ped*REAL(n)*(n+1))/(2.0*a*omega)
  END FUNCTION ef
!======================================================================
  FUNCTION hf (m, n) RESULT (h)
    INTEGER, INTENT(IN) :: n,m
    REAL  :: h
    h = 2.0*omega*REAL(m)/(REAL(n)*(n+1))
  END FUNCTION hf
!======================================================================
  FUNCTION of (n) RESULT (o)
    INTEGER, INTENT(IN) :: n
    REAL :: o
    o = 2.0*omega*SQRT(REAL(n**2)-1)/REAL(n)
  END FUNCTION of
!======================================================================
  FUNCTION pf(n) RESULT (p)
    INTEGER, INTENT(IN) :: n
    REAL :: p
    p = SQRT(REAL(n)*(n+1))
  END FUNCTION pf
!======================================================================
  FUNCTION qf (n) RESULT (q)
    INTEGER, INTENT(IN)  :: n
    REAL :: q
    q = a**2/pf(n)
  END FUNCTION qf
!======================================================================

  ! load and reload different data sets into NMI work space
  SUBROUTINE NMI_LoadData(ityp)

    INTEGER :: ityp      ! define data storage direction
                         ! >0 load something into NMI arrays
                         ! <0 use NMI arrays as source
    REAL    :: zrdt
    INTEGER :: ia1, ia2
#if (defined CRAY) || (defined SX)
    EXTERNAL util_reshape
#else
    INTRINSIC RESHAPE
#endif

    ia1 = nlev  *n2sp
    ia2 = nlevp1*n2sp

    SELECT CASE (ityp)
    CASE (NMI_LOAD_CLEAR)            ! clear NMI data
       svo_nmi(:) = 0.0
       sd_nmi (:) = 0.0
       stp_nmi(:) = 0.0

    CASE (NMI_LOAD_ECH_BUFA)         ! load ECHAM spectral arrays
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo_nmi, svo,ia1)
       CALL util_reshape(sd_nmi,  sd, ia1)
       CALL util_reshape(stp_nmi, stp,ia2)
#else
       svo_nmi(:) = RESHAPE(svo(:,:,:),(/ia1/))
       sd_nmi (:) = RESHAPE(sd (:,:,:),(/ia1/))
       stp_nmi(:) = RESHAPE(stp(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ECH_BUFA)        ! reload ECHAM spectral arrays
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo, svo_nmi,ia1)
       CALL util_reshape(sd,  sd_nmi, ia1)
       CALL util_reshape(stp, stp_nmi,ia2)
#else
       svo(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,nsp/))
       sd (:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,nsp/))
       stp(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,nsp/))
#endif

    CASE (NMI_LOAD_ECH_BUFB)         ! load data in second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo_dt, svo,ia1)
       CALL util_reshape(sd_dt,  sd, ia1)
       CALL util_reshape(stp_dt, stp,ia2)
#else
       svo_dt(:) = RESHAPE(svo(:,:,:),(/ia1/))
       sd_dt (:) = RESHAPE(sd (:,:,:),(/ia1/))
       stp_dt(:) = RESHAPE(stp(:,:,:),(/ia2/))
#endif   

    CASE (-NMI_LOAD_ECH_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo, svo_dt,ia1)
       CALL util_reshape(sd,  sd_dt, ia1)
       CALL util_reshape(stp, stp_dt,ia2)
#else
       svo(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,nsp/))
       sd (:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,nsp/))
       stp(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,nsp/))
#endif

    CASE (NMI_LOAD_ORD_BUFA)         ! load nudging data in first buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo_nmi, svobs3,ia1)
       CALL util_reshape(sd_nmi,  sdobs3,ia1)
       CALL util_reshape(stp_nmi, stobs3,ia2)
#else
       svo_nmi(:) = RESHAPE(svobs3(:,:,:),(/ia1/))
       sd_nmi (:) = RESHAPE(sdobs3(:,:,:),(/ia1/))
       stp_nmi(:) = RESHAPE(stobs3(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ORD_BUFA)        ! reload nudging spectral arrays from first buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svobs3, svo_nmi,ia1)
       CALL util_reshape(sdobs3, sd_nmi, ia1)
       CALL util_reshape(stobs3, stp_nmi,ia2)
#else
       svobs3(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,nsp/))
       sdobs3(:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,nsp/))
       stobs3(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,nsp/))
#endif

    CASE (NMI_LOAD_ORD_BUFB)         ! load nudging data in second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo_dt, svobs3,ia1)
       CALL util_reshape(sd_dt,  sdobs3,ia1)
       CALL util_reshape(stp_dt, stobs3,ia2)
#else
       svo_dt(:) = RESHAPE(svobs3(:,:,:),(/ia1/))
       sd_dt (:) = RESHAPE(sdobs3(:,:,:),(/ia1/))
       stp_dt(:) = RESHAPE(stobs3(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_ORD_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svobs3, svo_dt,ia1)
       CALL util_reshape(sdobs3, sd_dt, ia1)
       CALL util_reshape(stobs3, stp_dt,ia2)
#else
       svobs3(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,nsp/))
       sdobs3(:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,nsp/))
       stobs3(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,nsp/))
#endif

    CASE (NMI_LOAD_OIN_BUFB)         ! load nudging data in second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svo_dt, svobs,ia1)
       CALL util_reshape(sd_dt,  sdobs,ia1)
       CALL util_reshape(stp_dt, stobs,ia2)
#else
       svo_dt(:) = RESHAPE(svobs(:,:,:),(/ia1/))
       sd_dt (:) = RESHAPE(sdobs(:,:,:),(/ia1/))
       stp_dt(:) = RESHAPE(stobs(:,:,:),(/ia2/))
#endif

    CASE (-NMI_LOAD_OIN_BUFB)        ! reload ECHAM spectral arrays from second buffer
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svobs, svo_dt,ia1)
       CALL util_reshape(sdobs, sd_dt, ia1)
       CALL util_reshape(stobs, stp_dt,ia2)
#else
       svobs(:,:,:) = RESHAPE(svo_dt(:),(/nlev  ,2,nsp/))
       sdobs(:,:,:) = RESHAPE(sd_dt (:),(/nlev  ,2,nsp/))
       stobs(:,:,:) = RESHAPE(stp_dt(:),(/nlevp1,2,nsp/))
#endif

    CASE (NMI_COPY_AB)        ! copy NMI buffer  TEN -->> NMI
       svo_nmi(:) = svo_dt(:)
       sd_nmi (:) = sd_dt (:)
       stp_nmi(:) = stp_dt(:)

    CASE (NMI_COPY_BA)        ! copy NMI buffer NMI -->> TEN
       svo_dt(:) = svo_nmi(:)
       sd_dt (:) = sd_nmi (:)
       stp_dt(:) = stp_nmi(:)

    CASE (NMI_CALC_TEND)      ! calculate tendencies
       zrdt = 2.0 / twodt
       svo_nmi(:) = zrdt * (svo_nmi(:) - svo_dt(:))
       sd_nmi (:) = zrdt * (sd_nmi (:) - sd_dt (:))
       stp_nmi(:) = zrdt * (stp_nmi(:) - stp_dt(:))

    CASE (NMI_CALC_DIFF)       ! calculate difference
       svo_nmi(:) = (svo_nmi(:) - svo_dt(:))
       sd_nmi (:) = (sd_nmi (:) - sd_dt (:))
       stp_nmi(:) = (stp_nmi(:) - stp_dt(:))

    CASE (NMI_CALC_ADD)        ! add both NMI buffer, result in second buffer
       svo_dt(:) = svo_nmi(:) + svo_dt(:)
       sd_dt (:) = sd_nmi (:) + sd_dt (:)
       stp_dt(:) = stp_nmi(:) + stp_dt(:)

    case (NMI_STORE_FAST)      ! store filtered fields
       if (lfill_a) then
          svfast_accu(:,:,:) = svfast_accu(:,:,:) - svfast_a(:,:,:)
          sdfast_accu(:,:,:) = sdfast_accu(:,:,:) - sdfast_a(:,:,:)
          stfast_accu(:,:,:) = stfast_accu(:,:,:) - stfast_a(:,:,:)
          ifast_accu = ifast_accu + 1
       end if
       if (lfill_b) then
          svfast_a(:,:,:) = svfast_b(:,:,:)
          sdfast_a(:,:,:) = sdfast_b(:,:,:)
          stfast_a(:,:,:) = stfast_b(:,:,:)
       end if
#if (defined CRAY) || (defined SX)
       CALL util_reshape(svfast_b, svo_nmi,ia1)
       CALL util_reshape(sdfast_b, sd_nmi, ia1)
       CALL util_reshape(stfast_b, stp_nmi,ia2)
#else
       svfast_b(:,:,:) = RESHAPE(svo_nmi(:),(/nlev  ,2,nsp/))
       sdfast_b(:,:,:) = RESHAPE(sd_nmi (:),(/nlev  ,2,nsp/))
       stfast_b(:,:,:) = RESHAPE(stp_nmi(:),(/nlevp1,2,nsp/))
#endif
       if (lfill_a) then
          svfast_accu(:,:,:) = svfast_accu(:,:,:) + svfast_b(:,:,:)
          sdfast_accu(:,:,:) = sdfast_accu(:,:,:) + sdfast_b(:,:,:)
          stfast_accu(:,:,:) = stfast_accu(:,:,:) + stfast_b(:,:,:)
       end if
       if (lfill_b) lfill_a = .true.
       lfill_b = .true.

    CASE (NMI_CALC_DIAG)       ! make diagnostics
       CALL NMI_Diagnose

    CASE default
       CALL message('NMI_LoadData','illegal ityp selection')
       CALL finish('NMI_LoadData','Run aborted.')

    END SELECT

  END SUBROUTINE NMI_LoadData

  SUBROUTINE NMI_Diagnose
    ! calculate diagnostics of nmi-buffer
    ! norm in spectral space see details in mo_nudging

    INTEGER :: idxr, idxi ! fieldindex real and imaginary
    INTEGER :: iwz        ! zonal wavenumber
    INTEGER :: il,i1,i2,is
    REAL    :: sntp, snvo, sndi,delta

    INTEGER, PARAMETER :: BUFOPER_STP=1
    INTEGER, PARAMETER :: BUFOPER_DIV=2
    INTEGER, PARAMETER :: BUFOPER_VOR=3

    CALL message('NMI_Diagnose','Spectral space diagnostics:')
    CALL message('','   RMS[ abs(nofilter) - abs(filter) ] over meridional numbers')
    DO iwz = 0,nm
       WRITE(nmi_mess,*) ' ** zonal wave number ',iwz
       CALL message('',nmi_mess)
       i1 = nmp(iwz+1)*2*nlev
       i2 = nmp(iwz+1)*2*nlevp1
       DO il = 1,nlev
          sntp = 0.0
          snvo = 0.0
          sndi = 0.0
          DO is = 1,nnp(iwz+1)
             idxr = i1 + (is-1)*2*nlevp1 + il; idxi = idxr + nlevp1
             delta = diag_norm(BUFOPER_STP)-diag_norm(-BUFOPER_STP)
             sntp  = sntp + delta*delta             ! temperature and lnps
             idxr = i1 + (is-1)*2*nlev + il;   idxi = idxr + nlev
             delta = diag_norm(BUFOPER_DIV)-diag_norm(-BUFOPER_DIV)
             sndi  = sndi + delta*delta             ! divergence
             delta = diag_norm(BUFOPER_VOR)-diag_norm(-BUFOPER_VOR)
             snvo  = snvo + delta*delta             ! vorticity
          ENDDO
          sntp = SQRT(sntp/nnp(iwz+1))
          sndi = SQRT(sndi/nnp(iwz+1))
          snvo = SQRT(snvo/nnp(iwz+1))
          WRITE(nmi_mess,*) il,' temp ',sntp,' div ',sndi,' vor ',snvo
          CALL message('',nmi_mess)
       ENDDO
    ENDDO

  CONTAINS

    FUNCTION Diag_Norm(par)
      REAL  :: diag_norm
      INTEGER :: par
      SELECT CASE(par)
         ! buffer A operations
      CASE(BUFOPER_STP); diag_norm = stp_nmi(idxr)*stp_nmi(idxr) + stp_nmi(idxi)*stp_nmi(idxi)
      CASE(BUFOPER_DIV); diag_norm = sd_nmi (idxr)*sd_nmi (idxr) + sd_nmi (idxi)*sd_nmi (idxi)
      CASE(BUFOPER_VOR); diag_norm = svo_nmi(idxr)*svo_nmi(idxr) + svo_nmi(idxi)*svo_nmi(idxi)
         ! buffer B operations
      CASE(-BUFOPER_STP); diag_norm = stp_dt(idxr)*stp_dt(idxr) + stp_dt(idxi)*stp_dt(idxi)
      CASE(-BUFOPER_DIV); diag_norm = sd_dt (idxr)*sd_dt (idxr) + sd_dt (idxi)*sd_dt (idxi)
      CASE(-BUFOPER_VOR); diag_norm = svo_dt(idxr)*svo_dt(idxr) + svo_dt(idxi)*svo_dt(idxi)
      CASE default; diag_norm = 0.0
      END SELECT
    END FUNCTION Diag_Norm
  END SUBROUTINE NMI_Diagnose

END MODULE mo_nmi
