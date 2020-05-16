MODULE mo_dsst
  ! collect inputs: daily-mean surface wind speed (m/s),
  !                 daily-mean precipitation (mm/h), 
  !                 current net surface solar radiation (W/m**2)
  !         output: dsst (nlon,ngl)
  !
  ! Joshua Xiouhua Fu, IPRC/UH, Honolulu, Jun/2008
  !
  !
  ! transfer from atmospher to ocean
  !
  !   wind speed 10 meter
  !   solar, precipitation 
  !
  ! transfer from ocean to atmosphere
  !
  !   SST diurnal deviation in tropics
  !

  USE mo_exception,     ONLY: finish   
  USE mo_mpi,           ONLY: p_pe, p_io   !fu++
  USE mo_doctor,        ONLY: nout
  USE mo_kind,          ONLY: dp

  USE mo_memory_g3a,    ONLY: slmm

  USE mo_transpose,     ONLY: gather_gp, scatter_gp
  USE mo_decomposition, ONLY: dc=>local_decomposition, dcg=>global_decomposition
  USE mo_control,       ONLY: dtime, lctime, lwtime, nctime,  &
                           ngl, nlon, lcouple, imtt, ncbase, ntbase, lonudg

!fu++
  USE mo_couple,        ONLY: accu_time

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER, PUBLIC, SAVE :: dssto(:,:,:) !diurnal SST

  PUBLIC :: diurnal_init  !initializing space,call in control
  PUBLIC :: sst_diurnal   !calculate SST diurnal deviation,call in stepon

  REAL, POINTER :: slmfu1(:,:)
  REAL, POINTER :: gl_dsst(:,:,:)

  REAL,TARGET,ALLOCATABLE :: sstfu1(:,:,:)

  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fsolfx1(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: frainc1(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: frainl1(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fwind1(:,:)
  
  ! local fields for diurnal_collect

  ! pick up in VDIFF (m_vdiff.f90)
  REAL,ALLOCATABLE,PUBLIC :: wind10_diurnal(:,:) !accumulated
  ! pick up in CUDTDQ (cudtdq.f90)
  REAL, ALLOCATABLE,PUBLIC :: prconv_diurnal(:,:) !accumulated
  ! pick up in COND (m_cond.f90)
  REAL,ALLOCATABLE,PUBLIC :: prlarg_diurnal(:,:) !accumulated

  ! pick up in RADHEAT (m_radheat.f90)
  REAL, ALLOCATABLE, PUBLIC :: srads_diurnal(:,:) !spontanous value
  
  REAL, ALLOCATABLE, PUBLIC :: dsst(:,:)   !diurnal SST deviation   

  INTEGER  :: ioerr, i, j
  REAL     :: aa, bb, cc, dd, ee, ff
  REAL*4, ALLOCATABLE   :: dsst_fu(:,:)

!fu++  LOGICAL, PUBLIC :: ioflag =.FALSE.            !do computing on i/o PE

CONTAINS

  !================================================================================
  ! initialisation of SST-diurnal-cycle  

  SUBROUTINE diurnal_init


!fu++    ioflag = (p_parallel .AND. p_parallel_io) .OR. (.NOT.p_parallel)

    ALLOCATE (dssto(dc%nglon,dc%nglat,0:0))

    IF (p_pe == p_io) THEN
       ALLOCATE (slmfu1(nlon,ngl))
       ALLOCATE (dsst_fu(nlon,ngl))
       ALLOCATE (dsst(nlon,ngl))

       ALLOCATE (sstfu1(nlon,ngl,0:0))
       dsst = 0.0

       sstfu1(:,:,0) = dsst(:,:)

       gl_dsst  => sstfu1(:,:,0:0)

    END IF
!
!scatter dsst: SST diurnal cycle
!
    CALL scatter_gp (gl_dsst, dssto(:,:,0:0), dcg)

    CALL gather_gp (slmfu1,slmm,dcg)


    ! allocate ocean transfer data
    ALLOCATE (wind10_diurnal (dc%nglon, dc%nglat))
    ALLOCATE (prconv_diurnal (dc%nglon, dc%nglat))
    ALLOCATE (prlarg_diurnal (dc%nglon, dc%nglat))
    ALLOCATE ( srads_diurnal (dc%nglon, dc%nglat))

    wind10_diurnal  (:,:) = 0.0
    prconv_diurnal  (:,:) = 0.0
    prlarg_diurnal  (:,:) = 0.0
     srads_diurnal  (:,:) = 0.0

    IF (p_pe == p_io) THEN
      write(nout,*) 'ldsst=.true.,turn on SST diurnal cycle; fu++'
      IF(.NOT.ALLOCATED(fsolfx1)) ALLOCATE(fsolfx1(nlon,ngl))
      IF(.NOT.ALLOCATED(fwind1))  ALLOCATE(fwind1(nlon,ngl))
      IF(.NOT.ALLOCATED(frainc1))  ALLOCATE(frainc1(nlon,ngl))
      IF(.NOT.ALLOCATED(frainl1))  ALLOCATE(frainl1(nlon,ngl))
       fwind1 = 3.5                       !m/s
       frainc1 = 0.4/(3600.*1000)         !W/m**2
       frainl1 = 0.4/(3600.*1000)         !W/m**2

       Open(133, file='./sst_diurnal.dat',form='unformatted',&
           access='sequential')
    END IF

  END SUBROUTINE diurnal_init

  !================================================================================
  ! collect accumualte data (wind/rain) only for SST diurnal

  SUBROUTINE diurnal_collect

    USE mo_transpose,     ONLY: gather_gp

    REAL          :: quot
    REAL, POINTER :: collect_field(:,:)

    ! gather input arrays
    ! WIND10 10 meter wind speed [m/s]
    collect_field => fwind1(:,:)
    CALL gather_gp(collect_field,wind10_diurnal,dcg)

    ! PRCONV convective rain flux [W/m*m]
    collect_field => frainc1(:,:)
    CALL gather_gp(collect_field,prconv_diurnal,dcg)

    ! PRLARG large-scale rain flux [W/m*m] 
    collect_field => frainl1(:,:)
    CALL gather_gp(collect_field,prlarg_diurnal,dcg)

    IF (p_pe == p_io) THEN

       IF ( accu_time > 0.0 ) THEN
        quot = 1.0/accu_time
       ELSE
        CALL finish('diurnal_collect','wrong setting of accumulation time')
       END IF

       WRITE(nout,*) 'accum-time [sec]&factor[1/sec]',accu_time,quot

       fwind1(:,:)  = fwind1(:,:) *quot
       frainc1(:,:) = frainc1(:,:)*quot
       frainl1(:,:) = frainl1(:,:)*quot

    END IF

    ! zero exchange fields

    wind10_diurnal(:,:) = 0.0
    prconv_diurnal(:,:) = 0.0
    prlarg_diurnal(:,:) = 0.0

    RETURN
  END SUBROUTINE diurnal_collect
!
!================================================================================
! calculate SST diurnal-deviation

  SUBROUTINE sst_diurnal

    REAL :: frain(nlon,ngl)
    REAL, POINTER :: collect_field(:,:)

    IF (lctime) THEN
    write(nout,*) "collect input data for sst_diurnal!!"
      CALL diurnal_collect
    END IF

    ! SRADS solar surface flux [W/m*m]
    collect_field => fsolfx1(:,:)
    CALL gather_gp(collect_field,srads_diurnal,dcg)

   IF (p_pe == p_io) THEN

    !get total rain-fall

    frain = (frainc1+ frainl1)*3600.*1000.   !transfer to mm/h 

!    WRITE(nout,*) 'TOTAL RAINFALL at (160E,EQ):', frain(43,25)
   !calulate SST diurnal deviation
   !reference: Clayson and Curry, 1996, JGR
   !

      DO i=1,nlon
      DO j=1,ngl 

      dsst(i,j)=0.0

      IF ((fsolfx1(i,j).gt.0.0).and.(slmfu1(i,j).le.0.5)) THEN
!
! Day Time/Over Ocean
! 
      IF (fwind1(i,j).GE.2.0) THEN 
      aa=0.262
      bb=0.00265
      cc=0.028
      dd=-0.838
      ee=-0.00105
      ff=0.158
      ELSE
      aa=0.328
      bb=0.002
      cc=0.041
      dd=0.212
      ee=-0.000185
      ff=-0.329
      END IF
    
      dsst(i,j)=aa+bb*fsolfx1(i,j)+cc*frain(i,j)+ &
        dd*alog(fwind1(i,j))+ee*fsolfx1(i,j)*alog(fwind1(i,j))+ &
        ff*fwind1(i,j)
      IF (dsst(i,j).LE.0.0) dsst(i,j)=0.0

      ELSE

!night time/over land
      dsst(i,j)=0.

      END IF

      END DO
      END DO

! fu++ 
!    WRITE(nout,*) 'dsst:SST_diurnal_cycle ',dsst(24,24)
!    WRITE(nout,*) 'Total rain MAX,MIN',MAXVAL(frain),MINVAL(frain)
!    WRITE(nout,*) 'Total wind MAX,MIN',MAXVAL(fwind1),MINVAL(fwind1)
!    WRITE(nout,*) 'Total solar MAX,MIN',MAXVAL(fsolfx1),MINVAL(fsolfx1)
!    WRITE(nout,*) 'Total dsst MAX,MIN',MAXVAL(dsst),MINVAL(dsst)
!======================================================================c
! save diurnal cycle
!
    dsst_fu = dsst

    write(133) ((dsst_fu(i,j),i=1,nlon),j=1,ngl)

!define pointers
     sstfu1(:,:,0) = dsst(:,:)

    gl_dsst  => sstfu1(:,:,0:0)

    END IF
!
!scatter dsst: SST diurnal cycle
!
    CALL scatter_gp (gl_dsst,  dssto(:,:,0:0), dcg)


  END SUBROUTINE sst_diurnal

END MODULE mo_dsst

