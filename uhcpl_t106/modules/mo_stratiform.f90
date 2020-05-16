MODULE mo_stratiform
  ! collect inputs: convective precipitation (W/m**2)*24*3600*1000 (to mm/day)
  !                 temperature tendency due to convection (??) 
  !                 moisture tendency due to convection (??)
  !         output: dqmeso(nlon,ngl), dtmeso(nlon,ngl): moisture/temperature tendencies  
  !                                              due to stratiform rainfall
  ! Joshua Xiouhua Fu, IPRC/UH, Honolulu, Jun/2008
  !
  !
  ! INPUT will be collected from cudtdq.f90 
  !
  ! OUTPUT will be passed to cudtdq.f90 too 
  !
  !
  !

  USE mo_exception,     ONLY: finish   
  USE mo_mpi,           ONLY: p_pe, p_io   !fu++
  USE mo_doctor,        ONLY: nout
  USE mo_kind,          ONLY: dp

  USE mo_memory_g3a,    ONLY: slmm

  USE mo_transpose,     ONLY: gather_gp, scatter_gp
  USE mo_decomposition, ONLY: dc=>local_decomposition, dcg=>global_decomposition
  USE mo_control,       ONLY: dtime, lctime, lwtime, nctime,  &
                           ngl, nlon, nlev, imtt, ncbase, ntbase, lonudg

!fu++
! USE mo_couple,        ONLY: accu_time
!

  IMPLICIT NONE

  PRIVATE

  REAL(dp), POINTER, PUBLIC, SAVE :: dqmeso(:,:,:) !q tendency by meso-stratiform 
  REAL(dp), POINTER, PUBLIC, SAVE :: dtmeso(:,:,:) !t tendency by meso-stratiform 

  PUBLIC :: meso_init  !initializing space for meso-stratiform,call in control
  PUBLIC :: meso_heat  !calculate t,q tendencies by meso, call in cudtdq 

  REAL, POINTER :: slmfu2(:,:)
  REAL, POINTER :: gl_meso_dt(:,:,:)
  REAL, POINTER :: gl_meso_dq(:,:,:)

!REAL,TARGET,ALLOCATABLE :: sstfu1(:,:,:)

  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fprecip(:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fconvq(:,:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fconvt(:,:,:)
  REAL, ALLOCATABLE, SAVE, PUBLIC, TARGET :: fpress(:,:,:)
  
! local fields for meso_collect

! pick up in CUDTDQ (cudtdq.f90)
  REAL, ALLOCATABLE,PUBLIC :: prconv_meso(:,:)    !spontanous value 
  REAL, ALLOCATABLE,PUBLIC :: prconv_dq(:,:,:)    !spontanous value 
  REAL, ALLOCATABLE,PUBLIC :: prconv_dt(:,:,:)    !spontanous value 
! pick up in PHYSC (physc.f90)
  REAL, ALLOCATABLE,PUBLIC :: press_meso(:,:,:)    !spontanous value 
  
!pick up in COND (m_cond.f90)
!REAL,ALLOCATABLE,PUBLIC :: prlarg_diurnal(:,:) !accumulated
!pick up in RADHEAT (m_radheat.f90)
!REAL, ALLOCATABLE, PUBLIC :: srads_diurnal(:,:) !spontanous value
  
  REAL, ALLOCATABLE, PUBLIC, TARGET :: dqmes(:,:,:)   !q tendency   
  REAL, ALLOCATABLE, PUBLIC, TARGET :: dtmes(:,:,:)   !t tendency   

  INTEGER  :: ioerr, i, j, k, itop
  REAL     :: aa, alft, alfq, pi, pb, pt 
  REAL*4, ALLOCATABLE   :: fprecip_out(:,:)
  REAL*4, ALLOCATABLE   :: fconvq_out(:,:,:)
  REAL*4, ALLOCATABLE   :: fconvt_out(:,:,:)
  REAL*4, ALLOCATABLE   :: fpress_out(:,:,:)
  REAL*4, ALLOCATABLE   :: dtmes_out(:,:,:)
  REAL*4, ALLOCATABLE   :: dqmes_out(:,:,:)

!fu++  LOGICAL, PUBLIC :: ioflag =.FALSE.            !do computing on i/o PE

CONTAINS

!================================================================================
! initialisation of meso-init  

  SUBROUTINE meso_init


!fu++    ioflag = (p_parallel .AND. p_parallel_io) .OR. (.NOT.p_parallel)

    ALLOCATE (dqmeso(dc%nglon,nlev,dc%nglat))
    ALLOCATE (dtmeso(dc%nglon,nlev,dc%nglat))

    IF (p_pe == p_io) THEN
       ALLOCATE (slmfu2(nlon,ngl))
       ALLOCATE (dqmes(nlon,nlev,ngl))
       ALLOCATE (dtmes(nlon,nlev,ngl))
!
!      ALLOCATE (sstfu1(nlon,ngl,0:0))
!      dsst = 0.0
!      sstfu1(:,:,0) = dsst(:,:)
!      gl_dsst  => sstfu1(:,:,0:0)
!
    END IF
!
!scatter dsst: SST diurnal cycle
!
!    CALL scatter_gp (gl_dsst, dssto(:,:,0:0), dcg)

    CALL gather_gp (slmfu2,slmm,dcg)  !1:Land; 0:Sea


!allocate collected data
    ALLOCATE (prconv_meso(dc%nglon, dc%nglat))
    ALLOCATE (prconv_dq(dc%nglon, nlev, dc%nglat))
    ALLOCATE (prconv_dt(dc%nglon, nlev, dc%nglat))
    ALLOCATE (press_meso(dc%nglon, nlev, dc%nglat))

    prconv_meso(:,:) = 0.0
    prconv_dq(:,:,:) = 0.0
    prconv_dt(:,:,:) = 0.0
    press_meso(:,:,:) = 0.0

  IF (p_pe == p_io) THEN
   write(nout,*) 'FU++,nlon,ngl,nlev=',nlon,ngl,nlev

    IF(.NOT.ALLOCATED(fprecip)) ALLOCATE(fprecip(nlon,ngl))
    IF(.NOT.ALLOCATED(fconvq))  ALLOCATE(fconvq(nlon,nlev,ngl))
    IF(.NOT.ALLOCATED(fconvt))  ALLOCATE(fconvt(nlon,nlev,ngl))
    IF(.NOT.ALLOCATED(fpress))  ALLOCATE(fpress(nlon,nlev,ngl))

    IF(.NOT.ALLOCATED(fprecip_out)) ALLOCATE(fprecip_out(nlon,ngl))
    IF(.NOT.ALLOCATED(fconvq_out))  ALLOCATE(fconvq_out(nlon,ngl,nlev))
    IF(.NOT.ALLOCATED(fconvt_out))  ALLOCATE(fconvt_out(nlon,ngl,nlev))
    IF(.NOT.ALLOCATED(fpress_out))  ALLOCATE(fpress_out(nlon,ngl,nlev))
    IF(.NOT.ALLOCATED(dtmes_out))  ALLOCATE(dtmes_out(nlon,ngl,nlev))
    IF(.NOT.ALLOCATED(dqmes_out))  ALLOCATE(dqmes_out(nlon,ngl,nlev))

    Open(155, file='./conv_heat.dat',form='unformatted',&
        access='sequential')
  END IF

  END SUBROUTINE meso_init

!================================================================================
! collect accumualte data (wind/rain) only for SST diurnal
!
!SUBROUTINE meso_collect
!
!   USE mo_transpose,     ONLY: gather_gp
!
!    REAL          :: quot
!    REAL, POINTER :: collect_field(:,:)
!
!gather input arrays
!WIND10 10 meter wind speed [m/s]
!    collect_field => fwind1(:,:)
!    CALL gather_gp(collect_field,wind10_diurnal,dcg)
!
!PRCONV convective rain flux [W/m*m]
!    collect_field => frainc1(:,:)
!    CALL gather_gp(collect_field,prconv_diurnal,dcg)
!
!PRLARG large-scale rain flux [W/m*m] 
!    collect_field => frainl1(:,:)
!    CALL gather_gp(collect_field,prlarg_diurnal,dcg)
!
!    IF (p_pe == p_io) THEN
!
!       IF ( accu_time > 0.0 ) THEN
!        quot = 1.0/accu_time
!       ELSE
!        CALL finish('diurnal_collect','wrong setting of accumulation time')
!       END IF
!
!       WRITE(nout,*) 'accum-time [sec]&factor[1/sec]',accu_time,quot
!
!       fwind1(:,:)  = fwind1(:,:) *quot
!       frainc1(:,:) = frainc1(:,:)*quot
!       frainl1(:,:) = frainl1(:,:)*quot
!
!    END IF
!
! zero exchange fields
!
!    wind10_diurnal(:,:) = 0.0
!    prconv_diurnal(:,:) = 0.0
!    prlarg_diurnal(:,:) = 0.0
!
!    RETURN
! END SUBROUTINE meso_collect
!
!================================================================================
!collect prconv_meso, prconv_dq, prconv_dt 
!
  SUBROUTINE meso_heat

    REAL :: frain(nlon,ngl)
    REAL :: lstrat(nlon,ngl),qmax(nlon,ngl),tmax(nlon,ngl)
    REAL :: zq(nlev),zt(nlev)
    INTEGER :: ktop(nlon,ngl)

    REAL, POINTER :: collect_field2(:,:)
    REAL, POINTER :: collect_field(:,:,:)

!   IF (lctime) THEN
!    write(nout,*) "collect input data for sst_diurnal!!"
!      CALL diurnal_collect
!    END IF
!
!prconv_meso in CUDTDQ.F90:convective rainfall [W/m*m]
    collect_field2 => fprecip(:,:)
    CALL gather_gp(collect_field2,prconv_meso,dcg)
!prconv_dq in CUDTDQ.F90:convective moist tendency [??]
    collect_field => fconvq(:,:,:)
    CALL gather_gp(collect_field,prconv_dq,dcg)
!prconv_dt in CUDTDQ.F90:convective temperature tendency [??]
    collect_field => fconvt(:,:,:)
    CALL gather_gp(collect_field,prconv_dt,dcg)
!press_meso in PHYSC.F90:level pressures [??]
    collect_field => fpress(:,:,:)
    CALL gather_gp(collect_field,press_meso,dcg)

 IF (p_pe == p_io) THEN

    !get total rain-fall

    frain = fprecip*24.*3600.*1000. !transfer to mm/day,NOTYET 
    frain = frain*0.5/1000.   !(timestep=0.5hr; rhoh2o=1000.)

    do i=1,nlon
    do j=1,ngl
    fprecip_out(i,j)=frain(i,j)
    end do
    end do


!==============================================================!
!    WRITE(nout,*) 'TOTAL RAINFALL at (160E,EQ):', frain(43,25)
!calulate stratiform heating/moistening as a corrector 
!reference: Moncrieff and Liu 2006, JAS
!==============================================================!
!
!trigger mesoscale-forcing based on: rainfall and upper-troposphere
! convective heating (corresponding to cloud-top):prconv_dt 
!use vertically integrated prconv_dt/prconv_dq and their vertical
! structures to construct meso heating/moistening
!
      alft=0.0001
      alfq=0.0
      pi=3.1415926

      do i=1,nlon
      do j=1,ngl
      lstrat(i,j)=0.
      end do
      end do
     
      do i=1,nlon
      do j=15,35
!
!convection is higher than lev=12 (600 mb)
!
      if(fconvt(i,12,j).gt.0.0) then
      lstrat(i,j)=1.

!get maximum heating/moistening
      do k=1,17  !third bottom level (~950mb)
      zq(k)=fconvq(i,k,j)
      zt(k)=fconvt(i,k,j)
      end do
      qmax(i,j)=MAXVAL(-zq)
      tmax(i,j)=MAXVAL(zt)
   
!identfy the cloud top
      do k=1,12
      if (fconvt(i,k,j).gt.0) then
      ktop(i,j)=k-1
      goto 777
      end if
      end do
777   end if

      end do
      end do
 
      do k=1,nlev 
      do i=1,nlon
      do j=1,ngl 
      dqmes(i,k,j)=0.0
      dtmes(i,k,j)=0.0
      end do
      end do
      end do

      do i=1,nlon
      do j=1,ngl
!=====>convection is strong enough
      if(lstrat(i,j).gt.0.5) then
      pb=fpress(i,17,j)
      itop=ktop(i,j)
      pt=fpress(i,itop,j)
!      pstar=(pb-pt)/2.

      do k=itop,17
!     if(fpress(i,k,j).ge.pstar) then
      aa=2.*pi*(fpress(i,k,j)-pb)/(pb-pt)    !lower-level

      dtmes(i,k,j)=alft*tmax(i,j)*sin(aa)
      dqmes(i,k,j)=-alfq*qmax(i,j)*sin(aa)
      end do

      end if

      end do
      end do

    do k=1,ngl
    do j=1,nlev
    do i=1,nlon
    fconvq_out(i,k,j)=fconvq(i,j,k) 
    fconvt_out(i,k,j)=fconvt(i,j,k) 
    fpress_out(i,k,j)=fpress(i,j,k) 
    dtmes_out(i,k,j)=dtmes(i,j,k) 
    dqmes_out(i,k,j)=dqmes(i,j,k) 
    end do
    end do
    end do

!============================================================!
! save heaing/moistening profiles every 30 mins
!
!    write(155) ((fprecip_out(i,j),i=1,nlon),j=1,ngl)
!    do k=1,nlev
!    write(155) ((fconvq_out(i,j,k),i=1,nlon),j=1,ngl)
!    end do
!    do k=1,nlev
!    write(155) ((fconvt_out(i,j,k),i=1,nlon),j=1,ngl)
!    end do
!    do k=1,nlev
!    write(155) ((fpress_out(i,j,k),i=1,nlon),j=1,ngl)
!    end do
!    do k=1,nlev
!    write(155) ((dtmes_out(i,j,k),i=1,nlon),j=1,ngl)
!    end do
!    do k=1,nlev
!    write(155) ((dqmes_out(i,j,k),i=1,nlon),j=1,ngl)
!    end do
!
!===============================================================!
!      IF ((fsolfx1(i,j).gt.0.0).and.(slmfu1(i,j).le.0.5)) THEN
!
! Day Time/Over Ocean
! 
! fu++ 
!    WRITE(nout,*) 'dsst:SST_diurnal_cycle ',dsst(24,24)
!    WRITE(nout,*) 'Total rain MAX,MIN',MAXVAL(frain),MINVAL(frain)
!    WRITE(nout,*) 'Total wind MAX,MIN',MAXVAL(fwind1),MINVAL(fwind1)
!    WRITE(nout,*) 'Total solar MAX,MIN',MAXVAL(fsolfx1),MINVAL(fsolfx1)
!    WRITE(nout,*) 'Total dsst MAX,MIN',MAXVAL(dsst),MINVAL(dsst)
!======================================================================c
! save total heating/moistening 
! save convective heating/moistening
! save meso-scale heating/moistening
!
!    dsst_fu = dsst
!
!    write(133) ((dsst_fu(i,j),i=1,nlon),j=1,ngl)
!
!define pointers
!     sstfu1(:,:,0) = dsst(:,:)
!
    gl_meso_dt  => dtmes(:,:,:)
    gl_meso_dq  => dqmes(:,:,:)
!
    END IF
!
!scatter dtmes, dqmes 
!
    CALL scatter_gp (gl_meso_dt,  dtmeso(:,:,:), dcg)
    CALL scatter_gp (gl_meso_dq,  dqmeso(:,:,:), dcg)
!
!
END SUBROUTINE meso_heat

END MODULE mo_stratiform

