MODULE mo_memory_g3b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base, ONLY: create_list, delete_list, new_entry, get_entry,&
                            print_memory_table, print_memory_use, print_sinfo,&
                            get_info, memory_info, gptr, ass_entry
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: copy_g3b_to_g3a

  PUBLIC :: construct_g3b ! construct the g3b table
  PUBLIC :: destruct_g3b  ! destruct  the g3b table

  PUBLIC :: ass_entry
  PUBLIC :: new_entry
  PUBLIC :: get_entry
  PUBLIC :: get_info

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  PUBLIC :: memory_info            ! meta data

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: geosp(:,:)
  REAL(dp), POINTER, PUBLIC :: ts(:,:)
  REAL(dp), POINTER, PUBLIC :: ws(:,:)
  REAL(dp), POINTER, PUBLIC :: wl(:,:)
  REAL(dp), POINTER, PUBLIC :: sn(:,:)
  REAL(dp), POINTER, PUBLIC :: slm(:,:)
  REAL(dp), POINTER, PUBLIC :: az0(:,:)
  REAL(dp), POINTER, PUBLIC :: alb(:,:)
  REAL(dp), POINTER, PUBLIC :: ewov(:,:)
  REAL(dp), POINTER, PUBLIC :: nsov(:,:)
  REAL(dp), POINTER, PUBLIC :: nwov(:,:)
  REAL(dp), POINTER, PUBLIC :: neov(:,:)
  REAL(dp), POINTER, PUBLIC :: varor(:,:)
  REAL(dp), POINTER, PUBLIC :: forest(:,:)
  REAL(dp), POINTER, PUBLIC :: vgrat(:,:)
  REAL(dp), POINTER, PUBLIC :: vlt(:,:)
  REAL(dp), POINTER, PUBLIC :: wsmx(:,:)
  REAL(dp), POINTER, PUBLIC :: fao(:,:)
  REAL(dp), POINTER, PUBLIC :: td(:,:)
  REAL(dp), POINTER, PUBLIC :: tdcl(:,:)
  REAL(dp), POINTER, PUBLIC :: aps(:,:)
  REAL(dp), POINTER, PUBLIC :: aprl(:,:)
  REAL(dp), POINTER, PUBLIC :: aprc(:,:)
  REAL(dp), POINTER, PUBLIC :: aprs(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrgw(:,:)
  REAL(dp), POINTER, PUBLIC :: vdisgw(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcov(:,:)
  REAL(dp), POINTER, PUBLIC :: temp2(:,:)
  REAL(dp), POINTER, PUBLIC :: dew2(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10(:,:)
  REAL(dp), POINTER, PUBLIC :: u10(:,:)
  REAL(dp), POINTER, PUBLIC :: v10(:,:)
  REAL(dp), POINTER, PUBLIC :: srads(:,:)
  REAL(dp), POINTER, PUBLIC :: trads(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0(:,:)
  REAL(dp), POINTER, PUBLIC :: vdis(:,:)
  REAL(dp), POINTER, PUBLIC :: ustr(:,:)
  REAL(dp), POINTER, PUBLIC :: vstr(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfs(:,:)
  REAL(dp), POINTER, PUBLIC :: evap(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfl(:,:)
  REAL(dp), POINTER, PUBLIC :: wlm1(:,:)
  REAL(dp), POINTER, PUBLIC :: tsm1(:,:)
  REAL(dp), POINTER, PUBLIC :: tdm1(:,:)
  REAL(dp), POINTER, PUBLIC :: wsm1(:,:)
  REAL(dp), POINTER, PUBLIC :: wdm1(:,:)
  REAL(dp), POINTER, PUBLIC :: snm1(:,:)
  REAL(dp), POINTER, PUBLIC :: emter(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsol(:,:,:)
  REAL(dp), POINTER, PUBLIC :: runoff(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0u(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsu(:,:)
  REAL(dp), POINTER, PUBLIC :: albedo(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurf(:,:)
  REAL(dp), POINTER, PUBLIC :: tsn(:,:)
  REAL(dp), POINTER, PUBLIC :: td3(:,:)
  REAL(dp), POINTER, PUBLIC :: td4(:,:)
  REAL(dp), POINTER, PUBLIC :: td5(:,:)
  REAL(dp), POINTER, PUBLIC :: tsnm1(:,:)
  REAL(dp), POINTER, PUBLIC :: td3m1(:,:)
  REAL(dp), POINTER, PUBLIC :: td4m1(:,:)
  REAL(dp), POINTER, PUBLIC :: td5m1(:,:)
  REAL(dp), POINTER, PUBLIC :: tdclm1(:,:)
  REAL(dp), POINTER, PUBLIC :: seaice(:,:)
  REAL(dp), POINTER, PUBLIC :: siced(:,:)
  REAL(dp), POINTER, PUBLIC :: ustar3(:,:)
  REAL(dp), POINTER, PUBLIC :: teff(:,:)
  REAL(dp), POINTER, PUBLIC :: glac(:,:)
  REAL(dp), POINTER, PUBLIC :: aclc(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aclcac(:,:,:)
  REAL(dp), POINTER, PUBLIC :: snmel(:,:)
  REAL(dp), POINTER, PUBLIC :: runtoc(:,:)
  REAL(dp), POINTER, PUBLIC :: tslin(:,:)
  REAL(dp), POINTER, PUBLIC :: dsnac(:,:)
  REAL(dp), POINTER, PUBLIC :: t2max(:,:)
  REAL(dp), POINTER, PUBLIC :: t2min(:,:)
  REAL(dp), POINTER, PUBLIC :: tsmax(:,:)
  REAL(dp), POINTER, PUBLIC :: tsmin(:,:)
  REAL(dp), POINTER, PUBLIC :: wimax(:,:)
  REAL(dp), POINTER, PUBLIC :: topmax(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcv(:,:)
  REAL(dp), POINTER, PUBLIC :: qvi(:,:)
  REAL(dp), POINTER, PUBLIC :: alwcvi(:,:)
  REAL(dp), POINTER, PUBLIC :: runlnd(:,:)
  REAL(dp), POINTER, PUBLIC :: rgcgn(:,:)
  REAL(dp), POINTER, PUBLIC :: sodif(:,:)
  REAL(dp), POINTER, PUBLIC :: srafs(:,:)
  REAL(dp), POINTER, PUBLIC :: trafs(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0(:,:)
  REAL(dp), POINTER, PUBLIC :: emtef(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tclfs(:,:)
  REAL(dp), POINTER, PUBLIC :: sclfs(:,:)
  REAL(dp), POINTER, PUBLIC :: tclf0(:,:)
  REAL(dp), POINTER, PUBLIC :: sclf0(:,:)
  REAL(dp), POINTER, PUBLIC :: tke(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem1(:,:,:)
  REAL(dp), POINTER, PUBLIC :: auxil1(:,:)
  REAL(dp), POINTER, PUBLIC :: auxil2(:,:)
  REAL(dp), POINTER, PUBLIC :: drain(:,:)
  REAL(dp), POINTER, PUBLIC :: aprflux(:,:)     

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g3b
  TYPE (gptr), ALLOCATABLE, PUBLIC :: g3(:)

CONTAINS

  SUBROUTINE construct_g3b (lnlon, lnlpx, lnlev, lngl, nlon, nlev, ngl)

    USE mo_control,   ONLY: lmidatm
    USE mo_io_tables, ONLY: ng3xl, ng3xp
    USE mo_memory_g3a

    INTEGER, INTENT (in) :: lnlon, lnlpx, lnlev, lngl ! size on local PE
    INTEGER, INTENT (in) ::  nlon,         nlev,  ngl ! size of global arrays

    INTEGER :: nlp2, nlevp1, jx
    INTEGER :: lnlevp1
    INTEGER :: dim1(2), dim1p(2)
    INTEGER :: dim2(3), dim2p(3)
    INTEGER :: dim3(3), dim3p(3)
    INTEGER :: dim4(3), dim4p(3)
    INTEGER :: dimx(3), dimxp(3)
    CHARACTER (max_dim_name) :: dim1n(2), dim2n(3), dim3n(3), dim4n(3)
    CHARACTER (max_dim_name) :: yname

    ! construct the g3b table
    !
    ! all information specific to this table is set in this subroutine

    nlp2    = nlon  + 2
    nlevp1  = nlev  + 1
    lnlevp1 = lnlev + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    CALL create_list (g3b)

    ! assign pointers

    dim1p = (/ lnlpx, lngl  /)
    dim1  = (/  nlp2,  ngl  /)
    dim1n = (/ "nlp2","ngl "/)

    dim2p = (/ lnlpx,   lnlevp1, lngl    /)
    dim2  = (/  nlp2,    nlevp1,  ngl    /)
    dim2n = (/ "nlp2  ","nlevp1","ngl   "/)

    dim3p = (/ lnlpx, lnlev, lngl  /)
    dim3  = (/  nlp2,  nlev,  ngl  /)
    dim3n = (/ "nlp2","nlev","ngl "/)

    dim4p = (/ lnlpx,   2,   lngl  /)
    dim4  = (/  nlp2,   2,    ngl  /)
    dim4n = (/ "nlp2","n2  ","ngl "/)


    CALL ass_entry (g3b, 'GEOSP',  geosp,   g3a, 'GEOSPM',  dim1p,  dim1)
    CALL new_entry (g3b, 'TS',     ts,      dim1p,  dim1)
    CALL new_entry (g3b, 'WS',     ws,      dim1p,  dim1)
    CALL ass_entry (g3b, 'WL',     wl,      g3a, 'WLM',     dim1p,  dim1)
    CALL new_entry (g3b, 'SN',     sn,      dim1p,  dim1)
    CALL ass_entry (g3b, 'SLM',    slm,     g3a, 'SLMM',    dim1p,  dim1)
    ! only with fujitsu and nec compiler in call to vdiff
    CALL new_entry (g3b, 'AZ0',    az0,     dim1p,  dim1)
    CALL ass_entry (g3b, 'ALB',    alb,     g3a, 'ALBM',    dim1p,  dim1)
    CALL ass_entry (g3b, 'EWOV',   ewov,    g3a, 'EWOVM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'NSOV',   nsov,    g3a, 'NSOVM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'NWOV',   nwov,    g3a, 'NWOVM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'NEOV',   neov,    g3a, 'NEOVM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'VAROR',  varor,   g3a, 'VARORM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'FOREST', forest,  g3a, 'FORESTM', dim1p,  dim1)
    CALL ass_entry (g3b, 'VGRAT',  vgrat,   g3a, 'VGRATM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'VLT',    vlt,     g3a, 'VLTM',    dim1p,  dim1)
    CALL ass_entry (g3b, 'WSMX',   wsmx,    g3a, 'WSMXM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'FAO',    fao,     g3a, 'FAOM',    dim1p,  dim1)
    CALL new_entry (g3b, 'TD',     td,      dim1p,  dim1)
    CALL new_entry (g3b, 'TDCL',   tdcl,    dim1p,  dim1)
    CALL ass_entry (g3b, 'APS',    aps,     g3a, 'APSM',    dim1p,  dim1)
    CALL ass_entry (g3b, 'APRL',   aprl,    g3a, 'APRLM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'APRC',   aprc,    g3a, 'APRCM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'APRS',   aprs,    g3a, 'APRSM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'USTRGW', ustrgw,  g3a, 'USTRGWM', dim1p,  dim1)
    CALL ass_entry (g3b, 'VSTRGW', vstrgw,  g3a, 'VSTRGWM', dim1p,  dim1)
    CALL ass_entry (g3b, 'VDISGW', vdisgw,  g3a, 'VDISGWM', dim1p,  dim1)
    CALL new_entry (g3b, 'ACLCOV', aclcov,  dim1p,  dim1)
    CALL ass_entry (g3b, 'TEMP2',  temp2,   g3a, 'TEMP2M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'DEW2',   dew2,    g3a, 'DEW2M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'WIND10', wind10,  g3a, 'WIND10M', dim1p,  dim1)
    CALL ass_entry (g3b, 'U10',    u10,     g3a, 'U10M',    dim1p,  dim1)
    CALL ass_entry (g3b, 'V10',    v10,     g3a, 'V10M',    dim1p,  dim1)
    CALL ass_entry (g3b, 'SRADS',  srads,   g3a, 'SRADSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TRADS',  trads,   g3a, 'TRADSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SRAD0',  srad0,   g3a, 'SRAD0M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TRAD0',  trad0,   g3a, 'TRAD0M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'VDIS',   vdis,    g3a, 'VDISM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'USTR',   ustr,    g3a, 'USTRM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'VSTR',   vstr,    g3a, 'VSTRM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'AHFS',   ahfs,    g3a, 'AHFSM',   dim1p,  dim1)
    CALL new_entry (g3b, 'EVAP',   evap,    dim1p,  dim1)
    CALL new_entry (g3b, 'AHFL',   ahfl,    dim1p,  dim1)
    CALL ass_entry (g3b, 'WLM1',   wlm1,    g3a, 'WLM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'TSM1',   tsm1,    g3a, 'TSM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'TDM1',   tdm1,    g3a, 'TDM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'WSM1',   wsm1,    g3a, 'WSM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'WDM1',   wdm1,    g3a, 'WDM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'SNM1',   snm1,    g3a, 'SNM1M',   dim1p,  dim1)
    CALL ass_entry (g3b, 'EMTER',  emter,   g3a, 'EMTERM',  dim2p,  dim2)
    CALL ass_entry (g3b, 'TRSOL',  trsol,   g3a, 'TRSOLM',  dim2p,  dim2)
    CALL ass_entry (g3b, 'RUNOFF', runoff,  g3a, 'RUNOFFM', dim1p,  dim1)
    CALL ass_entry (g3b, 'SRAD0U', srad0u,  g3a, 'SRAD0UM', dim1p,  dim1)
    CALL ass_entry (g3b, 'SRADSU', sradsu,  g3a, 'SRADSUM', dim1p,  dim1)
    CALL ass_entry (g3b, 'TRADSU', tradsu,  g3a, 'TRADSUM', dim1p,  dim1)
    CALL ass_entry (g3b, 'ALBEDO', albedo,  g3a, 'ALBEDOM', dim1p,  dim1)
    CALL ass_entry (g3b, 'TSURF',  tsurf,   g3a, 'TSURFM',  dim1p,  dim1)
    CALL new_entry (g3b, 'TSN',    tsn,     dim1p,  dim1)
    CALL new_entry (g3b, 'TD3',    td3,     dim1p,  dim1)
    CALL new_entry (g3b, 'TD4',    td4,     dim1p,  dim1)
    CALL new_entry (g3b, 'TD5',    td5,     dim1p,  dim1)
    CALL ass_entry (g3b, 'TSNM1',  tsnm1,   g3a, 'TSNM1M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TD3M1',  td3m1,   g3a, 'TD3M1M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TD4M1',  td4m1,   g3a, 'TD4M1M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TD5M1',  td5m1,   g3a, 'TD5M1M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TDCLM1', tdclm1,  g3a, 'TDCLM1M', dim1p,  dim1)
    CALL ass_entry (g3b, 'SEAICE', seaice,  g3a, 'SEAICEM', dim1p,  dim1)
    CALL ass_entry (g3b, 'SICED',  siced,   g3a, 'SICEDM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'USTAR3', ustar3,  g3a, 'USTAR3M', dim1p,  dim1)
    CALL ass_entry (g3b, 'TEFF',   teff,    g3a, 'TEFFM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'GLAC',   glac,    g3a, 'GLACM',   dim1p,  dim1)
    CALL ass_entry (g3b, 'ACLC',   aclc,    g3a, 'ACLCM',   dim3p,  dim3)
    CALL ass_entry (g3b, 'ACLCAC', aclcac,  g3a, 'ACLCACM', dim3p,  dim3)
    CALL ass_entry (g3b, 'SNMEL',  snmel,   g3a, 'SNMELM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'RUNTOC', runtoc,  g3a, 'RUNTOCM', dim1p,  dim1)
    CALL ass_entry (g3b, 'TSLIN',  tslin,   g3a, 'TSLINM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'DSNAC',  dsnac,   g3a, 'DSNACM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'T2MAX',  t2max,   g3a, 'T2MAXM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'T2MIN',  t2min,   g3a, 'T2MINM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TSMAX',  tsmax,   g3a, 'TSMAXM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TSMIN',  tsmin,   g3a, 'TSMINM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'WIMAX',  wimax,   g3a, 'WIMAXM',  dim1p,  dim1)
    CALL new_entry (g3b, 'TOPMAX', topmax,  dim1p,  dim1)
    CALL ass_entry (g3b, 'ACLCV',  aclcv,   g3a, 'ACLCVM',  dim1p,  dim1)
    CALL new_entry (g3b, 'QVI',    qvi,     dim1p,  dim1)
    CALL new_entry (g3b, 'ALWCVI', alwcvi,  dim1p,  dim1)
    CALL ass_entry (g3b, 'RUNLND', runlnd,  g3a, 'RUNLNDM', dim1p,  dim1)
    CALL ass_entry (g3b, 'RGCGN',  rgcgn,   g3a, 'RGCGNM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SODIF',  sodif,   g3a, 'SODIFM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SRAFS',  srafs,   g3a, 'SRAFSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TRAFS',  trafs,   g3a, 'TRAFSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SRAF0',  sraf0,   g3a, 'SRAF0M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TRAF0',  traf0,   g3a, 'TRAF0M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'EMTEF',  emtef,   g3a, 'EMTEFM',  dim4p,  dim4)
    CALL ass_entry (g3b, 'TRSOF',  trsof,   g3a, 'TRSOFM',  dim4p,  dim4)
    CALL ass_entry (g3b, 'TCLFS',  tclfs,   g3a, 'TCLFSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SCLFS',  sclfs,   g3a, 'SCLFSM',  dim1p,  dim1)
    CALL ass_entry (g3b, 'TCLF0',  tclf0,   g3a, 'TCLF0M',  dim1p,  dim1)
    CALL ass_entry (g3b, 'SCLF0',  sclf0,   g3a, 'SCLF0M',  dim1p,  dim1)
    CALL new_entry (g3b, 'TKE',    tke,     dim3p,  dim3)
    ! only with fujitsu and nec compiler in call to vdiff
    CALL new_entry (g3b, 'TKEM1',  tkem1,   dim3p,  dim3)
    CALL ass_entry (g3b, 'AUXIL1', auxil1,  g3a, 'AUXIL1M', dim1p,  dim1)
    CALL ass_entry (g3b, 'AUXIL2', auxil2,  g3a, 'AUXIL2M', dim1p,  dim1)
    CALL ass_entry (g3b, 'DRAIN',  drain,   g3a, 'DRAINM',  dim1p,  dim1)

    IF (lmidatm) &
    CALL ass_entry (g3b, 'APRFLUX', aprflux, g3a, 'APRFLUXM', dim1p, dim1)
     
    IF (ng3xp > 0) THEN
       IF ( .NOT. ALLOCATED(g3))  ALLOCATE(g3(ng3xp))
       DO jx = 1, ng3xp
          WRITE(yname, '(a3,i2.2)' ) 'G3X',jx
          dimxp = (/ lnlpx, ng3xl(jx), lngl /)
          dimx  = (/  nlp2, ng3xl(jx),  ngl /)
          CALL new_entry (g3b, yname, g3(jx)%x, dimxp, dimx) 
       ENDDO
    ENDIF

  END SUBROUTINE construct_g3b

  SUBROUTINE destruct_g3b

    CALL delete_list (g3b)

  END SUBROUTINE destruct_g3b

  SUBROUTINE copy_g3b_to_g3a

    USE mo_linked_list,   ONLY: list_element
    USE mo_memory_g3a,    ONLY: g3a

    TYPE (list_element), POINTER :: g3b_list_element
    REAL, POINTER :: ptr4d(:,:,:,:)
    INTEGER       :: i
    CHARACTER (8) :: yname


    g3b_list_element => g3b%first_list_element

    DO WHILE (ASSOCIATED(g3b_list_element))
      IF (.NOT. g3b_list_element%field%info%assign) THEN
        yname = g3b_list_element%field%info%name
        IF (yname(1:3) == 'G3X') THEN
          yname(4:4) = 'M'
          yname(5:6) = g3b_list_element%field%info%name(4:5)
        ELSE
          DO i = 1, 8
            IF (yname(i:i) == ' ') EXIT
          END DO
          yname(i:i) = 'M'
        END IF
!        WRITE(*,*) 'copy_g3b_to_g3a : ', g3b_list_element%field%info%name, yname
        CALL get_entry (g3a, yname, ptr4d)
        ptr4d(:,:,:,:) = g3b_list_element%field%ptr(:,:,:,:)
      END IF
      g3b_list_element => g3b_list_element%next_list_element
    END DO

  END SUBROUTINE copy_g3b_to_g3a

END MODULE mo_memory_g3b
