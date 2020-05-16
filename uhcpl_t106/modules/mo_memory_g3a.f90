MODULE mo_memory_g3a

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: list
  USE mo_memory_base, ONLY: create_list, delete_list, new_entry, get_entry,&
                            print_memory_table, print_memory_use, print_sinfo,&
                            get_info, memory_info, gptr
  USE mo_netCDF,      ONLY: max_dim_name

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_g3a ! construct the g3a table
  PUBLIC :: destruct_g3a  ! destruct  the g3a table

  PUBLIC :: new_entry
  PUBLIC :: get_entry
  PUBLIC :: get_info

  PUBLIC :: print_memory_table     ! print information on sp table
  PUBLIC :: print_memory_use
  PUBLIC :: print_sinfo

  PUBLIC :: memory_info            ! meta data

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: geospm(:,:)
  REAL(dp), POINTER, PUBLIC :: tsm(:,:)
  REAL(dp), POINTER, PUBLIC :: wsm(:,:)
  REAL(dp), POINTER, PUBLIC :: wlm(:,:)
  REAL(dp), POINTER, PUBLIC :: snm(:,:)
  REAL(dp), POINTER, PUBLIC :: slmm(:,:)
  REAL(dp), POINTER, PUBLIC :: az0m(:,:)
  REAL(dp), POINTER, PUBLIC :: albm(:,:)
  REAL(dp), POINTER, PUBLIC :: ewovm(:,:)
  REAL(dp), POINTER, PUBLIC :: nsovm(:,:)
  REAL(dp), POINTER, PUBLIC :: nwovm(:,:)
  REAL(dp), POINTER, PUBLIC :: neovm(:,:)
  REAL(dp), POINTER, PUBLIC :: varorm(:,:)
  REAL(dp), POINTER, PUBLIC :: forestm(:,:)
  REAL(dp), POINTER, PUBLIC :: vgratm(:,:)
  REAL(dp), POINTER, PUBLIC :: vltm(:,:)
  REAL(dp), POINTER, PUBLIC :: wsmxm(:,:)
  REAL(dp), POINTER, PUBLIC :: faom(:,:)
  REAL(dp), POINTER, PUBLIC :: tdm(:,:)
  REAL(dp), POINTER, PUBLIC :: tdclm(:,:)
  REAL(dp), POINTER, PUBLIC :: apsm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprlm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprcm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprsm(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: vdisgwm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcovm(:,:)
  REAL(dp), POINTER, PUBLIC :: temp2m(:,:)
  REAL(dp), POINTER, PUBLIC :: dew2m(:,:)
  REAL(dp), POINTER, PUBLIC :: wind10m(:,:)
  REAL(dp), POINTER, PUBLIC :: u10m(:,:)
  REAL(dp), POINTER, PUBLIC :: v10m(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsm(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsm(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0m(:,:)
  REAL(dp), POINTER, PUBLIC :: trad0m(:,:)
  REAL(dp), POINTER, PUBLIC :: vdism(:,:)
  REAL(dp), POINTER, PUBLIC :: ustrm(:,:)
  REAL(dp), POINTER, PUBLIC :: vstrm(:,:)
  REAL(dp), POINTER, PUBLIC :: ahfsm(:,:)
  REAL(dp), POINTER, PUBLIC :: evapm(:,:)
  REAL(dp), POINTER, PUBLIC :: ahflm(:,:)
  REAL(dp), POINTER, PUBLIC :: wlm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: tsm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: tdm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: wsm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: wdm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: snm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: emterm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsolm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: runoffm(:,:)
  REAL(dp), POINTER, PUBLIC :: srad0um(:,:)
  REAL(dp), POINTER, PUBLIC :: sradsum(:,:)
  REAL(dp), POINTER, PUBLIC :: tradsum(:,:)
  REAL(dp), POINTER, PUBLIC :: albedom(:,:)
  REAL(dp), POINTER, PUBLIC :: tsurfm(:,:)
  REAL(dp), POINTER, PUBLIC :: tsnm(:,:)
  REAL(dp), POINTER, PUBLIC :: td3m(:,:)
  REAL(dp), POINTER, PUBLIC :: td4m(:,:)
  REAL(dp), POINTER, PUBLIC :: td5m(:,:)
  REAL(dp), POINTER, PUBLIC :: tsnm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: td3m1m(:,:)
  REAL(dp), POINTER, PUBLIC :: td4m1m(:,:)
  REAL(dp), POINTER, PUBLIC :: td5m1m(:,:)
  REAL(dp), POINTER, PUBLIC :: tdclm1m(:,:)
  REAL(dp), POINTER, PUBLIC :: seaicem(:,:)
  REAL(dp), POINTER, PUBLIC :: sicedm(:,:)
  REAL(dp), POINTER, PUBLIC :: ustar3m(:,:)
  REAL(dp), POINTER, PUBLIC :: teffm(:,:)
  REAL(dp), POINTER, PUBLIC :: glacm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: aclcacm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: snmelm(:,:)
  REAL(dp), POINTER, PUBLIC :: runtocm(:,:)
  REAL(dp), POINTER, PUBLIC :: tslinm(:,:)
  REAL(dp), POINTER, PUBLIC :: dsnacm(:,:)
  REAL(dp), POINTER, PUBLIC :: t2maxm(:,:)
  REAL(dp), POINTER, PUBLIC :: t2minm(:,:)
  REAL(dp), POINTER, PUBLIC :: tsmaxm(:,:)
  REAL(dp), POINTER, PUBLIC :: tsminm(:,:)
  REAL(dp), POINTER, PUBLIC :: wimaxm(:,:)
  REAL(dp), POINTER, PUBLIC :: topmaxm(:,:)
  REAL(dp), POINTER, PUBLIC :: aclcvm(:,:)
  REAL(dp), POINTER, PUBLIC :: qvim(:,:)
  REAL(dp), POINTER, PUBLIC :: alwcvim(:,:)
  REAL(dp), POINTER, PUBLIC :: runlndm(:,:)
  REAL(dp), POINTER, PUBLIC :: rgcgnm(:,:)
  REAL(dp), POINTER, PUBLIC :: sodifm(:,:)
  REAL(dp), POINTER, PUBLIC :: srafsm(:,:)
  REAL(dp), POINTER, PUBLIC :: trafsm(:,:)
  REAL(dp), POINTER, PUBLIC :: sraf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: traf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: emtefm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: trsofm(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tclfsm(:,:)
  REAL(dp), POINTER, PUBLIC :: sclfsm(:,:)
  REAL(dp), POINTER, PUBLIC :: tclf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: sclf0m(:,:)
  REAL(dp), POINTER, PUBLIC :: tkem(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tkem1m(:,:,:)
  REAL(dp), POINTER, PUBLIC :: auxil1m(:,:)
  REAL(dp), POINTER, PUBLIC :: auxil2m(:,:)
  REAL(dp), POINTER, PUBLIC :: drainm(:,:)
  REAL(dp), POINTER, PUBLIC :: aprfluxm(:,:)     ! for middle atmosphere only

  ! declaration of table with 3d-field entries

  TYPE (list), PUBLIC :: g3a
  TYPE (gptr), ALLOCATABLE, PUBLIC :: g3m(:)

CONTAINS

  SUBROUTINE construct_g3a (lnlon, lnlpx, lnlev, lngl, nlon, nlev, ngl)

    USE mo_control,   ONLY: lmidatm
    USE mo_io_tables, ONLY: ng3xl, ng3xp

    INTEGER, INTENT (in) :: lnlon, lnlpx, lnlev, lngl ! size on local PE
    INTEGER, INTENT (in) ::  nlon,         nlev,  ngl ! size of global arrays

    CHARACTER (max_dim_name) ::  dimxn(3)

    INTEGER :: nlp2, nlevp1, jx
    INTEGER :: lnlevp1
    INTEGER :: dim1(2), dim1p(2)
    INTEGER :: dim2(3), dim2p(3)
    INTEGER :: dim3(3), dim3p(3)
    INTEGER :: dim4(3), dim4p(3)
    INTEGER :: dimx(3), dimxp(3)
    CHARACTER (max_dim_name) :: dim1n(2), dim2n(3), dim3n(3), dim4n(3)
    CHARACTER (max_dim_name) :: yname

    ! construct the g3a table
    !
    ! all information specific to this table is set in this subroutine

    nlp2    = nlon  + 2
    nlevp1  = nlev  + 1
    lnlevp1 = lnlev + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields


    CALL create_list (g3a)

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


    CALL new_entry (g3a, 'GEOSPM',  geospm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSM',     tsm,     dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WSM',     wsm,     dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WLM',     wlm,     dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SNM',     snm,     dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SLMM',    slmm,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'AZ0M',    az0m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ALBM',    albm,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'EWOVM',   ewovm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'NSOVM',   nsovm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'NWOVM',   nwovm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'NEOVM',   neovm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VARORM',  varorm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'FORESTM', forestm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VGRATM',  vgratm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VLTM',    vltm,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WSMXM',   wsmxm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'FAOM',    faom,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TDM',     tdm,     dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TDCLM',   tdclm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'APSM',    apsm,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'APRLM',   aprlm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'APRCM',   aprcm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'APRSM',   aprsm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'USTRGWM', ustrgwm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VSTRGWM', vstrgwm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VDISGWM', vdisgwm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ACLCOVM', aclcovm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TEMP2M',  temp2m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'DEW2M',   dew2m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WIND10M', wind10m, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'U10M',    u10m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'V10M',    v10m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRADSM',  sradsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TRADSM',  tradsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRAD0M',  srad0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TRAD0M',  trad0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VDISM',   vdism,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'USTRM',   ustrm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'VSTRM',   vstrm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'AHFSM',   ahfsm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'EVAPM',   evapm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'AHFLM',   ahflm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WLM1M',   wlm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSM1M',   tsm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TDM1M',   tdm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WSM1M',   wsm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WDM1M',   wdm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SNM1M',   snm1m,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'EMTERM',  emterm,  dim2p,  dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (g3a, 'TRSOLM',  trsolm,  dim2p,  dim2, dimnames=dim2n, restart=.true.)
    CALL new_entry (g3a, 'RUNOFFM', runoffm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRAD0UM', srad0um, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRADSUM', sradsum, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TRADSUM', tradsum, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ALBEDOM', albedom, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSURFM',  tsurfm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSNM',    tsnm,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD3M',    td3m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD4M',    td4m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD5M',    td5m,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSNM1M',  tsnm1m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD3M1M',  td3m1m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD4M1M',  td4m1m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TD5M1M',  td5m1m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TDCLM1M', tdclm1m, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SEAICEM', seaicem, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SICEDM',  sicedm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'USTAR3M', ustar3m, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TEFFM',   teffm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'GLACM',   glacm,   dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ACLCM',   aclcm,   dim3p,  dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (g3a, 'ACLCACM', aclcacm, dim3p,  dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (g3a, 'SNMELM',  snmelm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'RUNTOCM', runtocm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSLINM',  tslinm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'DSNACM',  dsnacm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'T2MAXM',  t2maxm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'T2MINM',  t2minm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSMAXM',  tsmaxm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TSMINM',  tsminm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'WIMAXM',  wimaxm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TOPMAXM', topmaxm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ACLCVM',  aclcvm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'QVIM',    qvim,    dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'ALWCVIM', alwcvim, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'RUNLNDM', runlndm, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'RGCGNM',  rgcgnm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SODIFM',  sodifm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRAFSM',  srafsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TRAFSM',  trafsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SRAF0M',  sraf0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TRAF0M',  traf0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'EMTEFM',  emtefm,  dim4p,  dim4, dimnames=dim4n, restart=.true.)
    CALL new_entry (g3a, 'TRSOFM',  trsofm,  dim4p,  dim4, dimnames=dim4n, restart=.true.)
    CALL new_entry (g3a, 'TCLFSM',  tclfsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SCLFSM',  sclfsm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TCLF0M',  tclf0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'SCLF0M',  sclf0m,  dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'TKEM',    tkem,    dim3p,  dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (g3a, 'TKEM1M',  tkem1m,  dim3p,  dim3, dimnames=dim3n, restart=.true.)
    CALL new_entry (g3a, 'AUXIL1M', auxil1m, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'AUXIL2M', auxil2m, dim1p,  dim1, dimnames=dim1n, restart=.true.)
    CALL new_entry (g3a, 'DRAINM',  drainm,  dim1p,  dim1, dimnames=dim1n, restart=.true.)

    IF (lmidatm) &
    CALL new_entry (g3a, 'APRFLUXM', aprfluxm, dim1p,  dim1, dimnames=dim1n, restart=.true.)

    IF (ng3xp > 0) THEN
       IF ( .NOT. ALLOCATED(g3m))  ALLOCATE(g3m(ng3xp))
       DO jx = 1, ng3xp
          WRITE(yname,'(a5,i2.2)' ) 'ng3xl',jx
          dimxp = (/ lnlpx, ng3xl(jx), lngl /)
          dimx  = (/  nlp2, ng3xl(jx),  ngl /)
          dimxn(1) = "nlp2"
          dimxn(2) = yname
          dimxn(3) = "ngl"
          WRITE(yname,'(a9)' ) '         '
          WRITE(yname,'(a4,i2.2)' ) 'G3XM',jx
          CALL new_entry (g3a, yname, g3m(jx)%x, dimxp, dimx, dimnames=dimxn, restart=.true.) 
       ENDDO
    ENDIF

  END SUBROUTINE construct_g3a

  SUBROUTINE destruct_g3a

    CALL delete_list (g3a)

  END SUBROUTINE destruct_g3a

END MODULE mo_memory_g3a
