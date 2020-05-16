MODULE m_alloc_mods
  ! ------------------------------------------------------------------
  !
  ! module *mo_control* - control variables for model housekeeping.
  !
  ! u. schlese   dkrz-hamburg    dec-94
  !
  ! A. Rhodin    MPI-Hamburg     Jan-99:
  !   Subroutine m_control renamed to alloc_mods and moved from 
  !   module mo_control to module m_alloc_mods.
  ! ------------------------------------------------------------------

IMPLICIT NONE

CONTAINS

  SUBROUTINE alloc_mods

    USE mo_control,     ONLY: nlev, nlevp1
    USE mo_post,        ONLY: npplev
    USE mo_physc2,      ONLY: cevapcu
    USE mo_hdiff,       ONLY: diftcor
    USE mo_hyb,         ONLY: aktlrd, alpham, altrcp, ardprc, bb, ceta,      &
                              cetah, cpg, dela, delb, delpr, ralpha, ralphr, &
                              rddelb, rdelpr, rlnmar, rlnpr
    USE mo_semi_impl,   ONLY: vmax
    USE mo_truncation,  ONLY: ntrk, ntrm
    USE mo_stat_global, ONLY: dh, qh, th, voh, xh
    USE mo_stat_zonal,  ONLY: delph
    USE mo_rad2,        ONLY: cvdaed, cvdael, cvdaes, cvdaeu
    USE mo_diff,        ONLY: iq, ncdif

    IMPLICIT NONE

    LOGICAL, SAVE :: lnot_used = .TRUE.

    IF (lnot_used) THEN

       ! mo_physc2
       ALLOCATE (cevapcu(nlev))
       ! mo_post
       ALLOCATE (npplev(nlevp1,256))
       ! mo_hdiff
       ALLOCATE (diftcor(nlev))
       ! mo_hyp
       ALLOCATE (ralpha(nlev))
       ALLOCATE (rlnpr(nlev))
       ALLOCATE (dela(nlev))
       ALLOCATE (delb(nlev))
       ALLOCATE (rddelb(nlev))
       ALLOCATE (cpg(nlev))
       ALLOCATE (delpr(nlev))
       ALLOCATE (rdelpr(nlev))
       ALLOCATE (ralphr(nlev))
       ALLOCATE (alpham(nlev))
       ALLOCATE (ardprc(nlev))
       ALLOCATE (rlnmar(nlev))
       ALLOCATE (aktlrd(nlev))
       ALLOCATE (altrcp(nlev))
       ALLOCATE (ceta(nlev))
       ALLOCATE (cetah(nlevp1))
       ALLOCATE (bb(nlev,nlev))
       ! mo_stat_global
       ALLOCATE (dh(nlev))
       ALLOCATE (voh(nlev))
       ALLOCATE (qh(nlev))
       ALLOCATE (th(nlev))
       ALLOCATE (xh(nlev))
       ! mo_semi_impl
       ALLOCATE (vmax(nlev))
       ! mo_truncation
       ALLOCATE (ntrm(nlev))
!      ALLOCATE (ntrn(jpnlev)) ! because f90 namelist restriction
       ALLOCATE (ntrk(nlev))
       ! mo_stat_zonal
       ALLOCATE (delph(nlev))
       ! mo_rad2
       ALLOCATE (cvdaes(nlevp1))
       ALLOCATE (cvdael(nlevp1))
       ALLOCATE (cvdaeu(nlevp1))
       ALLOCATE (cvdaed(nlevp1))
       ! mo_diff
       ALLOCATE (ncdif(nlev))
       ALLOCATE (iq(nlev))

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE alloc_mods

END MODULE m_alloc_mods
