MODULE mo_test_trans
  !
  !+ $Id: mo_test_trans.f90,v 1.20 1999/11/30 09:22:11 m214030 Exp $
  !
  ! This module holds the routines to compare the decomposed fields
  ! on PE 1.. with the global fields calculated on PE 0 in debug mode
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  !

  USE mo_exception,     ONLY: finish
  USE mo_mpi,           ONLY: p_pe, p_io, p_send, p_recv
  USE mo_decomposition, ONLY: dcg => global_decomposition, &
                              dcl => local_decomposition,  &
                              debug_parallel, any_col_1d
  USE mo_transpose,     ONLY: gather_sp, gather_ls, gather_gp, &
                              gather_sa, tag_gather_gp
  USE mo_doctor,        ONLY: nerr

  IMPLICIT NONE

  PRIVATE

  ! Generic routines to compare fields 

  PUBLIC :: test_spectral   ! in spectral  space (generic)
  PUBLIC :: test_legendre   ! in Legendre  space (generic)
  PUBLIC :: test_gridpoint  ! in gridpoint space (generic)
  PUBLIC :: test_symasym    ! in Legendre  space (generic, [a]sym Fourier comp)
  PUBLIC :: test_zonmean    ! zonal mean in gridpoint space
  PUBLIC :: test_scalar     ! compare scalar value
  PUBLIC :: test_row        ! compare latitude row
  PUBLIC :: test_slrow      ! compare latitude row (SLT conform index)
  PUBLIC :: jsave, isave    ! latitude indices may be stored here

  ! Specific routines

  INTERFACE test_spectral
     MODULE PROCEDURE test_spectral2     ! (nlev, nsp for m=0)
     MODULE PROCEDURE test_spectral3     ! (nlev, nsp)
  END INTERFACE

  INTERFACE test_legendre
     MODULE PROCEDURE test_legendre0     ! (nlev,    nsp for m=0)
     MODULE PROCEDURE test_legendre3     ! (nlev, 2, nsp)
     MODULE PROCEDURE test_legendre4     ! (nmp1, nlevp1, nlat, nvar)
  END INTERFACE

  INTERFACE test_gridpoint
     MODULE PROCEDURE test_gridpoint2    ! (nlon,      nlat)
     MODULE PROCEDURE test_gridpoint3    ! (nlon, nlev,nlat)
     MODULE PROCEDURE test_gridpoint4
  END INTERFACE

  INTERFACE test_symasym
     MODULE PROCEDURE test_symasym2      ! (nlev[+1],          nhgl)
     MODULE PROCEDURE test_symasym4      ! (nlev[+1], 2, nmp1, nhgl)   
  END INTERFACE

  INTERFACE test_row
     MODULE PROCEDURE test_row0          ! scalar
     MODULE PROCEDURE test_row1          ! (nlon[+x])
     MODULE PROCEDURE test_row2          ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE test_row3          ! (nlon[+x], any, any)
  END INTERFACE

  INTERFACE test_slrow
     MODULE PROCEDURE test_slrow0        ! scalar
     MODULE PROCEDURE test_slrow1        ! (nlon[+x])
     MODULE PROCEDURE test_slrow2        ! (nlon[+x], nlev[+1])
     MODULE PROCEDURE test_slrow3        ! (nlon[+x], any, any)
  END INTERFACE

  ! module variables

                   ! latitude indices (2nd argument to test row)
                   ! may be stored here to be passed to subroutines
                   ! without changing the parameter list:
  INTEGER :: jsave ! local latitude                          N->S
  INTEGER :: isave ! local latitude, semilagrangian indexing S->N

CONTAINS

  !===========================================================================
  SUBROUTINE test_spectral3 (sp, name)
    REAL ,INTENT(in)              :: sp (:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(sp,1),SIZE(sp,2),SIZE(sp,3)))
       ENDIF
       CALL gather_sp(tmp,sp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sp)) THEN
             CALL finish('test_spectral3', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_spectral3
  !---------------------------------------------------------------------------
  SUBROUTINE test_spectral2 (sp, name)
    REAL              ,INTENT(in) :: sp (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(sp,1),SIZE(sp,2)))
       ENDIF
       CALL gather_sp(tmp,sp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sp)) THEN
             CALL finish('test_spectral2', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_spectral2
  !===========================================================================
  SUBROUTINE test_legendre4 (ls, name)
    REAL ,INTENT(in)              :: ls (:,:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2),SIZE(ls,3),SIZE(ls,4)))
       ENDIF
       CALL gather_ls (tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=ls)) THEN
             CALL finish('test_legendre4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre4
  !---------------------------------------------------------------------------
  SUBROUTINE test_legendre3 (ls, name)
    REAL ,INTENT(in)              :: ls (:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2),SIZE(ls,3)))
       ENDIF
       CALL gather_ls(tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
         IF (ANY (tmp/=ls)) THEN
           CALL finish('test_legendre3','decomposition test failed for '//name)
         ENDIF
         DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre3
  !---------------------------------------------------------------------------
  SUBROUTINE test_legendre0 (ls, name)
    REAL              ,INTENT(in) :: ls (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(ls,1),SIZE(ls,2)))
       ENDIF
       CALL gather_ls(tmp,ls,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=ls)) THEN
             CALL finish('test_legendre0', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_legendre0
  !===========================================================================
  SUBROUTINE test_symasym4 (sa, name)
    REAL ,INTENT(in)              :: sa (:,:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(sa,1),SIZE(sa,2),SIZE(sa,3),SIZE(sa,4)))
       ENDIF
       CALL gather_sa (tmp,sa,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sa)) THEN
             CALL finish('test_symasym4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_symasym4
  !---------------------------------------------------------------------------
  SUBROUTINE test_symasym2 (sa, name)
    REAL ,INTENT(in)              :: sa (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:)
    IF (any_col_1d) RETURN
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(sa,1),SIZE(sa,2)))
       ENDIF
       CALL gather_sa (tmp,sa,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp/=sa)) THEN
             CALL finish('test_symasym4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_symasym2
  !===========================================================================
  SUBROUTINE test_gridpoint4 (gp, name)
    REAL ,INTENT(in)              :: gp (:,:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:,:)
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2),SIZE(gp,3),SIZE(gp,4)))
         tmp = gp
       ENDIF
       CALL gather_gp (tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:,:,:)/=gp(:dcl%nlon,:,:,:))) THEN
             CALL finish('test_gridpoint4', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_gridpoint4
  !----------------------------------------------------------------------------
  SUBROUTINE test_gridpoint3 (gp, name)
    REAL ,INTENT(in)              :: gp (:,:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:,:)
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2),SIZE(gp,3)))
         tmp = gp
       ENDIF
       CALL gather_gp (tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:,:)/=gp(:dcl%nlon,:,:))) THEN
             CALL finish('test_gridpoint3', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_gridpoint3
  !----------------------------------------------------------------------------
  SUBROUTINE test_gridpoint2 (gp, name)
    REAL ,INTENT(in)              :: gp (:,:)
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL ,POINTER :: tmp (:,:)
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
         ALLOCATE (tmp(SIZE(gp,1),SIZE(gp,2)))
         tmp = gp
       ENDIF
       CALL gather_gp(tmp,gp,dcg,source=1)
       IF (p_pe==p_io) THEN
          IF (ANY (tmp(:dcl%nlon,:)/=gp(:dcl%nlon,:))) THEN
             CALL finish('test_gridpoint2', &
                  'decomposition test failed for '//name)
          ENDIF
          DEALLOCATE (tmp)
       ENDIF
    ENDIF
  END SUBROUTINE test_gridpoint2
  !===========================================================================
  SUBROUTINE test_zonmean (zm, name, abort)
    REAL ,INTENT(in)              :: zm (:,:) ! zonal mean (nlev,nlat)
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL ,POINTER :: tmp (:,:)
    INTEGER       :: i
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             ALLOCATE (tmp( SIZE(zm,1), dcg(i)%nglat))
             CALL p_recv (tmp, dcg(i)%pe, tag_gather_gp)
             IF (ANY(tmp(:,                :dcg(i)%nglh(1))   &
                  /=zm(:,dcg(i)%glats(1) :dcg(i)%glate(1)))) THEN
                  WRITE(nerr,*)'test_zonmean: decomposition test failed for '&
                    //name
                  IF (ab) CALL finish('test_zonmean',&
                    'decomposition test failed for '//name)
             ENDIF
             IF (ANY(tmp(:,dcg(i)%nglh(1)+1:)                 &
                  /=zm(:,dcg(i)%glats(2) :dcg(i)%glate(2)))) THEN
                  WRITE(nerr,*)'test_zonmean: decomposition test failed for '&
                    //name
                  IF (ab) CALL finish('test_zonmean',&
                    'decomposition test failed for '//name)
             ENDIF
             DEALLOCATE (tmp)
          END DO
       ELSE
          CALL p_send (zm, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_zonmean
  !===========================================================================
  SUBROUTINE test_scalar (sc, name)
    REAL ,INTENT(in)              :: sc   ! scalar
    CHARACTER (len=*) ,INTENT(in) :: name
    REAL    :: tmp
    INTEGER :: i
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             CALL p_recv (tmp, dcg(i)%pe, tag_gather_gp)
             IF (tmp/=sc) THEN
                WRITE (nerr,*)'test_scalar: decomposition test failed'
                WRITE (nerr,*)'test_scalar: PE,value=',p_io,sc
                WRITE (nerr,*)'test_scalar: PE,value=',dcg(i)%pe,tmp
                CALL finish('test_scalar','decomposition test failed for '//name)
             ENDIF
          END DO
       ELSE
          CALL p_send (sc, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_scalar
  !===========================================================================
  SUBROUTINE test_row1 (rw, j, name, abort)
    REAL              ,INTENT(in) :: rw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL ,POINTER :: tmp (:), buf(:)
    INTEGER       :: i, k, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE(tmp(dcl%nlon))
          tmp = 0.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon)) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   tmp(dcg(i)%glons(1):dcg(i)%glone(1)) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      tmp(dcg(i)%glons(2):dcg(i)%glone(2)) = buf
                   ELSE
                      tmp(:dcg(i)%glone(2) )=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:)
                      tmp( dcg(i)%glons(2):)=buf(:dcg(i)%nglon-dcg(i)%glone(2))
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(tmp/=rw(:dcl%nglon))) THEN
             DO k=1,dcl%nglon
                IF(tmp(k)/=rw(k)) WRITE(nerr,*) k, rw(k), tmp(k), rw(k)-tmp(k)
             END DO
             IF (ab) &
                  CALL finish('test_row1','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(tmp)
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon), p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_row1
  !---------------------------------------------------------------------------
  SUBROUTINE test_row2 (rw, j, name, abort)
    REAL              ,INTENT(in) :: rw (:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL ,POINTER :: tmp (:,:), buf(:,:)
    INTEGER       :: i, jg, k,l
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (tmp (dcl%nlon, SIZE(rw,2)))
          tmp = 0.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(rw,2))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   tmp(dcg(i)%glons(1):dcg(i)%glone(1),:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      tmp(dcg(i)%glons(2):dcg(i)%glone(2),:) = buf
                   ELSE
                      tmp(:dcg(i)%glone(2) ,:)=buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:)
                      tmp( dcg(i)%glons(2):,:)=buf(:dcg(i)%nglon-dcg(i)%glone(2),:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(tmp/=rw(:dcl%nglon,:))) THEN
             WRITE(nerr,*)'test_row2: FAIL: jg=',jg,name
             DO l=1,SIZE(tmp,2)
                DO k=1,SIZE(tmp,1)
                   IF(tmp(k,l)/=rw(k,l)) WRITE(nerr,*) &
                        'k,l, value(k,l) (1PE vs nPEs)',k,l,rw(k,l),tmp(k,l)
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row2','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(tmp)
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon,:), p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_row2
  !---------------------------------------------------------------------------
  SUBROUTINE test_row3 (rw, j, name, abort)
    REAL              ,INTENT(in) :: rw (:,:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: j        ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL ,POINTER :: tmp (:,:,:), buf(:,:,:)
    INTEGER       :: i, jg, k, l, m
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          jg = dcl%glat(j)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          ALLOCATE (tmp (dcl%nlon, SIZE(rw,2), SIZE(rw,3)))
          tmp = 0.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive rows only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                ALLOCATE(buf(dcg(i)%nglon, SIZE(rw,2), SIZE(rw,3))) 
                CALL p_recv (buf, dcg(i)%pe, tag_gather_gp)
                !------------------------------
                ! account for logitudinal shift
                !------------------------------
                IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1))) THEN
                   tmp(dcg(i)%glons(1):dcg(i)%glone(1),:,:) = buf
                ELSE
                   IF(dcg(i)%glats(2) <= dcg(i)%glate(2)) THEN
                      tmp(dcg(i)%glons(2):dcg(i)%glone(2),:,:) = buf
                   ELSE
                      tmp(:dcg(i)%glone(2) ,:,:) = &
                           buf( dcg(i)%nglon-dcg(i)%glone(2)+1:,:,:)
                      tmp( dcg(i)%glons(2):,:,:) = &
                           buf(:dcg(i)%nglon-dcg(i)%glone(2),:,:)
                   ENDIF
                ENDIF
                DEALLOCATE(buf)
             ENDIF
          END DO
          !--------
          ! compare
          !--------
          IF (ANY(tmp/=rw(:dcl%nglon,:,:))) THEN
             WRITE(nerr,*)'test_row3: FAIL: jg=',jg,name
             DO m=1,SIZE(tmp,3)
                DO l=1,SIZE(tmp,2)
                   DO k=1,SIZE(tmp,1)
                      IF(tmp(k,l,m)/=rw(k,l,m)) WRITE(nerr,*) &
                           'k,l,m, value(k,l,m) (1PE vs nPEs)',k,l,m,rw(k,l,m),tmp(k,l,m)
                   END DO
                END DO
             END DO
             IF (ab) &
                  CALL finish('test_row3','decomposition test failed for '//name)
          ENDIF
          DEALLOCATE(tmp)
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw(:dcl%nglon,:,:), p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_row3
  !---------------------------------------------------------------------------
  SUBROUTINE test_row0 (rw, jl, name, abort)
    REAL              ,INTENT(in) :: rw    ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    REAL          :: tmp
    INTEGER       :: i, jg
    LOGICAL       :: ab
    ab = .TRUE.; IF (PRESENT(abort)) ab = abort
    IF (debug_parallel>=0) THEN
       IF (p_pe==p_io) THEN
          jg = dcl%glat(jl)
          !---------------------------
          ! if pe == p_io receive rows
          !---------------------------
          tmp = 0.
          DO i=1,SIZE(dcg)
             IF(dcg(i)%pe==p_pe) CYCLE
             !-------------------------------------------------------
             ! receive row0 only from PE's with same global row index
             !-------------------------------------------------------
             IF((dcg(i)%glats(1) <= jg .AND. jg <= dcg(i)%glate(1)) .OR. &
                  (dcg(i)%glats(2) <= jg .AND. jg <= dcg(i)%glate(2))) THEN
                CALL p_recv (tmp, dcg(i)%pe, tag_gather_gp)
                !--------
                ! compare
                !--------
                IF (tmp/=rw) THEN
                   WRITE(nerr,*)'test_row0: FAIL: pe, jg=',dcg(i)%pe,jg,name
                   WRITE(nerr,*)'test_row0: pe, value',p_pe,rw
                   WRITE(nerr,*)'test_row0: pe, value',dcg(i)%pe,tmp
                   IF (ab) &
                        CALL finish('test_row0','decomposition test failed for '//name)
                ENDIF
             ENDIF
          END DO
       ELSE
          !-----------------------
          ! if pe /= p_io send row
          !-----------------------
          CALL p_send (rw, p_io, tag_gather_gp)
       ENDIF
    ENDIF
  END SUBROUTINE test_row0
  !============================================================================
  SUBROUTINE test_slrow3 (rw, jl, name, abort)
    REAL              ,INTENT(in) :: rw (:,:,:) ! slrow (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: jl         ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row3 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow3
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow2 (rw, jl, name, abort)
    REAL              ,INTENT(in) :: rw (:,:) ! row (nlon[+x], nlev[+1])
    INTEGER           ,INTENT(in) :: jl       ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row2 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow2
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow1 (rw, jl, name, abort)
    REAL              ,INTENT(in) :: rw   (:) ! row (nlon[+x])
    INTEGER           ,INTENT(in) :: jl       ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row1 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow1
  !---------------------------------------------------------------------------
  SUBROUTINE test_slrow0 (rw, jl, name, abort)
    REAL              ,INTENT(in) :: rw     ! scalar to compare
    INTEGER           ,INTENT(in) :: jl     ! slt row index
    CHARACTER (len=*) ,INTENT(in) :: name
    LOGICAL ,OPTIONAL ,INTENT(in) :: abort
    CALL test_row0 (rw, dcl%nglat-jl+1, name, abort)
  END SUBROUTINE test_slrow0
  !==============================================================================
END MODULE mo_test_trans
