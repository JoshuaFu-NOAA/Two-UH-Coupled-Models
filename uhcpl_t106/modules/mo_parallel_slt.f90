MODULE mo_parallel_slt

  ! Original source: T. Diehl, DKRZ, July 1999

  ! Syntax for variables: [n][pe]{snd,rcv}_{eq,po,ea,we}
  ! n indicates "number of ..."
  ! pe indicates processor
  ! snd = send, rcv = receive, eq = equator, po = pole, ea = east, we = west

  USE mo_decomposition, ONLY:gc => global_decomposition
  USE mo_decomposition, ONLY:dc => local_decomposition
  USE mo_transpose,     ONLY:indx
  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: nxpt_a(:,:,:)  ! nxpt_a(plev,platd,2); array 
  ! with current nxpt values

  ! send/receive "towards" equator
  INTEGER, ALLOCATABLE :: pe_snd_eq(:)   ! send lats to PEs with these id's
  INTEGER, ALLOCATABLE :: nsnd_eq(:)     ! # of lats to send to PEs
  INTEGER, ALLOCATABLE :: snd_eq(:,:,:)  ! my local lats to send
  INTEGER              :: npe_snd_eq     ! # of PEs to send to

  INTEGER, ALLOCATABLE :: nrcv_eq(:)     ! # of lats to receive from PEs
  INTEGER, ALLOCATABLE :: rcv_eq(:,:,:)  ! my local lats to receive
  INTEGER              :: npe_rcv_eq     ! # of PEs to receive from

  ! send/receive "towards" poles
  INTEGER, ALLOCATABLE :: pe_snd_po(:)
  INTEGER, ALLOCATABLE :: nsnd_po(:)
  INTEGER, ALLOCATABLE :: snd_po(:,:,:)
  INTEGER              :: npe_snd_po

  INTEGER, ALLOCATABLE :: nh(:,:,:)

  INTEGER, ALLOCATABLE :: nrcv_po(:)
  INTEGER, ALLOCATABLE :: rcv_po(:,:,:)
  INTEGER              :: npe_rcv_po
  INTEGER              :: nrcv_po_own
  INTEGER, ALLOCATABLE :: rcv_po_own(:,:)


  ! send to the east/receive from the west
  INTEGER, ALLOCATABLE :: blks_ea(:)

  ! send to the west/receive from the east
  INTEGER, ALLOCATABLE :: blks_we(:)

CONTAINS

  SUBROUTINE setup_overlap

    ! Set up communication scheme for the slt
    ! Current version does not work if second area of PEs is rotated !!
    ! Current version assumes static overlap in latitudinal direction !!
    ! This routine should be called after suslt in initialise

    USE mo_grid, ONLY: nxpt, jintmx, jstart, plato2

    ! local arrays
    INTEGER, ALLOCATABLE :: pe_a(:,:)
    INTEGER, ALLOCATABLE :: icountlats(:), ilat(:,:,:), znh(:,:,:)
    INTEGER, ALLOCATABLE :: ilat_own(:,:)
    LOGICAL, ALLOCATABLE :: lpe_tmp(:)

    INTEGER, ALLOCATABLE :: icount_jump(:,:)
    LOGICAL, ALLOCATABLE :: ljump(:)

    ! local scalars
    INTEGER :: row_pe, row_pe_tmp, lower_pe, sndpe_id, isndpe_id, jj
    INTEGER :: counter_lats, count2, sndpelats, mypelats, receiver, sender
    INTEGER :: nprocx, nprocy, pid_lt, pid_ln
    INTEGER :: row_end, row_stride, upper_pe, ih_rcv, ih_snd
    INTEGER :: blks_tmp, col_pe, src_ln, src, isrc
    INTEGER :: icountlats_own
    LOGICAL :: lfound

    ! Be careful about the handling of pid_lt=1: 
    ! - towards equator: only sends
    ! - Pole extensions of pole PEs filled in extyv/extys and not in overlap. 
    ! - Pole PEs are the only ones which might send latitudes from the southern
    !   overlap (on the southern hemisphere) and from the northern overlap (on
    !   the northern hemisphere).
    !
    ! If owner turns out to be the same as receiver: do not send/receive, just
    ! copy data across equator from northern area to southern area or vice versa 
    ! (in routine overlap).

    nprocy = dc%nproca
    nprocx = dc%nprocb
    pid_lt = dc%set_a
    pid_ln = dc%set_b

    ALLOCATE(pe_snd_eq(nprocy))
    ALLOCATE(nsnd_eq(nprocy))
    ALLOCATE(snd_eq(nxpt+jintmx,2,nprocy))

    ALLOCATE(nrcv_eq(nprocy))
    ALLOCATE(rcv_eq(nxpt+jintmx,2,nprocy))

    ALLOCATE(pe_snd_po(nprocy))
    ALLOCATE(nsnd_po(nprocy))
    ALLOCATE(snd_po(nxpt+jintmx,2,nprocy))

    ALLOCATE(nh(nxpt+jintmx,2,nprocy))

    ALLOCATE(nrcv_po(nprocy))
    ALLOCATE(rcv_po(nxpt+jintmx,2,nprocy))
    ALLOCATE(rcv_po_own(nxpt+jintmx,2))

    ALLOCATE(blks_ea(nprocx))
    ALLOCATE(blks_we(nprocx))

    ! local arrays
    ALLOCATE(pe_a(nxpt+jintmx,nprocy))
    ALLOCATE(icountlats(nprocy))
    ALLOCATE(ilat(nxpt+jintmx,2,nprocy))
    ALLOCATE(ilat_own(nxpt+jintmx,2))
    ALLOCATE(znh(nxpt+jintmx,2,nprocy))
    ALLOCATE(lpe_tmp(nprocy))

    ALLOCATE(icount_jump(2,nprocy))
    ALLOCATE(ljump(2))

    IF (pid_lt == 1) THEN
       mypelats = dc%nglat/2 + nxpt +jintmx
    ELSE
       mypelats = dc%nglat/2
    END IF

    !===========================================================================

    ! Set up scheme for receives of lats that are sent towards the equator,
    ! i.e. receives for southern overlap in southern hemisphere and 
    ! northern overlap in northern hemisphere (except for pid_lt=1).
    ! Determine owner of needed overlap lats and local lat numbers of owner.

    ! determine pe_a globally

    pe_a(:,:) = -999

    DO row_pe = 2,nprocy      ! determine receive info for these PEs
       row_pe_tmp = row_pe
       counter_lats = 0
       DO jj = 1,nxpt+jintmx  ! loop over all lats of an overlap region (for
          ! one hemisphere)

          ! Now loop over all PEs south of current PE (north of current PE for 
          ! northern hemisphere) until the PE is found which owns the current
          ! overlap latitude.
          DO lower_pe = row_pe_tmp - 1, 1, -1
             sndpe_id = dc%mapmesh(1,lower_pe)  ! choose any column, e.g. 1
             !             isndpe_id = indx(sndpe_id,gc)
             IF (lower_pe == 1 ) THEN
                sndpelats = gc(sndpe_id)%nglat/2 + nxpt +jintmx
             ELSE
                sndpelats = gc(sndpe_id)%nglat/2
             END IF
             ! if counter of lats for lower_pe <= total #lats for lower_pe,
             ! i.e. if needed lat is owned by lower_pe, then ...
             IF (counter_lats + 1 <= sndpelats) THEN

                ! row_pe rcvs from lower_pe
                receiver = row_pe
                sender = lower_pe
                pe_a(jj,receiver) = sender
                ! lat was found on this lower_pe; increment counter_lats and
                ! start looking for next lat on this same lower_pe
                counter_lats = counter_lats + 1
                EXIT
             END IF
             ! if lat was not found on this lower_pe, set counter_lats=0 and
             ! loop to next lower PE
             counter_lats = 0
          END DO
          ! for next lat, start searching at last lower_pe (and not at row_pe)
          row_pe_tmp    = lower_pe + 1
       END DO
    END DO

    ! local send values
    pe_snd_eq(:) = -999
    npe_snd_eq = 0

    DO row_pe = 1,nprocy
       DO jj = 1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN  ! If I am a sender...
             ! number of PEs to send to (inludes myself in cases of copying)
             npe_snd_eq = npe_snd_eq + 1
             ! dest id
             pe_snd_eq(npe_snd_eq) = dc%mapmesh(pid_ln,row_pe)
             EXIT
          END IF
       END DO
    END DO

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nsnd_eq(:) = 0
    snd_eq(:,:,:) = -999

    DO row_pe =1,nprocy
       DO jj=1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN
             ! icountlats counts the lats row_pe needs from me (per hemisph.)
             icountlats(row_pe) = icountlats(row_pe) + 1
             ! local lats (S)
             ilat(icountlats(row_pe),1,row_pe) = jstart + plato2 - &
                  &                                                icountlats(row_pe)
             ! local lats (N)        
             ilat(icountlats(row_pe),2,row_pe) = jstart + &
                  &                                                icountlats(row_pe) - 1 
          END IF
       END DO
    END DO

    nsnd_eq(:npe_snd_eq) = PACK(icountlats, icountlats /= 0)

    count2=0
    DO row_pe=1,nprocy
       IF (ilat(1,1,row_pe) /= -999) THEN
          count2 = count2 + 1
          snd_eq(:,:,count2) = ilat(:,:,row_pe)
       ENDIF
    ENDDO

    ! local receive values
    lpe_tmp(:) = .FALSE.
    npe_rcv_eq = 0

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       ! If I am a real receiver ...
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! number of PEs to receive from (does NOT include myself in cases
          ! of copying; these cases are handled in the send section only)
          lpe_tmp(pe_a(jj,row_pe)) = .TRUE.
       END IF
    END DO
    npe_rcv_eq = COUNT(lpe_tmp)

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nrcv_eq(:) = 0
    rcv_eq(:,:,:) = -999

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! icountlats counts the lats I need from pe_a(jj,row_pe) (per hemi.)
          sender = pe_a(jj,row_pe)
          icountlats(sender) = icountlats(sender) + 1
          ! local lats (S)
          ilat(icountlats(sender),1,sender) = nxpt + jintmx - jj + 1
          ! local lats (N)
          ilat(icountlats(sender),2,sender) = jstart + plato2 + jj - 1
       END IF
    END DO

    nrcv_eq(:) = icountlats(:)
    rcv_eq(:,:,:) = ilat(:,:,:)

    !============================================================================

    ! Set up scheme for receives of lats that are sent towards the poles and
    ! possibly across the equator,
    ! i.e. receives for northern overlap in southern hemisphere and 
    ! southern overlap in northern hemisphere.
    ! Determine owner of needed overlap lats and local lat numbers of owner.

    ! determine pe_a globally
    pe_a(:,:) = -999

    DO row_pe = 1,nprocy      ! determine receive info for these PEs
       row_pe_tmp = row_pe
       counter_lats = 0
       row_end = nprocy
       row_stride = 1
       upper_pe = row_pe_tmp + 1
       lfound = .FALSE.
       DO jj = 1,nxpt+jintmx  ! loop over all lats of an overlap region (for
          ! one hemisphere)

          ! Now loop over all PEs north of current PE (south of current PE for 
          ! northern hemisphere) until the PE is found which owns the current
          ! overlap latitude.
10        CONTINUE
          DO upper_pe = row_pe_tmp + 1, row_end, row_stride
             sndpe_id = dc%mapmesh(1,upper_pe)  ! choose any column, e.g. 1
             !             isndpe_id = indx(sndpe_id,gc)
             IF (upper_pe == 1 ) THEN
                sndpelats = gc(sndpe_id)%nglat/2 + nxpt +jintmx
             ELSE
                sndpelats = gc(sndpe_id)%nglat/2
             END IF
             ! if counter of lats for upper_pe <= total #lats for upper_pe,
             ! i.e. if needed lat is owned by upper_pe, then ...
             IF (counter_lats + 1 <= sndpelats) THEN

                ! row_pe rcvs from upper_pe
                receiver = row_pe
                sender = upper_pe
                pe_a(jj,receiver) = sender
                ! lat was found on this upper_pe; increment counter_lats and
                ! start looking for next lat on this same upper_pe
                counter_lats = counter_lats + 1
                lfound = .TRUE.
                EXIT
             END IF
             ! if lat was not found on this upper_pe, set counter_lats=0 and
             ! loop to next upper PE
             counter_lats = 0
          END DO
          ! for next lat, start searching at last upper_pe
          IF (lfound) THEN
             row_pe_tmp    = upper_pe - 1; lfound = .FALSE.; CYCLE
          END IF
          ! if not found on this hemisphere, jump to other hemisphere
          ! and repeat search for current jj
          IF (.NOT. lfound) THEN 
             row_stride = -1 
             row_end = 1
             row_pe_tmp = upper_pe - 2
             lfound = .FALSE.
             GO TO 10
          END IF
       END DO
    END DO

    ! local send values
    pe_snd_po(:) = -999
    npe_snd_po = 0

    DO row_pe = 1,nprocy
       DO jj = 1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN   ! If I am a sender...
             ! number of PEs to send to (inludes myself in cases of copying)
             npe_snd_po = npe_snd_po + 1
             ! dest index
             pe_snd_po(npe_snd_po) = dc%mapmesh(pid_ln,row_pe)
             EXIT
          END IF
       END DO
    END DO

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nsnd_po(:) = 0
    snd_po(:,:,:) = -999

    nh(:,:,:) = 0
    znh(:,:,:) = 0

    ljump(:) =.FALSE.
    icount_jump(:,:) = 0

    DO row_pe =1,nprocy
       DO jj=1,nxpt+jintmx
          IF (pe_a(jj,row_pe) == pid_lt) THEN
             ! icountlats counts the lats row_pe needs from me (per hemisph.)
             icountlats(row_pe) = icountlats(row_pe) + 1
             !================================================
             ! my local lats for northern overlap in southern hem. of receiver
             ih_rcv = 1 ! sending to this hemsiphere of receiver
             IF (pid_lt > row_pe) THEN   ! if I am north of row_pe on southern 
                ! hem., I might have to send from 
                ! both my hemispheres
                ih_snd = 1 ! my hemisphere to send from
                IF (icountlats(row_pe) > mypelats) THEN
                   ih_snd = 2 ! jump to my other hem.
                   ljump(ih_rcv) = .TRUE.
                END IF
             ELSE  ! if I am south of row_pe on southern hem. (or same), I only
                ! have to send from my northern hemisphere
                ih_snd = 2
             END IF
             IF (.NOT.(ljump(ih_rcv))) THEN
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + &
                     &                                                       icountlats(row_pe) - 1 
             ELSE
                icount_jump(ih_rcv,row_pe) = icount_jump(ih_rcv,row_pe) + 1
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + &
                     &                                                icount_jump(ih_rcv,row_pe) - 1
             END IF
             znh(icountlats(row_pe),ih_rcv,row_pe) = ih_snd
             !==================================================
             ! my local lats for southern overlap in northern hem. of receiver
             ih_rcv = 2
             IF (pid_lt > row_pe) THEN ! if I am south of row_pe on northern..
                ih_snd = 2
                IF (icountlats(row_pe) > mypelats) THEN
                   ih_snd = 1 ! jump to my other hem.
                   ljump(ih_rcv) = .TRUE.
                END IF
             ELSE ! if I am north of row_pe on northern hem. (or same), I only
                ! have to send from my southern hem.
                ih_snd = 1
             END IF
             IF (.NOT.(ljump(ih_rcv))) THEN             
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + plato2 - &
                     &                                                        icountlats(row_pe)
             ELSE
                icount_jump(ih_rcv,row_pe) = icount_jump(ih_rcv,row_pe) + 1
                ilat(icountlats(row_pe),ih_rcv,row_pe) = jstart + plato2 - &
                     &                                                   icount_jump(ih_rcv,row_pe)
             END IF
             znh(icountlats(row_pe),ih_rcv,row_pe) = ih_snd
             !==================================================
          END IF
       END DO
    END DO

    nsnd_po(:npe_snd_po) = PACK(icountlats, icountlats /= 0)

    count2=0
    DO row_pe=1,nprocy
       IF (ilat(1,1,row_pe) /= -999) THEN
          count2=count2+1
          snd_po(:,:,count2) = ilat(:,:,row_pe)
          nh(:,:,count2)     = znh(:,:,row_pe)
       ENDIF
    ENDDO

    ! local receive values
    lpe_tmp(:) = .FALSE.
    npe_rcv_po = 0

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       ! If I am a receiver ...
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! number of PEs to receive from (does NOT include myself in cases
          ! of copying; these cases are handled in the send section only)
          lpe_tmp(pe_a(jj,row_pe)) = .TRUE.
       END IF
    END DO
    npe_rcv_po = COUNT(lpe_tmp)

    icountlats(:) = 0
    ilat(:,:,:) = -999
    nrcv_po(:) = 0
    rcv_po(:,:,:) = -999

    icountlats_own = 0
    ilat_own(:,:) = -999
    nrcv_po_own = 0
    rcv_po_own(:,:) = -999

    row_pe = pid_lt
    DO jj = 1,nxpt+jintmx
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) /= pid_lt)) THEN
          ! icountlats counts the lats I need from pe_a(jj,row_pe) (per hemi.)
          sender = pe_a(jj,row_pe)
          icountlats(sender) = icountlats(sender) + 1
          ! local lats (S)
          ilat(icountlats(sender),1,sender) = jstart + plato2 + jj - 1
          ! local lats (N)
          ilat(icountlats(sender),2,sender) = nxpt + jintmx - jj + 1
       END IF
       IF ((pe_a(jj,row_pe) /= -999) .AND. (pe_a(jj,row_pe) == pid_lt)) THEN
          ! icountlats_own counts the lats I need from myself (per hemi.)
          icountlats_own = icountlats_own + 1
          ! local lats (S)
          ilat_own(icountlats_own,1) = jstart + plato2 + jj - 1
          ! local lats (N)
          ilat_own(icountlats_own,2) = nxpt + jintmx - jj + 1
       END IF
    END DO

    nrcv_po(:) = icountlats(:)
    rcv_po(:,:,:) = ilat(:,:,:)

    nrcv_po_own = icountlats_own
    rcv_po_own(:,:) = ilat_own(:,:)

    !============================================================================

    ! the major part of the setup for the longitudinal exchange must be 
    ! carried out dynamically, since the longitudinal overlap is dynamic 
    ! (the latitudinal overlap currently is static)

    ! blks_we(row_pe)/blks_ea(row_pe): 
    ! this many lons available from other PEs closer (than row_pe) to myself 
    ! in the west/east
    IF (nprocx > 1) THEN

       blks_we(:) = 0
       blks_tmp = 0
       col_pe = pid_ln
       DO src_ln = 1,nprocx - 1 ! this many PEs to the west of myself
          col_pe = col_pe - 1
          IF (col_pe < 1) col_pe = nprocx
          src = dc%mapmesh(col_pe,pid_lt)
          !          isrc = indx(src,gc)
          blks_we(col_pe) = blks_tmp
          blks_tmp = blks_tmp + gc(src)%nglon
       END DO

       blks_ea(:) = 0
       blks_tmp = 0
       col_pe = pid_ln
       DO src_ln = 1,nprocx - 1 ! this many PEs to the east of myself
          col_pe = col_pe + 1
          IF (col_pe > nprocx) col_pe = 1
          src = dc%mapmesh(col_pe,pid_lt)
          !          isrc = indx(src,gc)
          blks_ea(col_pe) = blks_tmp
          blks_tmp = blks_tmp + gc(src)%nglon
       END DO

    END IF

    !============================================================================


    ! deallocate local arrays
    DEALLOCATE(pe_a)
    DEALLOCATE(icountlats)
    DEALLOCATE(ilat)
    DEALLOCATE(ilat_own)
    DEALLOCATE(znh)
    DEALLOCATE(lpe_tmp)

    DEALLOCATE(icount_jump)
    DEALLOCATE(ljump)

  END SUBROUTINE setup_overlap

  SUBROUTINE overlap(ub,vb,fb)

    USE mo_grid, ONLY: plond, plev, platd, pcnst, nxpt, jintmx, istart, plon
    USE mo_mpi,  ONLY: p_isend, p_recv, p_probe, p_barrier, p_wait

    REAL ub(plond,plev,platd,2)
    REAL vb(plond,plev,platd,2)
    REAL fb(plond,plev,pcnst,platd,2)

    ! local arrays
    REAL :: buf_snd_eq(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_snd_eq)
    REAL :: buf_rcv_eq(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_rcv_eq)
    REAL :: buf_snd_po(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_snd_po)
    REAL :: buf_rcv_po(plond*plev*(nxpt+jintmx)*(2+pcnst)*2, npe_rcv_po)

    REAL :: buf_snd_ea(nxpt*plev*platd*(2+pcnst)*2,     dc%nprocb - 1)
    REAL :: buf_rcv_ea(nxpt*plev*platd*(2+pcnst)*2,     dc%nprocb - 1)
    REAL :: buf_snd_we((nxpt+1)*plev*platd*(2+pcnst)*2, dc%nprocb - 1)
    REAL :: buf_rcv_we((nxpt+1)*plev*platd*(2+pcnst)*2, dc%nprocb - 1)

    INTEGER :: tagtable_eq(dc%nproca), tagtable_po(dc%nproca)
    INTEGER :: tagtable_we(dc%nprocb), tagtable_ea(dc%nprocb)
    INTEGER :: nxpt_tmp(plev,platd,2)

    ! local scalars
    INTEGER :: nprocy, nprocx, jproc, imeslen, ih, ihrcv, n, lat, k, lon, m, i
    INTEGER :: tag_base_eq, tag_base_po, tagcount_eq, tagcount_po
    INTEGER :: tag_base_we, tag_base_ea, tagcount_we, tagcount_ea
    INTEGER :: tag, count, dest, idest, source, isource, src_lt, src_ln
    INTEGER :: pid_ln, pid_lt, col_pe, diff, iblkbeg, iblkend, total, elements
    REAL    :: rbuf

    !===========================================================================
    !  Update the overlap between processors.
    !  The arrays have been extended as if poles were present.
    !
    !  On entry:
    !    ub, vb and fb are extended arrays with pole extensions
    !    filled for processor zero and interior values
    !    set for all processors.
    !
    !  On exit:
    !    ub, vb and fb have been filled with the extended values
    !    from neighboring processors. The interior
    !    values are unchanged.
    !
    !=========================================================================

    nprocy = dc%nproca
    nprocx = dc%nprocb
    pid_ln = dc%set_b
    pid_lt = dc%set_a

    ! send towards equator

    tag_base_eq = 400+10+nprocx+2+nprocx   ! last change occurred in sltini...

    DO jproc = 1,npe_snd_eq
       imeslen = 0
       DO ih = 1,2
          DO n = 1,nsnd_eq(jproc) ! this many lats per my current hem.
             lat = snd_eq(n,ih,jproc)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   buf_snd_eq(imeslen,jproc) = ub(lon,k,lat,ih)
                   imeslen = imeslen +1 
                   buf_snd_eq(imeslen,jproc) = vb(lon,k,lat,ih)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nsnd_eq(jproc) ! this many lats per my current hem.
             lat = snd_eq(n,ih,jproc)
             DO m =1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      buf_snd_eq(imeslen,jproc) = fb(lon,k,m,lat,ih)
                   END DO
                END DO
             END DO
          END DO
       END DO

       tag = tag_base_eq + jproc
       count = imeslen
       dest = gc(pe_snd_eq(jproc))%pe
       CALL p_isend(buf_snd_eq(1,jproc),dest,tag,p_count=count)

    END DO

    ! send towards poles

    tag_base_po = tag_base_eq + nprocy

    DO jproc = 1,npe_snd_po
       imeslen = 0
       DO ihrcv = 1,2
          DO n = 1,nsnd_po(jproc) ! this many lats per hem. of receiver
             ! (note the difference to the send towards equator) 
             lat = snd_po(n,ihrcv,jproc)
             ih = nh(n,ihrcv,jproc)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   buf_snd_po(imeslen,jproc) = ub(lon,k,lat,ih)
                   imeslen = imeslen +1 
                   buf_snd_po(imeslen,jproc) = vb(lon,k,lat,ih)
                END DO
             END DO
          END DO
       END DO

       DO ihrcv = 1,2
          DO n = 1,nsnd_po(jproc) ! this many lats per hem. of receiver
             ! (note the difference to the send towards equator) 
             lat = snd_po(n,ihrcv,jproc)
             ih = nh(n,ihrcv,jproc)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      buf_snd_po(imeslen,jproc) = fb(lon,k,m,lat,ih)
                   END DO
                END DO
             END DO
          END DO
       END DO

       tag   = tag_base_po + jproc 
       count = imeslen
       dest  = gc(pe_snd_po(jproc))%pe

       IF (dest /= dc%pe) THEN
          CALL p_isend(buf_snd_po(1,jproc),dest,tag,p_count=count)
       ELSE ! copy across equator
          imeslen = 0
          DO ih =1,2
             DO n = 1,nrcv_po_own
                lat = rcv_po_own(n,ih)
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_snd_po(imeslen,jproc)
                      imeslen = imeslen + 1
                      vb(lon,k,lat,ih) = buf_snd_po(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO n = 1,nrcv_po_own
                lat = rcv_po_own(n,ih)
                DO m = 1,pcnst
                   DO k = 1,plev
                      DO lon = 1,plond
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_snd_po(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

       END IF
    END DO

    ! receive latitudes sent towards equator

    ! nprocy could be changed to npe_rcv_eq
    ! this would imply changes at other places...
    tagcount_eq = nprocy
    DO jproc = 1,nprocy
       tagtable_eq(jproc) = tag_base_eq + jproc
    END DO

    DO jproc = 1,npe_rcv_eq

       ! look for tags, not for sources; sources are not unique, could be 
       ! from pole sends as well.
       CALL p_probe(rbuf,tagcount_eq,tagtable_eq,source,tag,count)
       CALL p_recv(buf_rcv_eq(1,jproc),source,tag,p_count=count)

       isource = indx(source,gc)
       src_lt = gc(isource)%set_a

       imeslen = 0
       DO ih = 1,2
          DO n = 1,nrcv_eq(src_lt)
             lat = rcv_eq(n,ih,src_lt)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   ub(lon,k,lat,ih) = buf_rcv_eq(imeslen,jproc)
                   imeslen = imeslen +1 
                   vb(lon,k,lat,ih) = buf_rcv_eq(imeslen,jproc)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nrcv_eq(src_lt)
             lat = rcv_eq(n,ih,src_lt)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      fb(lon,k,m,lat,ih) = buf_rcv_eq(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    ! receive latitudes sent towards poles 

    ! nprocy could be changed to npe_rcv_po
    ! this would imply changes at other places...
    tagcount_po = nprocy
    DO jproc = 1,nprocy
       tagtable_po(jproc) = tag_base_po + jproc
    END DO

    DO jproc = 1,npe_rcv_po

       ! look for tags, not for sources; sources are not unique, could be 
       ! from pole sends as well.
       CALL p_probe(rbuf,tagcount_po,tagtable_po,source,tag,count)
       CALL p_recv(buf_rcv_po(1,jproc),source,tag,p_count=count)

       isource = indx(source,gc)
       src_lt = gc(isource)%set_a

       imeslen = 0
       DO ih = 1,2
          DO n = 1,nrcv_po(src_lt)
             lat = rcv_po(n,ih,src_lt)
             DO k = 1,plev
                DO lon = 1,plond
                   imeslen = imeslen + 1
                   ub(lon,k,lat,ih) = buf_rcv_po(imeslen,jproc)
                   imeslen = imeslen +1 
                   vb(lon,k,lat,ih) = buf_rcv_po(imeslen,jproc)
                END DO
             END DO
          END DO
       END DO

       DO ih = 1,2
          DO n = 1,nrcv_po(src_lt)
             lat = rcv_po(n,ih,src_lt)
             DO m = 1,pcnst
                DO k = 1,plev
                   DO lon = 1,plond
                      imeslen = imeslen + 1
                      fb(lon,k,m,lat,ih) = buf_rcv_po(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO
       END DO

    END DO

    !===================================================

    IF (nprocx > 1) THEN

       ! send towards east

       tag_base_ea = tag_base_po + nprocy

       col_pe = pid_ln
       diff = 0

       col:   DO jproc = 1,nprocx-1  ! max # of PEs to send to

          col_pe = col_pe + 1
          IF (col_pe > nprocx) col_pe = 1
          dest = gc(dc%mapmesh(col_pe,pid_lt))%pe
          tag = tag_base_ea + jproc

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = 1
                   iblkend = MIN(plon,nxpt_a(k,lat,ih) - diff)
                   DO i = iblkbeg, iblkend
                      lon = istart + plon - i
                      imeslen = imeslen + 1
                      buf_snd_ea(imeslen,jproc) = ub(lon,k,lat,ih)
                      imeslen = imeslen +1 
                      buf_snd_ea(imeslen,jproc) = vb(lon,k,lat,ih)
                   END DO
                END DO
             END DO
          END DO

          IF (imeslen == 0) EXIT col

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = 1
                      iblkend = MIN(plon,nxpt_a(k,lat,ih) - diff)
                      DO i = iblkbeg, iblkend
                         lon = istart + plon - i
                         imeslen = imeslen + 1
                         buf_snd_ea(imeslen,jproc) = fb(lon,k,m,lat,ih)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          count = imeslen
          CALL p_isend(buf_snd_ea(1,jproc),dest,tag,p_count=count)

          idest = indx(dest,gc)
          diff = diff + gc(idest)%nglon

       END DO col

       ! send towards west
       ! note the extra line sent towards the west (=to the eastern overlaps)

       tag_base_we = tag_base_ea + nprocx

       col_pe = pid_ln
       diff = 0

       col2:   DO jproc = 1,nprocx-1  ! max # of PEs to send to

          col_pe = col_pe - 1
          IF (col_pe < 1) col_pe = nprocx
          dest = gc(dc%mapmesh(col_pe,pid_lt))%pe
          tag = tag_base_we + jproc

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = 1
                   iblkend = MIN(plon,nxpt_a(k,lat,ih) + 1 - diff)
                   DO i = iblkbeg, iblkend
                      lon = istart - 1 + i
                      imeslen = imeslen + 1
                      buf_snd_we(imeslen,jproc) = ub(lon,k,lat,ih)
                      imeslen = imeslen +1 
                      buf_snd_we(imeslen,jproc) = vb(lon,k,lat,ih)
                   END DO
                END DO
             END DO
          END DO

          IF (imeslen == 0) EXIT col2

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = 1
                      iblkend = MIN(plon,nxpt_a(k,lat,ih) + 1 - diff)
                      DO i = iblkbeg, iblkend
                         lon = istart - 1 + i
                         imeslen = imeslen + 1
                         buf_snd_we(imeslen,jproc) = fb(lon,k,m,lat,ih)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          count = imeslen
          CALL p_isend(buf_snd_we(1,jproc),dest,tag,p_count=count)

          idest = indx(dest,gc)
          diff = diff + gc(idest)%nglon

       END DO col2

       ! receive longitudes sent towards east (i.e. receive my western overlap)

       ! nprocx could be changed to actual number of PEs to receive from
       ! this would imply changes at other places...
       tagcount_ea = nprocx - 1
       DO jproc = 1,nprocx - 1
          tagtable_ea(jproc) = tag_base_ea + jproc
       END DO

       total = (2+pcnst)*SUM(nxpt_a)
       elements = 0

       DO jproc = 1,nprocx - 1 ! max. # of PEs to receive from

          CALL p_probe(rbuf,tagcount_ea,tagtable_ea,source,tag,count)
          CALL p_recv(buf_rcv_ea(1,jproc),source,tag,p_count=count)

          isource = indx(source,gc)
          src_ln = gc(isource)%set_b

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = blks_we(src_ln) + 1
                   iblkend = iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) - &
                        &                                                          blks_we(src_ln)) - 1
                   DO i = iblkbeg, iblkend
                      lon = istart - i                
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_rcv_ea(imeslen,jproc)
                      imeslen = imeslen +1 
                      vb(lon,k,lat,ih) = buf_rcv_ea(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = blks_we(src_ln) + 1
                      iblkend = iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) - &
                           &                                                          blks_we(src_ln)) - 1
                      DO i = iblkbeg, iblkend
                         lon = istart - i                                   
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_rcv_ea(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          elements = elements + imeslen
          IF (elements == total) EXIT 

       END DO

       ! receive longitudes sent towards west (i.e. receive my eastern overlap)

       ! nprocx could be changed to actual number of PEs to receive from
       ! this would imply changes at other places...
       tagcount_we = nprocx - 1
       DO jproc = 1,nprocx - 1
          tagtable_we(jproc) = tag_base_we + jproc
       END DO

       nxpt_tmp(:,:,:) = nxpt_a(:,:,:) + 1
       total = (2+pcnst)*SUM(nxpt_tmp)
       elements = 0

       DO jproc = 1,nprocx - 1 ! max. # of PEs to receive from

          CALL p_probe(rbuf,tagcount_we,tagtable_we,source,tag,count)
          CALL p_recv(buf_rcv_we(1,jproc),source,tag,p_count=count)

          isource = indx(source,gc)
          src_ln = gc(isource)%set_b

          imeslen = 0
          DO ih = 1,2
             DO lat = 1,platd
                DO k = 1,plev
                   iblkbeg = blks_ea(src_ln) + 1
                   iblkend =iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) + 1&
                        &                                                        - blks_ea(src_ln)) - 1
                   DO i = iblkbeg, iblkend
                      lon = istart + plon - 1 + i            
                      imeslen = imeslen + 1
                      ub(lon,k,lat,ih) = buf_rcv_we(imeslen,jproc)
                      imeslen = imeslen +1 
                      vb(lon,k,lat,ih) = buf_rcv_we(imeslen,jproc)
                   END DO
                END DO
             END DO
          END DO

          DO ih = 1,2
             DO lat = 1,platd
                DO m = 1,pcnst
                   DO k = 1,plev
                      iblkbeg = blks_ea(src_ln) + 1
                      iblkend =iblkbeg + MIN(gc(isource)%nglon, nxpt_a(k,lat,ih) + 1&
                           &                                                        - blks_ea(src_ln)) - 1
                      DO i = iblkbeg, iblkend
                         lon = istart + plon - 1 + i                
                         imeslen = imeslen + 1
                         fb(lon,k,m,lat,ih) = buf_rcv_we(imeslen,jproc)
                      END DO
                   END DO
                END DO
             END DO
          END DO

          elements = elements + imeslen
          IF (elements == total) EXIT 

       END DO

       !========================================================================

       ! end of the multiple longitudinal processors if
    END IF

    CALL p_wait

    !    CALL p_barrier

  END SUBROUTINE OVERLAP

END MODULE mo_parallel_slt
