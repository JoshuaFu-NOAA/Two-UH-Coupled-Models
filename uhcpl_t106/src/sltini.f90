SUBROUTINE sltini(coslam,sinlam,dlam,ub,vb,fb,uxl,uxr,qxl,qxr)

  ! Description:
  !
  !...Prepare the extended arrays for use in the SLT routines
  !
  !   1)  Fill latitude extensions.
  !   2)  Fill longitude extensions.
  !   3)  Compute x-derivatives
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! Parallel version:  T. Diehl, DKRZ, July 1999
  ! 
  ! for more details see file AUTHORS
  !


  USE mo_grid,          ONLY: plon, plond, plev, platd, pcnst, nxpt, jintmx, &
                              plato2
  USE mo_mpi,           ONLY: p_probe,p_recv, p_send, p_isend, p_wait,&
                              p_barrier, p_communicator_d
  USE mo_slt,           ONLY: plimdr
  USE mo_parallel_slt,  ONLY: nxpt_a, overlap
  USE mo_decomposition, ONLY: gc => global_decomposition
  USE mo_decomposition, ONLY: dc => local_decomposition, debug_seriell
  USE mo_transpose,     ONLY: indx
  
  IMPLICIT NONE

  REAL :: coslam(plon), sinlam(plon)
  REAL :: dlam
  REAL :: ub(plond,plev,platd,2), vb(plond,plev,platd,2)
  REAL :: fb(plond,plev,pcnst,platd,2)
  REAL :: uxl(plond,plev,2,platd,2), uxr(plond,plev,2,platd,2)
  REAL :: qxl(plond,plev,pcnst,platd,2), qxr(plond,plev,pcnst,platd,2)

  INTEGER       :: itmp(plev,platd,2)
  INTEGER       :: nxpt_a_g(plev,platd,2,dc%nprocb)
  INTEGER       :: puvpts, pqpts, ih, j, k, jjfirst, jjlast, nei_nxpt
  INTEGER       :: tag_base, tag, count, dest, source, tag2, count2
  INTEGER       :: tagcount, tagtable(dc%nprocb), jproc, ibuf, isource, src_ln
  INTEGER       :: p_pe
  LOGICAL, SAVE :: first
  DATA first /.TRUE./

  ! external subroutines
  EXTERNAL extyv, extys, extx, cubxdr, limdx

  ! executable statements

  ibuf = 0
  puvpts = plond*plev
  pqpts = plond*plev*pcnst

  p_pe = dc%pe

  ! Fill latitude extensions beyond the southern- and northern-most
  ! latitudes in the global grid

  CALL extyv(plev,coslam,sinlam,ub)
  CALL extyv(plev,coslam,sinlam,vb)

  ! loop over fields now in extys

  CALL extys(pcnst,plev,fb)

  !=============================
  !either call p_barrier or use different message tags in following section
  
  CALL p_barrier (p_communicator_d)

  ! Begin of parallel extensions
  !======================================================================
  ! Synchronize the nxpt_a array for this timestep
  ! (The local update takes place in sltb1)

  IF (first .OR. dc%nprocb == 1) THEN
     !  Initialize the nxpt_a array
     nxpt_a(:,:,:) = nxpt
     first = .FALSE.
  ELSE
!     IF (.not. debug_seriell .or. dc%nprocb > 1) THEN
        
        !   Be sure the overlap extension is enough to support interpolation
        
        jjfirst = nxpt + jintmx + 1
        DO ih = 1,2
           DO j = 1,nxpt+jintmx
              DO k = 1,plev
                 itmp(k,j,ih) = nxpt_a(k,jjfirst,ih)
              END DO
           END DO
        END DO
        jjlast = plato2 + nxpt + jintmx
        DO ih = 1,2
           DO j = jjlast+1,platd
              DO k = 1,plev
                 itmp(k,j,ih) = nxpt_a(k,jjlast,ih)
              END DO
           END DO
        END DO

        DO ih = 1,2
           DO j = jjfirst,jjlast
              DO k = 1,plev
                 nei_nxpt = 0
                 IF (k > 1) nei_nxpt = MAX(nxpt_a(k-1,j,ih),nei_nxpt)
                 nei_nxpt = MAX(nxpt_a(k,j,ih),nei_nxpt)
                 IF (k <= plev-1) nei_nxpt = MAX(nxpt_a(k+1,j,ih),nei_nxpt)
                 IF (k <= plev-2) nei_nxpt = MAX(nxpt_a(k+2,j,ih),nei_nxpt)
                 nei_nxpt = MAX(nxpt_a(k,j-1,ih),nei_nxpt)
                 nei_nxpt = MAX(nxpt_a(k,j,ih),nei_nxpt)
                 nei_nxpt = MAX(nxpt_a(k,j+1,ih),nei_nxpt)
                 nei_nxpt = MAX(nxpt_a(k,j+2,ih),nei_nxpt)
                 itmp(k,j,ih) = nei_nxpt
              END DO
           END DO
        END DO
        DO ih = 1,2
           DO j = 1,platd
              DO k = 1,plev
                 nxpt_a(k,j,ih) = itmp(k,j,ih)
              END DO
           END DO
        END DO

        ! Synchronize the nxpt_a across the processors in the row

        ! On entry nxpt_a contains local data; on exit it contains the 
        ! vector max of the processor row

        tag_base = 400 + 10 + dc%nprocb  ! last change occurred in extys
        IF (dc%set_b /= 1) THEN
           tag = tag_base + 1
           count = SIZE(itmp)
           dest = gc(dc%mapmesh(1,dc%set_a))%pe
           ! send to column 1
           CALL p_isend(itmp,dest,tag,p_count=count)
           ! receive updated values from column 1
           source = gc(dc%mapmesh(1,dc%set_a))%pe
           tag2 = tag_base + 2
           count2 = SIZE(nxpt_a)
           CALL p_recv(nxpt_a,source,tag2,p_count=count2)
        ELSE
           ! collect
           nxpt_a_g(:,:,:,1) = itmp(:,:,:)
           tagcount = 1
           tagtable(1) = tag_base + 1
           DO jproc = 1,dc%nprocb - 1
              CALL p_probe(ibuf,tagcount,tagtable,source,tag,count)
              isource = indx(source,gc)
              src_ln = gc(isource)%set_b
              CALL p_recv(nxpt_a_g(1,1,1,src_ln),source,tag,p_count=count)
           END DO

           nxpt_a = MAXVAL(nxpt_a_g, DIM=4)

           ! distribute

           DO jproc = 1,dc%nprocb
              IF (gc(dc%mapmesh(jproc,dc%set_a))%pe /= dc%pe) THEN
                 dest = gc(dc%mapmesh(jproc,dc%set_a))%pe
                 tag = tag_base + 2
                 count = SIZE(nxpt_a)
                 CALL p_send(nxpt_a,dest,tag,p_count=count)
              END IF
           END DO
        END IF

!     END IF
  END IF


  CALL p_wait

  !  Update the overlap between processors
  !  The arrays have been extended as if poles were present
  
  CALL overlap(ub,vb,fb)

  ! end of parallel extensions
  !====================================================================

  ! Fill longitude extensions

  CALL extx(1,plev,ub(1,1,1,1))
  CALL extx(1,plev,ub(1,1,1,2))
  CALL extx(1,plev,vb(1,1,1,1))
  CALL extx(1,plev,vb(1,1,1,2))

  ! loop over fields now in extx     

  CALL extx(pcnst,plev,fb(1,1,1,1,1))
  CALL extx(pcnst,plev,fb(1,1,1,1,2))

  ! Compute x-derivatives

  DO ih = 1,2
     DO j = 1, platd
        CALL cubxdr(puvpts,2,puvpts-3,dlam,ub(1,1,j,ih), &
                    uxl(1,1,1,j,ih),uxr(1,1,1,j,ih))
        CALL cubxdr(puvpts,2,puvpts-3,dlam,vb(1,1,j,ih), &
                    uxl(1,1,2,j,ih),uxr(1,1,2,j,ih))
        CALL cubxdr(pqpts,2,pqpts-3,dlam,fb(1,1,1,j,ih), &
                    qxl(1,1,1,j,ih),qxr(1,1,1,j,ih))
        IF (plimdr) THEN
           CALL limdx(pqpts,2,pqpts-3,dlam,fb(1,1,1,j,ih), &
                      qxl(1,1,1,j,ih),qxr(1,1,1,j,ih))
        END IF
     END DO
  END DO

  CALL p_barrier (p_communicator_d)
 
END SUBROUTINE sltini
