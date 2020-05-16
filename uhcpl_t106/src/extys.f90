!+ fill latitude extensions of a scalar extended array
!+ $Id: extys.f90,v 1.20 2000/01/24 07:52:52 m214003 Exp $

SUBROUTINE extys(pkcnst,pkdim,fb)

  ! Description:
  !
  ! Fill latitude extensions of a scalar extended array.
  !
  ! Method:
  !
  ! This is done in 2 steps:
  !
  ! 1) Interpolate to the pole points; use the mean field value on the
  !    Gaussian latitude closest to the pole.
  ! 2) Add latitude lines beyond the poles.
  !
  ! Authors:
  !
  ! Original version:  J. Olson
  ! Standardized:      J. Rosinski, June 1992
  ! Reviewed:          D. Williamson, P. Rasch, August 1992
  ! f90 version:       L. Kornblueh, U. Schulzweida, May 1998
  ! Parallel version:  T. Diehl, July 1999
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_grid,          ONLY: plond, platd, plon, nxpt, jintmx, istart, &
                              istop, jn, js, plono2
  USE mo_decomposition, ONLY: dc=>local_decomposition
  USE mo_decomposition, ONLY: gc=>global_decomposition
  USE mo_transpose,     ONLY: indx
  USE mo_mpi,           ONLY: p_probe, p_recv, p_isend, p_wait

  IMPLICIT NONE

  !  Scalar arguments with intent(In):
  INTEGER, INTENT (IN) :: pkcnst, pkdim

  !  Array arguments with intent(InOut):
  REAL, INTENT (INOUT) :: fb(plond,pkdim,pkcnst,platd,2)

  ! pkcnst  Number of tracers
  ! pkdim   Vertical dimension
  ! fb      Output is same as on entry except with the pole latitude
  !         and extensions beyond it filled.

  ! Local arrays:
  REAL :: fb_buf(plon,2,pkdim), zbuf(2,pkdim)
  REAL :: fb_buf_g((plon+1)*2*pkdim,dc%nprocb)
  REAL :: fb_g(plon+1,2,pkdim,dc%nprocb)
  ! "plon+1" in fb_buf_g: need max. local #of lons
  REAL :: fb_g_one(dc%nlon,2,pkdim)
  REAL :: fbuf1(plon*pkdim*(nxpt+jintmx-1)*2,dc%nprocb)
  REAL :: fbuf2(plon*pkdim*(nxpt+jintmx-1)*2,dc%nprocb)
  INTEGER :: tagtable(dc%nprocb), pe_id2(dc%nprocb), pe_id(dc%nprocb)
  INTEGER :: nlons(dc%nprocb), ilon(plon,dc%nprocb)

  ! Local scalars: 
  REAL :: zave         ! accumulator for zonal averaging
  REAL :: fnlon        ! REAL(dc%nlon)
  REAL :: rbuf
  INTEGER :: i, ig, j, k
  INTEGER :: tag_base, tagcount, jproc, isource, src_ln
  INTEGER :: src, ig_mirror, node, inode, il, npes, idest
  INTEGER :: imeslen, jj, ii, ip, m
  INTEGER :: tagr1, mcountr1, sourcer1, tags1, mcounts1, dests1
  INTEGER :: tagr2, mcountr2, sourcer2, tags2, mcounts2, dests2
  INTEGER :: tagr3, mcountr3, sourcer3, tags3, mcounts3, dests3
  INTEGER :: len

  !  Intrinsic functions 
  INTRINSIC REAL, COUNT, PACK, SIZE


  !  Executable statements

  ! return, if I am not a pole PE
  ! This version currently works only if "nxpt+jintmx-1 <= plat/2";
  ! (nxpt+jintmx-1 = the lines beyond the pole line)
  ! otherwise, communication with non pole PEs will be necessary

  IF (dc%set_a /= 1) RETURN

  fnlon = REAL(dc%nlon)
  tag_base = 400
  fb_buf_g(:,:) = 0.0
  fb_buf(:,:,:) = 0.0
  fb_g(:,:,:,:) = 0.0

  ! loop over fields
  DO m=1,pkcnst

     !=======================================================================  
     ! Fill north and south pole line.

     IF (dc%set_a == 1) THEN  ! if I am a pole PE
  
     DO k = 1, pkdim
        ig = 0
        DO i = istart, istop
           ig = ig + 1
        
           IF (dc%set_b /= 1) THEN
              ! north
              fb_buf(ig,2,k) = fb(i,k,m,jn,2)
              ! south
              fb_buf(ig,1,k) = fb(i,k,m,js,1)
           ELSE
              ! north
              fb_g(ig,2,k,1) = fb(i,k,m,jn,2)
              ! south
              fb_g(ig,1,k,1) = fb(i,k,m,js,1)
           END IF
           
        END DO
     END DO
  
     IF (dc%nprocb > 1) THEN
        IF (dc%set_b /= 1) THEN
           tags1 = tag_base + m
           mcounts1 = SIZE(fb_buf)
           dests1 = gc(dc%mapmesh(1,dc%set_a))%pe
           CALL p_isend(fb_buf,dests1,tags1,p_count=mcounts1)
        ELSE
           tagcount = 1
           tagtable(1) = tag_base + m
           DO jproc = 1,dc%nprocb - 1
              CALL p_probe(rbuf,tagcount,tagtable,sourcer1,tagr1,mcountr1)
              isource = indx(sourcer1,gc)
              src_ln = gc(isource)%set_b 
              CALL p_recv(fb_buf_g(1,src_ln),sourcer1,tagr1,p_count=mcountr1)

              len=gc(isource)%nglon
              ig = 0
              DO k=1,pkdim
  	        fb_g(1:len,1,k,src_ln) = fb_buf_g(ig+1:ig+len,src_ln)
                ig =ig + len
	        fb_g(1:len,2,k,src_ln) = fb_buf_g(ig+1:ig+len,src_ln)
                ig =ig + len
              END DO

           END DO
        END IF
     END IF

     IF (dc%set_b == 1) THEN
        DO k = 1,pkdim
           ii = 0
           DO src_ln = 1,dc%nprocb
              src = dc%mapmesh(src_ln,dc%set_a)
              DO i = 1,gc(src)%nglon
                 ii = ii + 1
                 fb_g_one(ii,2,k) = fb_g(i,2,k,src_ln)
                 fb_g_one(ii,1,k) = fb_g(i,1,k,src_ln)
              END DO
           END DO
        END DO
        
        zbuf(:,:) = 0.0
        DO k = 1,pkdim
           DO i = 1,dc%nlon
              ! north
              zbuf(2,k) = zbuf(2,k) + fb_g_one(i,2,k)
              ! south
              zbuf(1,k) = zbuf(1,k) + fb_g_one(i,1,k)
           END DO
        END DO

        IF (dc%nprocb > 1) THEN
           DO jproc = 1,dc%nprocb
              IF (gc(dc%mapmesh(jproc,dc%set_a))%pe /= dc%pe) THEN
                 dests2 = gc(dc%mapmesh(jproc,dc%set_a))%pe
                 tags2 = tag_base + pkcnst + m
                 mcounts2 = SIZE(zbuf)
                 CALL p_isend(zbuf,dests2,tags2,p_count=mcounts2)
              END IF
           END DO
        END IF
     ELSE
	tagcount = 1
	tagtable(1) = tag_base + pkcnst + m
	CALL p_probe(rbuf,tagcount,tagtable,sourcer2,tagr2,mcountr2)
        CALL p_recv(zbuf,sourcer2,tagr2,p_count=mcountr2)
     END IF

     DO k = 1,pkdim
        ! north
        zave = zbuf(2,k)/fnlon
        DO i = istart, istop
           fb(i,k,m,jn+1,2) = zave
        END DO
        ! south
        zave = zbuf(1,k)/fnlon 
        DO i = istart, istop
           fb(i,k,m,js-1,1) = zave
        END DO
     END DO

     call p_wait

  END IF

  END DO

  !=======================================================================
  ! Fill lines beyond pole line
  ! Case only one pole PE
  IF (dc%nprocb == 1) THEN
     DO m=1,pkcnst  ! loop over fields
        ! Fill northern lines beyond pole line.
        DO j = jn + 2, platd
           DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i = istart, istart + plono2 - 1
                 fb(i,k,m,j,2) = fb(plono2+i,k,m,2*jn+2-j,2)
                 fb(plono2+i,k,m,j,2) = fb(i,k,m,2*jn+2-j,2)
              END DO
           END DO
        END DO

        ! Fill southern lines beyond pole line.   
        DO j = 1, js - 2
           DO k = 1, pkdim
!DIR$ IVDEP
!OCL NOVREC
              DO i = istart, istart + plono2 - 1
                 fb(i,k,m,j,1) = fb(plono2+i,k,m,2*js-2-j,1)
                 fb(plono2+i,k,m,j,1) = fb(i,k,m,2*js-2-j,1)
              END DO
           END DO
        END DO
     END DO

  ELSE  ! Case more than one pole PE

     pe_id2(:) = -999
     pe_id(:) = -999
     nlons(:) = 0
     ilon(:,:) = -999

     il = 0
     DO ig = dc%glons(1),dc%glone(1)
        il = il + 1
        IF ((ig >= 1) .AND. (ig <= dc%nlon/2)) THEN
           ig_mirror = dc%nlon/2 + ig
        ELSE
           ig_mirror = ig - dc%nlon/2
        END IF
        
        DO jproc = 1,dc%nprocb
           inode = dc%mapmesh(jproc,dc%set_a)
           node = gc(inode)%pe
           IF (node /= dc%pe) THEN
              IF ((ig_mirror >= gc(inode)%glons(1)) .AND. &
&                 (ig_mirror <= gc(inode)%glone(1))) THEN
                 pe_id2(jproc) = node
                 nlons(jproc) = nlons(jproc) + 1
                 ilon(nlons(jproc),jproc) = istart + il - 1
                 EXIT
              END IF
           END IF
        END DO
     END DO

     npes = COUNT(pe_id2 /= -999)
     pe_id(:npes) = PACK(pe_id2,pe_id2 /= -999)

     fbuf1(:,:) = 0.0
     fbuf2(:,:) = 0.0

     DO m=1,pkcnst ! loop over fields

        ! send part
        DO ip = 1,npes
           dests3 = pe_id(ip)
           idest = indx(dests3,gc)
           jproc = gc(idest)%set_b
           imeslen = 0

           ! north pole
           jj = 0
           DO j = jn+2,platd
              jj = jj+1
              DO k=1,pkdim
!DIR$ IVDEP
!OCL NOVREC
                 DO i=1,nlons(jproc)
                    ii = ilon(i,jproc)
                    imeslen = imeslen + 1
                    fbuf1(imeslen,ip) = fb(ii,k,m,2*jn+2-j,2)
                 END DO
              END DO
           END DO

           ! south pole
           jj = 0
           DO j = 1,js - 2
              jj = jj + 1
              DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
                 DO i=1,nlons(jproc)
                    ii = ilon(i,jproc)
                    imeslen = imeslen + 1
                    fbuf1(imeslen,ip) = fb(ii,k,m,2*js - 2 -j,1)
                 END DO
              END DO
           END DO
           
           tags3 = tag_base + 2*pkcnst + 10 + (dests3+1)*m
           mcounts3 = imeslen
           ! Note that mcount=0 if jn+2>platd or js-2<1
           CALL p_isend(fbuf1(1,ip),dests3,tags3,p_count=mcounts3)
        END DO

        ! receive part
        tagcount = npes
        DO ip = 1,npes
           tagtable(ip) = tag_base + 2*pkcnst + 10 + (dc%pe + 1)*m
        END DO
        
        DO ip = 1,npes
           CALL p_probe(rbuf,tagcount,tagtable,sourcer3,tagr3,mcountr3)
           CALL p_recv(fbuf2(1,ip),sourcer3,tagr3,p_count=mcountr3)
           
           isource = indx(sourcer3,gc)
           jproc = gc(isource)%set_b
           imeslen = 0
           
           ! north pole
           jj = 0
           DO j = jn+2,platd
              jj = jj + 1
              DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
                 DO i = 1,nlons(jproc)
                    ii = ilon(i,jproc)
                    imeslen = imeslen + 1
                    fb(ii,k,m,j,2) = fbuf2(imeslen,ip)
                 END DO
              END DO
           END DO

           ! south pole
           jj = 0
           DO j = 1,js-2
              jj = jj + 1
              DO k = 1,pkdim
!DIR$ IVDEP
!OCL NOVREC
                 DO i = 1,nlons(jproc)
                    ii = ilon(i,jproc)
                    imeslen = imeslen + 1
                    fb(ii,k,m,j,1) = fbuf2(imeslen,ip)
                 END DO
              END DO
           END DO

        END DO

call p_wait

     END DO

  END IF

  RETURN
END SUBROUTINE extys
