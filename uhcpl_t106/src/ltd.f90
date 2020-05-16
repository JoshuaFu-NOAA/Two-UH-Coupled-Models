!+ direct Legendre transforms
!+ $Id: ltd.f90,v 1.24 2000/02/14 13:21:37 m214003 Exp $

!OCL NOALIAS

#ifdef CRAY
#define dgemm sgemm
#endif

SUBROUTINE ltd

  ! Description:
  !
  ! Direct Legendre transforms.
  !
  ! Method:
  !
  ! This subroutine performs direct *Legendre transforms for
  ! the divergence equation,
  ! the temperature and surface pressure equations,
  ! the vorticity equation,
  ! the mean wind.
  !
  ! *ltd* is called from *scan1sl*
  !
  ! Results:
  ! *ltd* adds in the spectral arrays:-
  !      *ld*  the contribution of the current latitude line
  !      *ltp* the contribution of the current latitude line
  !      *lvo* the contribution of the current latitude line
  !      *lu0* the contribution of the current latitude line
  !
  ! Authors:
  !
  ! D. W. Dent, ECMWF, December 1984, original source
  ! U. Schlese, DKRZ, in 1991, and 1994, changed
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_memory_ls,     ONLY: ld, ltp, lu0, lvo
  USE mo_memory_f,      ONLY: fadl, fadm, far, fatp1, faul, fazl, fazm, fsdl, &
                              fsdm, fsr, fstp1, fsul, fszl, fszm
  USE mo_legendre,      ONLY: legmod
  USE mo_decomposition, ONLY: lc => local_decomposition, debug_seriell

  IMPLICIT NONE

  !  Local loop bounds

  INTEGER          :: nllev, nllevp1, nlmp1, nlnm0, lnsp
  INTEGER ,POINTER :: nlmp(:), nlnp(:)

  !  Local scalars: 

  INTEGER :: ims, ims2, ins, ins2, inu, is, iu, j, jh, jn, jhr, nhgl

  LOGICAL :: oldversion  ! .true. to get old results

  REAL, POINTER :: fdl (:,:,:)
  REAL, POINTER :: fdm (:,:,:)
  REAL, POINTER :: fr  (:,:,:)
  REAL, POINTER :: ftp1(:,:,:)
  REAL, POINTER :: fzl (:,:,:)
  REAL, POINTER :: fzm (:,:,:)
  REAL, POINTER :: ful (:)

  REAL, POINTER :: fdl2 (:,:,:,:)
  REAL, POINTER :: fdm2 (:,:,:,:)
  REAL, POINTER :: fr2  (:,:,:,:)
  REAL, POINTER :: ftp12(:,:,:,:)
  REAL, POINTER :: fzl2 (:,:,:,:)
  REAL, POINTER :: fzm2 (:,:,:,:)

  REAL :: pnmt(lc% nlat/2, lc% lnsp)
  REAL :: anmt(lc% nlat/2, lc% lnsp)
  REAL :: rnmt(lc% nlat/2, lc% lnsp)

  !  External subroutines 
  EXTERNAL dgemm


  !  Executable statements 

  oldversion = debug_seriell

!-- Set local loop bounds

  nllev   =  lc% nllev    ! number of levels
  nllevp1 =  lc% nllevp1  ! number of levels + 1
  nlmp1   =  lc% nlm      ! number of m wave numbers
  nlnm0   =  lc% nlnm0    ! number of coefficients for m=0
  lnsp    =  lc% lnsp     ! number of complex spectral coefficients on this pe
  nlmp    => lc% nlmp     ! displacement of the first point of columns
  nlnp    => lc% nlnp     ! number of points on each column
  nhgl    =  lc% nlat/2   ! global halv number of gaussian latitudes


    ! zero fields before summation

  lu0 = 0.

  IF (oldversion) THEN
    ld  = 0.
    ltp = 0.
    lvo = 0.
  ENDIF

  ! derive local wavenumber index
  ! calculate legendre coefficents for each latitude

  CALL legmod(pnmt=pnmt,anmt=anmt,rnmt=rnmt)

!-- 1. Legendre transforms

  ! Loop over hemispheres (1: north, 2: south)

  DO jh = 1, 2
    iu = 2 - jh

    IF (jh==1) THEN
      fdl2  => fsdl (:,:,:,:)
      fdm2  => fadm (:,:,:,:)
      fr2   => fsr  (:,:,:,:)
      fzl2  => fszl (:,:,:,:)
      fzm2  => fazm (:,:,:,:)
      ftp12 => fstp1(:,:,:,:)
    ELSE
      fdl2  => fadl (:,:,:,:)
      fdm2  => fsdm (:,:,:,:)
      fr2   => far  (:,:,:,:)
      fzl2  => fazl (:,:,:,:)
      fzm2  => fszm (:,:,:,:)
      ftp12 => fatp1(:,:,:,:)
    END IF

!-- 1.1 Transforms for d, vo, t and p

  ! Loop over latitudes

    DO jhr = 1, nhgl

      IF (oldversion) THEN

        fdl  => fdl2 (:,:,:,jhr)
        fdm  => fdm2 (:,:,:,jhr)
        fr   => fr2  (:,:,:,jhr)
        fzl  => fzl2 (:,:,:,jhr)
        fzm  => fzm2 (:,:,:,jhr)
        ftp1 => ftp12(:,:,:,jhr)

        DO j = 1, nlmp1

          ims  = nlmp(j) - iu   ! offset to local m columns  (spectral coef.)
          ins  = nlnp(j) + iu   ! column length

          DO jn = 2, ins, 2
            is  = ims  + jn
            ld(:,:,is)  = ld(:,:,is) + fdl(:,:,j)*pnmt(jhr,is) &
                                     - fdm(:,:,j)*anmt(jhr,is) &
                                     + fr (:,:,j)*rnmt(jhr,is)
            lvo(:,:,is) = lvo(:,:,is) + fzl (:,:,j)*pnmt(jhr,is) &
                                      - fzm (:,:,j)*anmt(jhr,is)
            ltp(:,:,is) = ltp(:,:,is) + ftp1(:,:,j)*pnmt(jhr,is)
          END DO
        END DO

      ENDIF

!-- 1.2 Transforms for the mean wind

      IF (jh==1) THEN
        ful  => fsul (:,    jhr)
      ELSE
        ful  => faul (:,    jhr)
      END IF

      inu = (nlnm0+iu)/2

      lu0(:,jh::2) = lu0(:,jh::2) + SPREAD(ful(:),2,inu) &
                                  * SPREAD(pnmt(jhr,jh:nlnm0:2),1,nllev)
    END DO ! jhr = 1, nhgl

    IF(.NOT.oldversion) THEN

      DO j = 1, nlmp1

        ims  = nlmp(j) - iu   ! offset to local m columns  (spectral coef.)
        ins  = nlnp(j) + iu   ! column length
        ins2 = ins/2
        ims2 = ims+2

        IF (nllev > 0 .AND. ins2 > 0) THEN
          CALL dgemm ('N','N',2*nllev,ins2,nhgl, 1.,fdl2 (1,1,j,1),2*nllev*nlmp1,pnmt(1,ims2),2*nhgl,0.,ld (1,1,ims2),4*nllev)
          CALL dgemm ('N','N',2*nllev,ins2,nhgl,-1.,fdm2 (1,1,j,1),2*nllev*nlmp1,anmt(1,ims2),2*nhgl,1.,ld (1,1,ims2),4*nllev)
          CALL dgemm ('N','N',2*nllev,ins2,nhgl,+1.,fr2  (1,1,j,1),2*nllev*nlmp1,rnmt(1,ims2),2*nhgl,1.,ld (1,1,ims2),4*nllev)
          CALL dgemm ('N','N',2*nllev,ins2,nhgl, 1.,fzl2 (1,1,j,1),2*nllev*nlmp1,pnmt(1,ims2),2*nhgl,0.,lvo(1,1,ims2),4*nllev)
          CALL dgemm ('N','N',2*nllev,ins2,nhgl,-1.,fzm2 (1,1,j,1),2*nllev*nlmp1,anmt(1,ims2),2*nhgl,1.,lvo(1,1,ims2),4*nllev)
        END IF

        IF (nllevp1 > 0 .AND. ins2 > 0) THEN
          CALL dgemm ('N','N',2*nllevp1,ins2,nhgl, 1.,ftp12(1,1,j,1),2*nllevp1*nlmp1,pnmt(1,ims2),2*nhgl,0.,ltp(1,1,ims2),4*nllevp1)
        ENDIF

      END DO ! j = 1, nlmp1

    END IF

  END DO

END SUBROUTINE ltd
