MODULE mo_transpose
!
!+ $Id: mo_transpose.f90,v 1.37 2002/12/19 09:51:02 m214003 Exp $
!
! this module holds the global distribution and transposition routines:
!   global distribution routines (global <-> pe)
!   transposition routines (pe <-> pe)
!
! Authors:
!
! A. Rhodin, MPI, August 1999, original source
!

!
! declare explicit shape arguments in routine (un)pack_fs_buf
! on vector machines
!
#if ((defined CRAY)&&(! defined _CRAYMPP)) || (defined SX) || (defined __uxp__)
#define EXPLICIT
#else
#undef EXPLICIT
#endif

  USE mo_exception,     ONLY: finish          ! abort in case of errors
  USE mo_doctor,        ONLY: nerr            ! stderr output unit
  USE mo_mpi,           ONLY: p_send,        &! MPI_send     routine
                              p_recv,        &! MPI_recv     routine
                              p_sendrecv,    &! MPI_sendrecv routine
                              p_bcast,       &! MPI_bcast    routine
                              p_pe,          &! this processor
                              p_nprocs,      &! number of processors
                              p_io            ! processor which performs I/O
  USE mo_decomposition, ONLY: pe_decomposed, &! decomposition table data type
                              debug_parallel,&! debug flag
                              dc=>local_decomposition

  IMPLICIT NONE
  PRIVATE
  !
  ! public routines:
  !
  PUBLIC :: indx       ! get index within decomposition table from processor id
  !
  !   scatter : global arrays -> local arrays
  !
  PUBLIC :: scatter_gp ! global grid point field -> local pe's
  PUBLIC :: scatter_ls ! global spectral   field -> Legendre space
  PUBLIC :: scatter_sp ! global spectral   field -> local pe's
  PUBLIC :: scatter_sa ! global Fourier sym/asym -> local pe's
  !
  !   gather  : global arrays <- local arrays
  !
  PUBLIC :: gather_gp  ! global grid point field <- local pe's 
  PUBLIC :: gather_gp3 ! global grid point field <- local pe's (nlon,nlev,nlat)
  PUBLIC :: gather_ls  ! global spectral   field <- Legendre space
  PUBLIC :: gather_sp  ! global spectral   field <- local pe's
  PUBLIC :: gather_sa  ! global Fourier sym/asym <- local pe's
  !
  ! transpositions:
  !
  PUBLIC :: tr_gp_fs   ! transpose grid point space <-> Fourier  space
  PUBLIC :: tr_fs_ls   ! transpose Fourier    space <-> Legendre space
  PUBLIC :: tr_ls_sp   ! transpose Legendre   space <-> spectral space
  !
  ! public constants
  !
  PUBLIC :: tag_gather_gp

  !
  ! interfaces (specific routines)
  !
  !  Gridpoint space
  !
  INTERFACE gather_gp
    MODULE PROCEDURE gather_gp432 ! gather gridp. field (nlon,nlev,ntrac,nlat)
                                  !                  or (nlon,nlev,nlat,1)
                                  !                  or (nlon,nlat,1,1)
    MODULE PROCEDURE gather_gp32  ! gather gridp. field (nlon,nlev,nlat)
                                  !                 or  (nlon,nlat,1)
    MODULE PROCEDURE gather_gp2   ! gather gridp. field (nlon,nlat)
  END INTERFACE

  INTERFACE scatter_gp
    MODULE PROCEDURE scatter_gp432! scatter gridp. field (nlon,nlev,ntrac,nlat)
                                  !                   or (nlon,nlev,nlat,1)
                                  !                   or (nlon,nlat,1,1)
    MODULE PROCEDURE scatter_gp32 ! scatter gridp. field (nlon,nlev,nlat)
                                  !                   or (nlon,nlat,1)
    MODULE PROCEDURE scatter_gp2  ! scatter gridp. field (nlon,nlat)
  END INTERFACE
  !
  !   Legendre space
  !
  INTERFACE scatter_ls
    MODULE PROCEDURE scatter_ls3  ! scatter full spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE scatter_ls0  ! scatter only m=0 wave number (nlev,   nnp1)
  END INTERFACE

  INTERFACE gather_ls
    MODULE PROCEDURE gather_ls3   ! gather full spectral field   (nlev, 2, nsp)
    MODULE PROCEDURE gather_ls0   ! gather only m=0 wave number  (nlev,   nnp1)
    MODULE PROCEDURE gather_ls4   !   spectral fields (nmp1*2,nlevp1,nlat,nvar)
  END INTERFACE
  !
  ! symetric and asymetric forier components
  ! decomposition in legendre space
  !
  INTERFACE scatter_sa
    MODULE PROCEDURE scatter_sa42 ! sym/asym components (nlev,2,nmp1,nhgl) 
                                  !                  or (nlev,nhgl,1,1)
    MODULE PROCEDURE scatter_sa2  ! sym/asym components (nlev,nhgl)
  END INTERFACE

  INTERFACE gather_sa
    MODULE PROCEDURE gather_sa42  ! sym/asym components (nlev,2,nmp1,nhgl) 
                                  !                  or (nlev,nhgl,1,1)
    MODULE PROCEDURE gather_sa2   ! sym/asym components (nlev,nhgl)
  END INTERFACE
  !
  ! spectral space
  !
  INTERFACE scatter_sp
    MODULE PROCEDURE scatter_sp4  ! scatter spectral field  (nlev, 2, nsp,1)
    MODULE PROCEDURE scatter_sp3  ! scatter spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE scatter_sp0  ! scatter only m=0 coeff. (nlev,   nnp1)
  END INTERFACE

  INTERFACE gather_sp
    MODULE PROCEDURE gather_sp4 ! scatter full spectral field  (nlev, 2, nsp,1)
    MODULE PROCEDURE gather_sp3 ! scatter full spectral field  (nlev, 2, nsp)
    MODULE PROCEDURE gather_sp0 ! scatter only m=0 wave number (nlev,   nnp1)
  END INTERFACE

  INTERFACE indx
    MODULE PROCEDURE indx0
    MODULE PROCEDURE indx2
  END INTERFACE
  !
  ! define tags
  !
  INTEGER, PARAMETER :: tag_scatter_gp   = 100
  INTEGER, PARAMETER :: tag_gather_gp    = 101
  INTEGER, PARAMETER :: tag_scatter_ls   = 110
  INTEGER, PARAMETER :: tag_gather_ls    = 111
  INTEGER, PARAMETER :: tag_scatter_sp   = 120
  INTEGER, PARAMETER :: tag_gather_sp    = 121
  INTEGER, PARAMETER :: tag_scatter_sa   = 130
  INTEGER, PARAMETER :: tag_gather_sa    = 131
  INTEGER, PARAMETER :: tag_tr_gp_fs     = 200
  INTEGER, PARAMETER :: tag_tr_fs_ls     = 210
  INTEGER, PARAMETER :: tag_tr_ls_sp     = 220
!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE scatter_gp432 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: i
    REAL    ,POINTER :: gl3(:,:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3th index
    !
    IF (p_pe == p_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    IF (size4 == 1) THEN
      gl3 => gl(:,:,:,1)
      CALL scatter_gp32 (gl3, lc(:,:,:,1), gl_dc)
    ELSE
      DO i=1,SIZE(lc,3)
        gl3 => gl(:,:,i,:)
        CALL scatter_gp32 (gl3, lc(:,:,i,:), gl_dc)
      END DO
    ENDIF
  END SUBROUTINE scatter_gp432
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp32 (gl, lc, gl_dc)
  !
  ! send global grid point field to pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
    INTEGER :: size3 ! size of 3rd dimension
    REAL    ,POINTER :: gl2(:,:)
    !
    ! call 2D scatter routine if 3rd dimension size is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) size3 = (SIZE(gl,3))
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      NULLIFY(gl2)
      IF (p_pe == p_io) gl2 => gl(:,:,1)
      CALL scatter_gp2 (gl2, lc(:,:,1), gl_dc)
    ELSE
      CALL scatter_gp3 (gl, lc, gl_dc)
    ENDIF
  END SUBROUTINE scatter_gp32
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp2 (gl, lc, gl_dc)
  !
  ! send global 2D grid point field to local pe's (nlon,nlat)
  !
  REAL                 ,POINTER     :: gl   (:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,1:nglat/2 ) = GL (glons(1):glone(1),glats(1):glate(1))
    !   LC (1:nglon,nglat/2+1:) = GL (glons(2):glone(2),glats(2):glate(2))
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:) ! buffer
    INTEGER              :: imype     ! index of this PE
    INTEGER              :: i         ! loop index
    INTEGER              :: pe        ! processor to communicate with
    INTEGER              :: nlon      ! global number of longitudes
    INTEGER              :: nglon     ! local
    imype = indx (p_pe, gl_dc)
    nlon = gl_dc(imype)% nlon
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:gl_dc(i)% nglh(1)) =                    &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,gl_dc(i)% nglat/2+1:) =                  &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
            !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:) = &
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,             &
                  gl_dc(i)% nglat/2+1:) =                          &
              gl (1: gl_dc(i)% glone(2),                           &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
          nglon = gl_dc(i)% nglon
          lc(       :nglon,:) = buf
          lc(nglon+1:     ,:) = 0.
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      nglon = gl_dc(imype)% nglon
      CALL p_recv (lc(:nglon,:), p_io, tag_scatter_gp)
      lc(nglon+1:,:) = 0.
    END IF
  END SUBROUTINE scatter_gp2
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp3 (gl, lc, gl_dc)
  !
  ! send global 3D grid point field to local pe's (nlon,nlev,nlat)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    ! Data structure:
    !
    !   global grid point field GL: (1:nlon (+2), 1:nlev, 1:nlat)
    !
    !   local  grid point field LC: (1:nglon(+2), 1:nlev, 1:nglat)
    !
    !   The actual size of the first index may or may be not larger than 
    !   NLON or NGLON
    !
    !   The local grid point field LC covers two distinct areas at opposite
    !   sides of the globe. The array elements correspond with the global
    !   field as follows:
    !
    !   LC (1:nglon,:,1:nglat/2 ) = GL (glons(1):glone(1),:,glats(1):glate(1))
    !   LC (1:nglon,:,nglat/2+1:) = GL (glons(2):glone(2),:,glats(2):glate(2))
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: nglon       ! local
    INTEGER              :: imype       ! index of this pe
    imype = indx (p_pe, gl_dc)
    nlon = gl_dc(1)% nlon
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        ALLOCATE (buf (gl_dc(i)% nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! pack first segment
        !
        buf(:,:,:gl_dc(i)% nglh(1)) =                   &
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1), : , &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1))
        !
        ! pack second segment
        !
        IF (gl_dc(i)% nglh(2)>0) THEN
          IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
            buf(:,:,gl_dc(i)% nglh(1)+1:) =                 &
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2), : , &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ELSE
            !
            ! pack second segment, split in longitudes
            !
            buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglh(1)+1:) = &
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
            buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,&
                  gl_dc(i)% nglh(1)+1:) = &
              gl (1: gl_dc(i)% glone(2), : ,          &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2))
          ENDIF
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_gp)
        ELSE
          nglon = gl_dc(i)% nglon
          lc(       :nglon,:,:) = buf
          lc(nglon+1:     ,:,:) = 0.
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      nglon = gl_dc(imype)% nglon
      CALL p_recv (lc(:nglon,:,:), p_io, tag_scatter_gp)
      lc(nglon+1:,:,:) = 0.
    END IF
  END SUBROUTINE scatter_gp3
!==============================================================================
  SUBROUTINE gather_gp432 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   ) 
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: i
    REAL    ,POINTER :: gl3(:,:,:)
    !
    ! call 2D gather routine if 4th dimension size is 1
    ! else loop over 3th index
    !
    IF (p_pe == p_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    IF (size4 == 1) THEN
      gl3 => gl(:,:,:,1)
      CALL gather_gp32 (gl3, lc(:,:,:,1), gl_dc, source)
    ELSE
      DO i=1,SIZE(lc,3)
        gl3 => gl(:,:,i,:)
        CALL gather_gp32 (gl3, lc(:,:,i,:), gl_dc, source)
      END DO
    ENDIF
  END SUBROUTINE gather_gp432
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp32 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlon,nlev,nlat) or (nlon,nlat,1)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size3 ! size of 3rd dimension
    REAL    ,POINTER :: gl2(:,:)
    !
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
    IF (p_pe == p_io) size3 = (SIZE(gl,3))
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      gl2 => gl(:,:,1)
      CALL gather_gp2 (gl2, lc(:,:,1), gl_dc, source)
    ELSE
      CALL gather_gp3 (gl, lc, gl_dc, source)
    ENDIF
  END SUBROUTINE gather_gp32
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp2 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nlon,nlat)
  !
  REAL                 ,POINTER     :: gl   (:,:)   ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:)   ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)   ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: nglon       ! local number of longitudes
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    imype = indx (p_pe, gl_dc)
    nlon  = gl_dc(1)% nlon
    !
    ! send if pe /= p_io
    !
    IF (p_pe /= p_io) THEN
      nglon = gl_dc(imype)% nglon
      CALL p_send (lc(:nglon,:), p_io, tag_gather_gp)
    ELSE
      DO i = 1, p_nprocs
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lc(:nglon,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon,                       &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2),                        &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
      !
      ! set elements with i>nlon to zero
      !
      gl (nlon+1:,:) = 0.
    ENDIF
  END SUBROUTINE gather_gp2
!------------------------------------------------------------------------------
  SUBROUTINE gather_gp3 (gl, lc, gl_dc, source)
  !
  ! receive global grid point field from local pe's (nlon,nlev,nlat)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    ! Data structure: As described in scatter_gp
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:) ! buffer
    INTEGER              :: i           ! loop index
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlon        ! global number of longitudes
    INTEGER              :: nglon       ! local number of longitudes
    INTEGER              :: imype       ! index of this pe
    INTEGER              :: src         ! source 
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    imype = indx (p_pe, gl_dc)
    nlon  = gl_dc(1)% nlon
    !
    ! send if pe /= p_io
    !
    IF (p_pe /= p_io) THEN
      nglon = gl_dc(imype)% nglon
      CALL p_send (lc(:nglon,:,:), p_io, tag_gather_gp)
    ELSE
      DO i = 1, p_nprocs
        pe    = gl_dc(i)% pe
        nglon = gl_dc(i)% nglon
        ALLOCATE (buf (nglon, SIZE(gl,2), gl_dc(i)% nglat))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_gp)
        ELSE
          buf = lc(:nglon,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack first segment
          !
          gl (gl_dc(i)% glons(1) : gl_dc(i)% glone(1),:,   &
              gl_dc(i)% glats(1) : gl_dc(i)% glate(1)) = &
            buf(:,:,:gl_dc(i)% nglh(1))
          !
          ! unpack second segment
          !
          IF (gl_dc(i)% nglh(2)>0) THEN
            IF (gl_dc(i)% glone(2)>gl_dc(i)% glons(2)) THEN
              gl (gl_dc(i)% glons(2) : gl_dc(i)% glone(2),:,   &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(:,:,gl_dc(i)% nglat/2+1:)
            ELSE
              !
              ! unpack second segment, split in longitudes
              !
              gl (gl_dc(i)% glons(2) : nlon, : ,                     &
                  gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) =       &
                buf(:nlon-gl_dc(i)% glons(2)+1,:,gl_dc(i)% nglat/2+1:)
              gl (1: gl_dc(i)% glone(2), : ,                      &
                     gl_dc(i)% glats(2) : gl_dc(i)% glate(2)) = &
                buf(gl_dc(i)% nglon-gl_dc(i)% glone(2)+1:,:,      &
                    gl_dc(i)% nglat/2+1:)
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
      !
      ! set elements with i>nlon to zero
      !
      gl (nlon+1:,:,:) = 0.
    ENDIF
  END SUBROUTINE gather_gp3
!==============================================================================
  SUBROUTINE scatter_ls3 (gl, lc, gl_dc)
  !
  ! send global spectral field to Legendre space (nlev, 2, nspec)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev    [+1]   ,2,1:nspec )
    !    local Legendre space:  LC (1:nllev | nllevp1,2,1:lnsp)
    !  
    !    The levels 1:NLLEVP1 in LC correspond to LLEVS:LLEVE in GL
    !
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LC covers NLM longitudinal wave numbers m=LM(i), i=1,NLM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,:,nlmp(i)+n+1)  
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:)  ! buffer (lev,2,nsp)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, nllevp1, lnsp, nlm
    INTEGER              :: ke, nk
    INTEGER              :: mp, np       ! displacement, no n waves per m
    INTEGER              :: mpgl         ! displacement on global processor
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1 = gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs   = gl_dc(i)% llevs      ! start level
        lleve   = gl_dc(i)% lleve      ! end level
        nlm     = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        lnsp    = gl_dc(i)% lnsp       ! number of spectral (complex) coeff.
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1     ! number of levels to send
        IF (nk * lnsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, 2, lnsp))
        !
        ! pack
        !
        mp=0
        DO im=1,nlm
          mp1  = gl_dc(i)% lm(im)+1
          np   = gl_dc(i)% nnp (mp1)
          mpgl = gl_dc(i)% nmp (mp1)
          buf (     :  ,:,mp  +1:mp  +np  ) = &
            gl(llevs:ke,:,mpgl+1:mpgl+np)
          mp = mp + np
        END DO
        IF (mp /= lnsp) THEN
          WRITE(nerr,*) 'scatter_ls: PE',pe,',mp/=lnsp:',mp,lnsp
          CALL finish('scatter_ls','mp/=lnsp')
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_ls)
        ELSE
          lc(:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:,:), p_io, tag_scatter_ls)
    END IF
  END SUBROUTINE scatter_ls3
!------------------------------------------------------------------------------
  SUBROUTINE scatter_ls0 (gl, lc, gl_dc)
  !
  ! send global spectral field to Legendre space (m=0 only) (nlev,nnp1)
  !
  REAL                 ,POINTER     :: gl   (:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition

    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev    [+1]   ,1:nnp1 )
    !    local Legendre space:  LC (1:nllev | nllevp1,1:nlnm0)
    !  
    !    The levels 1:NLLEVP1 in LC correspond to LLEVS:LLEVE in GL
    !
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,nmp(m+1)+n+1,:)
    !
    !    LC covers NLM longitudinal wave numbers m=LM(i), i=1,NLM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,nlmp(i)+n+1,:)  
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)    ! buffer (lev,nnp1)
    INTEGER              :: i            ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: llevs, lleve, nllevp1
    INTEGER              :: ke, nk       ! actuall last level, number of levs
    INTEGER              :: nlnm0        ! number of coefficiens for m=0 
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1  = gl_dc(i)% nllevp1  ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of coefficiens for m=0
        IF (nlnm0 == 0) CYCLE
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        ALLOCATE (buf (nk, nlnm0))
        !
        ! pack
        !
        buf (:,:) = gl(llevs:ke,:)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_ls)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:), p_io, tag_scatter_ls)
    END IF
  END SUBROUTINE scatter_ls0
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls3 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (nlev, 2, nsp)
  !
  REAL                 ,POINTER     :: gl    (:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc    (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source        ! -1=all;0=p_io;1=not p_io
    !
    ! Data structures: as described in scatter_ls
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:)  ! buffer (lev,2,nsp)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, nllevp1, lnsp, nlm
    INTEGER              :: ke, nk
    INTEGER              :: mp, np       ! displacement, no n waves per m
    INTEGER              :: mpgl         ! displacement on global processor
    INTEGER              :: src          ! local copy of 'source'
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1  = gl_dc(i)% nllevp1  ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        lnsp   = gl_dc(i)% lnsp       ! number of spectral (complex) coeff.
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, 2, lnsp))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          mp=0
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
            np   = gl_dc(i)% nnp (mp1)
            mpgl = gl_dc(i)% nmp (mp1)
            gl    (llevs:ke,:,mpgl+1:mpgl+np) = &
              buf (     :  ,:,mp  +1:mp  +np  )
            mp = mp + np
          END DO
          IF (mp /= lnsp) THEN
            WRITE(nerr,*) 'gather_ls: PE',pe,',mp/=lnsp:',mp,lnsp
            CALL finish('gather_ls','mp/=lnsp')
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls3
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls0 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (m=0 only) (nlev,nnp1)
  !
  REAL                 ,POINTER     :: gl    (:,:) ! global field
  REAL                 ,INTENT(in)  :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source      ! -1=all;0=p_io;1=not p_io
    !
    ! Data structures: as described in scatter_ls
    !
    REAL    ,ALLOCATABLE :: buf (:,:)    ! buffer (lev,2,nsp)
    INTEGER              :: i            ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: llevs, lleve, nllevp1, nlnm0
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1 = gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs   = gl_dc(i)% llevs      ! start level
        lleve   = gl_dc(i)% lleve      ! end level
        nlnm0   = gl_dc(i)% nlnm0      ! number of coefficients for m=0        
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1            ! number of levels to send
        IF (nk*nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk, nlnm0))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          gl    (llevs:ke,:) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) THEN
        CALL p_send (lc(:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls0
!------------------------------------------------------------------------------
  SUBROUTINE gather_ls4 (gl, lc, gl_dc, source)
  !
  ! receive global spectral field from Legendre space (nmp1*2,nlev,nlat,nvar)
  !
  REAL                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source          !
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlmp1*2,nlev,nlat,nvar)
    INTEGER              :: i, im          ! loop indices
    INTEGER              :: pe             ! processor to communicate with
    INTEGER              :: mp1            ! wavenumber m+1
    INTEGER              :: llevs, lleve, nllevp1, nlm, nlat, nvar
    INTEGER              :: ke, nk
    INTEGER              :: src            ! local copy of 'source'
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_ls
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1  = gl_dc(i)% nllevp1  ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m vavenumbers handled by pe
        nlat   = gl_dc(i)% nlat
        nvar   = SIZE(gl,4)
        ke = MIN (lleve,SIZE(gl,2))
        nk = ke - llevs + 1           ! number of levels to send
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlm*2,nk,nlat,nvar))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_ls)
        ELSE
          buf = lc(:,:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
            gl    (mp1*2-1,llevs:ke,:,:) = buf (im*2-1,:,:,:)
            gl    (mp1*2  ,llevs:ke,:,:) = buf (im*2  ,:,:,:)
          END DO
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:,:), p_io, tag_gather_ls)
      END IF
    END IF
  END SUBROUTINE gather_ls4
!==============================================================================
  SUBROUTINE scatter_sp4 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: size3 ! size of 3rd dimension
    REAL    ,POINTER :: gl3(:,:,:)
    REAL    ,POINTER :: gl2(:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3rd index, call 2D scatter routine if 3rd dim size is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) THEN
      size4 = (SIZE(gl,4))
      IF (size4/=1) CALL finish('scatter_sp4','size4/=1')
      size3 = (SIZE(gl,3))
    ENDIF
    CALL p_bcast (size3, p_io)
    IF (size3==1) THEN
      gl2 => gl(:,:,1,1)
      CALL scatter_sp0 (gl2, lc(:,:,1,1), gl_dc)
    ELSE
      gl3 => gl(:,:,:,1)
      CALL scatter_sp3 (gl3, lc(:,:,:,1), gl_dc)
    ENDIF
  END SUBROUTINE scatter_sp4
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sp3 (gl, lc, gl_dc)
  !
  ! send global spectral field to spectral space (nlev,2,nsp)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev ,2,1:nspec )
    !    local Legendre space:  LC (1:nlev ,2,1:snsp)
    !  
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LS covers NSM longitudinal wave numbers m=SM(i), i=1,NSM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LS(:,:,nsmp(i)+n+1-snn0(i))  
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:) ! buffer (lev,2,nsp)
    INTEGER              :: i, im       ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: mp1         ! wavenumber m+1
    INTEGER              :: snsp        ! number of spectral coefficients on PE
    INTEGER              :: nsm         ! number of m wavenumbers on pe
    INTEGER              :: mp, np      ! displacement, no n waves per m
    INTEGER              :: mpgl        ! displacement on global processor
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nsm     = gl_dc(i)% nsm        ! number of m vavenumbers handled by pe
        snsp    = gl_dc(i)% snsp       ! number of spectral (complex) coeff.
        nlev    = SIZE(gl,1)
        IF (snsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlev, 2, snsp))
        !
        ! pack
        !
        mp=0
        DO im=1,nsm
          mp1  = gl_dc(i)% sm(im)+1
          np   = gl_dc(i)% snnp(im)
          mpgl = gl_dc(i)% nmp (mp1) + gl_dc(i)% snn0(im)
          buf (:,:,mp  +1:mp  +np  ) = &
            gl(:,:,mpgl+1:mpgl+np)
          mp = mp + np
        END DO
        IF (mp /= snsp) THEN
          WRITE(nerr,*) 'scatter_sp: PE',pe,',mp/=snsp:',mp,snsp
          CALL finish('scatter_sp','mp/=snsp')
        ENDIF
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sp)
        ELSE
          lc(:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:,:), p_io, tag_scatter_sp)
    END IF
  END SUBROUTINE scatter_sp3
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sp0 (gl, lc, gl_dc)
  !
  ! send global spectral field to spectral space (m=0 only) (nlev, nnp1)
  !
  REAL                 ,POINTER     :: gl   (:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
    !
    !  Data structures:
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)   ! buffer (lev,nnp1)
    INTEGER              :: i           ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nsnm0       ! number of coefficiens for m=0 on PE 
    INTEGER              :: snn0        ! number of first n coefficient on PE
    INTEGER              :: nlev        ! number of levels
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nlev   = SIZE(gl,1)           ! number of levels
        nsnm0  = gl_dc(i)% nsnm0      ! number of coefficiens for m=0 on PE
        IF (nsnm0 == 0) CYCLE
        snn0   = gl_dc(i)% snn0(1)    ! number of first n coefficient on PE
        !
        ALLOCATE (buf (nlev, nsnm0))
        !
        ! pack
        !
        buf (:,:) = gl(:,1+snn0:nsnm0+snn0)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sp)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_recv (lc(:,:), p_io, tag_scatter_sp)
    END IF
  END SUBROUTINE scatter_sp0
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp4 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)    ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    INTEGER :: size4 ! size of 4rd dimension
    INTEGER :: size3 ! size of 3rd dimension
    REAL    ,POINTER :: gl3(:,:,:)
    REAL    ,POINTER :: gl2(:,:)
    !
    ! abort if 4th dimension size is 1
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
    IF (p_pe == p_io) THEN
      size4 = (SIZE(gl,4))
      IF (size4/=1) CALL finish('gather_sp4','size4/=1')
      size3 = (SIZE(gl,3))
    ENDIF
    CALL p_bcast (size3, p_io)
    IF (size3==1) THEN
      gl2 => gl(:,:,1,1)
      CALL gather_sp0 (gl2, lc(:,:,1,1), gl_dc, source)
    ELSE
      gl3 => gl(:,:,:,1)
      CALL gather_sp3 (gl3, lc(:,:,:,1), gl_dc, source)
    ENDIF
  END SUBROUTINE gather_sp4
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp3 (gl, lc, gl_dc, source)
  !
  ! gather global spectral field from spectral space (nlev,2,nsp)
  !
  REAL                 ,POINTER     :: gl   (:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)     ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    !
    !  Data structures:
    !
    !    global spectral field: GL (1:nlev ,2,1:nspec )
    !    local Legendre space:  LC (1:nlev ,2,1:snsp)
    !  
    !    GL covers the longitudinal wave numbers m=0..nm
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    GL(:,:,nmp(m+1)+n+1)
    !
    !    LC covers NSM longitudinal wave numbers m=SM(i), i=1,NSM.
    !    The coefficients corresponding to wave number (n,m) are stored in
    !    LC(:,:,nsmp(i)+n+1-snn0(i))  
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:) ! buffer (lev,2,nsp)
    INTEGER              :: i, im       ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: mp1         ! wavenumber m+1
    INTEGER              :: snsp        ! number of spectral coefficients on PE
    INTEGER              :: nsm         ! number of m wavenumbers on pe
    INTEGER              :: mp, np      ! displacement, no n waves per m
    INTEGER              :: mpgl        ! displacement on global processor
    INTEGER              :: src         ! source
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nsm     = gl_dc(i)% nsm        ! number of m vavenumbers handled by pe
        snsp    = gl_dc(i)% snsp       ! number of spectral (complex) coeff.
        nlev    = SIZE(gl,1)
        IF (snsp < 1) CYCLE
        !
        !
        ALLOCATE (buf (nlev, 2, snsp))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sp)
        ELSE
          buf = lc(:,:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack
          !
          mp=0
          DO im=1,nsm
            mp1  = gl_dc(i)% sm(im)+1
            np   = gl_dc(i)% snnp(im)
            mpgl = gl_dc(i)% nmp (mp1) + gl_dc(i)% snn0(im)
            gl    (:,:,mpgl+1:mpgl+np) = &
              buf (:,:,mp  +1:mp  +np  )
            mp = mp + np
          END DO
          IF (mp /= snsp) THEN
            WRITE(nerr,*) 'gather_sp: PE',pe,',mp/=snsp:',mp,snsp
            CALL finish('gather_sp','mp/=snsp')
          ENDIF
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) THEN
        CALL p_send (lc(:,:,:), p_io, tag_gather_sp)
      ENDIF
    END IF
  END SUBROUTINE gather_sp3
!------------------------------------------------------------------------------
  SUBROUTINE gather_sp0 (gl, lc, gl_dc, source)
  !
  ! gather global spectral field from spectral space (m=0 only) (nlev,nnp1)
  !
  REAL                 ,POINTER     :: gl   (:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                               ! -1=all;0=p_io;1=not p_io
    !
    !  Data structures:
    !
    !  local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)   ! buffer (lev,nnp1)
    INTEGER              :: i           ! loop indices
    INTEGER              :: pe          ! processor to communicate with
    INTEGER              :: nsnm0       ! number of coefficiens for m=0 on PE 
    INTEGER              :: snn0        ! number of first n coefficient on PE
    INTEGER              :: nlev        ! number of levels
    INTEGER              :: src         ! source
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nlev   = SIZE(gl,1)           ! number of levels
        nsnm0  = gl_dc(i)% nsnm0      ! number of coefficiens for m=0 on PE
        IF (nsnm0 == 0) CYCLE
        snn0   = gl_dc(i)% snn0(1)    ! number of first n coefficient on PE
        !
        ALLOCATE (buf (nlev, nsnm0))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sp)
        ELSE
          buf = lc(:,:)
        ENDIF
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          !
          ! unpack
          !
          gl(:,1+snn0:nsnm0+snn0) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc)>0) CALL p_send (lc(:,:), p_io, tag_gather_sp)
    END IF
  END SUBROUTINE gather_sp0
!==============================================================================
  SUBROUTINE gather_sa42 (gl, lc, gl_dc, source)
  !
  ! gather global grid point field from pe's (nlev,2,nmp1,nhgl) 
  !                                       or (nlev,nhgl,1,1)
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source         ! source to gather from
    !                                                 !-1=all;0=p_io;1=not p_io
    INTEGER :: size3 ! size of 3rd dimension * size of 4th dimension
    REAL    ,POINTER :: gl2(:,:)
    !
    ! call 2D gather routine if 3rd,4th index is 1
    ! else call 3D gather routine
    !
    IF (p_pe == p_io) size3 = SIZE(gl,3) * SIZE(gl,4)
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL gather_sa2 (gl2, lc(:,:,1,1), gl_dc, source)
    ELSE
      CALL gather_sa4 (gl, lc, gl_dc, source)
    ENDIF
  END SUBROUTINE gather_sa42
!------------------------------------------------------------------------------
  SUBROUTINE gather_sa4 (gl, lc, gl_dc, source)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],2,nmp1,nhgl), nlev,nmp1 split over PEs
  !
  REAL                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL                 ,INTENT(in)  :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source          !
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlev[+1],2,nmp1,nhgl)
    INTEGER              :: i, im          ! loop indices
    INTEGER              :: pe             ! processor to communicate with
    INTEGER              :: mp1            ! wavenumber m+1
    INTEGER              :: llevs, lleve, nllevp1, nlm, nhgl
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_sa
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1= gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m wavenumbers handled by pe
        nhgl   = gl_dc(i)% nlat/2.    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,2,nlm,nhgl))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sa)
        ELSE
          buf = lc(:,:,:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          DO im=1,nlm
            mp1  = gl_dc(i)% lm(im)+1
             gl    (llevs:ke,:,mp1,:) = buf (:,:,im,:)
          END DO
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_send (lc(:,:,:,:), p_io, tag_gather_sa)
      END IF
    END IF
  END SUBROUTINE gather_sa4
!------------------------------------------------------------------------------
  SUBROUTINE gather_sa2 (gl, lc, gl_dc, source)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],nhgl), nlev split over PEs, present only if nlnm0>0
  !
  REAL                 ,POINTER     :: gl    (:,:) ! global field
  REAL                 ,INTENT(in)  :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
  INTEGER, OPTIONAL    ,INTENT(in)  :: source      !
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)  ! buffer (nlev[+1],nhgl)
    INTEGER              :: i          ! loop indices
    INTEGER              :: pe         ! processor to communicate with
    INTEGER              :: llevs, lleve, nllevp1, nlnm0, nhgl
    INTEGER              :: ke, nk
    INTEGER              :: src
    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source
    !
    ! Data structures: as defined in scatter_sa
    !
    ! receive if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1= gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of n wavenumbers for m=0
        nhgl   = gl_dc(i)% nlat/2.    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1 .OR. nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,nhgl))
        !
        ! receive
        !
        IF (pe /= p_pe) THEN
          CALL p_recv( buf, pe, tag_gather_sa)
        ELSE
          buf = lc(:,:)
        ENDIF
        !
        ! unpack
        !
        IF(src==-1 .OR. (src==0 .AND. p_io==pe)      &
                   .OR. (src==1 .AND. p_io/=pe)) THEN
          gl    (llevs:ke,:) = buf (:,:)
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! send if (p_io /= p_pe)
      !
      i = indx (p_pe, gl_dc)
      nlnm0  = gl_dc(i)% nlnm0 ! number of n wavenumbers for m=0
      IF(SIZE(lc,1)>0 .AND. nlnm0>0 ) THEN
        CALL p_send (lc(:,:), p_io, tag_gather_sa)
      END IF
    END IF
  END SUBROUTINE gather_sa2
!==============================================================================
  SUBROUTINE scatter_sa42 (gl, lc, gl_dc)
  !
  ! scatter global grid point field from pe's (nlev,2,nmp1,nhgl) 
  !                                       or (nlev,nhgl,1,1)
  !
  REAL                 ,POINTER     :: gl   (:,:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc   (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc(0:)      ! global decomposition
    INTEGER :: size3 ! size of 3rd dimension * size of 4th dimension
    REAL    ,POINTER :: gl2(:,:)
    !
    ! call 2D scatter routine if 3rd,4th index is 1
    ! else call 3D scatter routine
    !
    IF (p_pe == p_io) size3 = SIZE(gl,3) * SIZE(gl,4)
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
      IF (p_pe == p_io) gl2 => gl(:,:,1,1)
      CALL scatter_sa2 (gl2, lc(:,:,1,1), gl_dc)
    ELSE
      CALL scatter_sa4 (gl, lc, gl_dc)
    ENDIF
  END SUBROUTINE scatter_sa42
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sa4 (gl, lc, gl_dc)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],2,nmp1,nhgl), nlev,nmp1 split over PEs
  !
  REAL                 ,POINTER     :: gl    (:,:,:,:) ! global field
  REAL                 ,INTENT(out) :: lc    (:,:,:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)       ! global decomposition
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:,:,:)  ! buffer (nlev[+1],2,nmp1,nhgl)
    INTEGER              :: i, im        ! loop indices
    INTEGER              :: pe           ! processor to communicate with
    INTEGER              :: mp1          ! wavenumber m+1
    INTEGER              :: llevs, lleve, nllevp1, nlm, nhgl
    INTEGER              :: ke, nk
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1= gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlm    = gl_dc(i)% nlm        ! number of m wavenumbers handled by pe
        nhgl   = gl_dc(i)% nlat/2.    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,2,nlm,nhgl))
        !
        ! pack
        !
        DO im=1,nlm
          mp1  = gl_dc(i)% lm(im)+1
          buf (:,:,im,:) = gl    (llevs:ke,:,mp1,:)
        END DO
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sa)
        ELSE
          lc(:,:,:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      IF(SIZE(lc,1)>0) THEN
        CALL p_recv (lc(:,:,:,:), p_io, tag_scatter_sa)
      END IF
    END IF
  END SUBROUTINE scatter_sa4
!------------------------------------------------------------------------------
  SUBROUTINE scatter_sa2 (gl, lc, gl_dc)
  !
  ! receive global fourier (symetric/asymetric part) from Legendre space 
  !   (nlev[+1],nhgl), nlev split over PEs, present only if nlnm0>0
  !
  REAL                 ,POINTER     :: gl    (:,:) ! global field
  REAL                 ,INTENT(out) :: lc    (:,:) ! local  field
  TYPE (pe_decomposed) ,INTENT(in)  :: gl_dc (:)   ! global decomposition
    !
    ! local variables
    !
    REAL    ,ALLOCATABLE :: buf (:,:)  ! buffer (nlev[+1],nhgl)
    INTEGER              :: i          ! loop indices
    INTEGER              :: pe         ! processor to communicate with
    INTEGER              :: llevs, lleve, nllevp1, nlnm0, nhgl
    INTEGER              :: ke, nk
    !
    ! send if pe = p_io
    !
    IF (p_pe == p_io) THEN
      DO i = 1, p_nprocs
        pe = gl_dc(i)% pe
        !
        ! short names
        !
        nllevp1= gl_dc(i)% nllevp1    ! number of levels handled by pe
        llevs  = gl_dc(i)% llevs      ! start level
        lleve  = gl_dc(i)% lleve      ! end level
        nlnm0  = gl_dc(i)% nlnm0      ! number of n wavenumbers for m=0
        nhgl   = gl_dc(i)% nlat/2.    ! half number of gaussian latitudes
        ke = MIN (lleve,SIZE(gl,1))
        nk = ke - llevs + 1           ! number of levels to receive
        IF (nk < 1 .OR. nlnm0 < 1) CYCLE
        !
        !
        ALLOCATE (buf (nk,nhgl))
        !
        ! pack
        !
        buf (:,:) = gl (llevs:ke,:)
        !
        ! send
        !
        IF (pe /= p_pe) THEN
          CALL p_send( buf, pe, tag_scatter_sa)
        ELSE
          lc(:,:) = buf
        ENDIF
        DEALLOCATE (buf)
      END DO
    ELSE
      !
      ! receive if (p_io /= p_pe)
      !
      i = indx (p_pe, gl_dc)
      nlnm0  = gl_dc(i)% nlnm0 ! number of n wavenumbers for m=0
      IF(SIZE(lc,1)>0 .AND. nlnm0>0 ) THEN
        CALL p_recv (lc(:,:), p_io, tag_scatter_sa)
      END IF
    END IF
  END SUBROUTINE scatter_sa2
!==============================================================================
  SUBROUTINE tr_gp_fs (gl_dc, sign, gp1, gp2, gp3, gp4, gp5, gp6, gp7,&
                       sf1, sf2, sf3, zm1, zm2, zm3, fs, fs0)
  !
  ! transpose
  !   sign= 1 : grid point space  -> Fourier space
  !   sign=-1 : grid point space <-  Fourier space
  ! 
  !
  TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
  INTEGER              ,INTENT(in)     :: sign             ! 1:gp>fs; -1:gp<fs
  REAL                 ,INTENT(inout)  :: gp1    (:,:,:)   ! gridpoint space 3d
  REAL                 ,INTENT(inout)  :: gp2    (:,:,:)   ! 
  REAL                 ,INTENT(inout)  :: gp3    (:,:,:)   ! 
  REAL                 ,INTENT(inout)  :: gp4    (:,:,:)   ! 
  REAL                 ,INTENT(inout)  :: gp5    (:,:,:)   ! 
  REAL                 ,INTENT(inout)  :: gp6    (:,:,:)   ! 
  REAL                 ,INTENT(inout)  :: gp7    (:,:,:)   ! 
  REAL ,OPTIONAL       ,INTENT(inout)  :: sf1    (:,:)     ! gridpoint space 2d
  REAL ,OPTIONAL       ,INTENT(inout)  :: sf2    (:,:)     ! gridpoint space 2d
  REAL ,OPTIONAL       ,INTENT(inout)  :: sf3    (:,:)     ! gridpoint space 2d
  REAL ,OPTIONAL       ,INTENT(inout)  :: zm1    (:,:)     ! zonal mean
  REAL ,OPTIONAL       ,INTENT(inout)  :: zm2    (:,:)     ! zonal mean
  REAL ,OPTIONAL       ,INTENT(inout)  :: zm3    (:,:)     ! zonal mean
  REAL                 ,INTENT(inout)  :: fs     (:,:,:,:) ! Fourier space
  REAL ,OPTIONAL       ,INTENT(inout)  :: fs0    (:,:,:)   ! zonal mean, Four.
    !
    ! Data structures:
    !
    ! local variables
    ! 
    INTEGER             :: i, j, k, l
    INTEGER             :: imype         ! index of this pe
    REAL   ,POINTER     :: buf (:,:,:,:) ! send buffer
    REAL   ,POINTER     :: bu0 (:,:,:)   ! send buffer zonal means
    REAL   ,POINTER     :: bufr(:,:,:,:) ! receive buffer
    REAL   ,POINTER     :: bu0r(:,:,:)   ! receive buffer zonal means
    LOGICAL             :: m0            ! receive zonal means
    LOGICAL             :: m0s           ! send zonal means
    INTEGER             :: seta          ! my set A
    INTEGER             :: nprocb        ! number of PEs in set A
    INTEGER,ALLOCATABLE :: idx_com (:)   ! PEs to communicate with 
    INTEGER             :: ks, ke, nk    ! vertical range of buffer
    INTEGER             :: nk0           ! vertical range of buffer bu0
    INTEGER             :: nglat, nglon  ! gridpoint space no. lons,lats
    INTEGER             :: glons(2)      ! first longitudes in gridspace
    INTEGER             :: glone(2)      ! last  longitudes in gridspace
    INTEGER             :: nglh(2)       ! number of lats in each domains
    INTEGER             :: nlon          ! global number of longitudes
    INTEGER             :: nvar          ! number of variables in fft buffer
    INTEGER             :: nva0          ! number of variables (zonal mean)
    INTEGER             :: nle0          ! number of levels (zonal mean)
    !
    ! loop 1: send if (dest_pe > src_pe); else recv
    !
    nvar   = SIZE (fs,4)
    imype  = indx (p_pe, gl_dc)
    seta   = gl_dc(imype)% set_a
    nprocb = gl_dc(imype)% nprocb
    IF (gl_dc(imype)% col_1d) RETURN
    ALLOCATE (idx_com (0:nprocb-1))
    CALL plan_b_sr (idx_com)
    nva0  = 0; IF (PRESENT(fs0)) nva0 = SIZE (fs0,3)
    nle0  = 0; IF (PRESENT(fs0)) nle0 = SIZE (fs0,1)

    i = imype         ! PE index (my)
    SELECT CASE (sign)
    CASE (1)
      !
      ! grid point -> Fourier
      !
      DO k = 0, nprocb-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
        CALL alloc_gp_fs_buf (gp_i=i, fs_i=j, buf=buf, bu0=bu0)
        CALL pack_gp_buf
        m0s = m0
        IF(i/=j) CALL alloc_gp_fs_buf (gp_i=l, fs_i=i, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_fs
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
      !
      ! zero latitudes > nlat
      !
      fs (gl_dc(imype)% nlon+1:,:,:,:) = 0.
    CASE (-1)
      !
      ! Fourier -> grid point
      !
      DO k = 0, nprocb-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
        CALL alloc_gp_fs_buf (gp_i=j, fs_i=i, buf=buf, bu0=bu0)
        CALL pack_fs_buf
        m0s = m0
        IF(i/=j) CALL alloc_gp_fs_buf (gp_i=i, fs_i=l, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_gp
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
    CASE default
      CALL finish ('tr_fs_buf','invalid SIGN parameter (not 1,-1)')
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE (idx_com)
  CONTAINS
!------------------------------------------------------------------------------
    SUBROUTINE alloc_gp_fs_buf (gp_i, fs_i, buf, bu0)
    INTEGER       :: gp_i, fs_i ! indices of grid point and Fourier space pe
    REAL ,POINTER :: buf(:,:,:,:), bu0(:,:,:)
    !
    ! derive bounds and allocate buffer
    !
      ks    = gl_dc(fs_i)% flevs
      ke    = gl_dc(fs_i)% fleve
      nk    = gl_dc(fs_i)% nflevp1
      nk0   = gl_dc(fs_i)% nflev
      nglat = gl_dc(gp_i)% nglat
      nglon = gl_dc(gp_i)% nglon
      glons = gl_dc(gp_i)% glons
      glone = gl_dc(gp_i)% glone
      nlon  = gl_dc(gp_i)% nlon
      nglh  = gl_dc(gp_i)% nglh
      m0    = PRESENT(fs0).AND.(fs_i==gp_i.OR.sign==-1)
      ALLOCATE (buf (nglon, nk, nglat, nvar) )
      IF (m0) ALLOCATE (bu0 (nk0, nglat, nva0))
    END SUBROUTINE alloc_gp_fs_buf
!------------------------------------------------------------------------------
    SUBROUTINE pack_gp_buf
    !
    ! pack message to send/recv buffer buf
    !
      !
      ! pack 2d arrays
      !
      IF (ke == gl_dc(imype)% nlev+1) THEN
        buf (:,nk,:,:) = 0.
        IF(PRESENT(sf1)) buf (:,nk,:,1) = sf1(:nglon,:)
        IF(PRESENT(sf2)) buf (:,nk,:,2) = sf2(:nglon,:)
        IF(PRESENT(sf3)) buf (:,nk,:,3) = sf3(:nglon,:)
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      ! pack 3d arrays
      !
      IF(nk > 0) THEN
        buf (:,:nk,:,1) = gp1 (:nglon,ks:ke,:)
        buf (:,:nk,:,2) = gp2 (:nglon,ks:ke,:)
        buf (:,:nk,:,3) = gp3 (:nglon,ks:ke,:)
        buf (:,:nk,:,4) = gp4 (:nglon,ks:ke,:)
        buf (:,:nk,:,5) = gp5 (:nglon,ks:ke,:)
        buf (:,:nk,:,6) = gp6 (:nglon,ks:ke,:)
        buf (:,:nk,:,7) = gp7 (:nglon,ks:ke,:)
      ENDIF
      !
      ! pack zonal mean
      !
      IF (m0) THEN
        IF(PRESENT(zm1)) bu0 (:,:,1) = zm1 (ks:ke,:)
        IF(PRESENT(zm2)) bu0 (:,:,2) = zm2 (ks:ke,:)
        IF(PRESENT(zm3)) bu0 (:,:,3) = zm3 (ks:ke,:)
      ENDIF
    END SUBROUTINE pack_gp_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_gp
    !
    ! unpack grid point space from send/recv buffer buf
    !
      !
      ! 2d arrays
      !
      IF (ke == gl_dc(imype)% nlev+1) THEN
        IF(PRESENT(sf1)) THEN
          sf1(:nglon,:) = bufr (:,nk,:,1); sf1(nglon+1:,:)=0.
        ENDIF
        IF(PRESENT(sf2)) THEN
          sf2(:nglon,:) = bufr (:,nk,:,2); sf2(nglon+1:,:)=0.
        ENDIF
        IF(PRESENT(sf3)) THEN
          sf3(:nglon,:) = bufr (:,nk,:,3); sf3(nglon+1:,:)=0.
        ENDIF
        ke = ke - 1
        nk = nk - 1
      ENDIF
      !
      ! unpack 3d arrays
      !
      IF(nk > 0) THEN
        gp1 (:nglon,ks:ke,:) = bufr (:,:nk,:,1); gp1 (nglon+1:,ks:ke,:)=0. 
        gp2 (:nglon,ks:ke,:) = bufr (:,:nk,:,2); gp2 (nglon+1:,ks:ke,:)=0.
        gp3 (:nglon,ks:ke,:) = bufr (:,:nk,:,3); gp3 (nglon+1:,ks:ke,:)=0.
        gp4 (:nglon,ks:ke,:) = bufr (:,:nk,:,4); gp4 (nglon+1:,ks:ke,:)=0.
        gp5 (:nglon,ks:ke,:) = bufr (:,:nk,:,5); gp5 (nglon+1:,ks:ke,:)=0.
        gp6 (:nglon,ks:ke,:) = bufr (:,:nk,:,6); gp6 (nglon+1:,ks:ke,:)=0.
        gp7 (:nglon,ks:ke,:) = bufr (:,:nk,:,7); gp7 (nglon+1:,ks:ke,:)=0.
      ENDIF
      !
      ! unpack zonal mean
      !
      IF (m0) THEN
        IF(PRESENT(zm1)) zm1 (ks:ke,:) = bu0r (:,:,1)
        IF(PRESENT(zm2)) zm2 (ks:ke,:) = bu0r (:,:,2)
        IF(PRESENT(zm3)) zm3 (ks:ke,:) = bu0r (:,:,3)
      ENDIF
    END SUBROUTINE unpack_buf_gp
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
    !
    ! unpack message to fourier buffer fs
    !
      !
      ! unpack first segment
      !
      fs(glons(1):glone(1),:,:nglat/2,:) = bufr(:,:,:nglat/2,:)
      !
      ! unpack second segment
      !
      IF (glone(2)>glons(2)) THEN
        fs(glons(2):glone(2),:,nglat/2+1:,:) = bufr(:,:,nglat/2+1:,:)
      ELSE
        !
        ! unpack second segment, split into longitudes
        !
        fs    (glons(2)        :nlon           ,:,nglat/2+1:,:) = &
          bufr(                :nlon-glons(2)+1,:,nglat/2+1:,:)

        fs    (1               :glone(2)       ,:,nglat/2+1:,:) = &
          bufr(nglon-glone(2)+1:               ,:,nglat/2+1:,:)
      ENDIF
      !
      ! unpack zonal mean
      !
      IF (m0) fs0 (:,:,:) = bu0r (:,:,:)
    END SUBROUTINE unpack_buf_fs
!------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf
    !
    ! pack fourier buffer fs to buffer
    !
      !
      ! pack first segment
      !
      buf(:,:,:nglh(1),:) = fs(glons(1):glone(1),:,:nglh(1),:)
      !
      ! pack second segment
      !
      IF(nglh(1)>0) THEN
        IF (glone(2)>glons(2)) THEN
          buf(:,:,nglh(1)+1:,:) = fs(glons(2):glone(2),:,nglh(1)+1:,:)
        ELSE
          !
          ! pack second segment, split into longitudes
          !
          buf  (                :nlon-glons(2)+1,:,nglh(1)+1:,:) = &
            fs (glons(2)        :nlon           ,:,nglh(1)+1:,:)
          buf  (nglon-glone(2)+1:               ,:,nglh(1)+1:,:) = &
            fs (1               :glone(2)       ,:,nglh(1)+1:,:)
        ENDIF
      ENDIF
      !
      ! pack zonal mean
      !
      IF (m0) bu0  (:,:,:) = fs0 (:,:,:)
    END SUBROUTINE pack_fs_buf
!------------------------------------------------------------------------------
    SUBROUTINE send_recv
    !
    ! send and receive buffer
    ! deallocate send buffer
    !
      IF(i/=j) THEN
        CALL p_sendrecv (buf,  gl_dc(j)% pe, &
                         bufr, gl_dc(l)% pe, tag_tr_ls_sp)
        IF (m0.AND.m0s) THEN
          CALL p_sendrecv (bu0,  gl_dc(j)% pe, &
                           bu0r, gl_dc(l)% pe, tag_tr_ls_sp)
        ELSE IF (m0s) THEN
          CALL p_send (bu0,  gl_dc(j)% pe, tag_tr_ls_sp)
        ELSE IF (m0) THEN
          CALL p_recv ( bu0r, gl_dc(l)% pe, tag_tr_ls_sp)
        ENDIF
        DEALLOCATE (buf)
        IF (m0s)  DEALLOCATE (bu0)
      ELSE
        bufr => buf
        IF (m0) bu0r => bu0
      ENDIF
    END SUBROUTINE send_recv
!------------------------------------------------------------------------------
    SUBROUTINE plan_b_sr (idx_com)
    !
    ! set up plan of PEs to communicate with
    !
    INTEGER :: idx_com (0:nprocb-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a /= seta) CYCLE
        idx_com (k) = i ! gl_dc(i)% pe
        IF (i == imype) n = k
        k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_b_sr
!------------------------------------------------------------------------------
  END SUBROUTINE tr_gp_fs
!==============================================================================
  SUBROUTINE tr_fs_ls (gl_dc, sign, fs, ls, fs0, ls0)
  !
  ! transpose
  !   sign= 1 : Fourier space  -> Legendre space
  !   sign=-1 : Fourier space <-  Legendre space
  !
  TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc  (:)       ! decomposition
  INTEGER              ,INTENT(in)     :: sign             ! 1:fs>ls; -1:gs<ls
  !
  ! Assumed shape array association:
  !
  REAL                 ,INTENT(inout)  :: fs   (:,:,:,:)   ! fs
  REAL                 ,INTENT(inout)  :: ls   (:,:,:,:)   ! ls
  REAL ,OPTIONAL       ,INTENT(inout)  :: fs0  (:,:,:)     ! fs, zonal means
  REAL ,OPTIONAL       ,INTENT(inout)  :: ls0  (:,:,:)     ! ls, zonal means
  !
  ! Array element sequence association to speed up indexing:
  !
    !
    ! local variables
    !
    INTEGER              :: i, j, k, l   ! loop indices
    INTEGER              :: imype        ! index of this pe
    INTEGER              :: nvar         ! number of variables
    INTEGER              :: nva0         ! number of variables (m=0 only)
    INTEGER              :: nlev         ! number of levels
    INTEGER              :: nlev0        ! number of levels (m=0 only)
    INTEGER              :: nflat        ! number of latitudes (2*nhgl)
    INTEGER              :: flats(2)     ! first latitude     in Fourier space
    INTEGER              :: flate(2)     ! last  latitude     in Fourier space
    INTEGER              :: nlm          ! number of m coeff. in Legendre space
    INTEGER ,POINTER     :: intr(:)      ! index array
    LOGICAL              :: m0           ! coeff. m=0 present in Legendre space
    LOGICAL              :: m0s
    INTEGER              :: setb         ! set B
    INTEGER              :: nproca       ! number of PEs in set A
    INTEGER              :: n2mp1        ! total number of coeff. from lgti 
    REAL    ,POINTER     :: buf(:,:,:,:) ! send buffer
    REAL    ,POINTER     :: bu0  (:,:,:) ! send buffer (zonal mean)
    REAL    ,POINTER     :: bufr(:,:,:,:)! receive buffer
    REAL    ,POINTER     :: bu0r  (:,:,:)! receive buffer (zonal mean)
    INTEGER ,ALLOCATABLE :: idx_com (:)  ! PEs to communicate with
    INTEGER              :: nb,nf,n2     ! explicit shapes
    !
    ! invariant local variables
    !
    nvar   = SIZE (fs,4)
    nva0   = 0; nlev0 = 0
    IF (PRESENT(fs0).AND.PRESENT(ls0)) THEN
      nva0 = SIZE (fs0,3); nlev0 = SIZE (fs0,1)
    ENDIF
    nlev   = SIZE (fs,2)
    imype  = indx (p_pe, gl_dc)
    setb   =  gl_dc(imype)% set_b
    nproca =  gl_dc(imype)% nproca
    n2mp1  = (gl_dc(imype)% nm + 1) * 2
    ALLOCATE (idx_com (0:nproca-1))
    CALL plan_a_sr (idx_com)
    SELECT CASE (sign)
    CASE (1)
      !
      ! Fourier -> Legendre space
      !
      i = imype         ! PE index (my)
      DO k = 0, nproca-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nproca-k,nproca)) ! PE index (recv)
        CALL alloc_fs_ls_buf (fs_i=i, ls_i=j, buf=buf, bu0=bu0)
        m0s = m0
        CALL pack_fs_buf
        IF(i/=j) CALL alloc_fs_ls_buf (fs_i=l, ls_i=i, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_ls
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
    CASE (-1)
      !
      ! Legendre -> Fourier
      !
      i = imype         ! PE index (my)
      DO k = 0, nproca-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nproca-k,nproca)) ! PE index (recv)
        CALL alloc_fs_ls_buf (fs_i=j, ls_i=i, buf=buf, bu0=bu0)
        m0s = m0
        CALL pack_ls_buf
        IF(i/=j) CALL alloc_fs_ls_buf (fs_i=i, ls_i=l, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_fs
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
      fs (n2mp1+1:,:,:,:) = 0. !zero coefficients not provided by inv. Legendre
    CASE default
      CALL finish ('tr_fs_ls','invalid SIGN parameter (not 1,-1)')
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE (idx_com)
  CONTAINS
!------------------------------------------------------------------------------
    SUBROUTINE send_recv
    !
    ! send and receive buffer
    ! deallocate send buffer
    !
      IF(i/=j) THEN
        CALL p_sendrecv (buf,  gl_dc(j)% pe, &
                         bufr, gl_dc(l)% pe, tag_tr_fs_ls)
        IF (m0.AND.m0s) THEN
          CALL p_sendrecv (bu0,  gl_dc(j)% pe, &
                           bu0r, gl_dc(l)% pe, tag_tr_fs_ls)
        ELSE IF (m0s) THEN
          CALL p_send (bu0,  gl_dc(j)% pe, tag_tr_fs_ls)
        ELSE IF (m0) THEN
          CALL p_recv ( bu0r, gl_dc(l)% pe, tag_tr_fs_ls)
        ENDIF
        DEALLOCATE (buf)
        IF (m0s)  DEALLOCATE (bu0)
      ELSE
        bufr => buf
        IF (m0) bu0r => bu0
      ENDIF
    END SUBROUTINE send_recv
!------------------------------------------------------------------------------
    SUBROUTINE alloc_fs_ls_buf (fs_i, ls_i, buf, bu0)
    INTEGER       :: fs_i, ls_i     ! indices of Fourier and Lagendre space pe
    REAL ,POINTER :: buf(:,:,:,:), bu0(:,:,:)
    !
    ! derive bounds and allocate buffer
    !
      nflat =  gl_dc(fs_i)% nflat
      flats =  gl_dc(fs_i)% flats
      flate =  gl_dc(fs_i)% flate
      nlm   =  gl_dc(ls_i)% nlm
      m0    =  gl_dc(ls_i)% nlnm0 > 0 .AND. PRESENT(fs0)
      intr  => gl_dc(ls_i)% intr
      ALLOCATE (buf(2*nlm,nlev,nflat,nvar))
      IF(m0) ALLOCATE (bu0 (nlev0,nflat,nva0))
      !
      ! determine sizes of buffer and fourier space variable
      ! for explicit shape dummy variables (to enforce vectorization)  
      !
      nb = 2*nlm
      nf = SIZE (fs,1)
      n2 = nlev * nflat * nvar
    END SUBROUTINE alloc_fs_ls_buf
!------------------------------------------------------------------------------
    SUBROUTINE pack_fs_buf
#ifdef EXPLICIT
      CALL pack_fs_buf_ex (fs(1,1,1,1), buf(1,1,1,1), intr, nf,nb,n2)
#else
      buf(:,:,:,:) = fs (intr,:,:,:)
#endif
      IF(m0) bu0 = fs0
    END SUBROUTINE pack_fs_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs
#ifdef EXPLICIT
      CALL unpack_buf_fs_ex (bufr(1,1,1,1), fs(1,1,1,1), intr, nb,nf,n2)
#else
      fs (intr,:,:,:) = bufr(:,:,:,:)
#endif
      IF(m0) fs0 = bu0r
    END SUBROUTINE unpack_buf_fs
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
      ls(:,:,flats(1):flate(1),:) = bufr (:,:,         :nflat/2,:)
      ls(:,:,flats(2):flate(2),:) = bufr (:,:,nflat/2+1:       ,:)
      IF(m0) THEN
        ls0 (:,flats(1):flate(1),:) = bu0r (:,         :nflat/2,:)
        ls0 (:,flats(2):flate(2),:) = bu0r (:,nflat/2+1:       ,:)
      ENDIF
    END SUBROUTINE unpack_buf_ls
!------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
      buf (:,:,         :nflat/2,:) = ls(:,:,flats(1):flate(1),:)
      buf (:,:,nflat/2+1:       ,:) = ls(:,:,flats(2):flate(2),:)
      IF(m0) THEN
        bu0 (:,         :nflat/2,:) = ls0 (:,flats(1):flate(1),:)
        bu0 (:,nflat/2+1:       ,:) = ls0 (:,flats(2):flate(2),:)
      ENDIF
    END SUBROUTINE pack_ls_buf
!------------------------------------------------------------------------------
    SUBROUTINE plan_a_sr (idx_com)
    !
    ! set up plan of PEs to communicate with
    !
    INTEGER :: idx_com (0:nproca-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_b /= setb) CYCLE
        idx_com (k) = i ! gl_dc(i)% pe
        IF (i == imype) n = k
        k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_a_sr
!------------------------------------------------------------------------------
  END SUBROUTINE tr_fs_ls
!==============================================================================
  SUBROUTINE tr_ls_sp (gl_dc, sign, ls1, sp1, ls2, sp2, ls3, sp3, ls0, sp0)
  !
  ! transpose
  !   sign= 1 : Legendre space  -> spectral space
  !   sign=-1 : Legendre space <-  spectral space
  !
  TYPE (pe_decomposed) ,INTENT(in)     :: gl_dc (:)     ! decomposition
  INTEGER              ,INTENT(in)     :: sign          ! 1:ls>sp; -1:ls<sp
  REAL                 ,INTENT(inout)  :: ls1   (:,:,:) ! Legendre space 
  REAL                 ,INTENT(inout)  :: sp1   (:,:,:) ! spectral space
  REAL                 ,INTENT(inout)  :: ls2   (:,:,:) ! Legendre space
  REAL                 ,INTENT(inout)  :: sp2   (:,:,:) ! spectral space
  REAL                 ,INTENT(inout)  :: ls3   (:,:,:) ! Legendre space
  REAL                 ,INTENT(inout)  :: sp3   (:,:,:) ! spectral space
  REAL ,OPTIONAL       ,INTENT(inout)  :: ls0   (:,:)   ! Legendre (m=0 only)
  REAL ,OPTIONAL       ,INTENT(inout)  :: sp0   (:,:)   ! spectral (m=0 only)
    !
    ! local variables
    !
    REAL    ,POINTER    :: buf (:,:,:,:) ! send buffer
    REAL    ,POINTER    :: bu0 (:,:)     ! send buffer (m=0 only)
    REAL    ,POINTER    :: bufr (:,:,:,:)! receive buffer
    REAL    ,POINTER    :: bu0r (:,:)    ! receive buffer (m=0 only)
    INTEGER             :: seta          ! this set A
    INTEGER             :: k, i, j, l    ! loop indices
    INTEGER             :: imype         ! decomposition table index of this pe
    INTEGER             :: nllevp1       ! number of levels in Legendre space
    INTEGER             :: llevs         ! first level in Legendre space
    INTEGER             :: lleve         ! last level in Legendre space
    INTEGER             :: nlnm0         ! number of coeff. with m=0 (Legendre)
    INTEGER             :: snsp          ! number of coefficients in sp. space
    INTEGER             :: ssps          ! first coefficients in spectral space
    INTEGER             :: sspe          ! last coefficients in spectral space
    INTEGER             :: nsnm0         ! number of coeff. with m=0 (spectral)
    INTEGER             :: snn0          ! first n for m=0 in spectral space
    LOGICAL             :: m0            ! transpose array with m=0 only (recv)
    LOGICAL             :: m0s           ! transpose array with m=0 only (send)
    INTEGER             :: ke, nk        ! actual last level, number of levels
    INTEGER             :: nprocb        ! number of PEs in set A
    INTEGER,ALLOCATABLE :: idx_com (:)   ! PEs to communicate with

    imype  = indx (p_pe, gl_dc)
    seta   = gl_dc(imype)% set_a
    nprocb = gl_dc(imype)% nprocb
    ALLOCATE (idx_com (0:nprocb-1))
    CALL plan_b_sr (idx_com)
    i = imype         ! PE index (my)
    SELECT CASE (sign)
    CASE (1)
      !
      ! Legendre space -> spectral space
      !
      DO k = 0, nprocb-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
        CALL alloc_ls_sp_buf (ls_i=i, sp_i=j, buf=buf, bu0=bu0)
        CALL pack_ls_buf
        m0s = m0
        IF(i/=j) CALL alloc_ls_sp_buf (ls_i=l, sp_i=i, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_sp
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
    CASE (-1)
      !
      ! Legendre space <- spectral space
      !
      DO k = 0, nprocb-1
        j = idx_com (           k        ) ! PE index (send)
        l = idx_com (MOD(nprocb-k,nprocb)) ! PE index (recv)
        CALL alloc_ls_sp_buf (ls_i=j, sp_i=i, buf=buf, bu0=bu0)
        CALL pack_sp_buf
        m0s = m0
        IF(i/=j) CALL alloc_ls_sp_buf (ls_i=i, sp_i=l, buf=bufr, bu0=bu0r)
        CALL send_recv
        CALL unpack_buf_ls
        DEALLOCATE (bufr)
        IF (m0)  DEALLOCATE (bu0r)
      END DO
    CASE default
      CALL finish ('tr_ls_sp','invalid SIGN parameter (not 1,-1)')
    END SELECT
    !
    ! frt at ECMWF doesnt yet deallocate
    !
    DEALLOCATE (idx_com)
  CONTAINS
!------------------------------------------------------------------------------
    SUBROUTINE alloc_ls_sp_buf (ls_i, sp_i, buf, bu0)
    INTEGER :: ls_i, sp_i              ! Legendre and spectral space pe indices
    REAL ,POINTER :: buf(:,:,:,:), bu0(:,:)
    !
    ! derive bounds and allocate buffer
    !
    nllevp1 = gl_dc(ls_i)% nllevp1! number of levels in Legendre space
    llevs   = gl_dc(ls_i)% llevs  ! first level in Legendre space
    lleve   = gl_dc(ls_i)% lleve  ! last level in Legendre space
    nlnm0   = gl_dc(ls_i)% nlnm0  ! number of coeff. with m=0 in Legendre space
    snsp    = gl_dc(sp_i)% snsp   ! number of coefficients in sp. space
    ssps    = gl_dc(sp_i)% ssps   ! first coefficients in spectral space
    sspe    = gl_dc(sp_i)% sspe   ! last coefficients in spectral space
    nsnm0   = gl_dc(sp_i)% nsnm0  ! number of coeff. with m=0 in spectral space
    m0      = PRESENT(ls0) .AND. PRESENT(sp0) .AND. nlnm0>0 .AND. nsnm0>0

    IF(m0) snn0 = gl_dc(sp_i)% snn0(1)

    ALLOCATE (buf (nllevp1, 2, snsp, 3))

    IF (m0) THEN
      ALLOCATE (bu0 (nllevp1, nsnm0))
    ENDIF

    END SUBROUTINE alloc_ls_sp_buf
!------------------------------------------------------------------------------
    SUBROUTINE pack_ls_buf
      buf (:SIZE(ls1,1),:,:,1)     = ls1 (:,:,ssps:sspe)
      buf (:SIZE(ls2,1),:,:,2)     = ls2 (:,:,ssps:sspe)
      buf (:SIZE(ls3,1),:,:,3)     = ls3 (:,:,ssps:sspe)
      IF (m0) bu0 (:SIZE(ls0,1),:) = ls0 (:,snn0+1:snn0+nsnm0)
    END SUBROUTINE pack_ls_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_sp
      ke = MIN(SIZE(sp1,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp1 (llevs:ke, :, :) = bufr (:nk,:,:,1)
      ke = MIN(SIZE(sp2,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp2 (llevs:ke, :, :) = bufr (:nk,:,:,2)
      ke = MIN(SIZE(sp3,1), lleve); nk = ke-llevs+1
      IF(nk>0) sp3 (llevs:ke, :, :) = bufr (:nk,:,:,3)
      IF (m0) THEN
        ke = MIN(SIZE(sp0,1), lleve); nk = ke-llevs+1
        IF(nk>0) sp0 (llevs:ke,:) = bu0r (:nk,:)
      ENDIF
    END SUBROUTINE unpack_buf_sp
!------------------------------------------------------------------------------
    SUBROUTINE pack_sp_buf
      ke = MIN(SIZE(sp1,1), lleve); nk = MAX(0, ke-llevs+1)
      buf (:nk,:,:,1) = sp1 (llevs:ke, :, :); buf(nk+1:,:,:,1) = 0.
      ke = MIN(SIZE(sp2,1), lleve); nk = MAX(0, ke-llevs+1)
      buf (:nk,:,:,2) = sp2 (llevs:ke, :, :); buf(nk+1:,:,:,2) = 0.
      ke = MIN(SIZE(sp3,1), lleve); nk = MAX(0, ke-llevs+1)
      buf (:nk,:,:,3) = sp3 (llevs:ke, :, :); buf(nk+1:,:,:,3) = 0.
      IF (m0) THEN
        ke = MIN(SIZE(sp0,1), lleve); nk = MAX(0, ke-llevs+1)
        bu0 (:nk,:) = sp0 (llevs:ke,:); bu0 (nk+1:,:) = 0.
      ENDIF
    END SUBROUTINE pack_sp_buf
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_ls
      ls1 (:,:,ssps:sspe) = bufr (:SIZE(ls1,1),:,:,1)
      ls2 (:,:,ssps:sspe) = bufr (:SIZE(ls2,1),:,:,2)
      ls3 (:,:,ssps:sspe) = bufr (:SIZE(ls3,1),:,:,3)
      IF (m0) THEN
        ls0 (:,snn0+1:snn0+nsnm0) = bu0r (:SIZE(ls0,1),:)
      ENDIF
    END SUBROUTINE unpack_buf_ls
!------------------------------------------------------------------------------
    SUBROUTINE send_recv
    !
    ! send and receive buffer
    ! deallocate send buffer
    !
      IF(i/=j) THEN
        CALL p_sendrecv (buf,  gl_dc(j)% pe, &
                         bufr, gl_dc(l)% pe, tag_tr_ls_sp)
        IF (m0.AND.m0s) THEN
          CALL p_sendrecv (bu0,  gl_dc(j)% pe, &
                           bu0r, gl_dc(l)% pe, tag_tr_ls_sp)
        ELSE IF (m0s) THEN
          CALL p_send (bu0,  gl_dc(j)% pe, tag_tr_ls_sp)
        ELSE IF (m0) THEN
          CALL p_recv ( bu0r, gl_dc(l)% pe, tag_tr_ls_sp)
        ENDIF
        DEALLOCATE (buf)
        IF (m0s)  DEALLOCATE (bu0)
      ELSE
        bufr => buf
        IF (m0) bu0r => bu0
      ENDIF
    END SUBROUTINE send_recv
!------------------------------------------------------------------------------
    SUBROUTINE plan_b_sr (idx_com)
    !
    ! set up plan of PEs to communicate with
    !
    INTEGER :: idx_com (0:nprocb-1)
      INTEGER :: i,k,n
      !
      ! get PE identification number
      !
      k = 0
      DO i = dc%spe, dc%epe
        IF (gl_dc(i)% set_a /= seta) CYCLE
        idx_com (k) = i ! gl_dc(i)% pe
        IF (i == imype) n = k
        k = k + 1
      END DO
      idx_com = CSHIFT (idx_com,n)
    END SUBROUTINE plan_b_sr
!------------------------------------------------------------------------------
  END SUBROUTINE tr_ls_sp
!==============================================================================
  FUNCTION indx0 (pe, gl_dc)
  INTEGER              ,INTENT(in) :: pe       ! processor id
  TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:) ! global decomposition
  INTEGER                          :: indx0     ! index
  !
  ! returns the index of a given PE in the global decomposition table
  !
    INTEGER :: i
    DO i = 1, SIZE(gl_dc)
      IF(gl_dc(i)% pe == pe) THEN
        indx0 = i
        RETURN
      ENDIF
    END DO
    WRITE (nerr,*) 'mo_transpose:indx - index not found in decomposition table'
    WRITE (nerr,*) '  reqired:',pe
    WRITE (nerr,*) '  found  :',gl_dc% pe
    CALL finish ('mo_transpose:indx','index not found in decomposition table')
  END FUNCTION indx0
!------------------------------------------------------------------------------
  FUNCTION indx2 (pe, gl_dc)
  INTEGER              ,INTENT(in) :: pe(:,:)       ! processor id
  TYPE (pe_decomposed) ,INTENT(in) :: gl_dc(:)      ! global decomposition
  INTEGER                          :: indx2(SIZE(pe,1),SIZE(pe,2))     ! index
  !
  ! returns the index of a given PE in the global decomposition table
  !
    INTEGER :: i
    indx2 = -1
    DO i = 1, SIZE(gl_dc)
      WHERE(gl_dc(i)% pe == pe)
        indx2 = i
      endwhere
    END DO
    IF(ANY(indx2==-1))THEN
      WRITE(nerr,*)'mo_transpose:indx - index not found in decomposition table'
      WRITE(nerr,*)'  reqired:',PACK(pe,mask=(indx2==-1))
      WRITE(nerr,*)'  found  :',gl_dc% pe
      CALL finish('mo_transpose:indx','index not found in decomposition table')
    END IF
  END FUNCTION indx2
!------------------------------------------------------------------------------
END MODULE mo_transpose
!==============================================================================
#ifdef EXPLICIT
    SUBROUTINE pack_fs_buf_ex (fs, buf, intr, nf,nb,n2)
    !
    ! explicit shape array argument version of pack_fs_buf 
    ! (enforces vectorization)
    !
    REAL    :: fs  (nf,n2)
    REAL    :: buf (nb,n2)
    INTEGER :: intr(nb)
    INTEGER :: nf,nb,n2
      INTEGER :: i
      DO i=1,nb
        buf(i,:) = fs (intr(i),:)
      END DO
    END SUBROUTINE pack_fs_buf_ex
!------------------------------------------------------------------------------
    SUBROUTINE unpack_buf_fs_ex (buf, fs, intr, nb,nf,n2)
    !
    ! explicit shape size array argument version of unpack_buf_fs
    ! (enforces vectorization)
    !
    REAL    :: buf (nb,n2)
    REAL    :: fs  (nf,n2)
    INTEGER :: intr(nb)
    INTEGER :: nf,nb,n2
      INTEGER :: i
      DO i=1,nb
        fs (intr(i),:) = buf(i,:)
      END DO
    END SUBROUTINE unpack_buf_fs_ex
#endif

