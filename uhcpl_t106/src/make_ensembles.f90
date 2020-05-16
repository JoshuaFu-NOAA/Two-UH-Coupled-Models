subroutine make_ensembles

  ! generates ensembles

  use mo_control, only: nens, nlp2,nlon, ngl, nlev
  use mo_kind,    only: dp

  ! list all potential perturbed fields
  use mo_memory_g1a, only: tm1, qm1, vom1, dm1, alpsm1
  use mo_memory_g2a, only: um1, vm1
  use mo_memory_g3a, only: tsm, tsurfm

  ! list the buffer numbers
  use mo_start_dataset, only: nhg1, nhg2, nhg3, nhgl1

  use mo_decomposition, ONLY: dc=>local_decomposition, &
                              dcg=>global_decomposition
  use mo_transpose, only: gather_gp, scatter_gp
  use mo_io,  only: no_ens, restart
  use mo_mpi, only: p_pe, p_io
  use mo_doctor, only: nout

  implicit none

  !fu++ add pointers and targets for global domains
  REAL, POINTER :: um1fu(:,:,:)
  REAL, POINTER :: vm1fu(:,:,:)
  REAL, POINTER :: tm1fu(:,:,:)
  REAL, POINTER :: qm1fu(:,:,:)
  REAL, POINTER :: um1funew(:,:,:)
  REAL, POINTER :: vm1funew(:,:,:)
  REAL, POINTER :: tm1funew(:,:,:)
  REAL, POINTER :: qm1funew(:,:,:)
  !fu++ for perurbations
  REAL, TARGET :: up1(nlp2,nlev,ngl)
  REAL, TARGET :: vp1(nlp2,nlev,ngl)
  REAL, TARGET :: tp1(nlp2,nlev,ngl)
  REAL, TARGET :: qp1(nlp2,nlev,ngl)
  !fu++ for output/input
  REAL*4 :: up2(nlon,ngl,nlev)
  REAL*4 :: vp2(nlon,ngl,nlev)
  REAL*4 :: tp2(nlon,ngl,nlev)
  REAL*4 :: qp2(nlon,ngl,nlev)

  !fu++
  LOGICAL :: Ishift
  INTEGER :: i, i1, j, jr, k,k1
  CHARACTER*3 fid(20)  !maximum 999 assembles 

  external posts1
  data fid/'101','102','103','104','105','106','107',&
           '108','109','110','111','112','113','114',&
           '115','116','117','118','119','120'/


  ! EXCUTABLE FROM HERE
 
  ! generate ensembles
  IF (p_pe == p_io) THEN
     write(nout,*) 'RESTART FILE G1 ',trim(restart(nhg1)%nc_file_name)
     write(nout,*) 'RESTART FILE G2 ',trim(restart(nhg2)%nc_file_name)
     write(nout,*) 'RESTART FILE G3 ',trim(restart(nhg3)%nc_file_name)
     write(nout,*) 'RESTART FILE GL ',trim(restart(nhgl1)%nc_file_name)
  end IF
 
  !fu++   gather g1,g2,g3 from local domains to global domain
  !allocate new global spaces to variables: um1, vm1, tm1,qm1
!	ALLOCATE (um1(dc%nglon,dc%nlev,dc%nglat)
!	ALLOCATE (vm1(dc%nglon,dc%nlev,dc%nglat)
!	ALLOCATE (tm1(dc%nglon,dc%nlev,dc%nglat)
!	ALLOCATE (qm1(dc%nglon,dc%nlev,dc%nglat)

	ALLOCATE (um1fu(nlp2,nlev,ngl))
	ALLOCATE (vm1fu(nlp2,nlev,ngl))
	ALLOCATE (tm1fu(nlp2,nlev,ngl))
	ALLOCATE (qm1fu(nlp2,nlev,ngl))

	ALLOCATE (um1funew(nlp2,nlev,ngl))
	ALLOCATE (vm1funew(nlp2,nlev,ngl))
	ALLOCATE (tm1funew(nlp2,nlev,ngl))
	ALLOCATE (qm1funew(nlp2,nlev,ngl))
!    end if
        
        call gather_gp (um1fu,um1,dcg)
        call gather_gp (vm1fu,vm1,dcg)
        call gather_gp (tm1fu,tm1,dcg)
        call gather_gp (qm1fu,qm1,dcg)
	
        IF (p_pe == p_io) THEN
	do k=1,ngl
	do j=1,nlev
!	if(mod(j,2) == 1) then
!	jr=(j+1)/2
!	else
!	jr=ngl+1-j/2
!	end if 
        do i=1,nlon
	up2(i,k,j)=um1fu(i,j,k)
	vp2(i,k,j)=vm1fu(i,j,k)
	tp2(i,k,j)=tm1fu(i,j,k)
	qp2(i,k,j)=qm1fu(i,j,k)
	end do
	end do
	end do

!==	up2(:,:,:) = um1fu(:,:,:)
!==	vp2(:,:,:) = vm1fu(:,:,:)
!==	tp2(:,:,:) = tm1fu(:,:,:)
!==	qp2(:,:,:) = qm1fu(:,:,:)
!!
!  fu++ don't write out 06/12/2008
!	open(110,file='../DATAIO/model.out',&
!              form='unformatted',access='sequential')
!        do k=1,nlev
!	write(110) ((up2(i,j,k),i=1,nlon),j=1,ngl)
!	end do
!	do k=1,nlev
!	write(110) ((vp2(i,j,k),i=1,nlon),j=1,ngl)
!	end do
!	do k=1,nlev
!	write(110) ((tp2(i,j,k),i=1,nlon),j=1,ngl)
!	end do
!	do k=1,nlev
!	write(110) ((qp2(i,j,k),i=1,nlon),j=1,ngl)
!        end do
!	close(110)	
	END IF
!!fu++ =================================end gathering
!     if(p_pe == p_io)  write(nout,*) 'done 1 write data'

  do no_ens=1,nens
     !read in perturbation fields in 3-D global domain
!        Ishift=.true.
     IF (p_pe == p_io) THEN
	open(110,file='./perturb.'//fid(no_ens),&
              form='unformatted',access='sequential')
        do k=1,nlev
	read(110) ((up2(i,j,k),i=1,nlon),j=1,ngl)
	end do
	do k=1,nlev
	read(110) ((vp2(i,j,k),i=1,nlon),j=1,ngl)
	end do
	do k=1,nlev
	read(110) ((tp2(i,j,k),i=1,nlon),j=1,ngl)
	end do
	do k=1,nlev
	read(110) ((qp2(i,j,k),i=1,nlon),j=1,ngl)
        end do
	close(110)	
     !fu++ add perturbations in 3-D global domain
	do k=1,nlev
	do j=1,ngl
        do i=1,nlp2
	i1=i
	if(i.gt.nlon) i1=nlon
	k1=nlev-k+1
	up1(i,k,j)=up2(i1,j,k1)+um1fu(i,k,j)
	vp1(i,k,j)=vp2(i1,j,k1)+vm1fu(i,k,j)
	tp1(i,k,j)=tp2(i1,j,k1)+tm1fu(i,k,j)
	qp1(i,k,j)=qp2(i1,j,k1)+qm1fu(i,k,j)
	end do
	end do
	end do
	um1funew => up1(:,:,:)
	vm1funew => vp1(:,:,:)
	tm1funew => tp1(:,:,:)
	qm1funew => qp1(:,:,:)
     END IF

     if(p_pe == p_io)  write(nout,*) 'done 2 read new+perturbation'
     !scatter into local domain
!       IF (Ishift) THEN
       call scatter_gp (um1funew, um1,dcg)
       call scatter_gp (vm1funew, vm1,dcg)
       call scatter_gp (tm1funew, tm1,dcg)
       call scatter_gp (qm1funew, qm1,dcg)
!       END IF
          
     if(p_pe == p_io)  write(nout,*) 'done 3 scatter'
     ! write new set of restartfiles
     IF (p_pe == p_io) THEN
        write(nout,*) 'RESTART FILE G1 ',trim(restart(nhg1)%nc_file_name)
        write(nout,*) 'RESTART FILE G2 ',trim(restart(nhg2)%nc_file_name)
        write(nout,*) 'RESTART FILE G3 ',trim(restart(nhg3)%nc_file_name)
        write(nout,*) 'RESTART FILE GL ',trim(restart(nhgl1)%nc_file_name)
     end IF

     call posts1

     if(p_pe == p_io)  write(nout,*) 'done 4 new restart files'
  end do

!	if(p_pe == p_io) then
!	DEALLOCATE (um1fu)
!	DEALLOCATE (vm1fu)
!	DEALLOCATE (tm1fu)
!	DEALLOCATE (qm1fu)
!	DEALLOCATE (um1funew)
!	DEALLOCATE (vm1funew)
!	DEALLOCATE (tm1funew)
!	DEALLOCATE (qm1funew)
!	end if

  no_ens = 0

end subroutine make_ensembles
