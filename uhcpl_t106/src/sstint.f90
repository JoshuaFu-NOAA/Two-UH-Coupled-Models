subroutine sstint(ia,ja,sstfu)

  !c=======merge oceanic modeling SST with atmospheric model
  !c          SST forcing
  !c-----------------------------------------------------------
  !c  atmospheric grids are on the Gaussian grids
  !c-----------------------------------------------------------
  !c  oceanic grids are on the .5x.5 system
  !c===========================================================
  USE mo_landsea,     ONLY:itopo,bzone,agx,agy,sgx,sgy
  USE mo_doctor,      ONLY:nout

  USE mo_couple,      ONLY: restart_tmp

  parameter(im=720, jm=121)	
  common /var0/ u2(im,jm,2),v2(im,jm,2),h2(im,jm,2),&
       u1(im,jm,2),v1(im,jm,2),h1(im,jm,2),tmp(im,jm,2)
  common /fulat/jbs,jbe       !j=1 is North Pole

  !	parameter(ia=96,ja=48,io=720,jo=121)
  !	integer aland,bzone
  !	common /csbc1/aland(ia,ja),itopo(io,jo),
  !     1      iland(ia,ja),agx(ia),agy(ja),sgx(io),sgy(jo)

  real  sstfu(ia,ja),ssto(im,jm)
  real  hlp(im,ja),help(ia,ja)
  real*4 sta(ia,ja),sto(im,jm)

  data ii2,jj2/2,2/

  !c===========================================================
  !c interpolating ssto(im,jm)=>ssta
  !c        in bzone=1 area; 
  !c   then smothing ssta in bzone=0 area.
  !c===========================================================

  sstfu(:,:)=0.

  ssto(:,:)=tmp(:,:,2)
!ik  ssto(:,:)=restart_tmp(:,:,2)

  !print*,'jbs,jbe=',jbs,jbe
  !print*,'ssto=',(ssto(i,31),i=1,im)
  !print*,'itopo=',(itopo(i,31),i=1,im)
  !print*,'agx=',agx
  !print*,'sgx=',sgx

  !c===>interpolating in y direction

  hlp(:,:)=0.

  !c-------------------------------------------
  do j=jbs,jbe,-1
     do i=1,im
        !c--->find the oceanic grids in y direction
        do jj=jj2,jm
           !fu++      if(itopo(i,jj).gt.0) then  !!!!IN ocn model domain
	   if(agy(j).ge.sgy(jj-1).and.agy(j).lt.sgy(jj)) then
              !c	   jj2=jj
              aa1=ssto(i,jj-1)
              aa2=ssto(i,jj)
              hlp(i,j)=aa1+(aa2-aa1)*(agy(j)-sgy(jj-1))/(sgy(jj)-sgy(jj-1))
	   end if
           !fu++	   end if

        end do
	
     end do
  end do

  !c===>interplating in x direction
  do i=1,ia
     do j=jbs,jbe,-1
        !c--find the oceanic grids in x direction
        if(i.eq.1.and.bzone(1,j).eq.1) then  !in ocean 
	   sstfu(1,j)=hlp(1,j)      !start from 0E in longitude
	   goto 77
        end if
	   
        if(bzone(i,j).ne.1) goto 77

        !c--------------------------------------- 
        do ii=ii2,im
              
           if(agx(i).gt.sgx(im)) then
              aa1=hlp(im,j)
              aa2=hlp(1,j)
              sx=360.+sgx(1)
              sstfu(i,j)=aa1+(aa2-aa1)*(agx(i)-sgx(im))/(sx-sgx(im))
              goto 77
           end if

           ia1=ii-1
           ib=ii
           
           if(agx(i).gt.sgx(ia1).and.agx(i).le.sgx(ib)) then
              !c	   ii2=ib
	
              aa1=hlp(ia1,j)
              aa2=hlp(ib,j)
              sstfu(i,j)=aa1+(aa2-aa1)*(agx(i)-sgx(ia1))/(sgx(ib)-sgx(ia1))
              goto 77
           end if
        end do
	
77      continue
     end do
  end do

  !c------------------------------------------------------------------------
  !c  smothing at bzone=0 area
  !c========================================================================
  !c	do i=1,ia
  !c	   do j=1,ja
  !c	   help(i,j)=ssta(i,j)
  !c	   end do
  !c	end do
  !c	call smothbz(bzone,help,ssta,ia,ja,0.5)
  !c------------------------------------------------------------------------
  !c	do i=1,ia
  !c	do j=1,ja
  !c	sta(i,j)=sstfu(i,j)
  !c	end do
  !c	end do
  !c	do i=1,io
  !c	do j=1,jo
  !c	sto(i,j)=ssto(i,j)
  !c	end do
  !c	end do
  !c	open(115,file='ssta.dat',form='unformatted',
  !c     1    access='sequential')
  !c	write(115) sta
  !c	close(115)
  
  !c	open(116,file='ssto.dat',form='unformatted',
  !c     1    access='sequential')
  !c	write(116) sto
  !c	close(116)
  !c------------------------------------------------------------------c

  WRITE(nout,*)'atmospheric inputs along the equator:'
  WRITE(nout,*)'SST',(sstfu(i,25),i=1,ia,3)

  !c--------------------------------------------------------------------

  return

end subroutine sstint

