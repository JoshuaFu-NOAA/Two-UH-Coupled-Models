	subroutine getsst(sstfu,ssto,ia,ja)
!c=======merge oceanic modeling SST with atmospheric model
!c          SST forcing
!c-----------------------------------------------------------
!c  atmospheric grids are on the Gaussian grids
!c-----------------------------------------------------------
!c  oceanic grids are on the .5x.5 system
!c===========================================================
        USE mo_landsea,    ONLY: bzone, itopo, agx, agy, &
                                 sgx, sgy
        USE mo_doctor,     ONLY: nout
	parameter(io=720,jo=121)
!	integer aland,bzone
!	common /csbc1/aland(ia,ja),itopo(io,jo),
!     1      iland(ia,ja),agx(ia),agy(ja),sgx(io),sgy(jo)
	common /fulat/jbs, jbe

        dimension sstfu(ia,ja),ssto(io,jo)
	dimension hlp(io,ja),help(ia,ja)
	real*4 sta(ia,ja),sto(io,jo)

!	data jbs,jbe/34,15/      !j=1 is Noth Pole
        data ii2,jj2/2,2/
!c===========================================================
!c interpolating ssto(io,jo)=>ssta
!c        in bzone=1 area; 
!c   then smothing ssta in bzone=0 area.
!c===========================================================
!c===>interpolating in y direction
	do j=1,ja
	   do i=1,io
	   hlp(i,j)=0.
	   end do
        end do
!c-------------------------------------------
	do j=jbs,jbe,-1
	   do i=1,io
!c--->find the oceanic grids in y direction
	   do jj=jj2,jo
!fu+	   if(itopo(i,jj-1).gt.0.and.itopo(i,jj).gt.0) then 
!fu+3/25/2002	   if(itopo(i,jj-1).gt.0.or.itopo(i,jj).gt.0) then 
!cccc !!!!IN ocn model domain (change to all regions)
	   if(agy(j).ge.sgy(jj-1).and.agy(j).lt.sgy(jj)) then
!c	   jj2=jj
	   aa1=ssto(i,jj-1)
	   aa2=ssto(i,jj)
	   hlp(i,j)=aa1+(aa2-aa1)*(agy(j)-sgy(jj-1))/(sgy(jj)-sgy(jj-1))
!c-----------avoid zero
	   if(aa1.eq.0.) then
	   hlp(i,j)=aa2
	   goto 551
	   end if
	   if(aa2.eq.0.) then
	   hlp(i,j)=aa1
	   goto 551
	   end if
551        continue
!c------------------------------------
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
	   do ii=ii2,io
	   
	   if(agx(i).gt.sgx(io)) then
	   aa1=hlp(io,j)
	   aa2=hlp(1,j)
	   sx=360.+sgx(1)
	   sstfu(i,j)=aa1+(aa2-aa1)*(agx(i)-sgx(io))/(sx-sgx(io))

!c-----------avoid zero
	   if(aa1.eq.0.) then
	   sstfu(i,j)=aa2
	   end if
	   if(aa2.eq.0.) then
	   sstfu(i,j)=aa1
	   end if
!c------------------------------------
	   goto 77
	   end if

	   ia1=ii-1
	   ib=ii

	   if(agx(i).gt.sgx(ia1).and.agx(i).le.sgx(ib)) then
!c	   ii2=ib
	
	   aa1=hlp(ia1,j)
	   aa2=hlp(ib,j)
	   sstfu(i,j)=aa1+(aa2-aa1)*(agx(i)-sgx(ia1))/(sgx(ib)-sgx(ia1))
!c-----------avoid zero
	   if(aa1.eq.0.) then
	   sstfu(i,j)=aa2
	   end if
	   if(aa2.eq.0.) then
	   sstfu(i,j)=aa1
	   end if
!c------------------------------------
	   goto 77
	   end if
 	   end do
	
77         continue
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
	write(nout,*)'new-step atmospheric inputs at (140W,EQ):'
	write(nout,*)'SST',sstfu(65,25)
!c--------------------------------------------------------------------
	return
	end subroutine getsst

	subroutine smothbz(bzone,t,p,m,n,s)
	real t(m,n),p(m,n)
        integer bzone(m,n)

!c-----------------------------------------------------------------------
	c1=1.-s*s
	c2=s*(1.-s)/2.
	c3=s*s/4.

	do i=1,m
	  do j=1,n
	  p(i,j)=t(i,j)
	  end do
	 end do	
	
	do 10 j=1,n
	do 10 i=1,m
	if(bzone(i,j).eq.0) then
	ia=max(i-1,1)
	ib=min(i+1,m)
	if(i.eq.1) then
	ia=m
	ib=2
	end if
	if(i.eq.m) then
	ia=m-1
	ib=1
	end if
	ja=max(j-1,1)
	jb=min(j+1,n)
	p(i,j)=t(i,j)*c1+c2*(t(ia,j)+t(ib,j)+t(i,ja)&
     		+t(i,jb)-4.*t(i,j))+c3*(t(ia,ja)+t(ia,jb)&
     		+t(ib,ja)+t(ib,jb))
	end if
10	continue

	return
	end subroutine smothbz

