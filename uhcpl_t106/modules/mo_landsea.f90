MODULE mo_landsea

       !mo_landsea.f90
       !xiouhua Fu (sep 2000/fu++ change to T106 on dec 2007)

       USE mo_control,     ONLY:ngl, nlon
       USE mo_doctor,      ONLY:nout,nerr 

       INTEGER, ALLOCATABLE, SAVE :: aland(:,:)	
       INTEGER, SAVE :: itopo(720,121)
       INTEGER, ALLOCATABLE, SAVE :: iland(:,:)	
       INTEGER, ALLOCATABLE, TARGET, SAVE :: bzone(:,:)	
       REAL, ALLOCATABLE, SAVE :: agx(:)
       REAL, ALLOCATABLE, SAVE :: agy(:)
       REAL, SAVE :: sgx(720)
       REAL, SAVE :: sgy(121)
       REAL, ALLOCATABLE, TARGET, SAVE :: sstatm(:,:)
       CONTAINS
	
       subroutine amask(slmm,vlat,vlon)
!c===get the land-sea mask on the atmospheric grid
!c     and Land-Sea mask for ocean model and bzone
!c--------------------------------------------------------
!c  Atmospheric grid (j=1) is from North Pole
!c--------------------------------------------------------
	parameter (imo=720,jmo=121)
!c	integer*4 imask
!c	common /csbc1/aland(ima,jma),itopo(imo,jmo),
!c     1         iland(ima,jma),agx(ima),agy(jma),sgx(imo),sgy(jmo)
        common /fulat/jbs,jbe

!fu++	integer aslm(nlon,ngl)
        integer ieas(ngl),iwes(ngl)
        REAL vlat(ngl),vlon(nlon)
	dimension isea(nlon,ngl),hlp(nlon,jmo),slmm(:,:), &
          topo(imo,jmo),jland(nlon,ngl)
	REAL*4 atopo(imo,jmo)
        REAL*4 aslmm(nlon,ngl)
 
	character*110  filein

	data dxo/0.5/,dyo/0.5/
	data ibs/2/
!	data jbs/34/,jbe/15/,ibs/2/

	data nind,npac,natl/1,1,1/ !active three ocean

!	vlon(1)=0.0
!	do i=2,nlon
!	vlon(i)=vlon(i-1)+3.75
!	end do

	ALLOCATE(aland(1:nlon,1:ngl))
	ALLOCATE(iland(1:nlon,1:ngl))
	ALLOCATE(bzone(1:nlon,1:ngl))
	ALLOCATE(agx(1:nlon))
	ALLOCATE(agy(1:ngl))

        !get the latitude & longitude values
	agx = vlon
        agy = vlat
        !
        !try to find jbs and jbe from vlat(ngl)
	do j=1,ngl
	if(vlat(j).lt.30.0) then
        jbe=j-2
	write(nout,*)'jbe=,vlat=',jbe,vlat(jbe)
        goto 61 	
	end if
	end do
61      continue

        do j=1,ngl
	if(vlat(j).lt.-30.0) then
        jbs=j+1
	write(nout,*)'jbs=,vlat=',jbs,vlat(jbs)
	goto 62 
        end if
        end do 
62      continue

	do i=1,nlon
	do j=1,ngl
        aslmm(i,j)=slmm(i,j)
	end do
	end do
!====================================================!
	open(330,file='slmm.dat',&
           &form='unformatted',access='sequential')
	write(330) ((aslmm(i,j),i=1,nlon),j=1,ngl)
	close(330)
!====================================================!

	
	do i=1,nlon
           do j=1,ngl
           
	   if(slmm(i,j).ge.1.) then       !over land
	   aland(i,j)=0
	   iland(i,j)=0
	   else
	   aland(i,j)=1
	   iland(i,j)=1
	   end if
           jland(i,j)=0
	   end do
	end do
!c===adjust the original 'aland' to a modified version
!c======for T30 version of AGCM
!c-------Pursian Gulf
!	iland(11,18)=0
!	iland(11,19)=0
!	iland(12,20)=0
!	iland(13,21)=0
!	iland(14,21)=0
!c-------Marinetime continent
!	iland(28,22)=0
!	iland(29,25)=0
!	iland(30,25)=0
!	iland(33,23)=0
!	iland(34,23)=0
!	iland(34,22)=0
!	iland(34,21)=0
!	iland(34,20)=0
!c--------Australian
!	iland(34,26)=0
!	iland(35,25)=0
!	iland(35,26)=0
!	iland(36,26)=0
!	iland(36,27)=0
!	iland(37,26)=0
!	iland(37,27)=0
!	iland(38,27)=0
!	iland(39,27)=0
!	iland(38,28)=0
!	iland(39,28)=0
!	iland(38,29)=0
!c--------Central American
!	iland(71,18)=0
!	iland(72,19)=0
!	iland(74,20)=0
!	iland(75,21)=0
!	iland(76,21)=0
!	iland(76,22)=0
!c--------end of adjustment-------------	
!c===view the structure of 'iland'
!
!	call view(iland,nlon,ngl,1,nlon,1,ngl,1)
!
!c===get 'bzone' for atmospheric model
!c===get modified active ocean domain (nind,npac,natl)
	if(min(nind,npac,natl).eq.1) then
	do j=1,ngl
	 do i=1,nlon
	 jland(i,j)=iland(i,j)
	 end do
	end do
	goto 777
	end if
!c=====active Indian ocean ONLY
       if(nind.eq.1) then
	do j=52,109
	do i=18,108
	jland(i,j)=iland(i,j)
	end do
	end do
        end if
!c===============================================c
!c	do j=27,33
!c	do i=29,37
!c	jland(i,j)=iland(i,j)
!c	end do
!c	end do
!c       end if
!c	iland(75,22)=0
!c	iland(76,22)=0
!c====active Pacific ocean ONLY
	if(npac.eq.1) then
	do j=52,109
	is=204
	do i=is,268
	if(iland(i,j).eq.1.and.iland(i+1,j).ne.1) then
	ieas(j)=i
	goto 33
	end if
        end do
33      end do

	do j=52,109
	do i=204,ieas(j)
	jland(i,j)=iland(i,j)
	end do
	end do

	do j=52,109
	do i=108,204
	jland(i,j)=iland(i,j)
	end do
	end do
	end if
!c======active Atlantic ocean ONLY
!c	if(natl.eq.1) then
!c	do j=16,33
!c	is=52
!c	do i=is,82
!c	if(iland(i,j).eq.1.and.iland(i+1,j).ne.1) then
!c	ieas(j)=i+1
!c	goto 34
!c	end if
!c	end do
!c34      end do
!c
!c	do j=16,33
!c	do i=ieas(j),102
!c	ii=i
!c	if(i.gt.96) ii=i-96
!c	jland(ii,j)=iland(ii,j)
!c	end do
!c	end do
!c	end if
!c
!c===get the active ocean domain on atmospheric grid
777	do i=1,nlon
           do j=1,ngl
	   bzone(i,j)=2
	   isea(i,j)=0
           if(j.ge.(jbe+2).and.j.le.(jbs-2).and.jland(i,j).gt.0) then
	   bzone(i,j)=1
	   isea(i,j)=1
	   end if
	   end do
	end do
!c--get the transient zone between the active ocean and 'land mass'
!c--primarily over the land in this case

	do i=1,nlon
	   do j=jbe+1,jbs-1
	   ia=i-1
	   ib=i+1
	   ja=j-1
	   jb=j+1
	   if(i.eq.1) ia=nlon 
	   if(i.eq.nlon) ib=1
	   
	   if(bzone(i,j).eq.2) then
	   if(bzone(ia,j).eq.1.or.bzone(ib,j).eq.1 &
         .or.bzone(i,ja).eq.1.or.bzone(i,jb).eq.1)then
!c     1    .or.bzone(ia,ja).eq.1.or.bzone(ia,jb).eq.1
!c     1    .or.bzone(ib,ja).eq.1.or.bzone(ib,jb).eq.1)then
	   bzone(i,j)=0
	   end if
	   end if

	   end do
	end do	
 
	   	   
!	call view(bzone,nlon,ngl,1,nlon,1,ngl,1)
!
!c===get 'topo' for intermediate ocean model
!c==linearly interpolating isea(nlon,ngl) to (imo,jmo)
!	agx(1)=0.
!	do i=2,nlon
!	agx(i)=agx(i-1)+dx
!	end do
!c======================================================c

	sgx(1)=0.
	do i=2,imo
	sgx(i)=sgx(i-1)+dxo
	end do

	sgy(1)=-30.
	do j=2,jmo
	sgy(j)=sgy(j-1)+dyo
	end do
!c=======================================================c
	write(nout,*)'agx=',agx
	write(nout,*)'sgx=',sgx
	write(nout,*)'sgy=',sgy
!c=======================================================c
!c===>interplating in y direction========================c
	do j=1,jmo
	   do i=1,nlon
!c--find the atmospheric grids in y direction============c
            jb1=jbs
	   do jj=jb1,jbe,-1
	   if(sgy(j).ge.agy(jj).and.sgy(j).lt.agy(jj-1)) then
	   jb1=jj
	   aa1=float(isea(i,jj))
	   aa2=float(isea(i,jj-1))
	   hlp(i,j)=aa1+(aa2-aa1)*(sgy(j)-agy(jj))/(agy(jj-1)-agy(jj))
	   goto 78
	   end if
 	   end do
	
78	   end do
	end do

!c===>interplating in x direction<========================c
	do i=1,imo
	   do j=1,jmo
!c--find the atmospheric grids in x direction<============c
	   if(i.eq.1) then
	   topo(1,j)=hlp(1,j)      !start from 0
	   end if

	   do ii=ibs,nlon
	   ia=ii-1
	   ib=ii

	   if(sgx(i).lt.agx(1)) then
	   aa1=hlp(nlon,j)
	   aa2=hlp(1,j)
	   ax=-360.+agx(nlon)
	   topo(i,j)=aa1+(aa2-aa1)*(sgx(i)-ax)/(agx(1)-ax)
	   goto 77
	   end if
	  
	   if(sgx(i).gt.agx(nlon)) then
	   aa1=hlp(nlon,j)
	   aa2=hlp(1,j)
	   ax=360.+agx(1)
	   topo(i,j)=aa1+(aa2-aa1)*(sgx(i)-agx(nlon))/(ax-agx(nlon))
	   goto 77
	   end if

	   if(sgx(i).gt.agx(ia).and.sgx(i).le.agx(ib)) then
	   ibs=ib
	
	   aa1=hlp(ia,j)
	   aa2=hlp(ib,j)
	   topo(i,j)=aa1+(aa2-aa1)*(sgx(i)-agx(ia))/(agx(ib)-agx(ia))
	   goto 77
	   end if
 	   end do
	
77	end do
	end do

!c=====>view topo<============================
	do i=1,imo
	do j=1,jmo
	itopo(i,j)=0
	if(topo(i,j).ge.0.5) itopo(i,j)=1
	end do
        end do
!c------------------------------------------------------>
!c------------->ajustment of itopo----------------------->
	do j=86,91
	do i=241,244
	itopo(i,j)=0
	end do
	end do

	do j=81,86
	do i=243,248
	itopo(i,j)=0
	end do
	end do

	do j=74,81
	do i=247,253
	itopo(i,j)=0
	end do
	end do

	do j=68,73
	do i=252,258
	itopo(i,j)=0
	end do
	end do

	do j=62,67
	do i=257,263
	itopo(i,j)=0
	end do
	end do

	do j=56,61
	do i=262,268
	itopo(i,j)=0
	end do
	end do

	do j=50,55
	do i=267,273
	itopo(i,j)=0
	end do
	end do
!c========================================================>
!c              Australian & Newzland
!c
	do j=40,49
	do i=285,289
	itopo(i,j)=0
	end do
	end do

!c=======================================================>
!c    close Central America
!c
        do j=78,79
        do i=559,566
        itopo(i,j)=0
        end do
        end do
!c====================end of adjustment===================>

	do i=1,imo
	do j=1,jmo
	atopo(i,j)=float(itopo(i,j))
	end do
	end do

!c===================================================>
!c===save 'itopo'
        open(10,file='itop.dat',form='unformatted')
	write(10) itopo
	close(10)
        open(10,file='itop.look',form='unformatted',&
        access='sequential')
	write(10) ((atopo(i,j),i=1,imo),j=1,jmo)
	close(10)
!c======================================================================
!c
!c	call viewtop(itopo,1,imo,jmo,1,-1)

        return	
	end subroutine amask

	subroutine view(fland,nlon,ngl,is,ie,js,je,jv)
	integer fland(nlon,ngl)
	integer*1 gland(nlon,ngl)

	!change length
	gland = fland
        !
	im=nlon/2
	write(nout,110)
110     format(/, 1x,'Atmospheric land/sea mask follows:'/)

	do j=js,je,jv
	write(nout,'(i2,1x,96i1)')j,gland(is:ie,j)
	end do

!c======save for print-out
	write(12,111)
111     format(//////, 1x,'Atmospheric land/sea mask follows:'////)

	write(12,121)((gland(i,j),i=is,im),j=js,je,jv)
121     format(48i1)
	write(12,111)
	write(12,123)((gland(i,j),i=im,ie),j=js,je,jv)
123     format(49i1)


	return
	end subroutine view
  
	subroutine viewtop(top,is,ie,js,je,jv)
	integer top
	parameter(ima=720,jma=121)
	dimension top(ima,jma)
        integer*1 atop(ima,jma)

	!change length
	atop=top
        !
	write(nout,110)
110     format(/, 1x,'Oceanic land/sea mask follows:'/)

	write(nout,120)((atop(i,j),i=is,ie,2),j=js,je,jv)
120     format(1x,90i1)

!c======save for print-out
	write(13,111)
111     format(//////, 1x,'Oceanic land/sea mask follows:'////)

	write(13,121)((atop(i,j),i=is,60),j=js,je,jv)
121     format(1x,60i1)

	write(14,111)
	write(14,121)((atop(i,j),i=61,120),j=js,je,jv)

	write(15,111)
	write(15,121)((atop(i,j),i=121,ie),j=js,je,jv)

	return
	end subroutine viewtop

	END MODULE mo_landsea
