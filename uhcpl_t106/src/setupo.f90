!c===================================================================c
!c******************* O C E A N (0-360,30S-30N)+++++++++*************c
!c  set up the grid,top,constant,land-sea mask,coffi,datard,initialo
!c*******************************************************************c

subroutine  setupo(mtt) 

  USE mo_doctor,    ONLY:nout

!qbao
! USE mo_couple,    ONLY: restart_itt
!qbao

  include 'pram1.inc'
  include 'com1.inc'

  !c
  !c------------- coupling with AGCM------------------------------------
  !c	common /scouple/tua(im,jm),tva(im,jm),wfua(im,jm),qsola(im,jm),
  !c     1    hfluxa(im,jm)
  common /fudate/runday,dta,dto,day1,tseg,numcpl,&
       numatm,numocn

  !c-------------------------------------------------------------------c

  data idxcpl/1/

  !cfu

  data cdat/0.5E-3/               !Wu et al(????)

  !c
  !c===================================================================c
  !c===> Set up grids
  !c

!ik not used
!ik  iloop=1        !SST-save-handler
  ncouple=numocn
  ndays=int(runday)
  dt=dto

  call grids2

  !c
  !c===> Define geometry & topogrphy. set land to epsln
  !c
  !c    Read in topography, sept. 23, 1991
  !c

  call topog

  !c===> Get some constants
  !c

  call constant

  write(nout,*) 'tspd=',tspd,'ncouple=',ncouple

  !c
  !c===> Initialize variables
  !c

  call getzero

  !c	call zerot(tmpc)
  !c
  !c===> Get sea-land mask (Rightnow,works only for Pacifc !!)
  !c

  call sea(topo,eqdpth,nbounx,nbouny)

  !c
  !c===> Define coefficients
  !c

  call coffi

  !c
  !c===>Read input data
  !c

  call datard

!qbao
! itt = restart_itt
  call restart_setupo
!qbao

  !c
  !c===>Set initial conditons
  !c	write(nout,*) 'calling initial', ' itt=',itt

!IK  if(mtt.lt.100) then  
  IF (itt == 0) then

     ! initial with zero or use data from UNIT 513
     call initialo       !first month

  else

     call initial(mtt)

  end if

  !c
  !c===> Compute dt based on cfl criterian for "c" grid
  !c     and coriolis parameter
  !c	call getdt
  !c
  !c===>Out put parameters
  !c	call pramout
  !c=====read in atmospheric initial condition for ocean
  !c	call initiala
  !c	call uatouo(tu,au)
  !c	call uatouo(tv,av)
  !c	call uatouo(qsh,sh)    !by ocean itself
  !c	call uatouo(qlh,hlat)
  !c	call uatouo(qfu,aq)
  !c	call uatouo(wfu,wda)

  return

end subroutine setupo


!c***** c m e s h ******************************************************
!c
      subroutine cmesh (dx,apb,amb,wid,n,ndir)
!c
!c  cosine mesh
!c  apb=a+b = coarse resolution
!c  amb=a-b=fine resolution
!c  wid = width between coarse & fine resolution
!c  n = number of variable spaces calculated to span wid
!c  ndir = (1,-1) = spacing (increases,decreases) with increasing index
!c

      USE mo_doctor,        ONLY:nout

      dimension dx(501)
      data pi/3.14159/
      a=.5*(apb+amb)
      b=.5*(apb-amb)
      c=wid/a
      n=c+.1
      rn=1./float(n)
      do 20 ii=1,n
      if (ndir.eq.1) i=n+1-ii
      if (ndir.eq.-1)  i=ii
   20 dx(i)=a+b*cos(pi*rn*(ii-.5))
      write(nout,10) apb,amb,wid,c
   10 format (1x,'MAX RES=',e11.4,'  MIN RES=',e11.4,'  WIDTH=',e11.4,&
      '  NUMBER OF VARIABLE SPACES=',e14.7)
      return
      end subroutine cmesh
!c
!c***** c t o m g p ****************************************************
!c
      subroutine ctomgp (mgp,c,ca,ia)
!c
!c     convert lat,lon or depth coordinate to closest model grid
!c     point i,j,or k.
!c
!c     c  = coordinate to convert....same units as ca
!c     ca = coordinate array corresponding to model grid points
!c     must be monotonically increasing
!c     ia = dimension of ca
!c     mgp= closest model grid point to c
!c

      USE mo_doctor,        ONLY:nout

      dimension ca(ia)
!c
      if (c.lt.ca(1).or.c.gt.ca(ia)) then
      write (nout,9) c,ca(1),ca(ia),ia
      if (c.lt.ca(1))  mgp=1
      if (c.gt.ca(ia)) mgp=ia
9     format (1x,'COORDINATE OUT OF GRID RANGE !!!',/,20x,&
        'COORDINATE =',g14.7,' GRID GOES FROM ',g14.7,' TO ',&
         g14.7,' WITH ',i5,'POINTS')

      return
      else
      do 10 i=2,ia
      if (c.gt.ca(i)) go to 10
      mgp=i
      if (ca(i)-c.gt.c-ca(i-1)) mgp=i-1
      go to 11
10    continue
11    continue
      endif
      return
      end subroutine ctomgp
!c
!c***** e r r **********************************************************
!c
      subroutine err (i,j,a)
      character*8 a
      write(6,10) i,j,a
   10 format (2i10,a8)
      if (i.eq.-12345678) return
      stop
      end subroutine err
!c
!c***** g r i d s 2 ****************************************************
!c
      subroutine grids2
      USE mo_doctor,        ONLY:nout

      parameter (im=720,jm=121,mx=im*jm,nx=1,ny=1,nharms=6)
      common /grid/ xh(im),yh(jm),xz(im),yz(jm),dxh(im),dyh(jm),&
                   dxz(im),dyz(jm),wface,sface
      degcm=111.324e5
!c
!c     define grid spacing in degrees.
!c
      call mesh (dxh,dyh)
      xh(1)=wface+.5*dxh(1)
      yh(1)=sface+.5*dyh(1)
      xz(1)=wface
      yz(1)=sface
      dxz(1)=dxh(1)
      dyz(1)=dyh(1)
      do 200 i=2,im
      dxz(i)=.5*(dxh(i)+dxh(i-1))
      xh(i)=xh(i-1)+dxz(i)
200   xz(i)=xz(i-1)+dxh(i-1)
      do 300 j=2,jm
      dyz(j)=.5*(dyh(j)+dyh(j-1))
      yh(j)=yh(j-1)+dyz(j)
300   yz(j)=yz(j-1)+dyh(j-1)
      write (nout,301) xh
      write (nout,302) yh
301   format (///,' XH(I),I=1,IM',(/,1x,10f9.3))
302   format (///,' YH(J),J=1,JM',(/,1x,10f9.3))
!c
!c     convert dxh dxz dyh dyz to cm.
!c
      do 325 i=1,im
      dxh(i)=dxh(i)*degcm
      dxz(i)=dxz(i)*degcm
325   continue
      do 350 j=1,jm
      dyh(j)=dyh(j)*degcm
      dyz(j)=dyz(j)*degcm
350   continue
      return
      end subroutine grids2
!c
!c***** m e s h ********************************************************
!c
      subroutine mesh (dx,dy)
      USE mo_doctor,        ONLY:nout

      parameter (im=720,jm=121,mx=im*jm,nx=1,ny=1,nharms=6)
      dimension dx(im),dy(jm)
      dimension xmax(nx),xmin(nx),xwid(nx),id(nx)
      dimension ymax(ny),ymin(ny),ywid(ny),jd(ny)
      data xmax/0.5/
      data xmin/0.5/
      data xwid/360.0/
      data ymax/0.5/
      data ymin/0.5/
      data ywid/61.0/
      data id/1/
      data jd/1/
      m=1
      do 5 l=1,nx
      call cmesh (dx(m),xmax(l),xmin(l),xwid(l),n,id(l))
      m=m+n
    5 continue
      m=m-1
!c      write(6,101) m
!c      write (6,1)
!c      write(6,30) dx
      sum=0.
      do 10 i=1,m
   10 sum=sum+dx(i)
      write(nout,20) sum
      m=1
      do 8 l=1,ny
      call cmesh (dy(m),ymax(l),ymin(l),ywid(l),n,jd(l))
      m=m+n
    8 continue
!c      write(6,102) m
      m=m-1
!c      write(6,2)
!c      write(6,30) dy
      sum=0.
      do 15 j=1,m
   15 sum=sum+dy(j)
!c      write(6,20) sum
    1 format (/,10x,'DXT=')
    2 format (/,10x,'DYT=')
   20 format (1x,e14.7)
   30 format (1x,8e14.6)
  101 format(/,10x,'NO OF GRIDSPACES IN X = ',i4)
  102 format(/,10x,'NO OF GRIDSPACES IN Y = ',i4)
      return
      end subroutine mesh
!c
	subroutine sea(topo,eqdpth,nbounx,nbouny)
	parameter(im=720,jm=121)
	real topo(im,jm)
	logical bx1,bx2,by1,by2
	common /var4/ bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)
	integer*1 nbx1(im,jm),nbx2(im,jm),nby1(im,jm),nby2(im,jm)

	do 7 i=1,im
	do 7 j=1,jm
	nbx1(i,j)=0
	nbx2(i,j)=0
	nby1(i,j)=0
	nby2(i,j)=0
7       continue

	do 5 j=1,jm
	do 5 i=1,im
	if(nbounx.eq.1) then
	il1=i-1+im*((im-(i-1))/im)
	ir1=i+1-im*(i/im)
	else
	il1=i-1
	ir1=i+1
	if(i.eq.1) il1=1
	if(i.eq.im) ir1=im
	end if

	if(nbouny.eq.1) then
	jl1=j-1+jm*((jm-(j-1))/jm)
	jr1=j+1-jm*(j/jm)
	else
	jl1=j-1
	jr1=j+1
	if(j.eq.1) jl1=1
	if(j.eq.jm) jr1=jm
	end if

	if(topo(i,j).eq.eqdpth.and.topo(il1,j).ne.eqdpth)then
	bx1(i,j)=.true.
	nbx1(i,j)=1
	end if

	if(topo(i,j).eq.eqdpth.and.topo(ir1,j).ne.eqdpth)then
	bx2(i,j)=.true.
	nbx2(i,j)=1
	end if

	if(topo(i,j).eq.eqdpth.and.topo(i,jl1).ne.eqdpth)then
	by1(i,j)=.true.
	nby1(i,j)=1
	end if

	if(topo(i,j).eq.eqdpth.and.topo(i,jr1).ne.eqdpth)then
	by2(i,j)=.true.
	nby2(i,j)=1
	end if
5	continue

!	open(10,file='sea.index')
!	write(10,*)'bx1'
!	write(10,*)'bx2'
!	write(10,*)'by1'
!	write(10,*)'by2'
!	close(10)

	return
	end subroutine sea


	subroutine getzero
	include 'pram1.inc'
	include 'com1.inc'

        do 10 i=1,ijm
         ui(i)=0.
         vi(i)=0.
         hlni(i)=0.
         u1i(i)=0.
         v1i(i)=0.
         h1lni(i)=0.
         tmpi(i)=0
10     continue

	do 15 i=1,im
	do 15 j=1,jm
	bx1(i,j)=.false.
	bx2(i,j)=.false.
	by1(i,j)=.false.
	by2(i,j)=.false.
15	continue

	do 30 j=1,jm
	do 30 i=1,im
	tmpc(i,j)=0.
	um(i,j)=0.
	vm(i,j)=0.
	hm(i,j)=0.
	tmpm(i,j)=0.
	wem(i,j)=0.
	tem(i,j)=0.
	xahm(i,j)=0.
	yahm(i,j)=0.
	zahm(i,j)=0.
	tfcm(i,j)=0.
	tdfm(i,j)=0.
	tndm(i,j)=0.
	u1m(i,j)=0.
	v1m(i,j)=0.
	h1m(i,j)=0.
	dltsstm(i,j)=0.
	htcm(i,j)=0.
	tfc1m(i,j)=0.
	tfc2m(i,j)=0.
	tfc3m(i,j)=0.
	tfc4m(i,j)=0.
	tfc5m(i,j)=0
	tfc6m(i,j)=0
	help(i,j)=0.
30	continue

        do 31 j=1,jm
        do 31 i=1,im
         b1(i,j)=epsln
         b2(i,j)=epsln
         tse(i,j)=epsln
         we(i,j)=epsln
         hstar(i,j)=epsln
31      continue
	return
	end subroutine getzero


	subroutine coffi
	include 'pram1.inc'
	include 'com1.inc'
!c	logical be(im,jm)
        do 10 j=1,jm
        jm1=j-1
        if (j.eq.1) jm1=j
        do 10 i=1,im
        im1=i-1
        if (i.eq.1)im1=i
        if (be(i,j)) then
         c1(i,j)=1./dxz(i)
         c2(i,j)=1./dyz(j)
        else
         c1(i,j)=0.
         c2(i,j)=0.
        end if
!c
!c     nf = 0 ==> sphere, otherwise f-plane
!c
        if (nf.eq.0) then
         c3(i,j)=cos(yz(j)*degrad)
         c4(i,j)=1./(4.*cos(yh(j)*degrad))
         c5(i,j)=1./(cos(yh(j)*degrad)*dxz(i))
         c6(i,j)=1./dyz(j)
         c7(i,j)=cos(yh(j)*degrad)
         c8(i,j)=1./(cos(yh(j)*degrad)*dxh(i))
         c9(i,j)=1./(cos(yh(j)*degrad)*dyh(j))
         cf(i,j)=2.*omega*sin(yz(j)*degrad)*c3(i,j)
         ff(j)=2.*omega*sin((yz(j)+0.5)*degrad)
        else
         c3(i,j)=cos(theta*degrad)
         c4(i,j)=1./(4.*cos(theta*degrad))
         c5(i,j)=1./(cos(theta*degrad)*dxz(i)) 
         c6(i,j)=1./dyz(j)
         c7(i,j)=cos(theta*degrad)
         c8(i,j)=1./(cos(theta*degrad)*dxh(i))
         c9(i,j)=1./(cos(theta*degrad)*dyh(j))
         cf(i,j)=2.*omega*sin(theta*degrad)*c3(i,j)
         ff(j)=2.*omega*sin(theta*degrad)
        end if
10      continue

!c
!c===> Build area array for integrals
!c
         do 15 j=1,jm
         do 15 i=1,im
!c
!c        nf = 0 ==> sphere, otherwise f-plane
!c
         if (nf.eq.0) then
          area(i,j)=dxh(i)*dyh(j)*cos(yh(j)*degrad)
         else
          area(i,j)=dxh(i)*dyh(j)*cos(theta*degrad)
         end if
         if (topo(i,j).eq.epsln) area(i,j)=0.
15      continue
	return
	end subroutine coffi

	subroutine datard
	include 'pram1.inc'
	include 'com1.inc'
	character*100 filein
        real*4 ctx4(im,jm,nharms),stx4(im,jm,nharms),&
          cty4(im,jm,nharms),sty4(im,jm,nharms),&
          cta4(im,jm,nharms),sta4(im,jm,nharms),&
          csst4(im,jm,nharms),ssst4(im,jm,nharms),&
          ccld4(im,jm,nharms),scld4(im,jm,nharms),&
          cqin4(im,jm,nharms),sqin4(im,jm,nharms),&
          cqout4(im,jm,nharms),sqout4(im,jm,nharms),&
          txm4(im,jm),tym4(im,jm),sstm4(im,jm),qinm4(im,jm),&
          tairm4(im,jm),qoutm4(im,jm),cldm4(im,jm)
        real*4  cqlw4(jm)
!c
!c==============================================================
!c  read in mean+6 harmonics of surface winds,sst,tair,cloud
!c      ,net heat flux,solar radiation,from the climatological
!c      monthly mean of Oberhuber(1988)
!c------------------------------------------------------------- 
!c===>Read observed winds in 2X1 grid

	filein='tx.harmon' !ecmwf
	open(10,file=filein,form='unformatted')
	read(10) txm4,ctx4,stx4
	close(10)
!c
	filein='ty.harmon'
	open(10,file=filein,form='unformatted')
	read(10) tym4,cty4,sty4
	close(10)
!c-------------------------------------------------------------
!c===>Read observed SST in 2X1 grid

	filein='sst.harmon'
	open(10,file=filein,form='unformatted')
	read(10) sstm4,csst4,ssst4
	close(10)
!c-------------------------------------------------------------
!c===>Read observed AIR TEMPERATURE in 2X1 grid

	filein='tair.harmon'
	open(10,file=filein,form='unformatted')
	read(10) tairm4,cta4,sta4
	close(10)
!c-------------------------------------------------------------
!c===>Read observed Cloud-Fraction in 2X1 grid

	filein='cld.harmon'
	open(10,file=filein,form='unformatted')
	read(10) cldm4,ccld4,scld4
	close(10)
!c-------------------------------------------------------------
!c===>Read observed NET HEAT FLUX in 2X1 grid

	filein='net.harmon'
	open(10,file=filein,form='unformatted')
	read(10) qinm4,cqin4,sqin4
	close(10)
!c-------------------------------------------------------------
!c===>Read in heat-flux correction term in 2X1 grid
!c	filein='../DATAIO/flx.mean'
!c	open(10,file=filein,form='unformatted',&
!c        access='sequential')
!c	read(10) ((flx(i,j),i=1,im),j=1,jm) 
!c	close(10)
!c-------------------------------------------------------------
!c===>Read observed SOLAR RADIATION in 2X1 grid

	filein='sol.harmon'
	open(10,file=filein,form='unformatted')
	read(10) qoutm4,cqout4,sqout4
	close(10)

	filein='cqlw.dat'
	open(10,file=filein,form='unformatted')
!         access='sequential')
	read(10) cqlw4
	close(10)

!	rewind 65
!	read(65,*) cqlw4
!c==========================================================

	do j=1,jm
	cqlw(j)=cqlw4(j)	
	do i=1,im
	txm(i,j)=txm4(i,j)	
	tym(i,j)=tym4(i,j)	
	tairm(i,j)=tairm4(i,j)	
	sstm(i,j)=sstm4(i,j)	
	qinm(i,j)=qinm4(i,j)	
	qoutm(i,j)=qoutm4(i,j)	
	cldm(i,j)=cldm4(i,j)	
	  do nhm=1,nharms
	  ctx(i,j,nhm)=ctx4(i,j,nhm)
	  stx(i,j,nhm)=stx4(i,j,nhm)
	  cty(i,j,nhm)=cty4(i,j,nhm)
	  sty(i,j,nhm)=sty4(i,j,nhm)
	  cta(i,j,nhm)=cta4(i,j,nhm)
	  sta(i,j,nhm)=sta4(i,j,nhm)
	  csst(i,j,nhm)=csst4(i,j,nhm)
	  ssst(i,j,nhm)=ssst4(i,j,nhm)
	  cqin(i,j,nhm)=cqin4(i,j,nhm)
	  sqin(i,j,nhm)=sqin4(i,j,nhm)
	  cqout(i,j,nhm)=cqout4(i,j,nhm)
	  sqout(i,j,nhm)=sqout4(i,j,nhm)
	  ccld(i,j,nhm)=ccld(i,j,nhm)
	  scld(i,j,nhm)=scld(i,j,nhm)
	  end do
	end do
	end do
!c
!c===>Convert qinm and qoutm to units of (erg/cm**2/s)
!c
        do j=1,jm
        do i=1,im

        qinm(i,j)=qinm(i,j)*1000.
        qoutm(i,j)=qoutm(i,j)*1000.

	  do k=1,nharms

          cqin(i,j,k)=cqin(i,j,k)*1000.
          sqin(i,j,k)=sqin(i,j,k)*1000.
          cqout(i,j,k)=cqout(i,j,k)*1000.
          sqout(i,j,k)=sqout(i,j,k)*1000.

	  end do
	end do
	end do
!c===============END READING OF OCEANIC FORCING DATA========
	return
	end subroutine datard

	subroutine initialo
        USE mo_doctor,        ONLY:nout
	
	include 'pram1.inc'
	include 'com1.inc'
!c	logical be(im,jm)
	real*4 u24(im,jm,2),v24(im,jm,2),h24(im,jm,2),&
         tmp4(im,jm,2),u14(im,jm,2),v14(im,jm,2),&
         h14(im,jm,2)
	integer*4 itt4

        if (itt.eq.0) then
         do 10 k=1,2
         do 10 j=1,jm
         do 10 i=1,im
	  u2(i,j,k)=0.
	  v2(i,j,k)=0.
	  u1(i,j,k)=0.
	  v1(i,j,k)=0.
	  if(be(i,j))then
	   h2(i,j,k)=eqdpth
	   h1(i,j,k)=eqdpth1
           tmp(i,j,k)=tmpmax
           endif
10       continue
        else
         rewind 513
         read(513) itt4,u24,v24,h24,tmp4,u14,v14,h14 
	 write(nout,*) 'OCEAN STEPS FROM UNIT 513 ITT!!!!!=',itt4
         rewind 513
         itt=0			!!! reset itt
        end if

	do i=1,im
	do j=1,jm
	  do k=1,2

	  u2(i,j,k)=u24(i,j,k)
	  v2(i,j,k)=v24(i,j,k)
	  h2(i,j,k)=h24(i,j,k)
!c	  h2(i,j,k)=eqdpth
	  u1(i,j,k)=u14(i,j,k)
	  v1(i,j,k)=v14(i,j,k)
	  h1(i,j,k)=h14(i,j,k)
!c	  h1(i,j,k)=eqdpth1
	  tmp(i,j,k)=tmp4(i,j,k)

	  if(.NOT.be(i,j)) then
	  h2(i,j,k)=eqdpth
	  h1(i,j,k)=eqdpth1
	  end if

	  end do
	end do
	end do
!c
!c===>Impose cyclic boundary condition
!c
        if (nbouny.eq.1) then
        do 16 k=1,2
        do 16 i=1,im
         u2(i,1,k)=u2(i,jm-1,k)
         u2(i,jm,k)=u2(i,2,k)
         v2(i,1,k)=v2(i,jm-1,k)
         v2(i,jm,k)=v2(i,2,k)
         h2(i,1,k)=h2(i,jm-1,k)
         h2(i,jm,k)=h2(i,2,k)
         u1(i,1,k)=u1(i,jm-1,k)
         u1(i,jm,k)=u1(i,2,k)
         v1(i,1,k)=v1(i,jm-1,k)
         v1(i,jm,k)=v1(i,2,k)
         h1(i,1,k)=h1(i,jm-1,k)
         h1(i,jm,k)=h1(i,2,k)
         tmp(i,1,k)=tmp(i,jm-1,k)
         tmp(i,jm,k)=tmp(i,2,k)
16     continue
       end if
!c       if (nbounx.eq.1) then
!c        do 17 k=1,2
!c        do 17 j=1,jm
!c         u2(1,j,k)=u2(im-1,j,k)
!c         u2(im,j,k)=u2(2,j,k)
!c         v2(1,j,k)=v2(im-1,j,k)
!c         v2(im,j,k)=v2(2,j,k)
!c         h2(1,j,k)=h2(im-1,j,k)
!c         h2(im,j,k)=h2(2,j,k)
!c         u1(1,j,k)=u1(im-1,j,k)
!c         u1(im,j,k)=u1(2,j,k)
!c         v1(1,j,k)=v1(im-1,j,k)
!c         v1(im,j,k)=v1(2,j,k)
!c         h1(1,j,k)=h1(im-1,j,k)
!c         h1(im,j,k)=h1(2,j,k)
!c         tmp(1,j,k)=tmp(im-1,j,k)
!c         tmp(im,j,k)=tmp(2,j,k)
!c17      continue
!c        end if
	return
	end subroutine initialo


subroutine initial(mtt)
  USE mo_doctor,        ONLY:nout

!qbao
! USE mo_couple,    ONLY: restart_u1, restart_u2, restart_itt,&
!      restart_v1, restart_v2, restart_h1, restart_h2, restart_tmp
!qbao

  include 'pram1.inc'
  include 'com1.inc'

  !c	logical be(im,jm)

  if (itt.eq.0) then
     do  k=1,2
        do  j=1,jm
           do  i=1,im
              u2(i,j,k)=0.
              v2(i,j,k)=0.
              u1(i,j,k)=0.
              v1(i,j,k)=0.
              if(be(i,j))then
                 h2(i,j,k)=eqdpth
                 h1(i,j,k)=eqdpth1
                 tmp(i,j,k)=tmpmax
              endif
           end do
        end do
     end do
  else

     !ik
     !ik         rewind 51
     !ik         read(51) itt,u2,v2,h2,tmp,u1,v1,h1 
     !ik         rewind 51
     !ik read the restart data

!qbao
     call restart_initial 
!    u1  = restart_u1
!    u2  = restart_u2
!    v1  = restart_v1
!    v2  = restart_v2
!    h1  = restart_h1
!    h2  = restart_h2
!    tmp = restart_tmp
!    itt = restart_itt
!qbao

     write(nout,*) 'OCEAN STEPS IN INITIAL ITT!!!!!=',itt

  end if

  !c=====this is used to keep the correctness of the S-N
  !c=====boundary of the observed SST

  if(mtt.eq.1000) then  !start from Jan1
     itt=0
  end if

  !c
  !c===>Impose cyclic boundary condition
  !c

  if (nbouny.eq.1) then
     do  k=1,2
        do  i=1,im
           u2(i,1,k)=u2(i,jm-1,k)
           u2(i,jm,k)=u2(i,2,k)
           v2(i,1,k)=v2(i,jm-1,k)
           v2(i,jm,k)=v2(i,2,k)
           h2(i,1,k)=h2(i,jm-1,k)
           h2(i,jm,k)=h2(i,2,k)
           u1(i,1,k)=u1(i,jm-1,k)
           u1(i,jm,k)=u1(i,2,k)
           v1(i,1,k)=v1(i,jm-1,k)
           v1(i,jm,k)=v1(i,2,k)
           h1(i,1,k)=h1(i,jm-1,k)
           h1(i,jm,k)=h1(i,2,k)
           tmp(i,1,k)=tmp(i,jm-1,k)
           tmp(i,jm,k)=tmp(i,2,k)
        end do
     end do
  end if

  !c       if (nbounx.eq.1) then
  !c        do 17 k=1,2
  !c        do 17 j=1,jm
  !c         u2(1,j,k)=u2(im-1,j,k)
  !c         u2(im,j,k)=u2(2,j,k)
  !c         v2(1,j,k)=v2(im-1,j,k)
  !c         v2(im,j,k)=v2(2,j,k)
  !c         h2(1,j,k)=h2(im-1,j,k)
  !c         h2(im,j,k)=h2(2,j,k)
  !c         u1(1,j,k)=u1(im-1,j,k)
  !c         u1(im,j,k)=u1(2,j,k)
  !c         v1(1,j,k)=v1(im-1,j,k)
  !c         v1(im,j,k)=v1(2,j,k)
  !c         h1(1,j,k)=h1(im-1,j,k)
  !c         h1(im,j,k)=h1(2,j,k)
  !c         tmp(1,j,k)=tmp(im-1,j,k)
  !c         tmp(im,j,k)=tmp(2,j,k)
  !c17      continue
  !c        end if

  return
end subroutine initial

!c
!c=========================topography
	subroutine topog
!c============coupled to AGCM
        USE mo_landsea,     ONLY: itopo

	include 'pram1.inc'
	include 'com1.inc'

!c-----set north;south as coastal boundary
	do i=1,im
	itopo(i,1)=0
	itopo(i,121)=0
	end do
!c==================================================================
	do i=1,im
	   do j=1,jm
	   topo(i,j)=float(itopo(i,j))
	   end do
	end do
!c------------------------------------------------------------
!c        if (nbounx.eq.1) then !periodic x direction 
!c        do 35 j=1,jm
!c         topo(1,j)=topo(im-1,j)
!c         topo(im,j)=topo(2,j)
!c35      continue
!c        end if
        if (nbouny.eq.1) then  !periodic y direction
        do 36 i=1,im
         topo(i,1)=topo(i,jm-1)
         topo(i,jm)=topo(i,2)
36      continue
        end if
!c----------------------------------------------------------
	do 37 j=1,jm
	do 37 i=1,im
	be(i,j)=.false.
	if(topo(i,j).gt.epsln) be(i,j)=.true.
	if(be(i,j)) topo(i,j)=eqdpth
!c	boun(i,j)=0.
37	continue

!c----get boun(im,jm) for Te smothing
!c	do 10 i=1,im
!c	do 12 j=25,37
!c	if(be(i,j)) then  !in the ocean
!c	boun(i,j)=1.0          !5s
!c	end if
!c12      continue
!c10      continue

!c	open(10,file='boun.index')
!c	write(10,*)'boun'
!c	write(10,100)((i-10*(i/10)),i=1,im)
!c	do j=jm,1,-1
!c	write(10,200) j, (boun(i,j),i=1,im),j
!c	end do
!c	close(10)

!c100	format(4x,180i1)
!c200     format(1x,i2,1x,180i1,1x,i2)
	return
	end subroutine topog

	subroutine constant
	include 'pram1.inc'
	include 'com1.inc'
        data epsln,grav,nergy,itt,tfr,tfnz,tfnm,wface,sface&
           /1.e-10,980.6,30,1,2.,0.5,0.5,-0.25,-30.25/   !modify here
!c	                   -------

        data ansamp,nsmth,alpha,maps,nrstrt,nf,netr,tmn,tmpmax,thmexp, &
           rhoa,cd,cw,wm,tn,r,gamma1,gamma2,nbuoy,nbounx,nbouny,&
            forcef,forcev,deltat,nseason,nsmthstr,hepsln,nbeta,beta,&
            nqin &
            /30.,1,.5,0,1,0,2,14.5,30.,2.5e-4,1.2e-3,1.5e-3,4.2e7,1.25, &
!c                                                               ----
             0.2,0.55,2.857e-2,4.e-4,1,1,0,1.,1.,0.,2,0,5.0e-5,1, &
!c                                      ---       --    ------
             2.28e-11,2/
!c                     ---

	data tcount,tadd0,nladv,cdfv,cdfv1,wdh,theta &
             /0.,0.,1,1.,1.,1000.,45/

        degrad=pi/180.
        omega=pi/43082.
        degcm=111.324e5
        earthr=degcm/degrad
	ufcoef=sqrt(rhoa*cdt0*10000.)
!c        q0unit=4.18*1000./30.4/8.64*(1.-albedo)
!cfu(7/14/97) solar harmonics is used
        q0unit=(1.-albedo)
	if(ntxy.eq.0) then
	 h0=1000.
	 crmt=1.
	 txyk=0.
         txyk1=0.
	else
	 crmt=0.
	 txyk=0.e-5
         txyk1=1.e-5
	end if
!c
	wm0=wm
	wm30=0.3			!!!
	yscale=-alog(wm30/wm0)/60.      !0.5 degree in y direction
	do j=1,jm
	wmy1(j)=wm0*exp(-abs(j-61)*yscale)
!c	wmy1(j)=1.25
	end do
!c
!c
!c    set newtonian cooling and rayleigh damping time scales
!c
        tfrr=1./(tfr*86400.)
        tfnrz=1./(tfnz*86400.)
        tfnrm=1./(tfnm*86400.)
!c
!c    friction layers along northern and southern boundary to
!c    kill costal kelvin waves
!c
        do 10 j=1,jm
        do 10 i=1,im
        dysouth=abs(yh(j)-yh(2))
        dynorth=abs(yh(j)-yh(jm-1))
        dxwest=abs(xh(i)-200.)   !160W
        dxeast=abs(xh(i)-xh(im-1))
	fe1(i,j)=0.
        fr(i,j)=tfrr*(exp(-dysouth/1.5)+exp(-dynorth/1.5))
!c     1         +1./150./86400.
        fn(i,j)=tfnrm*(exp(-dysouth/1.5)+exp(-dynorth/1.5))
        fe(i,j)=tfrr*(exp(-dysouth/1.5)+exp(-dynorth/1.5))
!c     1        +1./150./86400.
	if(xh(i).ge.200.and.xh(i).le.290) then
        fe1(i,j)=tfnrm
        goto 333	
	end if
	if(xh(i).ge.180.and.xh(i).le.200) then
        fe1(i,j)=tfnrm*exp(-dxwest/1.5)
	end if
333     visc(i,j)=2.0e7	!!! cm2/s(from 1.e8 to 1.e7)
        dffs(i,j)=2.0e7
!c------which will affect the short period instability (McCreary 1992)

10      continue
        tspd=86400./dt
	ntspdo=int(tspd)              !total steps per day
!ik not used        ittend=tspd*ndays+.5
        nergy=tspd*nergy+.5
        ansamp=tspd*ansamp       ! ping - changed nsamp to real (x/day)
        insamp=int(ansamp)
!c        daysave=ndays-360      !! save last year's 12 monthly mean
	daysave=0   !! save all monthly mean

        if (insamp.eq.0) insamp =1
!ctim...............................................
!c
        insmpir=tspd*5		!!! 5 days
!c       insmpir=insamp
        ittwrt=tspd*1.		!!! 1 days
        ittatm=itt
!c     set forcing parameters
!c
!c     timep= forcing period (1 yr)
!c     phlag= phase lag
!c
!c     since the winds start at Jan. 15, phlag = 15 days to
!c     set the forcing starting from Jan 1.
!c            i.e., phlag=timep/24.
!c
        timep=360*86400.0
        phlag=timep/24.
!c--------------------------------------------------------------
!c      phlag=timep/4.		!!! March 15
!c           =timep/4.*3.	!!! Sep.  15
!c--------------------------------------------------------------
!c      imthlck=3		!!! March 15
!c             =9		!!! Sep.  15
!c--------------------------------------------------------------
!c	lcktim=1	as usual
!c	lcktim=0	phase locking for a particular moment
!c--------------------------------------------------------------
	lcktim=1	!!!
	imthlck=3	!!! acting only when lcktim=0 and nqin=2
	iyear=int(360.*24.*3600./dt+0.001) !!! number of steps/year
!c----ten-year total steps
	nitt=3600*ntspdo
!c
	return
        end subroutine constant
!c------------------------------------------------------------------c
!c     end of the IOMset-up code (3/17/99)
!c------------------------------------------------------------------c

