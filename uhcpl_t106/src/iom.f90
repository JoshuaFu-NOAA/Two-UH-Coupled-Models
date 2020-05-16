!c===================================================================c
!c******************* O C E A N (0-360,30S-30N)+++++++++*************c
!c===================================================================c
!c this version ocean model is the 2-layers ocean model (M & Y 1990) 
!c  coupled  with atmospheric model just over Pacific. 
!c*******************************************************************c
!c========extend to global tropics(Samuel Fu,03/09/99)
!c*******************************************************************c
!c   1, Gaspar (1988) entrainment (include new solar penetration)
!c   2, use atmospheric model Qa and Ta (latent,sensible,longwave)
!c   3,in 'PHYSICS', vertical momentum & heat mixing
!c      Tr=14.,amu(i,j)=1.23
!c   4, mixed-layer friction umum + shear instability (Price 1973) 
!c-------------------------------------------------------------------c
!c    Samuel Fu (01/08/1998,coupled to AGCM)
!c-------------------------------------------------------------------c
!c===>uncoupled mode (set ncpl=0 in pram1.inc)
!c====================================================================

subroutine  iom(tcp,loop) 

  USE mo_doctor,         ONLY: nout

!qbao
! USE mo_couple, ONLY: restart_u1, restart_u2, restart_v1, restart_v2, &
!      restart_h1, restart_h2, restart_tmp, restart_itt
!qbao

  !c******************************************************************** 
  include 'pram1.inc'
  include 'com1.inc'
  external funcwe
  !c
  !c------------- coupling with AGCM-----------------------------------
  common /scouple/tua(im,jm),tva(im,jm),wfua(im,jm),qsola(im,jm),&
       hfluxa(im,jm),sst_nudg(im,jm)
  common /fudate/runday,dta,dto,day1,tseg,numcpl,&
       numatm,numocn
  dimension wfu(im,jm),tcp(im,jm)

  !c=====================================================================
  !c    try to save daily data
  !c=========add anomaly wind forcing over western+central Pacific
  !c===from 180-140W

  icouple=itt

  !c	icouple=int(itt/tspd)
  !c	if(loop.eq.1.and.nano.eq.1.and.
  !c     1     icouple.lt.300)then
  !c	do j=1,jm
  !c	ry=abs(float(j-31))/20.
  !c	do i=1,im
  !c	if(i.ge.91.and.i.le.111) rx=0.
  !c	if(i.lt.91) rx=float(91-i)/3.
  !c	if(i.gt.111) rx=float(i-111)/3.
  !c	uana=6.*exp(-(rx**2+ry**2)) !first 10 mths
  !c	uano(i,j)=rhoa*cdt0*uana*6.*10000.
  !c	end do
  !c	end do
  !c	do i=1,im
  !c	do j=1,jm
  !c	tua(i,j)=tua(i,j)+uano(i,j)
  !c	end do
  !c	end do
  !c	end if

  !c---------------------------------------------------------------------
  !c                          START TIME INTEGRATION
  !c1111111111111111111111111111111111111111111111111111111111111111111111

  call bound0(h2,1,nbounx,nbouny)
  call bound0(h1,1,nbounx,nbouny)
  call bound0(tmp,1,nbounx,nbouny)

2000 continue
  itt=itt+1

  call getforce(wfu) !windspeed after adjustment

  !c	print*,'u,v,sst=',tu(136,31),tv(136,31),tmp(136,31,n)
  !c	print*,'U,tair,cld=',wfu(136,31),tair(136,31),cld(136,31)
  !c	print*,'qn,qo,qlw,ql,qs,qp=',tfc1(136,31),tfc2(136,31),
  !c     1    tfc3(136,31),tfc4(136,31),tfc5(136,31),tfc6(136,31)
  !c
  !c==> Zonal advection 
  !c

  do j=1,jm
     do i=1,im
        u2(i,j,np)  = u2(i,j,n)
        v2(i,j,np)  = v2(i,j,n)
        h2(i,j,np)  = h2(i,j,n)
        u1(i,j,np)  = u1(i,j,n)
        v1(i,j,np)  = v1(i,j,n)
        h1(i,j,np)  = h1(i,j,n)
        tmp(i,j,np) = tmp(i,j,n)
     enddo
  enddo

  !c	print*,'itt,tmp(n),u2,u1',itt,tmp(76,22,n),
  !c    1     u2(76,22,n),u1(76,22,n)
  do j=1,jm
     do i=1,im
	if(be(i,j))then
           uh(i)=u2(i,j,np)     !!! zonal advect velocity
           u1h(i)=u1(i,j,np)    !!! zonal advect velocity
           ui(i)=u2(i,j,n)
           vi(i)=v2(i,j,n)
           if(h2(i,j,n).le.0.) then
              call out_put
           end if
           hlni(i)=alog(h2(i,j,n))
           u1i(i)=u1(i,j,n)
           v1i(i)=v1(i,j,n)
           h1lni(i)=alog(h1(i,j,n))
           tmpi(i)=tmp(i,j,n)
	endif
     enddo
	
     do i=1,im
	if(be(i,j))then
           if(uh(i).ge.0)nsign=1
           if(uh(i).lt.0)nsign=-1
           call mux(ux1,ux2,ux3,uh,i,j,ijm,dt,nbounx)
           call advecx(ada,ui,i,j,ux1,ux2,ux3,nsign)
           adu(i)=ada
           call advecx(ada,vi,i,j,ux1,ux2,ux3,nsign)
           adv(i)=ada
           call advecx(ada,hlni,i,j,ux1,ux2,ux3,nsign)
           adh(i)=ada

           if(u1h(i).ge.0)nsign=1
           if(u1h(i).lt.0)nsign=-1
           call mux(ux1,ux2,ux3,u1h,i,j,ijm,dt,nbounx)
           call advecx(ada,u1i,i,j,ux1,ux2,ux3,nsign)
           adu1(i)=ada
           call advecx(ada,v1i,i,j,ux1,ux2,ux3,nsign)
           adv1(i)=ada
           call advecx(ada,h1lni,i,j,ux1,ux2,ux3,nsign)
           adh1(i)=ada
           call advecx(ada,tmpi,i,j,ux1,ux2,ux3,nsign)
           adt(i)=ada
	endif
     enddo

     do i=1,im
	if(be(i,j))then
           u2(i,j,n)=adu(i)
           v2(i,j,n)=adv(i)
           h2(i,j,n)=adh(i)    !!! Here h is actually ln(h)
           u1(i,j,n)=adu1(i)
           v1(i,j,n)=adv1(i)
           h1(i,j,n)=adh1(i)  !!! Here h1 is actually ln(h1)
           tmp(i,j,n)=adt(i)
           !c ---get zonal advection 
           xah(i,j)=(adt(i)-tmp(i,j,np))/dt
           
	endif
     enddo
     
  enddo

  !c	write(6,*)'done zonal advection'
  !c
  !c==> Meridional advection 
  !c
  do i=1,im
     do j=1,jm
	if(be(i,j))then
           uh(j)=v2(i,j,np)   !!! meridional advect velocity
           u1h(j)=v1(i,j,np) !!! meridional advect velocity
           ui(j)=u2(i,j,n)
           vi(j)=v2(i,j,n)
           hlni(j)=h2(i,j,n)  
           u1i(j)=u1(i,j,n)
           v1i(j)=v1(i,j,n)
           h1lni(j)=h1(i,j,n)
           tmpi(j)=tmp(i,j,n)
	endif
     enddo
	
     do j=1,jm
	if(be(i,j))then
           if(uh(j).ge.0)nsign=1
           if(uh(j).lt.0)nsign=-1
           call muy(uy1,uy2,uy3,uh,i,j,ijm,dt,nbouny)
           call advecy(ada,ui,i,j,uy1,uy2,uy3,nsign)
           adu(j)=ada
           call advecy(ada,vi,i,j,uy1,uy2,uy3,nsign)
           adv(j)=ada
           call advecy(ada,hlni,i,j,uy1,uy2,uy3,nsign)
           adh(j)=ada

           if(u1h(j).ge.0)nsign=1
           if(u1h(j).lt.0)nsign=-1
           call muy(uy1,uy2,uy3,u1h,i,j,ijm,dt,nbouny)
           call advecy(ada,u1i,i,j,uy1,uy2,uy3,nsign)
           adu1(j)=ada
           call advecy(ada,v1i,i,j,uy1,uy2,uy3,nsign)
           adv1(j)=ada
           call advecy(ada,h1lni,i,j,uy1,uy2,uy3,nsign)
           adh1(j)=ada
           call advecy(ada,tmpi,i,j,uy1,uy2,uy3,nsign)
           adt(j)=ada
	endif
     enddo

     do j=1,jm
	if(be(i,j))then
           adhxy(i,j)=(adh(j)-alog(h2(i,j,np)))/dt 
           !c                      !!! advect tendency of d(ln(h2))/dt

           adh1xy(i,j)=(adh1(j)-alog(h1(i,j,np)))/dt 
           !c                      !!! advect tendency of d(ln(h1))/dt

           !c ---get meridional advection 
           yah(i,j)=(adt(j)-tmp(i,j,n))/dt

           u2(i,j,n)=adu(j)
           v2(i,j,n)=adv(j)
           u1(i,j,n)=adu1(j)
           v1(i,j,n)=adv1(j)
           tmp(i,j,n)=adt(j)
	endif
     enddo

  enddo

  call bound0(tmp,1,nbounx,nbouny)

  !c	print*,'itt,tmp(n),u2,u1',itt,tmp(76,22,n),u2(76,22,n),u1(76,22,n)
  !c	print*,'itt,au,av',itt,tu(76,22),tv(76,22)

  !c	write(6,*)'done meridional advection' !!h and h1 are not updated
  !c
  !c===> Smoth boundary region
  !c
  if(mod(itt,132).eq.0) then
     call bduv(1,nbounx,nbouny)
     call bdht(tmp,1,nbounx,nbouny)
  end if
  
  call bound0(tmp,1,nbounx,nbouny)

  do j=1,jm
     do i=1,im
	tmp(i,j,np)=tmp(i,j,n)
     enddo
  enddo
  !c	write(6,*)'done bound'
  !c
  !c===> Geostrophic adjust process
  !c
  !c------------------------------------------------------------------c
  call adjust                  !!! temperature is not updated	
  call bound0(h1,2,nbounx,nbouny)
  call bound0(h2,2,nbounx,nbouny)

  if(mod(itt,132).eq.0) then    !smothing boundary
     call bduv(2,nbounx,nbouny)
     call bdht(h1,2,nbounx,nbouny)
     call bdht(h2,2,nbounx,nbouny)
  end if

  call bound0(h1,2,nbounx,nbouny)
  call bound0(h2,2,nbounx,nbouny)
  
  !c	write(6,*)'done adjust'
  !c......
  !c......
  do j=1,jm
     do i=1,im
        if(be(i,j)) then
           h(i,j)=h1(i,j,np)+h2(i,j,np)
           us(i,j)=u1(i,j,np)-u2(i,j,np)
           vs(i,j)=v1(i,j,np)-v2(i,j,np)

           !c---------------------------------------------------------

           u(i,j)=u1(i,j,np)-us(i,j)*h2(i,j,np)/h(i,j)
           v(i,j)=v1(i,j,np)-vs(i,j)*h2(i,j,np)/h(i,j)

           !c**********************************************************
           !cFU(10/20/97)
           !c**** using wang's subsurface temp as entrained temp
           !fu	wjh=upw(i,j)/h(i,j)
           !fu	whh=900.+500.*tanh((wjh-z00)/h00)
           !c------from 5S-5N
           !c	if(j.le.36.and.j.ge.26) then
           !c	tse(i,j)=2.-1.9*tanh((h(i,j)-5500.)/4000.)
           !c	whh=100.+(wdh-100.)*(tanh((float(i)-56.)/5.)+1.)/2.
           !c	end if
           !c fu++3/31/2003	tse(i,j)=0.5
           !c	tee=14.
           !	whh=wdh
           !        fuaa=wgm*whh*(tmp(i,j,np)-tmn)&
           !           /h2(i,j,np)
           !	tseh(i,j)=fuaa**2
           !c-----------shear-instability adjustment------------------ 
           !	ssp=(us(i,j)*us(i,j)+vs(i,j)*vs(i,j))/100./100. !m
           !	rup=.23*9.8*tse(i,j)*h1(i,j,np)
           !	rdn=1025.*ssp*100.
           !	ri(i,j)=rup/rdn
           !	if(ri(i,j).le.rict) then
           !	tse(i,j)=tse(i,j)*ri(i,j)/rict
           !	end if
           !	if(tseh(i,j).lt.0.001) tseh(i,j)=0.001
           !	if(tseh(i,j).gt.2.) tseh(i,j)=2.

           tse(i,j)=0.5   ! 5/5/2003
           hfu=h(i,j)
           if(hfu.gt.fumax) hfu=fumax
           abc=tta1+(tta2-tta1)*((fumax-hfu)/(fumax-fumin))**2.5
           tseh(i,j)=abc*hfu

           !c==========================================================
	
	endif
     enddo
  enddo

  !c****smoth tse along the boundaries
  !c	call  smoth2(be,boun,tse,tse,im,jm,0.5)
  !c============================================================
  !c	call out_put

  call entrain(wfu)  !scalar windspeed

  !c	call bound1(we)
  !c	call bound1(wet)
  !c	call filter0(we,be,nbounx)
  !c	call filter0(wet,be,nbounx)
  !c	print*,'itt,tmp,u,u1,BP',itt,tmp(76,22,np),u2(76,22,n),u1(76,22,n)

  call physics(sst_nudg)

  !c	print*,'itt,tmp,u,u1,AP',itt,tmp(76,22,np),u2(76,22,n),u1(76,22,n)
  !c	print*,'done physics!'

  call bound0(h1,1,nbounx,nbouny)
  call bound0(h2,1,nbounx,nbouny)
  call bound0(tmp,1,nbounx,nbouny)

  if(mod(itt,132).eq.0) then       !smothing every two-day
     call bduv(1,nbounx,nbouny)
     call bdht(h1,1,nbounx,nbouny)
     call bdht(h2,1,nbounx,nbouny)
     call bdht(tmp,1,nbounx,nbouny)
  end if
  
  call bound0(h1,1,nbounx,nbouny)
  call bound0(h2,1,nbounx,nbouny)
  call bound0(tmp,1,nbounx,nbouny)

  !c===============================================================	
  if(mod(itt,62).eq.0) then      !smothing every ten days ?
     call filterm(u2,be,n,nbounx)   !do shapiro smothing
     call filterm(v2,be,n,nbounx)
     call filterm(h2,be,n,nbounx)
     call filterm(u1,be,n,nbounx)
     call filterm(v1,be,n,nbounx)
     call filterm(h1,be,n,nbounx)
     call filterm(tmp,be,n,nbounx)
  endif

  call bound0(h1,1,nbounx,nbouny)
  call bound0(h2,1,nbounx,nbouny)
  call bound0(tmp,1,nbounx,nbouny)
  !c	write(6,*)'done filter'

  do j=1,jm
     do i=1,im
	tcp(i,j)=sst(i,j)    !observed sst outside
	if(be(i,j))then
           !cFu(2/13/97)

           if(tmp(i,j,n).le.18) tmp(i,j,n)=18
           if(tmp(i,j,n).ge.34.5) tmp(i,j,n)=34.5
           tcp(i,j)=tmp(i,j,n)                       !model sst

           !c-----------------------------------------------------------------

           dltsst(i,j)=tmp(i,j,n)-sst(i,j)
	else
           tmp(i,j,1)=sst(i,j)   !over land
           tmp(i,j,2)=sst(i,j)   !over land
           dltsst(i,j)=0.
	endif
     enddo
  enddo

  !c
  !c===> Get daily mean of MLT for atmos model
  !c
  !c	if(ncpl.eq.1)then
  !c 	 call tempm(tmpc,tmp,1,tspd)
  !c	endif
  !c
  !c===> For monitor
  !c
  !c	if(mod(itt,6).eq.0)then
  !c9999	continue
  !c	write(6,*)'itt=',itt
  !c	i0=40
  !c	j0=31
  !c	k0=n
  !c	write(6,*)u1(i0,j0,k0),v1(i0,j0,k0),h1(i0,j0,k0)/100.
  !c	write(6,*)u(i0,j0,k0),v(i0,j0,k0),h(i0,j0,k0)/100.
  !c	write(6,*)te(i0,j0),tmp(i0,j0,k0),tmpc(i0,j0)
  !c	write(6,*)we(i0,j0),hstar(i0,j0)/100.,h1(i0,j0,2)/100.
  !c	write(6,*)' Surface Zonal Wind from Atmosphric Model(128E-280E)'
  !c	write(6,999)(tu(ilook,j0),ilook=5,81,4)
  !c	write(6,*)' surface meridional wind'
  !c	write(6,999)(tv(ilook,j0),ilook=20,80,10)
  !c	write(6,*)' Mixed Layer Temperature (128E-280E)'
  !c	write(6,999)(tmp(ilook,j0,1),ilook=5,81,4)
  !c	write(6,*)' Daily Mean Mixed-Layer Temperature (128E-280E)'
  !c	write(6,999)(tmpc(ilook,j0),ilook=5,81,4)
  !c	write(6,*)' entrainment rate'
  !c	write(6,'(1x,7e10.3)')(we(ilook,j0),ilook=20,80,10)
  !c	endif
  !c999	format(1x,10f6.2)
  !c	write(6,*)qnet(i0,j0)*cw*h1(i0,j0,k0)/1000.,tfc1(i0,j0),tfc6(i0,j0)
  !c	write(6,'(4e12.3)')qnet(i0,j0),zah(i0,j0),tnd(i0,j0),tdf(i0,j0)
  !c	j0=3
  !c	write(6,*)u1(i0,j0,k0),v1(i0,j0,k0),h1(i0,j0,k0)/100.
  !c	write(6,*)we(i0,j0),te(i0,j0),tmp(i0,j0,k0)
  !c	call cycle
  
  !c
  !c******************************************************************* 
  !c
  !c     time smoothing for all predicted variables
  !c        (from original Tim Li's model)
  !c      if (mod(itt,nsmth).ne.0) go to 4800
  !c      call robert(uh,bu,im,jm,nm,n,np,0.1)
  !c      call robert(vh,bv,im,jm,nm,n,np,0.1)
  !c      call robert(e,be,im,jm,nm,n,np,0.1)
  !c      call robert(uh1,bu,im,jm,nm,n,np,0.1)
  !c      call robert(vh1,bv,im,jm,nm,n,np,0.1)
  !c      call robert(e1,be,im,jm,nm,n,np,0.1)
  !c      call robert(tmp,be,im,jm,nm,n,np,0.1)
  !c4800  continue
  !c
  !c===>Write restart file at the end of every year
  !c
  !ik        if (lwtime.and.loop.eq.numocn)then
  !ik         write(51) itt,u2,v2,h2,tmp,u1,v1,h1  !only final year!
  !ik         write(nout,*) 'save new 51 for next year'
  !ik	 rewind 51
  !ik	endif
  !ik store restart fields for later use
  
  ! store restart data for next iteration
  IF (loop == numocn) THEN
!qbao
     call restart_ocn
!    restart_itt = itt
!    restart_u1  = u1
!    restart_u2  = u2
!    restart_v1  = v1
!    restart_v2  = v2
!    restart_h1  = h1
!    restart_h2  = h2
!    restart_tmp = tmp
!qbao
  END IF

  !c------------------------------------------------------------------
  !c
  !c********************** end numerical integration *******************
  !c
  !c------------------------------------------------------------------
  !c
  !c====> Calculate domain mean energy
  !c
  !c	if(mod(nday,nergy).eq.0)then
  !c	 call energy
  !c	endif
  !c
  !c===> Detect stability
  !c
  do j=1,jm
     do i=1,im
	timc1=abs(u1(i,j,n))
	timc2=abs(h1(i,j,n))
	timc3=abs(tmp(i,j,n))
	if(timc1.gt.1000..or.timc2.gt.50000.0.or.timc3&
             .gt.40.) then
           write(nout,*)'step,day,current,h1,sst,i,j'
           write(nout,*) itt,day,timc1,timc2,timc3,i,j
           write(nout,*)'u,v,wind=',tua(i,j),tva(i,j),wfua(i,j)
           write(nout,*)'surface heat fluxes:'
           write(nout,*)'qsol=',qout(i,j),'qlw=',qlw(i,j),'qnet=',qnet(i,j)
           write(nout,*)'qin=',qin(i,j),'qlh=',qlh(i,j),'qsh=',qsh(i,j)
           write(nout,*)'sst balance:'
           write(nout,*)'tfc=',tfc(i,j),'zah=',zah(i,j),'tnd=',tnd(i,j)
           write(nout,*)'tdf=',tdf(i,j),'tmp(np)=',tmp(i,j,np)

           !c	call out_put
           stop
	end if
     end do
  end do
  !c
  !c=============================================================
  !cFU===>   save the daily-mean results

!ik daysave is zero and day is not set in the coupled mode, cut it
!ik  if(day.ge.daysave) then

  !c=====>get two-dimension fields.
  do i=1,im
     do j=1,jm
        u22(i,j)=u2(i,j,n)
        v22(i,j)=v2(i,j,n)
        h22(i,j)=h2(i,j,n)
        tmp2(i,j)=tmp(i,j,n)
        u12(i,j)=u1(i,j,n)
        v12(i,j)=v1(i,j,n)
        h12(i,j)=h1(i,j,n)
     end do
  end do
  !c=====>
  call daym(numocn,wfu,loop)

!ik  end if

  !c------------------------------------------------------------
  !c   couple with atmospheric model every 12 steps(one day)
  !c------------------------------------------------------------
  !c	call tempm(tmpc,tmp,n,ncouple)
  !c-----------------------------------------------------------

  if(loop.eq.numocn) then

     icouple=icouple+1

     !c------------end of atmospheric wind and PBL tair qair readings
!ik     WRITE(nout,*)'icouple,iloop=',icouple,iloop
     WRITE(nout,*)'icouple,loop=',icouple,loop
     WRITE(nout,*)'nday=',nday,'ts(70E)=',tmp(73,61,n),&
          'ts(90W)=',tmp(273,61,n),'ts(6W)=',tmp(356,61,n),&
          'upw=',upw(252,61)
     WRITE(nout,*)'h1=',h1(252,61,n)/100.,'h=',h(252,61)/100.&
          ,'u1,u2=',u1(252,61,n),u2(252,61,n),&
          'we=',we(252,61)*100000.,'be=',be(252,61) 
     !c-----------------------------------------------------------------c

  end if

  !c******************************************************************
  !c?????????????????????????????????????????????????????????????????
  !c
  !c          END OF THE TIME INTEGRATION CYCLE
  !c
  !cqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq
  !c	if(itt.lt.ittend) go to 2000
  !c	print*,'ocean new 51 file!!!!'
  !c       write(51) itt,u2,v2,h2,tmp,u1,v1,h1
  
  return	
  
end subroutine iom
!------- END IOM


!c
!c------------------------------------------------------
!c
      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB,II,JJ)
!c      DIMENSION XB1(1),XB2(1)
      EXTERNAL FX
      NBB=NB
      NB=0
      X=X1
      DX=(X2-X1)/N
      FP=FX(X,II,JJ)
      DO 11 I=1,N
        X=X+DX
        FC=FX(X,II,JJ)
        IF(FC*FP.LT.0.) THEN
          NB=NB+1
!c          XB1(NB)=X-DX
!c          XB2(NB)=X
          XB1=X-DX
          XB2=X
        ENDIF
        FP=FC
        IF(NBB.EQ.NB)RETURN
11    CONTINUE
      RETURN
      END
!c
!c**** function funcwe ***************************************
!c
        function funcwe(x,i,j)
	parameter(im=720,jm=121)
        common /forc/ tx(im,jm),ty(im,jm),tu(im,jm),tv(im,jm),&
                     sst(im,jm),qnet(im,jm),qpet(im,jm),qout(im,jm),&
                     qin(im,jm),ef1(im,jm),ef2(im,jm),uf(im,jm),&
                      wmy(im,jm),wmy1(jm)

        common /scalar/ grav,omega,itt,ndays,tfr,tfnz,tfnm,epsln,dt,&
                       degrad,nergy,ansamp,alpha,nrstrt,maps,nsmth,nf,&
                       netr,tmn,tmpmax,thmexp,rhoa,cd,cw,wm,tn,r,&
                       gamma1,gamma2,nbuoy,nbounx,nbouny,forcef,forcev,&
                       deltat,nseason,nsmthstr,hepsln,nbeta,beta,&
                       ittend,tspd,txyk1,txyk,cdfv,ufcoef,&
                       q0unit,earthr,daysave,ittwrt,ittatm,phlag,&
                       lcktim,imthlck,day,nday,iyear,tcount,tadd0,&
                       nladv,cdfv1,wdh,theta,timep,nqin,ncouple&
                      ,iloop
!c
	data hr/1700./,rr/0.35/

         wndstr=2.*wmy(i,j)*uf(i,j)**3
         bckdss=2.*hepsln*x

 	pp=exp(-x/hr)-2.*hr*(1.-exp(-x/hr))&
          /x	
        htflux=thmexp*grav*(qin(i,j)+rr*qout(i,j)*pp)/cw
!c         htflux=thmexp*grav*(qin(i,j)+wef1*qout(i,j))/cw
!c         petwe=-(x-2./gamma2)*thmexp*grav/cw*qout(i,j)*(1.-r)
!c	 petwe=petwe*(1.-exp(-gamma2*x))
	petwe=0.
         funcwe=wndstr-bckdss-0.5*x*((1.-tn)*abs(htflux)+(1.+tn)*htflux)&
              +petwe

        return
        end


	subroutine monthm(p,pm,dt,tadd,tcount,tadd0,ifile)

	parameter(im=720,jm=121)
	dimension p(im,jm),pm(im,jm)

	data tmonth/2592000./
	
!c-->    Summarize monthly mean field
!c
        do 20 j=1,jm
	do 20 i=1,im
	 pm(i,j)=pm(i,j)+tadd*p(i,j)
20	continue
!c
!c-->    Write out month mean field if tcount > tmonth
!c
        if(tcount.ge.tmonth)then
	 do 30 j=1,jm
	 do 30 i=1,im
	  pm(i,j)=pm(i,j)/tmonth
30	 continue
	 write(ifile)pm
!c
!c-->    Save residue of pm 
!c
         tadd0=tcount-tmonth
         do 40 j=1,jm
	 do 40 i=1,im
          pm(i,j)=tadd0*p(i,j)
40	 continue
	endif	
	return
	end

	subroutine daym(numocn,wfu,loop)

!qbao
!         USE mo_couple, ONLY: output_unit
!qbao

	include 'pram1.inc'
	include 'com1.inc'
        dimension wfu(im,jm)
!c-->    accumulate 23 ocean variables 
!c
!ik	if(iloop.eq.1) then
	if(loop.eq.1) then

	do i=1,im
	do j=1,jm

	f9m(i,j)=0.

	if(be(i,j)) then
	f1m(i,j)=0.
	f2m(i,j)=0.
	f3m(i,j)=0.
	f4m(i,j)=0.
	f5m(i,j)=0.
	f6m(i,j)=0.
	f7m(i,j)=0.
	f8m(i,j)=0.
	f10m(i,j)=0.
	f11m(i,j)=0.
	f12m(i,j)=0.
	f13m(i,j)=0.
	f14m(i,j)=0.
	f15m(i,j)=0.
	f16m(i,j)=0.
	f17m(i,j)=0.
	f18m(i,j)=0.
	f19m(i,j)=0.
	f20m(i,j)=0.
	f21m(i,j)=0.
	f22m(i,j)=0.
	f23m(i,j)=0.
	goto 777
	end if

	f1m(i,j)=-9999.
	f2m(i,j)=-9999.
	f3m(i,j)=-9999.
	f4m(i,j)=-9999.
	f5m(i,j)=-9999.
	f6m(i,j)=-9999.
	f7m(i,j)=-9999.
	f8m(i,j)=-9999.
!	f9m(i,j)=-9999.
	f10m(i,j)=-9999.
	f11m(i,j)=-9999.
	f12m(i,j)=-9999.
	f13m(i,j)=-9999.
	f14m(i,j)=-9999.
	f15m(i,j)=-9999.
	f16m(i,j)=-9999.
	f17m(i,j)=-9999.
	f18m(i,j)=-9999.
	f19m(i,j)=-9999.
	f20m(i,j)=-9999.
	f21m(i,j)=-9999.
	f22m(i,j)=-9999.
	f23m(i,j)=-9999.
777     continue	
	end do
	end do

	end if

        do j=1,jm
	do i=1,im

	f9m(i,j)=f9m(i,j)+tmp2(i,j)/float(numocn*1)
	
	if(be(i,j)) then
	f1m(i,j)=f1m(i,j)+tseh(i,j)/float(numocn*1)
	f2m(i,j)=f2m(i,j)+flx(i,j)/float(numocn*1)
	f3m(i,j)=f3m(i,j)+u12(i,j)/float(numocn*1)
	f4m(i,j)=f4m(i,j)+v12(i,j)/float(numocn*1)
	f5m(i,j)=f5m(i,j)+u22(i,j)/float(numocn*1)
	f6m(i,j)=f6m(i,j)+v22(i,j)/float(numocn*1)
	f7m(i,j)=f7m(i,j)+h(i,j)/float(numocn*1)
	f8m(i,j)=f8m(i,j)+h12(i,j)/float(numocn*1)

	f10m(i,j)=f10m(i,j)+we(i,j)/float(numocn*1)

	f11m(i,j)=f11m(i,j)+tfc1(i,j)/float(numocn*1)
	f12m(i,j)=f12m(i,j)+tfc2(i,j)/float(numocn*1)
	f13m(i,j)=f13m(i,j)+tx(i,j)/float(numocn*1)
	f14m(i,j)=f14m(i,j)+ty(i,j)/float(numocn*1)
	f15m(i,j)=f15m(i,j)+wfu(i,j)/float(numocn*1)
	f16m(i,j)=f16m(i,j)+tfc6(i,j)/float(numocn*1)
	f17m(i,j)=f17m(i,j)+xah(i,j)/float(numocn*1)
	f18m(i,j)=f18m(i,j)+yah(i,j)/float(numocn*1)
	f19m(i,j)=f19m(i,j)+zah(i,j)/float(numocn*1)
	f20m(i,j)=f20m(i,j)+tfc(i,j)/float(numocn*1)
	f21m(i,j)=f21m(i,j)+tdf(i,j)/float(numocn*1)
	f22m(i,j)=f22m(i,j)+tnd(i,j)/float(numocn*1)
	f23m(i,j)=f23m(i,j)+dltsst(i,j)/float(numocn*1)
	end if
	end do
	end do
!c
!c-->    Write out daily mean field
!c
!ik save every ocean cycle the data
        if(loop.eq.numocn)then !save every-numocn-steps
	
!c	write(output_unit) f2m,f3m,f4m,f5m,f6m,f7m,f8m,f9m,f10m,f11m,
!c     $   f12m,f16m,f17m,f18m,f19m,f20m,f21m,f22m,f23m
!write(output_unit) ((f1m(i,j),i=1,im),j=1,jm)
!write(output_unit) ((f2m(i,j),i=1,im),j=1,jm)
!
!fu++============================================================c
        write(117) ((f3m(i,j),i=1,im),j=1,jm)  !u1
        write(117) ((f4m(i,j),i=1,im),j=1,jm)  !v1
        write(117) ((f5m(i,j),i=1,im),j=1,jm)  !u2
        write(117) ((f6m(i,j),i=1,im),j=1,jm)  !v2
        write(117) ((f7m(i,j),i=1,im),j=1,jm)  !h
        write(117) ((f8m(i,j),i=1,im),j=1,jm)  !h1
        write(117) ((f9m(i,j),i=1,im),j=1,jm)  !sst
        write(117) ((f10m(i,j),i=1,im),j=1,jm) !we
!cfu++============================================================c
!write(117) ((f11m(i,j),i=1,im),j=1,jm)
!write(117) ((f12m(i,j),i=1,im),j=1,jm)
!c----------------------------------------------------------------c
        write(117) ((f13m(i,j),i=1,im),j=1,jm) !tx
        write(117) ((f14m(i,j),i=1,im),j=1,jm) !ty
        write(117) ((f15m(i,j),i=1,im),j=1,jm) !wfu
!c================================================================c
!write(117) ((f16m(i,j),i=1,im),j=1,jm)
!c
!c----------------------------------------------------------------c
        write(117) ((f17m(i,j),i=1,im),j=1,jm)  !xadv 
        write(117) ((f18m(i,j),i=1,im),j=1,jm)  !yadv
        write(117) ((f19m(i,j),i=1,im),j=1,jm)  !zadv
        write(117) ((f20m(i,j),i=1,im),j=1,jm)  !flux

!c----------------------------------------------------------------c
!write(output_unit) ((f21m(i,j),i=1,im),j=1,jm)
!write(output_unit) ((f22m(i,j),i=1,im),j=1,jm)
!write(output_unit) ((f23m(i,j),i=1,im),j=1,jm)

         endif	
!ik	iloop=iloop+1
!ik wrong for 365 day years
!ik	if(iloop.gt.(numocn*5)) iloop=1
	return
	end

	
	subroutine heatm(p,pm,dt,tadd,tcount,tadd0,ifile)
	parameter(im=720,jm=121)
	dimension p(im,jm),pm(im,jm)
	data tmonth/2592000./
	
!c-->    Summarize monthly mean field
!c
        do 20 j=1,jm
	do 20 i=1,im
	 pm(i,j)=pm(i,j)+tadd*p(i,j)
20	continue
!c
!c-->    Write out month mean field if tcount > tmonth
!c
        if(tcount.ge.tmonth)then
	 do 30 j=1,jm
	 do 30 i=1,im
	  pm(i,j)=pm(i,j)/tmonth
30	 continue
	 write(ifile)pm
!c
!c-->    Save residue of pm 
!c
         tadd0=tcount-tmonth
         do 40 j=1,jm
	 do 40 i=1,im
          pm(i,j)=tadd0*p(i,j)
40	 continue
	endif	
	return
	end

	subroutine get2d(p3,p2,np)
	parameter(im=720,jm=121)
	dimension p3(im,jm,2),p2(im,jm)
	do 10 j=1,jm
	do 10 i=1,im
	p2(i,j)=p3(i,j,np)
10	continue
	return
	end 

	subroutine mux(ux1,ux2,ux3,u,i,j,ijm,dt,nbounx)
	dimension u(ijm)
	parameter(im=720,jm=121)
        common /ci/ c1(im,jm),c2(im,jm),c3(im,jm),c4(im,jm),&
                   c5(im,jm),c6(im,jm),c7(im,jm),c8(im,jm),c9(im,jm)
        common /grid/ xh(im),yh(jm),xz(im),yz(jm),dxh(im),dyh(jm),&
                     dxz(im),dyz(jm),wface,sface
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	dx=dxh(i)*c7(i,j)
	dudx=(u(ip)-u(iq))/2./dx
	dudx2=(u(ip)-2.*u(i)+u(iq))/dx/dx
	
	x1=u(i)*dt/dx
	x2=-u(i)*dt*dt*dudx/dx
	x3=u(i)*dt*dt*dt/dx*(dudx**2+u(i)/2.*dudx2)
	ux1=x1
	ux2=x1+x2
	ux3=x1+x2+x3
	return
	end
	
	subroutine muy(uy1,uy2,uy3,u,i,j,ijm,dt,nbouny)
	dimension u(ijm)
	parameter(im=720,jm=121)
        common /ci/ c1(im,jm),c2(im,jm),c3(im,jm),c4(im,jm),&
           c5(im,jm),c6(im,jm),c7(im,jm),c8(im,jm),c9(im,jm)
        common /grid/ xh(im),yh(jm),xz(im),yz(jm),dxh(im),dyh(jm),&
           dxz(im),dyz(jm),wface,sface

	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1
	dy2=dyh(j)**2
	dudy=(u(jp)-u(jq))/2./dyh(j)
	dudy2=(u(jp)-2.*u(j)+u(jq))/dy2
	
	y1=u(j)*dt/dyh(j)
	y2=-u(j)*dt*dt*dudy/dyh(j)
	y3=u(j)*dt*dt*dt/dyh(j)*(dudy**2+u(j)/2.*dudy2)
	uy1=y1
	uy2=y1+y2
	uy3=y1+y2+y3
	return
	end

	subroutine advecx(ada,a,i,j,ux1,ux2,ux3,n0)
	include 'pram1.inc'
	include 'com1.inc'
	dimension a(ijm)

!c	parameter(im=720,jm=121,ijm=im+jm)
!c        common /ci/ c1(im,jm),c2(im,jm),c3(im,jm),c4(im,jm),
!c     1              c5(im,jm),c6(im,jm),c7(im,jm),c8(im,jm),c9(im,jm)

!c        common /grid/ xh(im),yh(jm),xz(im),yz(jm),dxh(im),dyh(jm),
!c     1                dxz(im),dyz(jm),wface,sface
!c	common /var4/i1(jm),i2(jm),j1(im,2),j2(im,2)

	il1=i-1
	il2=i-2
	ir1=i+1
	ir2=i+2
	if(nbounx.eq.1.and.il1.lt.1) il1=im+il1
	if(nbounx.eq.1.and.il2.lt.1) il2=im+il2
	if(nbounx.eq.1.and.ir1.gt.im) ir1=ir1-im
	if(nbounx.eq.1.and.ir2.gt.im) ir2=ir2-im
	if(nbounx.eq.0.and.il1.lt.1) il1=1
	if(nbounx.eq.0.and.il2.lt.1) il2=1
	if(nbounx.eq.0.and.ir1.gt.im) ir1=im
	if(nbounx.eq.0.and.ir2.gt.im) ir2=im
	if(bx1(i,j).or.bx2(i,j))then
	 if(n0.eq.1)then
	  ada=a(i)-ux1*(a(i)-a(il1))
	 else
	  ada=a(i)-ux1*(a(ir1)-a(i))
	 endif
	else
	 a1=-ux3/2.*(a(ir1)-a(il1))
	 a2=ux3**2/2.*(a(ir1)-2.*a(i)+a(il1))
	 rn0=float(n0)*ux3*(1-ux3**2)/6.
	if(n0.eq.1) then
	 a3=rn0*(a(ir1)-3.*a(i)+3.*a(il1)-a(il2))
	else
	 a3=rn0*(a(il1)-3.*a(i)+3.*a(ir1)-a(ir2))
	end if
	 ada=a(i)+a1+a2+a3
	endif
	return
	end

	subroutine advecy(ada,a,i,j,uy1,uy2,uy3,n0)
	include 'pram1.inc'
	include 'com1.inc'
	dimension a(ijm)

!c	parameter(im=720,jm=121,ijm=im+jm)
!c        common /ci/ c1(im,jm),c2(im,jm),c3(im,jm),c4(im,jm),
!c     1              c5(im,jm),c6(im,jm),c7(im,jm),c8(im,jm),c9(im,jm)

!c        common /grid/ xh(im),yh(jm),xz(im),yz(jm),dxh(im),dyh(jm),
!c     1                dxz(im),dyz(jm),wface,sface
!c	common /var4/ i1(jm),i2(jm),j1(im,2),j2(im,2)

	jl1=j-1
	jl2=j-2
	jr1=j+1
	jr2=j+2
	if(nbouny.eq.1.and.jl1.lt.1) jl1=jm+jl1
	if(nbouny.eq.1.and.jl2.lt.1) jl2=jm+jl2
	if(nbouny.eq.1.and.jr1.gt.jm) jr1=jr1-jm
	if(nbouny.eq.1.and.jr2.gt.jm) jr2=jr2-jm
	if(nbouny.eq.0.and.jl1.lt.1) jl1=1
	if(nbouny.eq.0.and.jl2.lt.1) jl2=1
	if(nbouny.eq.0.and.jr1.gt.jm) jr1=jm
	if(nbouny.eq.0.and.jr2.gt.jm) jr2=jm
	if(by1(i,j).or.by2(i,j))then
	 if(n0.eq.1)then
	  ada=a(j)-uy1*(a(j)-a(jl1))
	 else
	  ada=a(j)-uy1*(a(jr1)-a(j))
	 endif
	else
	 a1=-uy3/2.*(a(jr1)-a(jl1))
	 a2=uy3**2/2.*(a(jr1)-2.*a(j)+a(jl1))
	 rn0=float(n0)*uy3*(1-uy3**2)/6.
	if(n0.eq.1) then
	 a3=rn0*(a(jr1)-3.*a(j)+3.*a(jl1)-a(jl2))
	else
	 a3=rn0*(a(jl1)-3.*a(j)+3.*a(jr1)-a(jr2))
	end if
	 ada=a(j)+a1+a2+a3
	endif
	return
	end
	

	subroutine bduv(m,nbounx,nbouny)
	include 'pram1.inc'
	common /var0/u2(im,jm,2),v2(im,jm,2),h2(im,jm,2),&
            u1(im,jm,2),v1(im,jm,2),h1(im,jm,2),tmp(im,jm,2)
	logical bx1,bx2,by1,by2
	common /var4/bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)
	logical be
	common /logi/be(im,jm)

	do j=1,jm
	do i=1,im
	if(bx1(i,j)) then
	 do i0=i,i+5
	  r=0.25*float(6+i-i0)/6.
	  ii=i0
	  if(nbounx.eq.1.and.ii.gt.im) ii=ii-im
	  if(nbounx.eq.0.and.ii.gt.im) ii=im
	if(be(ii,j)) then
	ip=ii+1
	iq=ii-1
	  if(nbounx.eq.1.and.ii.eq.im) ip=1
	  if(nbounx.eq.1.and.ii.eq.1) iq=im
	  if(nbounx.eq.0.and.ip.gt.im) ip=im
	  if(nbounx.eq.0.and.iq.lt.1) iq=1

	  u2(ii,j,m)=u2(ii,j,m)+r*(u2(iq,j,m)-2.*u2(ii,j,m)+u2(ip,j,m))
	  v2(ii,j,m)=v2(ii,j,m)+r*(v2(iq,j,m)-2.*v2(ii,j,m)+v2(ip,j,m))
	  u1(ii,j,m)=u1(ii,j,m)+r*(u1(iq,j,m)-2.*u1(ii,j,m)+u1(ip,j,m))
	  v1(ii,j,m)=v1(ii,j,m)+r*(v1(iq,j,m)-2.*v1(ii,j,m)+v1(ip,j,m))
	else
	 goto 1
	end if
	 enddo
	endif

1	continue
	if(bx2(i,j)) then
	 do i0=i,i-5,-1
	  r=0.25*float(6+i0-i)/6.
	  ii=i0
	  if(nbounx.eq.1.and.ii.lt.1) ii=im+ii
	  if(nbounx.eq.0.and.ii.lt.1) ii=1
	if(be(ii,j)) then
	ip=ii+1
	iq=ii-1
	  if(nbounx.eq.1.and.ii.eq.im) ip=1
	  if(nbounx.eq.1.and.ii.eq.1) iq=im
	  if(nbounx.eq.0.and.ip.gt.im) ip=im
	  if(nbounx.eq.0.and.iq.lt.1) iq=1
	  u2(ii,j,m)=u2(ii,j,m)+r*(u2(iq,j,m)-2.*u2(ii,j,m)+u2(ip,j,m))
	  v2(ii,j,m)=v2(ii,j,m)+r*(v2(iq,j,m)-2.*v2(ii,j,m)+v2(ip,j,m))
	  u1(ii,j,m)=u1(ii,j,m)+r*(u1(iq,j,m)-2.*u1(ii,j,m)+u1(ip,j,m))
	  v1(ii,j,m)=v1(ii,j,m)+r*(v1(iq,j,m)-2.*v1(ii,j,m)+v1(ip,j,m))
	else
	goto 2
	end if
	 enddo
	endif
2       continue
	enddo
	enddo

	jbd=1
	if(jbd.eq.1)then
	do i=1,im
	do j=1,jm

	if(by1(i,j)) then
	 do j0=j,j+5
	  r=0.25*float(6+j-j0)/6.
	  jj=j0
	  if(nbouny.eq.1.and.jj.gt.jm) jj=jj-jm
	  if(nbouny.eq.0.and.jj.gt.jm) jj=jm
	if(be(i,jj)) then
	jp=jj+1
	jq=jj-1
	  if(nbouny.eq.1.and.jj.eq.jm) jp=1
	  if(nbouny.eq.1.and.jj.eq.1) jq=jm
	  if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	  if(nbouny.eq.0.and.jq.lt.1) jq=1
	   u2(i,jj,m)=u2(i,jj,m)+r*(u2(i,jq,m)-2.*u2(i,jj,m)+u2(i,jp,m))
	   v2(i,jj,m)=v2(i,jj,m)+r*(v2(i,jq,m)-2.*v2(i,jj,m)+v2(i,jp,m))
	   u1(i,jj,m)=u1(i,jj,m)+r*(u1(i,jq,m)-2.*u1(i,jj,m)+u1(i,jp,m))
	   v1(i,jj,m)=v1(i,jj,m)+r*(v1(i,jq,m)-2.*v1(i,jj,m)+v1(i,jp,m))
	else
	 goto 3
	end if
	  enddo
	endif
3       continue

	if(by2(i,j)) then
	 do j0=j,j-5,-1
	  r=0.25*float(6+j0-j)/6.
	  jj=j0
	  if(nbouny.eq.1.and.jj.lt.1) jj=jm+jj
	  if(nbouny.eq.0.and.jj.lt.1) jj=1
	if(be(i,jj)) then
	jp=jj+1
	jq=jj-1
	  if(nbouny.eq.1.and.jj.eq.jm) jp=1
	  if(nbouny.eq.1.and.jj.eq.1) jq=jm
	  if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	  if(nbouny.eq.0.and.jq.lt.1) jq=1
	   u2(i,jj,m)=u2(i,jj,m)+r*(u2(i,jq,m)-2.*u2(i,jj,m)+u2(i,jp,m))
	   v2(i,jj,m)=v2(i,jj,m)+r*(v2(i,jq,m)-2.*v2(i,jj,m)+v2(i,jp,m))
	   u1(i,jj,m)=u1(i,jj,m)+r*(u1(i,jq,m)-2.*u1(i,jj,m)+u1(i,jp,m))
	   v1(i,jj,m)=v1(i,jj,m)+r*(v1(i,jq,m)-2.*v1(i,jj,m)+v1(i,jp,m))
	else
	 goto 4
	end if
	  enddo
	endif
4       continue
	enddo
	enddo

	endif

	return
	end


	subroutine bdht(h,m,nbounx,nbouny)
	include 'pram1.inc'
	dimension h(im,jm,2)
	logical bx1,bx2,by1,by2
	common /var4/bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)
	logical be
	common /logi/be(im,jm)

	do j=1,jm
	do i=1,im
	if(bx1(i,j)) then
	 do i0=i,i+5
	  r=0.25*float(6+i-i0)/6.
	  ii=i0
	  if(nbounx.eq.1.and.ii.gt.im) ii=ii-im
	  if(nbounx.eq.0.and.ii.gt.im) ii=im
	if(be(ii,j)) then
	ip=ii+1
	iq=ii-1
	  if(nbounx.eq.1.and.ii.eq.im) ip=1
	  if(nbounx.eq.1.and.ii.eq.1) iq=im
	  if(nbounx.eq.0.and.ip.gt.im) ip=im
	  if(nbounx.eq.0.and.iq.lt.1) iq=1

	  h(ii,j,m)=h(ii,j,m)+r*(h(iq,j,m)-2.*h(ii,j,m)+h(ip,j,m))
	else
	 goto 1
	end if
	 enddo
	endif

1	continue
	if(bx2(i,j)) then
	 do i0=i,i-5,-1
	  r=0.25*float(6+i0-i)/6.
	  ii=i0
	  if(nbounx.eq.1.and.ii.lt.1) ii=im+ii
	  if(nbounx.eq.0.and.ii.lt.1) ii=1
	if(be(ii,j)) then
	ip=ii+1
	iq=ii-1
	  if(nbounx.eq.1.and.ii.eq.im) ip=1
	  if(nbounx.eq.1.and.ii.eq.1) iq=im
	  if(nbounx.eq.0.and.ip.gt.im) ip=im
	  if(nbounx.eq.0.and.iq.lt.1) iq=1
	  h(ii,j,m)=h(ii,j,m)+r*(h(iq,j,m)-2.*h(ii,j,m)+h(ip,j,m))
	else
	goto 2
	end if
	 enddo
	endif
2       continue
	enddo
	enddo

	jbd=1
	if(jbd.eq.1)then
	do i=1,im
	do j=1,jm

	if(by1(i,j)) then
	 do j0=j,j+5
	  r=0.25*float(6+j-j0)/6.
	  jj=j0
	  if(nbouny.eq.1.and.jj.gt.jm) jj=jj-jm
	  if(nbouny.eq.0.and.jj.gt.jm) jj=jm
	if(be(i,jj)) then
	jp=jj+1
	jq=jj-1
	  if(nbouny.eq.1.and.jj.eq.jm) jp=1
	  if(nbouny.eq.1.and.jj.eq.1) jq=jm
	  if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	  if(nbouny.eq.0.and.jq.lt.1) jq=1
	   h(i,jj,m)=h(i,jj,m)+r*(h(i,jq,m)-2.*h(i,jj,m)+h(i,jp,m))
	else
	 goto 3
	end if
	  enddo
	endif
3       continue

	if(by2(i,j)) then
	 do j0=j,j-5,-1
	  r=0.25*float(6+j0-j)/6.
	  jj=j0
	  if(nbouny.eq.1.and.jj.lt.1) jj=jm+jj
	  if(nbouny.eq.0.and.jj.lt.1) jj=1
	if(be(i,jj)) then
	jp=jj+1
	jq=jj-1
	  if(nbouny.eq.1.and.jj.eq.jm) jp=1
	  if(nbouny.eq.1.and.jj.eq.1) jq=jm
	  if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	  if(nbouny.eq.0.and.jq.lt.1) jq=1
	   h(i,jj,m)=h(i,jj,m)+r*(h(i,jq,m)-2.*h(i,jj,m)+h(i,jp,m))
	else
	 goto 4
	end if
	  enddo
	endif
4       continue
	enddo
	enddo

	endif

	return
	end


	subroutine bound0(u,m,nbounx,nbouny)
	include 'pram1.inc'
	dimension u(im,jm,2)
	logical bx1,bx2,by1,by2
	common /var4/bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)

	do j=1,jm
	do i=1,im
	if(bx1(i,j))then
	i0=i-1
	if(nbounx.eq.1.and.i.eq.1) i0=im
	u(i0,j,m)=u(i,j,m)
	endif

	if(bx2(i,j))then
	i0=i+1
	if(nbounx.eq.1.and.i.eq.im) i0=1
	u(i0,j,m)=u(i,j,m)
	endif

	if(by1(i,j))then
	j0=j-1
	if(nbouny.eq.1.and.j.eq.1) j0=jm
	u(i,j0,m)=u(i,j,m)
	endif

	if(by2(i,j))then
	j0=j+1
	if(nbouny.eq.1.and.j.eq.jm) j0=1
	u(i,j0,m)=u(i,j,m)
	endif
	enddo
	enddo
	return
	end
	

	subroutine bound1(u,nbounx,nbouny)
	include 'pram1.inc'
	dimension u(im,jm)
	logical bx1,bx2,by1,by2
	common /var4/bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)

	do j=1,jm
	do i=1,im
	if(bx1(i,j))then
	i0=i-1
	if(nbounx.eq.1.and.i.eq.1) i0=im
	u(i0,j)=u(i,j)
	endif

	if(bx2(i,j))then
	i0=i+1
	if(nbounx.eq.1.and.i.eq.im) i0=1
	u(i0,j)=u(i,j)
	endif

	if(by1(i,j))then
	j0=j-1
	if(nbouny.eq.1.and.j.eq.1) j0=jm
	u(i,j0)=u(i,j)
	endif

	if(by2(i,j))then
	j0=j+1
	if(nbouny.eq.1.and.j.eq.jm) j0=1
	u(i,j0)=u(i,j)
	endif
	enddo
	enddo
	return
	end
	

	subroutine adjust
	include 'pram1.inc'
!c	logical be(im,jm)
	include 'com1.inc'
	real divuv(im,jm),divuv1(im,jm)
!c
!c===>Update h2,h1 first
!c
	do j=1,jm
	do i=1,im
	upw(i,j)=0.

	if(be(i,j))then
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1

	 div=-(u2(ip,j,n)-u2(iq,j,n))/2.*c8(i,j)&
          -(v2(i,jp,n)*c3(i,jp)-v2(i,jq,n)*c3(i,jq))/2.*c9(i,j)
	 hln=alog(h2(i,j,np))
	 h2(i,j,np)=hln+dt*(adhxy(i,j)+div)
	 h2(i,j,np)=exp(h2(i,j,np))
	 if(h2(i,j,np).le.h2min)h2(i,j,np)=h2min

	 div=-(u1(ip,j,n)-u1(iq,j,n))/2.*c8(i,j)&
          -(v1(i,jp,n)*c3(i,jp)-v1(i,jq,n)*c3(i,jq))/2.*c9(i,j)
	 h1ln=alog(h1(i,j,np))
	 h1(i,j,np)=h1ln+dt*(adh1xy(i,j)+div)
	 h1(i,j,np)=exp(h1(i,j,np))
!c======get the upwelling as a parameter for the
!c entrainment depth (Fu 8/26/99)
	upw(i,j)=-div*h1(i,j,np)
!c===>
	 if(h1(i,j,np).le.h1min)h1(i,j,np)=h1min
	 if(h1(i,j,np).ge.h1max)h1(i,j,np)=h1max
	
	end if
	end do
	end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	call bound0(h1,2,nbounx,nbouny)
	call bound0(h2,2,nbounx,nbouny)
!c
!c===>Now update u,v,u1,v1 with using new h, h1
!c
	do j=1,jm
	do i=1,im
	if(be(i,j))then
!cfu (tmn=14.)
	b1(i,j)=thmexp*grav*(tmp(i,j,n)-tmn)
!c	b2(i,j)=thmexp*grav*((tmp(i,j,n)+tmn)/2.-tmn)
        b2(i,j)=b1(i,j)/2.
!c        b(i,j)=b2(i,j)*(1.+h1(i,j,np)/h(i,j,np))
!c	 b(i,j)=(b1(i,j)*h1(i,j,np)+b2(i,j)*(h(i,j,np)-h1(i,j,np)))
!c     1    /h(i,j,np)

!cZF using more precise buoyancy expression
!c        b(i,j)=2./3.*b2(i,j)*
!c     1  (1+h1(i,j,np)/h(i,j,np)+(h1(i,j,np)/h(i,j,np))**2)
	endif
	enddo
	enddo

	call bound1(b1,nbounx,nbouny)
	call bound1(b2,nbounx,nbouny)
!c	call bound1(b,nbounx,nbouny)
	do j=1,jm
	do i=1,im
	if(be(i,j))then
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1

	bhl=b2(iq,j)*(h1(iq,j,np)+h2(iq,j,np))
	bhr=b2(ip,j)*(h1(ip,j,np)+h2(ip,j,np))

	ha12=h1(i,j,np)+h2(i,j,np)/2.
        pgx=(-(bhr-bhl)+ha12*(b2(ip,j)-b2(iq,j)))*c5(i,j)/2.
	
	rx=dt*pgx+u2(i,j,n)+dt*ff(j)/2.*v2(i,j,n)
	bhl=b2(i,jq)*(h1(i,jq,np)+h2(i,jq,np))
	bhr=b2(i,jp)*(h1(i,jp,np)+h2(i,jp,np))
        pgy=(-(bhr-bhl)+ha12*(b2(i,jp)-b2(i,jq)))*c6(i,j)/2.
	ry=dt*pgy+v2(i,j,n)-dt*ff(j)/2.*u2(i,j,n)
	delta=1+(dt/2.*ff(j))**2
	u2(i,j,np)=(rx+dt/2.*ff(j)*ry)/delta
	v2(i,j,np)=(ry-dt/2.*ff(j)*rx)/delta

	endif
	enddo
	enddo

!c
!c===> Now update u1,v1
!c

	do j=1,jm
	do i=1,im
	if(be(i,j))then
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1

	bhl=b1(iq,j)*(2.*h1(iq,j,np)+h2(iq,j,np))
	bhr=b1(ip,j)*(2.*h1(ip,j,np)+h2(ip,j,np))
        pgx=(-(bhr-bhl)+h1(i,j,np)*(b1(ip,j)-b1(iq,j)))*c5(i,j)/4.
	rx=dt*pgx+u1(i,j,n)+dt*ff(j)/2.*v1(i,j,n)
	bhl=b1(i,jq)*(2.*h1(i,jq,np)+h2(i,jq,np))
	bhr=b1(i,jp)*(2.*h1(i,jp,np)+h2(i,jp,np))
        pgy=(-(bhr-bhl)+h1(i,j,np)*(b1(i,jp)-b1(i,jq)))*c6(i,j)/4.
	ry=dt*pgy+v1(i,j,n)-dt*ff(j)/2.*u1(i,j,n)
	delta=1+(dt/2.*ff(j))**2
	u1(i,j,np)=(rx+dt/2.*ff(j)*ry)/delta
	v1(i,j,np)=(ry-dt/2.*ff(j)*rx)/delta
	tmp(i,j,np)=tmp(i,j,n)

	endif
	enddo
	enddo
!c
!c divergence damping
!c
	ndiv=0
	if(ndiv.eq.1)then
	alphad=1.e8
	do j=1,jm
	do i=1,im
	if(be(i,j))then
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1
        divuv(i,j)=(u2(ip,j,np)-u2(iq,j,np))/2.*c8(i,j)+&
         (v2(i,jp,np)*c3(i,jp)-v2(i,jq,np)*c3(i,jq))/2.*c9(i,j)
	divuv1(i,j)=(u1(ip,j,np)-u1(iq,j,np))/2.*c8(i,j)+&
         (v1(i,jp,np)*c3(i,jp)-v1(i,jq,np)*c3(i,jq))/2.*c9(i,j)
	endif
	enddo
	enddo

	do j=1,jm
	do i=1,im
	if(be(i,j))then
	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1
	u2(i,j,np)=u2(i,j,np)+&
         dt*alphad*(divuv(ip,j)-divuv(iq,j))/2.*c8(i,j)
	v2(i,j,np)=v2(i,j,np)+&
         dt*alphad*(divuv(i,jp)-divuv(i,jq))/2.*c9(i,j)
	u1(i,j,np)=u1(i,j,np)+&
         dt*alphad*(divuv1(ip,j)-divuv1(iq,j))/2.*c8(i,j)
	v1(i,j,np)=v1(i,j,np)+&
         dt*alphad*(divuv1(i,jp)-divuv1(i,jq))/2.*c9(i,j)
	endif
	enddo
	enddo
	endif
	return
	end

	subroutine physics(sst_nudg)

        USE mo_control,        ONLY: lonudg
	include 'pram1.inc'
	include 'com1.inc'
        real sst_nudg(im,jm)
!c	logical be(im,jm)
!c--=7(86)
	umum=25.          !shear mixing
	umu=0.          !ML friction
	hk=0.003

	do j=1,jm

!c----within 2-degree latitude (Chang 1993)
!c	anml=float(j-31)*float(j-31)/16.           !31=equator

	do i=1,im
	if(be(i,j))then
        call diff(diss,h1,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
	damp=-fe(i,j)*(h1(i,j,np)-eqdpth1)
	h1(i,j,n)=h1(i,j,np)+dt*(diss+we(i,j)+damp)
!c	h1(i,j,n)=eqdpth1+(h1(i,j,np)-eqdpth1+dt*diss)/(1+dt*fe(i,j))
!c	h1(i,j,n)=(h1(i,j,np)+dt*(diss+fe(i,j)*eqdpth1))/(1+dt*fe(i,j))
	if(h1(i,j,n).le.h1min)h1(i,j,n)=h1min
	if(h1(i,j,n).ge.h1max)h1(i,j,n)=h1max

        call diff(diss,h2,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
	damp=-fe(i,j)*(h2(i,j,np)-eqdpth)
	h2(i,j,n)=h2(i,j,np)+dt*(diss-we(i,j)+damp)
!c	h2(i,j,n)=eqdpth+(h2(i,j,np)-eqdpth+dt*diss)/(1+dt*fe(i,j))
!c	h2(i,j,n)=(h2(i,j,np)+dt*(diss+fe(i,j)*eqdpth))/(1+dt*fe(i,j))
	if(h2(i,j,n).le.h2min)h2(i,j,n)=h2min

!cffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
	ww=we(i,j)
	if(ww.ge.0.) ww=0.
        call diff(diss,u2,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
        rfu=-fr(i,j)*u2(i,j,np)
!c        fin=-txyk*u2(i,j,np)*sqrt(u2(i,j,np)**2+v2(i,j,np)**2)
!cfu
	fin=-0.*u2(i,j,np)/h2(i,j,np)
	rhand=fin/h2(i,j,np)+rfu+diss-ww*us(i,j)/h2(i,j,np)
	u2(i,j,n)=u2(i,j,np)+dt*rhand
!c	rhand=fin/h2(i,j,np)+diss
!c	u2(i,j,n)=(u2(i,j,np)+dt*rhand)/(1.+fr(i,j)*dt)

        call diff(diss,v2,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
        rfv=-fr(i,j)*v2(i,j,np)
!c        fin=-txyk*v2(i,j,np)*sqrt(u2(i,j,np)**2+v2(i,j,np)**2)
!cfu
	fin=-0.*v2(i,j,np)/h2(i,j,np)
	rhand=fin/h2(i,j,np)+rfv+diss-ww*vs(i,j)/h2(i,j,np)
	v2(i,j,n)=v2(i,j,np)+dt*rhand
!c	rhand=fin/h2(i,j,np)+diss
!c	v2(i,j,n)=(v2(i,j,np)+dt*rhand)/(1.+fr(i,j)*dt)
!c------------------------------------------------------------------
	um1=umum
	um2=umu
!c	um1=1.+umum*exp(-anml)
!c	um2=1.+umu*exp(-anml)
!c------------------------------------------------------------------
        call diff(diss,u1,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
        rfu=-fr(i,j)*u1(i,j,np)
!c        fin=-txyk1*us(i,j)*sqrt(us(i,j)**2+vs(i,j)**2)
        fin=-um1*us(i,j)/h1(i,j,np)-&
         um2*u1(i,j,np)/h1(i,j,np)           !fu
	ww=we(i,j)
	if(ww.lt.0)ww=0
	rhand=(tx(i,j)+fin)/h1(i,j,np)+rfu+diss-ww*us(i,j)/h1(i,j,np)
	u1(i,j,n)=u1(i,j,np)+dt*rhand
!c	rhand=(tx(i,j)+fin)/h1(i,j,np)+diss-ww*us(i,j)/h1(i,j,np)
!c	u1(i,j,n)=(u1(i,j,np)+dt*rhand)/(1.+fr(i,j)*dt)

        call diff(diss,v1,i,j,np,nbounx,nbouny)
        diss=cdfv*visc(i,j)*diss
        rfv=-fr(i,j)*v1(i,j,np)
!c        fin=-txyk1*vs(i,j)*sqrt(us(i,j)**2+vs(i,j)**2)
        fin=-um1*vs(i,j)/h1(i,j,np)-&
          um2*v1(i,j,np)/h1(i,j,np)         !fu
	rhand=(ty(i,j)+fin)/h1(i,j,np)+rfv+diss-ww*vs(i,j)/h1(i,j,np)
	v1(i,j,n)=v1(i,j,np)+dt*rhand
!c	rhand=(ty(i,j)+fin)/h1(i,j,np)+diss-ww*vs(i,j)/h1(i,j,np)
!c	v1(i,j,n)=(v1(i,j,np)+dt*rhand)/(1.+fr(i,j)*dt)

        call diff(diss,tmp,i,j,np,nbounx,nbouny)
        tdf(i,j)=dffs(i,j)*diss
!c	tdf(i,j)=tdf(i,j)-hk*tse(i,j)/h1(i,j,np) !mixing
!c     1      /h1(i,j,np)

        if(lonudg.and.(sst_nudg(i,j).gt.10.)) then                                          !fu
!20090820_ctl	tnd(i,j)=-(tmp(i,j,np)-sst_nudg(i,j))/(3.*24.*3600.) !three-day relaxation
	tnd(i,j)=-(tmp(i,j,np)-sst_nudg(i,j))/(.2*24.*3600.) !5-hrs relaxation
!c=====================================================>
!c(1/7/2000)fu++relax tmp=>sst in the eastern coast!!
!c with half-day relaxing time
!fu+3/25/2002	tnd(i,j)=tnd(i,j)-fe1(i,j)*(tmp(i,j,np)-sst(i,j))
!======================================================>
!fu+1/16/2003    heat-flux-correction
!      flx(i,j)=-(tmp(i,j,np)-sst(i,j))/(2.*24.*3600.)
!c
        else
!relaxing model-ocean-boundary SST to observed climatological sst
	tnd(i,j)=-fn(i,j)*(tmp(i,j,np)-sst(i,j))
        end if

	ww=wet(i,j)
	if(ww.lt.0)  ww=0.
	zah(i,j)=-ww/h1(i,j,np)-hk*tseh(i,j)/h1(i,j,np)&
           /h1(i,j,np)

	tfc(i,j)=qnet(i,j)/cw/h1(i,j,np)
!c	tfc(i,j)=qnet(i,j)
!c  /1. !too warm   /2. !too cold
!	fuflx=flx(i,j)/1.25
!	if(flx(i,j).ge.1.30e-6) fuflx=1.30e-6/1.25
	rhand=tfc(i,j)+zah(i,j)+tnd(i,j)+tdf(i,j)
	tmp(i,j,n)=tmp(i,j,np)+dt*rhand
	
	endif
	enddo
	enddo
	return
	end

	subroutine diff(diss,h,i,j,n,nbounx,nbouny)
	parameter(im=720,jm=121)
	dimension h(im,jm,2)
        common /ci/ c1(im,jm),c2(im,jm),c3(im,jm),c4(im,jm),&
           c5(im,jm),c6(im,jm),c7(im,jm),c8(im,jm),c9(im,jm)

	ip=i+1
	if(nbounx.eq.1.and.ip.gt.im) ip=1
	if(nbounx.eq.0.and.ip.gt.im) ip=im
	iq=i-1
	if(nbounx.eq.1.and.iq.lt.1) iq=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	jp=j+1
	if(nbouny.eq.1.and.jp.gt.jm) jp=1
	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	jq=j-1
	if(nbouny.eq.1.and.jq.lt.1) jq=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1
        hx1=(h(i,j,n)-h(iq,j,n))*c5(i,j)
        hx2=(h(ip,j,n)-h(i,j,n))*c5(ip,j)
        hy1=(h(i,j,n)-h(i,jq,n))*c2(i,j)*c7(i,j)
        hy2=(h(i,jp,n)-h(i,j,n))*c2(i,jp)*c7(i,jp)
        diss=(hx2-hx1)*c8(i,j)+(hy2-hy1)*c9(i,j)
	return
	end

	subroutine filter0(u,be,nbounx)
	parameter(im=720,jm=121,fw=0.1)
	real help(im,jm),u(im,jm)
	logical be(im,jm)

	do 5 j=1,jm
	do 5 i=1,im
	 help(i,j)=u(i,j)
5	continue
	call shapiro(help,be,nbounx)
	do 6 j=1,jm
	do 6 i=1,im
	if(be(i,j)) then
         u(i,j)=(1.0-fw)*u(i,j)+fw*help(i,j)
	end if
6	continue
	return
	end

	subroutine filterm(u,be,n,nbounx)
	parameter(im=720,jm=121,fw=0.1)
	real help(im,jm),u(im,jm,2)
	logical be(im,jm)
	do 5 j=1,jm
	do 5 i=1,im
	 help(i,j)=u(i,j,n)
5	continue
	call shapiro(help,be,nbounx)
	do 6 j=1,jm
	do 6 i=1,im
	if(be(i,j)) then
         u(i,j,n)=(1.0-fw)*u(i,j,n)+fw*help(i,j)
	end if
6	continue
	return
	end

	subroutine shapiro(u,be,nbounx)
	parameter(im=720,jm=121)
	real u(im,jm),a(im,jm)
	logical be(im,jm)
	logical bx1,bx2,by1,by2
	common /var4/bx1(im,jm),bx2(im,jm),by1(im,jm),by2(im,jm)
	r=1./16.
	do 3 j=1,jm
	do 3 i=1,im
	a(i,j)=u(i,j)
3	continue

	do 5 j=1,jm
	do 5 i=1,im

	if(be(i,j))then
	ip=i+1
	ip2=i+2
	iq=i-1
	iq2=i-2
	if(nbounx.eq.1.and.ip.gt.im) ip=ip-im
	if(nbounx.eq.1.and.ip2.gt.im) ip2=ip2-im

	if(nbounx.eq.1.and.iq.lt.1) iq=iq+im
	if(nbounx.eq.1.and.iq2.lt.1) iq2=iq2+im

	if(nbounx.eq.0.and.ip.gt.im) ip=im
	if(nbounx.eq.0.and.ip2.gt.im) ip2=im
	if(nbounx.eq.0.and.iq.lt.1) iq=1
	if(nbounx.eq.0.and.iq2.lt.1) iq2=1

	if(bx1(i,j).or.bx2(i,j)) then
        a(i,j)=0.25*(u(iq,j)+2.*u(i,j)+u(ip,j))
	else
	a(i,j)=r*(-u(iq2,j)+4.*u(iq,j)+10.*u(i,j)+4.*u(ip,j)-u(ip2,j))
	 end if
	end if

5	continue

	call bound1(a,nbounx,nbouny)

	do 6 i=1,im
	do 6 j=1,jm
	if(be(i,j)) then

	jp=j+1
	jp2=j+2
	jq=j-1
	jq2=j-2
	if(nbouny.eq.1.and.jp.gt.jm) jp=jp-jm
	if(nbouny.eq.1.and.jp2.gt.jm) jp2=jp2-jm

	if(nbouny.eq.1.and.jq.lt.1) jq=jq+jm
	if(nbouny.eq.1.and.jq2.lt.1) jq2=jq2+jm

	if(nbouny.eq.0.and.jp.gt.jm) jp=jm
	if(nbouny.eq.0.and.jp2.gt.jm) jp2=jm
	if(nbouny.eq.0.and.jq.lt.1) jq=1
	if(nbouny.eq.0.and.jq2.lt.1) jq2=1
	if(by1(i,j).or.by2(i,j)) then
	u(i,j)=0.25*(a(i,jq)+2.*a(i,j)+a(i,jp))
	else
	u(i,j)=r*(-a(i,jq2)+4.*a(i,jq)+10.*a(i,j)+4.*a(i,jp)-a(i,jp2))
	 end if
	end if	
6	continue
	return
	end

!c============
	subroutine entrain(wfu)
        USE mo_doctor,   ONLY:nout

	include 'pram1.inc'
	include 'com1.inc'
	external funcwe
	dimension wfu(im,jm)
!c	logical be(im,jm)
	data hr/1700./,rr/0.35/

        do 10 j=1,jm

	anml=float(j-31)*float(j-31)/16.

        do 10 i=1,im
!c	wei=6500./h1(i,j,np)*exp(-anml)
!c	if(wei.lt.1.0) wei=1.0
	wei=1.

	wmy(i,j)=wmy1(j)*wei
        if(be(i,j)) then
!c        vvvv=sqrt(tu(i,j)**2+tv(i,j)**2)
!c       if(vvvv.lt.3.) vvvv=3.	!!!  m/s(for transient)
!ctim(fu)
	vvvv=wfu(i,j)
        uf(i,j)=ufcoef*vvvv
        wndstr=2.*wmy(i,j)*uf(i,j)**3
        bckdss=2.*hepsln*h1(i,j,np)
!c        htflux=thmexp*grav*(qin(i,j)+ef*qout(i,j))/cw
!cfu Gaspar(1988)
 	pp=exp(-h1(i,j,np)/hr)-2.*hr*(1.-exp(-h1(i,j,np)/hr))&
          /h1(i,j,np)	

        htflux=thmexp*grav*(qin(i,j)+rr*qout(i,j)*pp)/cw
        turbvl=amax1(3.,3.*uf(i,j))
!cppp
!c        petwe=-(h1(i,j,np)-2./gamma2)*thmexp*grav/cw*qout(i,j)*(1.-r)
	petwe=0.
        wup=wndstr-bckdss-&
           0.5*h1(i,j,np)*((1.-tn)*abs(htflux)+(1.+tn)*htflux)+&
           petwe
!c
!c      term of wme*(us**2+vs**2) is not concerned now !
!c
       if(wup.le.0.) then
        hgus1=h1(i,j,np)
        hgus2=0.1
        nh12=int(h1(i,j,np)/100.)+1
        nbh12=1
!c	 write(6,*)'i=',i,' j=',j
        call ZBRAK(funcwe,hgus1,hgus2,nh12,hrs1,hrs2,nbh12,i,j)
         hrs=0.5*(hrs1+hrs2)
         if(hrs.lt.h1min) hrs=h1min		!!!
         if(nbh12.eq.0) hrs=h1(i,j,np)
         hstar(i,j)=hrs

!c         h1u=wndstr+2.*thmexp*grav/cw/gamma2*qout(i,j)*(1.-r)
!c	 h1d=0.5*((1.-tn)*abs(htflux)+(1.+tn)*htflux)+
!c     1       thmexp*grav/cw*qout(i,j)*(1.-r)
!c         hstar(i,j)=h1u/h1d
       endif
!c
!cfu
	ddt=tse(i,j)
        wdown=h1(i,j,np)*thmexp*grav*ddt
         if (wdown.le.0.) then
          write(nout,*)'wdown<=0 !!!  temp < tmn ?'
          write(nout,*) itt,day,wdown,tmp(i,j,np),tse(i,j),i,j
	  call out_put
          stop
         end if
!c
         if(wup.gt.0) then
            we(i,j)=wup/wdown
	    wet(i,j)=we(i,j)*tseh(i,j)
!c	if((tmp(i,j,np)-te(i,j)).lt.0.) then
!c	    wet(i,j)=wet(i,j)/2. !when tmp(i,j,np) <te(i,j)+0.0
!c	end if
!ctim
	    if(we(i,j).gt.welmt) we(i,j)=welmt	      !!! cm/s
!c
         else
!c
!c	 note: arbitary adjustment time of 2*dt2
!c		we<0 only affects h1 next step
!c
            we(i,j)=(hstar(i,j)-h1(i,j,np))/2./dt
	    welmtn=-welmt/5.
	    if(we(i,j).lt.welmtn) we(i,j)=welmtn	      !!! cm/s
	    wet(i,j)=0.
         end if
         end if
10     continue
	return
	end

!c================
	subroutine ijmp(im1,jm1,ip1,jp1,i,j,im,jm)
	im1=i-1
	jm1=j-1
	ip1=i+1
	jp1=j+1
	if(i.eq.1)im1=1
	if(j.eq.1)jm1=1
	if(i.eq.im)ip1=im
	if(j.eq.jm)jp1=jm
	return
	end

!c=================
	subroutine robert(u,bu,im,jm,nm,n,np,s)
	dimension u(im,jm,3)
	logical bu(im,jm)
        do j=1,jm
        do i=1,im
	if(bu(i,j))then
          u(i,j,n)=u(i,j,n)+0.5*s*(u(i,j,np)-2.*u(i,j,n)+u(i,j,nm))
	endif
        enddo
        enddo
	return
	end

!c=================
	subroutine out_put
	include 'pram1.inc'
	include 'com1.inc'
!c        write(134) ((u2(l,m,n),l=1,im),m=1,jm)
!c        write(135) ((v2(l,m,n),l=1,im),m=1,jm)
!c        write(136) ((h2(l,m,n),l=1,im),m=1,jm)
!c        write(137) ((tmp(l,m,n),l=1,im),m=1,jm)
!c        write(138) ((we(l,m),l=1,im),m=1,jm)
!c        write(139) ((te(l,m),l=1,im),m=1,jm)
!c        write(140) ((xah(l,m),l=1,im),m=1,jm)
!c        write(141) ((yah(l,m),l=1,im),m=1,jm)
!c        write(142) ((zah(l,m),l=1,im),m=1,jm)
!c        write(143) ((tfc(l,m),l=1,im),m=1,jm)
!c        write(144) ((tu(l,m),l=1,im),m=1,jm)
!c        write(145) ((tv(l,m),l=1,im),m=1,jm)
!c        write(146) ((u1(l,m,n),l=1,im),m=1,jm)
!c        write(147) ((v1(l,m,n),l=1,im),m=1,jm)
!c        write(148) ((h1(l,m,n),l=1,im),m=1,jm)
	return
	end

!c=============================
	subroutine cycle
	include 'pram1.inc'
	include 'com1.inc'
!c
!c
!c******************************************************************
!c
!c    impose cyclic boundary condition
!c
      if (nbouny.eq.1) then
       do 4670 i=1,im
         u2(i,1,n)=u2(i,jm-1,n)
         u2(i,jm,n)=u2(i,2,n)
         v2(i,1,n)=v2(i,jm-1,n)
         v2(i,jm,n)=v2(i,2,n)
         h2(i,1,n)=h2(i,jm-1,n)
         h2(i,jm,n)=h2(i,2,n)
         u1(i,1,n)=u1(i,jm-1,n)
         u1(i,jm,n)=u1(i,2,n)
         v1(i,1,n)=v1(i,jm-1,n)
         v1(i,jm,n)=v1(i,2,n)
         h1(i,1,n)=h1(i,jm-1,n)
         h1(i,jm,n)=h1(i,2,n)
         tmp(i,1,n)=tmp(i,jm-1,n)
         tmp(i,jm,n)=tmp(i,2,n)
 4670  continue
       end if
       if (nbounx.eq.1) then
       do 4680 j=1,jm
         u2(1,j,n)=u2(im-1,j,n)
         u2(im,j,n)=u2(2,j,n)
         v2(1,j,n)=v2(im-1,j,n)
         v2(im,j,n)=v2(2,j,n)
         h2(1,j,n)=h2(im-1,j,n)
         h2(im,j,n)=h2(2,j,n)
         u1(1,j,n)=u1(im-1,j,n)
         u1(im,j,n)=u1(2,j,n)
         v1(1,j,n)=v1(im-1,j,n)
         v1(im,j,n)=v1(2,j,n)
         h1(1,j,n)=h1(im-1,j,n)
         h1(im,j,n)=h1(2,j,n)
         tmp(1,j,n)=tmp(im-1,j,n)
         tmp(im,j,n)=tmp(2,j,n)
 4680   continue
        end if
	return
	end

!c==============
	subroutine zerot(t)
	parameter(im=720,jm=121)
	real t(im,jm)
	do j=1,jm
	do i=1,im
	t(i,j)=0.
	enddo
	enddo
	return
	end

!c==============
	subroutine tempm(tmpc,tmp,m,nspd)
	parameter(im=720,jm=121)
	real tmpc(im,jm),tmp(im,jm,2)
	do j=1,jm
	do i=1,im
	tmpc(i,j)=tmpc(i,j)+tmp(i,j,m)/float(nspd)
	enddo
	enddo
	return
	end

!c===============
	subroutine getdt
	include 'pram1.inc'
	include 'com1.inc'
        if (dt.ne.0.) go to 800
        cmin=amin1(c3(1,1),c3(1,jm))
        hmax=0.
        dxmin=dxh(1)
        dymin=dyh(1)
        do 700 j=1,jm
        dymin=amin1(dyh(j),dymin)
        do 700 i=1,im
        if (j.eq.1) dxmin=amin1(dxh(i),dxmin)
        hmax=amax1(topo(i,j),hmax)
700     continue
!c
!c===> Redgrav=4.17 cm/s2
!c
!c===> This corresponds to C=250 (cm/s)
!c
        sws=sqrt(redgrav*hmax)
        dsq=1./(cmin*dxmin)**2+1./dymin**2
        cfl=1./(2.*sqrt(redgrav*hmax*dsq))
        if (nf.eq.0) then
        cor1=abs(2.*omega*sin(yz(1)*degrad))
        cor2=abs(2.*omega*sin(yz(jm)*degrad))
        else
        cor1=abs(2.*omega*sin(theta*degrad))
        cor2=abs(2.*omega*sin(theta*degrad))
        end if
        cor=1./amax1(cor1,cor2)
        percnt=.8
        dt=percnt*amin1(cfl,cor)
        write (6,9000) dt,cfl,cor,cmin,dxmin,dymin,hmax,sws
9000    format(1x,'dt,cfl,cor,cmin,dxmin,dymin,hmax,sws='&
             /1x,8e11.4)
        write(6,9001) nbuoy,nbounx,nbouny
9001    format (1x,'nbuoy=',i4,1x,'nbounx=',i4,1x,'nbouny=',i4,1x)
800     continue
	return
	end

!c=====================
!	subroutine pramout
!	include 'pram1.inc'
!	include 'com1.inc'
!        open(11,file='parameter')
!        write(11,*) '--------------------------------------------------'
!        write(11,*) '******          Model Parameters            ******'
!        write(11,*) '--------------------------------------------------'
!        write(11,*) ' '
!        write(11,*) '1.  Tr  Reference temperature           : ',tmn
!        write(11,*) '2.  Cd  Drad coeff.                     : ',cdt0
!        write(11,*) '3.  Cl  Moisture transfer coeff.        : ',cd
!        write(11,*) '4.  R   Solar rad. penetration coeff.   : ',r
!        write(11,*) '5.  gam Solar rad. attenuation          : ',gamma2
!        write(11,*) '6.  h0  Depth of constant shear         : ',h0
!        write(11,*) '7.  dh  Thickness of entrainment layer  : ',wdh
!        write(11,*) '8.  ms  Wind stirring coefficient       : ',wm,wm30
!        write(11,*) '9.  mb  Convective mixing coeff.        : ',1-tn
!        write(11,*) '10. mu  Horizontal laplasion diffusion  : ',cdfv*&
!            visc(1,1)
!        write(11,*) '11. rh  damping coeff. in h-eq.         : ',feinde
!        write(11,*) '12. r   rayleigh friction coeff.        : ',frinde
!        write(11,*) '--------------------------------------------------'
!        write(11,*) ' Models :'
!        write(11,*) 'Model A : ',modela, '  Model B : ',modelb
!        write(11,*) 'Model C : ',modelc, '  Model D : ',modeld
!        write(11,*) 'Model E : ',modele, '  Model F : ',modelf 
!        write(11,*) '--------------------------------------------------'
!        write(11,*) ' Remark :'
!        write(11,*) '1. B.C Standard : '
!        write(11,*) '   tfn zonal newtonian damping time scale  : ',tfnz
!        write(11,*) '   tfn medir newtonian damping time scale  : ',tfnm
!        write(11,*) '   tfr rayleigh friction time scale  : ',tfr
!        write(11,*) '   cdf diffusion term in h equation  : ',cdfe
!        write(11,*) '   cd1 diffusion term in h1 equation : ',cdfe1
!        write(11,*) '2. Heat flux calculation method      : ',nqin
!        write(11,*) '3. Heat flux correction              : ',fninde
!        write(11,*) '4. Friction at ML base               : ',etxy
!        close(11)
!	return
!	end
!c------------------------------------------------------------------c
!c-------read in wind, calculate surface heat fluxes here!----------c
!c------------------------------------------------------------------c

	subroutine getforce(wfu)
        USE mo_doctor,     ONLY:nout
	include 'pram1.inc'
	include 'com1.inc'
	dimension wfu(im,jm),tcp(im,jm)
!c------------- coupling with AGCM------------------------------------
	common /scouple/tua(im,jm),tva(im,jm),wfua(im,jm),qsola(im,jm),&
         hfluxa(im,jm),sst_nudg(im,jm)
!c===================================================================
	data hr/1700./,rr/0.35/
!c  turn on annual cycle after ten-year
!c	mitt=1
!c	if(itt.gt.nitt) mitt=itt-nitt
	mitt=itt
!c----------------------------------------------------------------
        day=float(mitt)/tspd
        yr=day/360
        rmth=yr*12.
        mth=int(rmth)
	nday=day+0.01

        if (day.le.90.and.nsmthstr.eq.1)then
         fgrad=day/90.
        else
         fgrad=1.
        endif
	if(itt.eq.1)then
	write(nout,*)'itt=',itt,' tspd=',tspd,' fgrad=',fgrad
	write(nout,*)'nday=',nday,' nsmthstr=',nsmthstr,' dt=',dt
	endif
!c
!c=== Set monthly mean climatological winds in terms of Fourier series
!c
      if(ncpl.eq.0) then    !!! uncoupled mode (see pram1.inc)
        do 20 j=1,jm
        do 20 i=1,im
         txa(i,j)=0.
         tya(i,j)=0.
        do 21 nhm=1,nharms
         arg=2.*pi*float(nhm)*(dt*float(itt)*lcktim-phlag)/timep
         txa(i,j)=txa(i,j)+ctx(i,j,nhm)*cos(arg)+stx(i,j,nhm)*sin(arg)
         tya(i,j)=tya(i,j)+cty(i,j,nhm)*cos(arg)+sty(i,j,nhm)*sin(arg) 
21      continue
20      continue
!c
!c===> Start forcing gradually, Sept. 16, 1991
!c
        do 23 j=1,jm
        do 23 i=1,im
         tu(i,j)=txm(i,j)+forcev*txa(i,j)*fgrad
         tv(i,j)=tym(i,j)+forcev*tya(i,j)*fgrad
23      continue
!c
!c================Set  Q0, Qsw==================================
!c     use observed solar radiation
!c
        if(nqin.eq.2) then  !!! use observed heat flux
         do 25 j=1,jm
         do 25 i=1,im
!c         qina(i,j)=0.
	  clda(i,j)=0.
          qouta(i,j)=0.
          do 24 nhm=1,nharms
           arg=2.*pi*float(nhm)*(dt*float(itt)*lcktim-phlag)/timep
!c           qina(i,j)=qina(i,j)+cqin(i,j,nhm)*cos(arg)
!c     1              +sqin(i,j,nhm)*sin(arg) 
           clda(i,j)=clda(i,j)+ccld(i,j,nhm)*cos(arg)&
                   +scld(i,j,nhm)*sin(arg) 

           qouta(i,j)=qouta(i,j)+cqout(i,j,nhm)*cos(arg)&
                    +sqout(i,j,nhm)*sin(arg) 
24        continue
25       continue
!c
!c===> Start forcing gradually, Sept. 16, 1991
!c
	
        do 26 j=1,jm
        do 26 i=1,im
!c        qin(i,j)=qinm(i,j)+forcef*qina(i,j)*fgrad
         cld(i,j)=cldm(i,j)+forcef*clda(i,j)*fgrad
!c	 cld(i,j)=cld(i,j)*0.1  !maximum is 8
!c==========================================================
!c  or from 8 to 10
         cld(i,j)=cld(i,j)*10./8.*0.1
!c----------------------------------------------------------
         qout(i,j)=qoutm(i,j)+forcef*qouta(i,j)*fgrad
26      continue
       end if
!c
!c===>Set  SST
!c
        do 30 j=1,jm
        do 30 i=1,im
         ssta(i,j)=0.
	 taira(i,j)=0.

!cfu SST and AT data both start on Jan 15, so phase lag is -0.5 month
	
        do 31 nhm=1,nharms
         arg=2.*pi*float(nhm)*(dt*float(mitt)*lcktim-phlag)/timep
         ssta(i,j)=ssta(i,j)+csst(i,j,nhm)*cos(arg)&
                +ssst(i,j,nhm)*sin(arg) 
         taira(i,j)=taira(i,j)+cta(i,j,nhm)*cos(arg)&
                +sta(i,j,nhm)*sin(arg) 
	
31      continue
         sst(i,j)=sstm(i,j)+forcef*ssta(i,j)*fgrad
         tair(i,j)=tairm(i,j)+forcef*taira(i,j)*fgrad
30      continue
!c
!c===TIM=====> Start forcing gradually, Sept. 16, 1991
!c	imth=mod(mth,12)	!!!
!c	imth1=imth
!c	imth2=imth+1
!c	if(imth.eq.0) imth1=12
!c	cldt=rmth-mth
!c
!c	if(lcktim.eq.0) then
!c	imth1=imthlck
!c	imth2=imthlck
!c	cldt=1.
!c	end if
!ctim  !calculate heat flux
!c        do 32 j=1,jm
!c         if(nqin.ge.2.and.forcef.eq.1.)then  
!c         q0s(j)=0.
!cfu     qo,cl,h data start on Jan 15, so phase lag is -0.5 month
!c            (11/6/97 try)
!c        do 130 nhm=1,nharms
!c         arg=2.*pi*float(nhm)*(dt*float(mitt)*lcktim-phlag)/timep
!c         q0s(j)=q0s(j)+csol(j,nhm)*cos(arg)
!c     1            +ssol(j,nhm)*sin(arg) 
!c130      continue
!c       q0s(j)=(1.-cldt)*q0sol(j,imth1)+cldt*q0sol(j,imth2)
!c	q0s(j)=soli(j)+q0s(j)               !plus mean
!c	 endif
!c	do 32 i=1,im      !Coads cloud	
!c        cld(i,j)=(1.-cldt)*cloud(i,j,imth1)+cldt*cloud(i,j,imth2)
!c32      continue

!c-----------------------------------------------------------------
        do 35 j=1,jm
        do 35 i=1,im
         tcp(i,j)=sst(i,j)
         if(be(i,j)) tcp(i,j)=tmp(i,j,n) !model sst
35      continue

!cc=================================================================
!c===>Calculate wind stress, Latent/sensible heat from u,v,tmp
!c
	do 40 j=1,jm
	do 40 i=1,im
	 wfu(i,j)=4.                  !min
	if(be(i,j)) then
         uuvv=sqrt(tu(i,j)**2+tv(i,j)**2)
         if(uuvv.lt.4.) uuvv=4.        !Umin(m/s)   	
         wfu(i,j)=uuvv         

         tx(i,j)=rhoa*cdt0*tu(i,j)*uuvv*10000.
         ty(i,j)=rhoa*cdt0*tv(i,j)*uuvv*10000.

	 if(nqin.ge.2) then
	  tsts=tcp(i,j)
	  eses=6.11*10**(7.5*tsts/(tsts+237.3))
	  qsqs=0.622*eses/(1013.-0.378*eses)
!cZF===> Use COADS surface air temperature
  	  qaqa=(0.972*(1.32+tair(i,j))/1.03-8.92)*0.001
!c	  qaqa=(0.972*tsts-8.92)*0.001
	  wvp=1013./(1.+0.622/qaqa)
	  tabs=tair(i,j)+273.15
	  tsa=tsts-tair(i,j)
!c	  tata=1.03*tsts-1.32    !use model SST
!c	  tabs=tata+273.15
!c	  tsa=tsts-tata
!c	  wvp=1013./(0.378+0.622/qaqa)
	  if(wvp.le.0) wvp=0
!c
!c---fu (adjust cloud fraction if tsts>29.5)
!c 	  if(tsts.gt.29.5) then
!c          clde=0.125*(tsts-29.5)
!c	  cld(i,j)=cld(i,j)+clde
!c	  if(cld(i,j).ge.1.0) cld(i,j)=1.0
!c	  end if
!c==========================================================>	
!c
!cEND
!c===> Qsw: Solar radiation
!c not  with observed solar radiation
!c
!c	  qout(i,j)=q0s(j)*(1.-(asol(j)+0.38*cld(i,j))*
!c     1              cld(i,j))*q0unit*1000.
!cfu use q0s harmonics(q0unit has been changed accordingly)
!c===> Qlw: Longwave radiation
!c
	  qlw(i,j)=0.97*5.673e-8*tabs*tabs*tabs*(tabs*(0.39-&
               0.05*sqrt(wvp))*(1.-cqlw(j)*cld(i,j)*cld(i,j))&
               +4.*tsa)*1000.
!c
!c===> Ql: Latent heat flux
!c
	  qlh(i,j)=1.2*cd*uuvv*2.5*1000000.*(qsqs-qaqa)*1000.
!c         (directly from atmospheric model 10/5/97)
	  if(qlh(i,j).lt.0.) qlh(i,j)=0.
!c
!c===> Qs: Sensible heat flux
!c
	  qsh(i,j)=1.2*cd*uuvv*1004.*tsa*1000.
	  if(qsh(i,j).lt.0.) qsh(i,j)=0.
	endif
	endif
40	continue

	do 50 j=1,jm
	do 50 i=1,im
!c
!c===> Q0: Total heat flux
!c
	  qin(i,j)=qout(i,j)-qlw(i,j)-qlh(i,j)-qsh(i,j)
!c
!c===> Q_h: Penetrated heat flux
!c 
          qpet(i,j)=qout(i,j)*rr*exp(-h1(i,j,n)/hr)
!c
!c===> Net heat flux input in MLT equation
!c
          qnet(i,j)=qin(i,j)-qpet(i,j)
          qnet(i,j)=qnet(i,j)*(1.-exp(yh(1)-yh(j))-exp(yh(j)-yh(jm)))
!c
!c===> Thermal forcing components for out put, just to make sure
!c     your ocean is not going to be boiling or be icing :)
!c
          tfc1(i,j)=qin(i,j)/1000.           !! total heat flux
          tfc2(i,j)=qout(i,j)/1000.          !! solar radiation
          tfc3(i,j)=qlw(i,j)/1000.           !! long wave radiation
          tfc4(i,j)=qlh(i,j)/1000.           !! latent heat flux
          tfc5(i,j)=qsh(i,j)/1000.           !! sensible heat flux
          tfc6(i,j)=qpet(i,j)/1000.          !!penetrate solar radiation
50	  continue
	goto 77      
	end if
!c===================================================================
!c     COUPLE to AGCM
!c-------------------------------------------------------------------
	if(ncpl.eq.1) then
!cc==================================================================
!c
!fu++ using observed ECMWF surface wind 
!        do 720 j=1,jm
!        do 720 i=1,im
!         txa(i,j)=0.
!         tya(i,j)=0.
!         do 721 nhm=1,nharms
!         arg=2.*pi*float(nhm)*(dt*float(itt)*lcktim-phlag)/timep
!         txa(i,j)=txa(i,j)+ctx(i,j,nhm)*cos(arg)+stx(i,j,nhm)*sin(arg)
!         tya(i,j)=tya(i,j)+cty(i,j,nhm)*cos(arg)+sty(i,j,nhm)*sin(arg) 
!721      continue
!720      continue
!c
!c===> Start forcing gradually, Sept. 16, 1991
!c
!        do 723 j=1,jm
!        do 723 i=1,im
!         tu(i,j)=txm(i,j)+forcev*txa(i,j)*fgrad
!         tv(i,j)=tym(i,j)+forcev*tya(i,j)*fgrad
!723      continue
!c
!c==============fu++========================================================
	do 121 j=1,jm
	do 121 i=1,im
         ssta(i,j)=0.
!cfu SST and AT data both start on Jan 15, so phase lag is -0.5 month
	
        do 131 nhm=1,nharms
         arg=2.*pi*float(nhm)*(dt*float(mitt)*lcktim-phlag)/timep
         ssta(i,j)=ssta(i,j)+csst(i,j,nhm)*cos(arg)&
                 +ssst(i,j,nhm)*sin(arg) 
	
131      continue
         sst(i,j)=sstm(i,j)+ssta(i,j)
121     continue

	do 140 j=1,jm
	do 140 i=1,im
!	wfu(i,j)=4.
	if(be(i,j)) then
!c===>Set  SST
!c
!         ssta(i,j)=0.
!cfu SST and AT data both start on Jan 15, so phase lag is -0.5 month
!        do 131 nhm=1,nharms
!         arg=2.*pi*float(nhm)*(dt*float(mitt)*lcktim-phlag)/timep
!         ssta(i,j)=ssta(i,j)+csst(i,j,nhm)*cos(arg)&
!                 +ssst(i,j,nhm)*sin(arg) 
!131      continue
!         sst(i,j)=sstm(i,j)+ssta(i,j)

!c=======>AGCM gives wind stress directly
	  tx(i,j)=tua(i,j)
	  ty(i,j)=tva(i,j)
	  wfu(i,j)=wfua(i,j)
!c===>Calculate wind stress from ECMWF surface wind 
!	uuvv=sqrt(tu(i,j)*tu(i,j)+tv(i,j)*tv(i,j))
!	if(uuvv.lt.4.) uuvv=4.
!        tx(i,j)=rhoa*cdt0*tu(i,j)*uuvv*10000.
!        ty(i,j)=rhoa*cdt0*tv(i,j)*uuvv*10000.
!c==========================================================>	
!c===>  solar radiation

	  qout(i,j)=qsola(i,j)*1000.
!c===> Q0: Total heat flux
!c
	  qin(i,j)=hfluxa(i,j)*1000.
!c
!c===> Q_h: Penetrated heat flux
!c 
          qpet(i,j)=qout(i,j)*rr*exp(-h1(i,j,n)/hr)
!c
!c===> Net heat flux input in MLT equation
!c
          qnet(i,j)=qin(i,j)-qpet(i,j)
          qnet(i,j)=qnet(i,j)*(1.-exp(yh(1)-yh(j))-exp(yh(j)-yh(jm)))
!c
!c===> Thermal forcing components for out put, just to make sure
!c     your ocean is not going to be boiling or be icing :)
!c
          tfc1(i,j)=qin(i,j)/1000.           !! total heat flux
          tfc2(i,j)=qout(i,j)/1000.          !! solar radiation
!c          tfc3(i,j)=qlw(i,j)/1000.           !! long wave radiation
!c          tfc4(i,j)=qlh(i,j)/1000.           !! latent heat flux
!c          tfc5(i,j)=qsh(i,j)/1000.           !! sensible heat flux
          tfc6(i,j)=qpet(i,j)/1000.          !!penetrate solar radiation
	
	end if
140	  continue

	end if

!c================END====================================================

77	return
	end

!c----------------------------------------------------------------------
!c                             END OF THE IOM CODE
!c======================================================================
!c
!c===>	On line document of parameterizition of Te	
!c       Pray for answer which one is the best : )
!c
!c      do 2150 j=1,jm
!c      do 2150 i=1,im
!c      if(be(i,j)) then
!c
!c**1** computer relation between e and te  (Seager, 1988)
!c
!c      hee=e(i,j,n)*0.001
!c      h20=q0+q1*hee+q2*hee**2+q3*hee**3+q4*hee**4+q5*hee**5
!c     1        +q6*hee**6
!c      if(hee.lt.-6.3) h20=0.2
!c      if(hee.gt.10000.) h20=19000.
!c      te(i,j)=p0+p1*h20+p2*h20**2+p3*h20**3
!c     1        +p4*h20**4+p5*h20**5+p6*h20**6
!c     1        +p7*h20**7+p8*h20**8+p9*h20**9
!c      if (h20.lt.0.3) te(i,j)=17.00
!c      if (h20.gt.4.) te(i,j)=29.80
!c
!c
!c**2** using the observed subsurface temp as entrained temp
!c         (Chang Ping, 1992)
!c
!c      te(i,j)=(subtemp(i,j)+
!c     1         deltat*dtdz(i,j)*(h(i,j)-amh(i,j)))
!c
!c
!c**3** Zebiak and Cane (1987)
!c
!c	if(nseason.eq.0.or.deltat.eq.0.) then
!c	tdlt=0.
!c	else
!c	hdlt=h(i,j)-amh(i,j)
!c	if(hdlt.ge.0.) then
!c	tdlt=t1o*(tanh((amh(i,j)+hdlt)/b1o)-tanh(amh(i,j)/b1o))
!c	else
!c	tdlt=t2o*(tanh((amh(i,j)-hdlt)/b2o)-tanh(amh(i,j)/b2o))
!c	end if
!c	if(tdlt.gt.5.) tdlt=5.
!c	if(tdlt.lt.-5.) tdlt=-5.
!c	end if
!c	tsub=subtemp(i,j)+tdlt
!c	te(i,j)=gmo*tsub+(1.-gmo)*tmp(i,j,n)
!c	if(te(i,j).gt.tmp(i,j,n)) te(i,j)=tmp(i,j,n)
!c
!c**4** Wang, 1992
!c
!c      te(i,j)=tmp(i,j,n)-wgm*wdh*(tmp(i,j,n)-tmn)/(h(i,j,n)-h1(i,j,n))
!c      if((tmp(i,j,n)-te(i,j)).ge.3.0) te(i,j)=tmp(i,j,n)-3.0
!c      if(te(i,j).lt.tmn) te(i,j)=tmn
!c      end if
!c2150  continue
!c
!c     define buoyancy b(i,j)=thmexp*grav*(tmp-tmn)
!c     where thmexp is thermal expansion coef., tmn is
!c     reference temperature in abyss layer.    
!c
!c
!c
!c-------------------------------------------------------------
!c     calculate entrainment velocity we(i,j)
!c
!c     netr = 0, no entrainment     ==> we(i,j)=0
!c     netr = 1, const. mixed layer ==> we(i,j)=d(hu)/dx+d(hv)/dy
!c     netr = 2, Kraus and Turner mixed model ==>
!c     netr = 3, Csanady mixed model ==>
!c     netr = 4, Cane and Zebiak model ==>
!c     
!c     in Kraus and Turner model
!c
!c     we=2*wm*uf**3/b/h+tn*(qout-ef1*qin)/rhoo/cw/(tmp-tmn)
!c    
!c     in Csanady model
!c
!c     we=0.12*uf**2/ri/(1+0.12/ri)
!c     
!c     in Cane and Zebiak
!c
!c     we=h1*(du1/dx+dv1/dy)
!c
!c     where
!c   
!c---------------------------------------------------------------------
!c     in Kraus and Turner model
!c     ef1=1+e1-2*e2/h; 
!c     e1=r*exp(-gamma1*h)+(1-r)*exp(-gamma2*h);
!c     e2=r*(1-exp(-gamma1*h))/gamma1+(1-r)*(1-exp(-gamma2*h))/gamma2;
!c     
!c     uf = friction velocity:
!c     uf=(rhoa*cd/rhoo)**0.5*wsp, or
!c     uf=((tx**2+ty**2)**0.5/rhoo)**0.5
!c
!c     b = bouyancy:  b = thmexp*(tmp-tmn)
!c
!c     qout=sensible + latent + longwave fluxes (eng/cm**2=1.e-3(w/m**2))
!c     qin  = insolation flux (eng/cm**2=g/s**3=1e-3(kg/s**3=w/m**2)
!c     wsp = wind speed (cm/s)
!c     tx and ty = wind stress in x and y direction (dyn/cm**2)
!c  
!c----------------------------------------------------------------------
!c     constants:
!c
!c     thmexp: thermal expansion coefficient: 2.0e-4 (1/k)
!c     rhoo: density of ocean water: 1.0 (g/cm**3)
!c     rhoa: density of air: 1.2e-3 (g/cm**3)
!c     cd: drag coefficient: 1.3e-3
!c     cw: heat capacity of water: 4.2e7 (erg/g/k)
!c     wm: wind mixing entrainment calibration coef.: 1.25
!c     tn: thermal mixing entrainment calibration coef.: 0.3
!c     r: solar radiation penetration coefficient: 0.58
!c     gamma1: e-folding depth of solar radiation: 1/35 (1/cm)
!c     gamma2: e-folding depth of solar radiation: 1/2300 (1/cm)    
!c     tmn: temperature in abyssal layer: 14.0 (c)
!c===>namelist /input/ itt,ndays,wface,sface,tfr,tfnz,tfnm,nergy,nsamp
!c    1,nsmth,maps,nrstrt,grav,epsln,alpha,eqdpth
!c
!c
!c    nergy = approx number of days between energy integrals
!c    ndays = approx number of days to integrate
!c    ansamp= number of days between instantaneous samples
!c    nsmth = number of time steps between time smoothing
!c    itt   = (0,n)=start from (ic,restart)
!c    maps  = (0,1)=(no maps,maps) at end of integration
!c    nrstrt= (0,1)=(do not write,write) restart at end
!c    dt    = time step in sec.if dt=0 then dt is calculated based on
!c            cfl & coriolis restrictions
!c    wface = longitude (in degrees) of the western edge of the
!c            first free surface box
!c    sface = latitude (in degrees) of the southern edge of the
!c            first free surface box
!c      tfr = time scale (days) for rayleigh friction
!c      tfn = time scale (days) for newtonian damping
!c    alpha = coefficient of timewise smoothing
!c    eqdpth= mixed depth depth in centimeters
!c    degnon= degree of nonlinearity (measured by the expected
!c            height responce at the source compared with the equivalent
!c            depth)
!c    pforc = period of forcing in days.
!c    nf    = 0 model is on sphere, otherwise model is on f-plane
!c    netr  = set for different entrainment
!c    nbuoy = (0,others) if 0, temperature is a passive tracer
!c    nbounx= (0,1) if 1, set cyclic boundary condition in x-dir.
!c    nbouny= (0,1) if 1, set cyclic boundary condition in y-dir.
!c    forcef= between 0 and 1 which contrals the amplitude of forcing.
!c            if 0, only annual mean forcing is used.
!c_______________________________________________________________________

	subroutine smoth2(be,boun,t,p,m,n,s)
!c---------just smoth the coastline region!!!
	logical be(m,n)
	real t(m,n),p(m,n),boun(m,n),ss(720,121)
	c1=1.-s*s
	c2=s*(1.-s)/2.
	c3=s*s/4.
	do 5 j=1,n
	do 5 i=1,m
	ss(i,j)=0.
	p(i,j)=t(i,j)
5       continue

	do 10 j=1,n
	do 10 i=1,m

	if(boun(i,j).eq.1.0) then
	
	ia=max(i-1,1)
	ib=min(i+1,m)
	ja=max(j-1,1)
	jb=min(j+1,n)

!c--consider the land boundary
	nn=0
	do 20 l=ia,ib
	do 20 k=ja,jb
	ss(l,k)=0.
	if(be(l,k)) then          !ocean domain
	nn=nn+1
	ss(l,k)=1.               !weight every point
	end if
20      continue
	if(nn.eq.9) goto 99

	p(i,j)=0	
	do 30 l=ia,ib
	do 30 k=ja,jb
	p(i,j)=p(i,j)+t(l,k)*ss(l,k)/float(nn)
30      continue
	goto 98
99	p(i,j)=t(i,j)*c1+c2*(t(ia,j)+t(ib,j)+t(i,ja)&
     		+t(i,jb)-4.*t(i,j))+c3*(t(ia,ja)+t(ia,jb)&
     		+t(ib,ja)+t(ib,jb))
98      continue
	end if

10	continue

	return
	end

!c------------------------------------------------------------------c
!c     end of the IOM code (3/17/99)
!c------------------------------------------------------------------c

