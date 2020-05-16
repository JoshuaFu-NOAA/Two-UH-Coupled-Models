	program stratiform

!   test the impacts of stratiform rainfall on
!   heating/moistening profiles
! 
	parameter (kk=10)

	REAL p(kk),sheat(kk)

	data p/950,900,850,800,700,600,500,400,&
             300,250/ 

	alf=0.5
	conv=0.0003           !K/day

	ps=950             !mb
	pt=250              !mb
	pi=3.1415926

	pstar=(ps-pt)/2.

	print*,'pstar=',pstar

	do k=1,kk

	if((p(k).ge.pstar).and.(p(k).le.ps)) then
	aa=2.*pi*(p(k)-ps)/(ps-pt)  !lower level
	else
	aa=2.*pi*(p(k)-ps)/(ps-pt)  !upper level
	endif
	
	sheat(k)=alf*conv*sin(aa)

	print*,'k,p,heat=',k,p(k),sheat(k)

	end do

	print*,'maxheat=',MAXVAL(sheat)

	STOP
	END
	

