!c====this code is designed to call the IOM(as a
!c==subroutines, then accumulate the inputs for the
!c==ocean model
!c== nrday:days in a month
!c== imday:total days in a specific month

subroutine ocean(sstfu,taux,tauy,sflux,qsol,wind,ia,ja,fu_sst)

  USE mo_landsea,  ONLY: aland, bzone, itopo, iland,&
       agx, agy, sgx, sgy
  USE mo_couple,   ONLY: restart_sst    !fu++

  parameter(io=720,jo=121)

  !c-------------------------------------------------------
  !c  atmospheric grids are on the Gaussian grids
  !c-------------------------------------------------------
  !c  oceanic grids are on the 2x1 system
  !c========================================================
  !	integer aland,bzone	
  !	common /csbc1/aland(ia,ja),itopo(io,jo),&
  !           iland(ia,ja),agx(ia),agy(ja),sgx(io),sgy(jo)
  common /scouple/tx(io,jo),ty(io,jo),wfu(io,jo),qsl(io,jo),&
       flux(io,jo),sst_nudg(io,jo)

  common /fudate/runday,dta,dto,day1,tseg,numcpl,numatm,numocn
  dimension sst(io,jo),sstfu(ia,ja),&
       taux(ia,ja),tauy(ia,ja),sflux(ia,ja),qsol(ia,ja),&
       ssti(io,jo),wind(ia,ja),fu_sst(ia,ja)

  !c=== prepare the oceanic inputs (taux,tauy,qsol,sflux,)
  !c=== from atmospheric model
  !c*********************************************************
  !c    BECAUSE THE ATMOSPHERIC MODEL STARTS FIRST!!!!!
  !c*********************************************************
  !c--get ocean grid inputs from atmospheric model

  call winda(taux,tauy,qsol,sflux,wind,ia,ja,fu_sst) 

  ssti(1:io,1:jo)=0.0
  
  !c--run IOM
  do loop = 1,numocn

     call iom(sst,loop) 

     ! accumulates SST field

     ssti(1:io,1:jo) = ssti(1:io,1:jo) + sst(1:io,1:jo)

  end do

  ! average SST field
  ssti(1:io,1:jo) = ssti(1:io,1:jo)/float(numocn)

  !c========================================================
  !c     get sstfu after one-day ocean integration
  !c-------------------------------------------------------
  call getsst(sstfu,ssti,ia,ja)

  restart_sst  = sstfu    !fu++

  return	
end subroutine ocean

