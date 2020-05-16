!c===================================================================c 
!c  modification for pgi compiler
!c  Qing Bao . Oct. 20, 2005 
!c===================================================================c

subroutine restart_setupo
  
   USE mo_couple, ONLY:  restart_itt
        common /scalar/ grav,omega,itt,ndays,tfr,tfnz,tfnm,epsln,dt,&
                       degrad,nergy,ansamp,alpha,nrstrt,maps,nsmth,nf,&
                       netr,tmn,tmpmax,thmexp,rhoa,cd,cw,wm,tn,r,&
                       gamma1,gamma2,nbuoy,nbounx,nbouny,forcef,forcev,&
                       deltat,nseason,nsmthstr,hepsln,nbeta,beta,&
                       ittend,tspd,txyk1,txyk,cdfv,ufcoef,&
                       q0unit,earthr,daysave,ittwrt,ittatm,phlag,&
                       lcktim,imthlck,day,nday,iyear,tcount,tadd0,&
                       nladv,cdfv1,wdh,theta,timep,nqin,ncouple,&
                       iloop
      
   !c******************************************************************** 
!  include 'pram1.inc'
!  include 'com1.inc'
   !c******************************************************************** 

     itt = restart_itt
end subroutine restart_setupo
